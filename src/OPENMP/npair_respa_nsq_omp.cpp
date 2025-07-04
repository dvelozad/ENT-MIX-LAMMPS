// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair_respa_nsq_omp.h"
#include "npair_omp.h"
#include "omp_compat.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int NEWTON, int TRI>
NPairRespaNsqOmp<NEWTON, TRI>::NPairRespaNsqOmp(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   multiple respa lists
   Newtoff
     N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
     pair added to list if atoms i and j are both owned and i < j
     pair added if j is ghost (also stored by proc owning j)
  Newton
     N^2 / 2 search for neighbor pairs with full Newton's 3rd law
     pair added to list if atoms i and j are both owned and i < j
     if j is ghost only me or other proc adds pair
     decision based on itag,jtag tests
     use itag/jtag comparision to eliminate half the interactions
     itag = jtag is possible for long cutoffs that include images of self
  Newton + Triclinic:
     for triclinic, must use delta to eliminate half the I/J interactions
     cannot use I/J exact coord comparision as for orthog
     b/c transforming orthog -> lambda -> orthog for ghost atoms
     with an added PBC offset can shift all 3 coords by epsilon
------------------------------------------------------------------------- */

template<int NEWTON, int TRI>
void NPairRespaNsqOmp<NEWTON, TRI>::build(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int bitmask = (includegroup) ? group->bitmask[includegroup] : 0;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == Atom::TEMPLATE) ? 1 : 0;
  const double delta = 0.01 * force->angstrom;

  NPAIR_OMP_INIT;

  const int respamiddle = list->respamiddle;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list)
#endif
  NPAIR_OMP_SETUP(nlocal);

  int i, j, n, itype, jtype, n_inner, n_middle, imol, iatom;
  tagint itag, jtag, tagprev;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *neighptr, *neighptr_inner, *neighptr_middle;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nall = atom->nlocal + atom->nghost;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int *ilist_inner = list->ilist_inner;
  int *numneigh_inner = list->numneigh_inner;
  int **firstneigh_inner = list->firstneigh_inner;

  int *ilist_middle, *numneigh_middle, **firstneigh_middle;
  if (respamiddle) {
    ilist_middle = list->ilist_middle;
    numneigh_middle = list->numneigh_middle;
    firstneigh_middle = list->firstneigh_middle;
  }

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  MyPage<int> &ipage_inner = list->ipage_inner[tid];
  ipage.reset();
  ipage_inner.reset();

  MyPage<int> *ipage_middle;
  if (respamiddle) {
    ipage_middle = list->ipage_middle + tid;
    ipage_middle->reset();
  }

  int which = 0;
  int minchange = 0;

  for (i = ifrom; i < ito; i++) {

    n = n_inner = 0;
    neighptr = ipage.vget();
    neighptr_inner = ipage_inner.vget();
    if (respamiddle) {
      n_middle = 0;
      neighptr_middle = ipage_middle->vget();
    }

    itype = type[i];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    // loop over remaining atoms, owned and ghost

    for (j = i + 1; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;

      if (NEWTON) {
        if (j >= nlocal) {
          jtag = tag[j];
          if (itag > jtag) {
            if ((itag + jtag) % 2 == 0) continue;
          } else if (itag < jtag) {
            if ((itag + jtag) % 2 == 1) continue;
          } else if (TRI) {
            if (fabs(x[j][2] - ztmp) > delta) {
              if (x[j][2] < ztmp) continue;
            } else if (fabs(x[j][1] - ytmp) > delta) {
              if (x[j][1] < ytmp) continue;
            } else {
              if (x[j][0] < xtmp) continue;
            }
          } else {
            if (x[j][2] < ztmp) continue;
            if (x[j][2] == ztmp) {
              if (x[j][1] < ytmp) continue;
              if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
            }
          }
        }
      }

      jtype = type[j];
      if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq <= cutneighsq[itype][jtype]) {
        if (molecular != Atom::ATOMIC) {
          if (!moltemplate)
            which = find_special(special[i], nspecial[i], tag[j]);
          else if ((imol >= 0) && onemols[imol]->special)
            which = find_special(onemols[imol]->special[iatom], onemols[imol]->nspecial[iatom],
                                 tag[j] - tagprev);
          else
            which = 0;
          if (which == 0)
            neighptr[n++] = j;
          else if ((minchange = domain->minimum_image_check(delx, dely, delz)))
            neighptr[n++] = j;
          else if (which > 0)
            neighptr[n++] = j ^ (which << SBBITS);
        } else
          neighptr[n++] = j;

        if (rsq < cut_inner_sq) {
          if (which == 0)
            neighptr_inner[n_inner++] = j;
          else if (minchange)
            neighptr_inner[n_inner++] = j;
          else if (which > 0)
            neighptr_inner[n_inner++] = j ^ (which << SBBITS);
        }

        if (respamiddle && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
          if (which == 0)
            neighptr_middle[n_middle++] = j;
          else if (minchange)
            neighptr_middle[n_middle++] = j;
          else if (which > 0)
            neighptr_middle[n_middle++] = j ^ (which << SBBITS);
        }
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

    ilist_inner[i] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    ipage.vgot(n_inner);
    if (ipage_inner.status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[i] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      ipage_middle->vgot(n_middle);
      if (ipage_middle->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
    }
  }
  NPAIR_OMP_CLOSE;
  list->inum = nlocal;
  list->inum_inner = nlocal;
  if (respamiddle) list->inum_middle = nlocal;
}

namespace LAMMPS_NS {
template class NPairRespaNsqOmp<0,0>;
template class NPairRespaNsqOmp<1,0>;
template class NPairRespaNsqOmp<1,1>;
}
