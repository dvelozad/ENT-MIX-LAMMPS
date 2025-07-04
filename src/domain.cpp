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

/* ----------------------------------------------------------------------
   Contributing author (triclinic) : Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "domain.h"
#include "style_region.h"   // IWYU pragma: keep

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_deform.h"
#include "force.h"
#include "kspace.h"
#include "lattice.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "output.h"
#include "region.h"
#include "thermo.h"
#include "universe.h"
#include "update.h"

#include <cstring>
#include <cmath>

using namespace LAMMPS_NS;
using namespace MathExtra;

static constexpr double BIG =   1.0e20;
static constexpr double SMALL = 1.0e-4;
static constexpr double BONDSTRETCH = 1.1;

/* ----------------------------------------------------------------------
   one instance per region style in style_region.h
------------------------------------------------------------------------- */

template <typename T> static Region *region_creator(LAMMPS *lmp, int narg, char ** arg)
{
  return new T(lmp, narg, arg);
}

/* ----------------------------------------------------------------------
   default is periodic
------------------------------------------------------------------------- */

Domain::Domain(LAMMPS *lmp) : Pointers(lmp)
{
  box_exist = 0;
  box_change = 0;
  deform_flag = deform_vremap = deform_groupbit = 0;

  dimension = 3;
  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  boundary[0][0] = boundary[0][1] = 0;
  boundary[1][0] = boundary[1][1] = 0;
  boundary[2][0] = boundary[2][1] = 0;

  minxlo = minxhi = 0.0;
  minylo = minyhi = 0.0;
  minzlo = minzhi = 0.0;

  triclinic = triclinic_general = 0;

  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  xy = xz = yz = 0.0;

  h[3] = h[4] = h[5] = 0.0;
  h_inv[3] = h_inv[4] = h_inv[5] = 0.0;
  h_rate[0] = h_rate[1] = h_rate[2] =
    h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;

  prd_lamda[0] = prd_lamda[1] = prd_lamda[2] = 1.0;
  prd_half_lamda[0] = prd_half_lamda[1] = prd_half_lamda[2] = 0.5;
  boxlo_lamda[0] = boxlo_lamda[1] = boxlo_lamda[2] = 0.0;
  boxhi_lamda[0] = boxhi_lamda[1] = boxhi_lamda[2] = 1.0;

  lattice = nullptr;
  auto args = new char*[2];
  args[0] = (char *) "none";
  args[1] = (char *) "1.0";
  set_lattice(2,args);
  delete[] args;

  copymode = 0;

  region_map = new RegionCreatorMap();

#define REGION_CLASS
#define RegionStyle(key,Class) \
  (*region_map)[#key] = &region_creator<Class>;
#include "style_region.h"   // IWYU pragma: keep

#undef RegionStyle
#undef REGION_CLASS
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  if (copymode) return;

  for (auto &reg : regions) delete reg;
  regions.clear();
  delete lattice;
  delete region_map;
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  // set box_change flags if box size/shape/sub-domains ever change
  // due to shrink-wrapping or fixes that change box size/shape/sub-domains

  box_change_size = box_change_shape = box_change_domain = 0;

  // flags for detecting, if multiple fixes try to change the
  // same box size or shape parameter

  int box_change_x=0, box_change_y=0, box_change_z=0;
  int box_change_yz=0, box_change_xz=0, box_change_xy=0;
  const auto &fixes = modify->get_fix_list();

  if (nonperiodic == 2) box_change_size = 1;
  for (const auto &fix : fixes) {
    if (fix->box_change & Fix::BOX_CHANGE_SIZE)   box_change_size = 1;
    if (fix->box_change & Fix::BOX_CHANGE_SHAPE)  box_change_shape = 1;
    if (fix->box_change & Fix::BOX_CHANGE_DOMAIN) box_change_domain = 1;
    if (fix->box_change & Fix::BOX_CHANGE_X)      box_change_x++;
    if (fix->box_change & Fix::BOX_CHANGE_Y)      box_change_y++;
    if (fix->box_change & Fix::BOX_CHANGE_Z)      box_change_z++;
    if (fix->box_change & Fix::BOX_CHANGE_YZ)     box_change_yz++;
    if (fix->box_change & Fix::BOX_CHANGE_XZ)     box_change_xz++;
    if (fix->box_change & Fix::BOX_CHANGE_XY)     box_change_xy++;
  }

  std::string mesg = "Must not have multiple fixes change box parameter ";

#define CHECK_BOX_FIX_ERROR(par)                                        \
  if (box_change_ ## par > 1) error->all(FLERR,(mesg + #par))

  CHECK_BOX_FIX_ERROR(x);
  CHECK_BOX_FIX_ERROR(y);
  CHECK_BOX_FIX_ERROR(z);
  CHECK_BOX_FIX_ERROR(yz);
  CHECK_BOX_FIX_ERROR(xz);
  CHECK_BOX_FIX_ERROR(xy);
#undef CHECK_BOX_FIX_ERROR

  box_change = 0;
  if (box_change_size || box_change_shape || box_change_domain) box_change = 1;

  // check for fix deform

  deform_flag = deform_vremap = deform_groupbit = 0;
  for (const auto &fix : fixes)
    if (utils::strmatch(fix->style,"^deform")) {
      deform_flag = 1;
      if ((dynamic_cast<FixDeform *>(fix))->remapflag == Domain::V_REMAP) {
        deform_vremap = 1;
        deform_groupbit = fix->groupbit;
      }
    }

  // region inits

  for (auto &reg : regions) reg->init();
}

/* ----------------------------------------------------------------------
   set initial global box
   assumes boxlo/hi and triclinic tilts are already set
   expandflag = 1 if need to expand box in shrink-wrapped dims
   not invoked by read_restart since box is already expanded
   if don't prevent further expansion, restarted triclinic box
     with unchanged tilt factors can become a box with atoms outside the box
------------------------------------------------------------------------- */

void Domain::set_initial_box(int expandflag)
{
  // error checks for orthogonal and triclinic domains

  if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
    error->one(FLERR,"Box bounds are invalid or missing");

  if (dimension == 2 && (xz != 0.0 || yz != 0.0))
    error->all(FLERR,"Cannot skew triclinic box in z for 2d simulation");

  // check on triclinic tilt factors

  if (triclinic) {
    if ((fabs(xy/(boxhi[1]-boxlo[1])) > 0.5 && yperiodic) ||
        ((fabs(xz)+fabs(yz))/(boxhi[2]-boxlo[2]) > 0.5 && zperiodic)) {
      if (comm->me == 0)
        error->warning(FLERR,"Triclinic box skew is large. LAMMPS will run inefficiently.");
    }
  }

  // set small based on box size and SMALL
  // this works for any unit system

  small[0] = SMALL * (boxhi[0] - boxlo[0]);
  small[1] = SMALL * (boxhi[1] - boxlo[1]);
  small[2] = SMALL * (boxhi[2] - boxlo[2]);

  // if expandflag, adjust box lo/hi for shrink-wrapped dims

  if (!expandflag) return;

  if (boundary[0][0] == 2) boxlo[0] -= small[0];
  else if (boundary[0][0] == 3) minxlo = boxlo[0];
  if (boundary[0][1] == 2) boxhi[0] += small[0];
  else if (boundary[0][1] == 3) minxhi = boxhi[0];

  if (boundary[1][0] == 2) boxlo[1] -= small[1];
  else if (boundary[1][0] == 3) minylo = boxlo[1];
  if (boundary[1][1] == 2) boxhi[1] += small[1];
  else if (boundary[1][1] == 3) minyhi = boxhi[1];

  if (boundary[2][0] == 2) boxlo[2] -= small[2];
  else if (boundary[2][0] == 3) minzlo = boxlo[2];
  if (boundary[2][1] == 2) boxhi[2] += small[2];
  else if (boundary[2][1] == 3) minzhi = boxhi[2];
}

/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_global_box()
{
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0/h[0];
  h_inv[1] = 1.0/h[1];
  h_inv[2] = 1.0/h[2];

  prd_half[0] = xprd_half = 0.5*xprd;
  prd_half[1] = yprd_half = 0.5*yprd;
  prd_half[2] = zprd_half = 0.5*zprd;

  if (triclinic) {
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);

    boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
    boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
    boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
    boxlo_bound[2] = boxlo[2];

    boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
    boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
    boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
    boxhi_bound[2] = boxhi[2];
  }

  // update general triclinic box if defined
  // reset general tri ABC edge vectors from restricted tri box

  if (triclinic_general) {
    double aprime[3],bprime[3],cprime[3];

    // A'B'C' = edge vectors of restricted triclinic box

    aprime[0] = boxhi[0] - boxlo[0];
    aprime[1] = aprime[2] = 0.0;
    bprime[0] = xy;
    bprime[1] = boxhi[1] - boxlo[1];
    bprime[2] = 0.0;
    cprime[0] = xz;
    cprime[1] = yz;
    cprime[2] = boxhi[2] - boxlo[2];

    // transform restricted A'B'C' to general triclinic ABC

    MathExtra::matvec(rotate_r2g,aprime,avec);
    MathExtra::matvec(rotate_r2g,bprime,bvec);
    MathExtra::matvec(rotate_r2g,cprime,cvec);
  }
}

/* ----------------------------------------------------------------------
   set lamda box params
   assumes global box is defined and proc assignment has been made
   uses comm->xyz_split or comm->mysplit
     to define subbox boundaries in consistent manner
------------------------------------------------------------------------- */

void Domain::set_lamda_box()
{
  if (comm->layout != Comm::LAYOUT_TILED) {
    int *myloc = comm->myloc;
    double *xsplit = comm->xsplit;
    double *ysplit = comm->ysplit;
    double *zsplit = comm->zsplit;

    sublo_lamda[0] = xsplit[myloc[0]];
    subhi_lamda[0] = xsplit[myloc[0]+1];
    sublo_lamda[1] = ysplit[myloc[1]];
    subhi_lamda[1] = ysplit[myloc[1]+1];
    sublo_lamda[2] = zsplit[myloc[2]];
    subhi_lamda[2] = zsplit[myloc[2]+1];

  } else {
    double (*mysplit)[2] = comm->mysplit;

    sublo_lamda[0] = mysplit[0][0];
    subhi_lamda[0] = mysplit[0][1];
    sublo_lamda[1] = mysplit[1][0];
    subhi_lamda[1] = mysplit[1][1];
    sublo_lamda[2] = mysplit[2][0];
    subhi_lamda[2] = mysplit[2][1];
  }
}

/* ----------------------------------------------------------------------
   set local subbox params for orthogonal boxes
   assumes global box is defined and proc assignment has been made
   uses comm->xyz_split or comm->mysplit
     to define subbox boundaries in consistent manner
   ensure subhi[max] = boxhi
------------------------------------------------------------------------- */

void Domain::set_local_box()
{
  if (triclinic) return;

  if (comm->layout != Comm::LAYOUT_TILED) {
    int *myloc = comm->myloc;
    int *procgrid = comm->procgrid;
    double *xsplit = comm->xsplit;
    double *ysplit = comm->ysplit;
    double *zsplit = comm->zsplit;

    sublo[0] = boxlo[0] + xprd*xsplit[myloc[0]];
    if (myloc[0] < procgrid[0]-1) subhi[0] = boxlo[0] + xprd*xsplit[myloc[0]+1];
    else subhi[0] = boxhi[0];

    sublo[1] = boxlo[1] + yprd*ysplit[myloc[1]];
    if (myloc[1] < procgrid[1]-1) subhi[1] = boxlo[1] + yprd*ysplit[myloc[1]+1];
    else subhi[1] = boxhi[1];

    sublo[2] = boxlo[2] + zprd*zsplit[myloc[2]];
    if (myloc[2] < procgrid[2]-1) subhi[2] = boxlo[2] + zprd*zsplit[myloc[2]+1];
    else subhi[2] = boxhi[2];

  } else {
    double (*mysplit)[2] = comm->mysplit;

    sublo[0] = boxlo[0] + xprd*mysplit[0][0];
    if (mysplit[0][1] < 1.0) subhi[0] = boxlo[0] + xprd*mysplit[0][1];
    else subhi[0] = boxhi[0];

    sublo[1] = boxlo[1] + yprd*mysplit[1][0];
    if (mysplit[1][1] < 1.0) subhi[1] = boxlo[1] + yprd*mysplit[1][1];
    else subhi[1] = boxhi[1];

    sublo[2] = boxlo[2] + zprd*mysplit[2][0];
    if (mysplit[2][1] < 1.0) subhi[2] = boxlo[2] + zprd*mysplit[2][1];
    else subhi[2] = boxhi[2];
  }
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxlo/hi
   for triclinic, atoms must be in lamda coords (0-1) before reset_box is called
------------------------------------------------------------------------- */

void Domain::reset_box()
{
  // perform shrink-wrapping

  // nothing to do for empty systems

  if (atom->natoms == 0) return;

  // compute extent of atoms on this proc
  // for triclinic, this is done in lamda space

  if (nonperiodic == 2) {
    double extent[3][2],all[3][2];

    extent[2][0] = extent[1][0] = extent[0][0] = BIG;
    extent[2][1] = extent[1][1] = extent[0][1] = -BIG;

    double **x = atom->x;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      extent[0][0] = MIN(extent[0][0],x[i][0]);
      extent[0][1] = MAX(extent[0][1],x[i][0]);
      extent[1][0] = MIN(extent[1][0],x[i][1]);
      extent[1][1] = MAX(extent[1][1],x[i][1]);
      extent[2][0] = MIN(extent[2][0],x[i][2]);
      extent[2][1] = MAX(extent[2][1],x[i][2]);
    }

    // compute extent across all procs
    // flip sign of MIN to do it in one Allreduce MAX

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);

    // for triclinic, convert back to box coords before changing box

    if (triclinic) lamda2x(atom->nlocal);

    // in shrink-wrapped dims, set box by atom extent
    // if minimum set, enforce min box size settings
    // for triclinic, convert lamda extent to box coords, then set box lo/hi
    // decided NOT to do the next comment - don't want to sneakily change tilt
    // for triclinic, adjust tilt factors if 2nd dim is shrink-wrapped,
    //   so that displacement in 1st dim stays the same

    if (triclinic == 0) {
      if (xperiodic == 0) {
        if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - small[0];
        else if (boundary[0][0] == 3)
          boxlo[0] = MIN(-all[0][0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = all[0][1] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(all[0][1]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - small[1];
        else if (boundary[1][0] == 3)
          boxlo[1] = MIN(-all[1][0]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = all[1][1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(all[1][1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
      }
      if (zperiodic == 0) {
        if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - small[2];
        else if (boundary[2][0] == 3)
          boxlo[2] = MIN(-all[2][0]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = all[2][1] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(all[2][1]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
      }

    } else {
      double lo[3],hi[3];
      if (xperiodic == 0) {
        lo[0] = -all[0][0]; lo[1] = 0.0; lo[2] = 0.0;
        lamda2x(lo,lo);
        hi[0] = all[0][1]; hi[1] = 0.0; hi[2] = 0.0;
        lamda2x(hi,hi);
        if (boundary[0][0] == 2) boxlo[0] = lo[0] - small[0];
        else if (boundary[0][0] == 3) boxlo[0] = MIN(lo[0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = hi[0] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(hi[0]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        lo[0] = 0.0; lo[1] = -all[1][0]; lo[2] = 0.0;
        lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = all[1][1]; hi[2] = 0.0;
        lamda2x(hi,hi);
        if (boundary[1][0] == 2) boxlo[1] = lo[1] - small[1];
        else if (boundary[1][0] == 3) boxlo[1] = MIN(lo[1]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = hi[1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(hi[1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
        //xy *= (boxhi[1]-boxlo[1]) / yprd;
      }
      if (zperiodic == 0) {
        lo[0] = 0.0; lo[1] = 0.0; lo[2] = -all[2][0];
        lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = 0.0; hi[2] = all[2][1];
        lamda2x(hi,hi);
        if (boundary[2][0] == 2) boxlo[2] = lo[2] - small[2];
        else if (boundary[2][0] == 3) boxlo[2] = MIN(lo[2]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = hi[2] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(hi[2]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
        //xz *= (boxhi[2]-boxlo[2]) / xprd;
        //yz *= (boxhi[2]-boxlo[2]) / yprd;
      }
    }
  }

  // reset box whether shrink-wrapping or not

  set_global_box();
  set_local_box();

  // if shrink-wrapped & kspace is defined (i.e. using MSM), call setup()
  // also call init() (to test for compatibility) ?

  if (nonperiodic == 2 && force->kspace) {
    //force->kspace->init();
    force->kspace->setup();
  }

  // if shrink-wrapped & triclinic, re-convert to lamda coords for new box
  // re-invoke pbc() b/c x2lamda result can be outside [0,1] due to roundoff

  if (nonperiodic == 2 && triclinic) {
    x2lamda(atom->nlocal);
    pbc();
  }
}

/* ----------------------------------------------------------------------
   define and store a general triclinic simulation box
   3 edge vectors of box = avec/bvec/cvec caller
   origin of edge vectors = origin_caller = lower left corner of box
   create mapping to restricted triclinic box
   set boxlo[3], boxhi[3] and 3 tilt factors
   create rotation matrices for general <--> restricted transformations
------------------------------------------------------------------------- */

void Domain::define_general_triclinic(double *avec_caller, double *bvec_caller,
                                      double *cvec_caller, double *origin_caller)
{
  if (triclinic || triclinic_general)
    error->all(FLERR,"General triclinic box edge vectors are already set");

  triclinic = triclinic_general = 1;

  avec[0] = avec_caller[0];
  avec[1] = avec_caller[1];
  avec[2] = avec_caller[2];

  bvec[0] = bvec_caller[0];
  bvec[1] = bvec_caller[1];
  bvec[2] = bvec_caller[2];

  cvec[0] = cvec_caller[0];
  cvec[1] = cvec_caller[1];
  cvec[2] = cvec_caller[2];

  // error check on cvec for 2d systems

  if (dimension == 2 && (cvec[0] != 0.0 || cvec[1] != 0.0))
    error->all(FLERR,"General triclinic box edge vector C invalid for 2d system");

  // rotate_g2r = rotation matrix from general to restricted triclnic
  // rotate_r2g = rotation matrix from restricted to general triclnic

  double aprime[3],bprime[3],cprime[3];
  general_to_restricted_rotation(avec,bvec,cvec,rotate_g2r,aprime,bprime,cprime);
  MathExtra::transpose3(rotate_g2r,rotate_r2g);

  // set restricted triclinic boxlo, boxhi, and tilt factors

  boxlo[0] = origin_caller[0];
  boxlo[1] = origin_caller[1];
  boxlo[2] = origin_caller[2];

  boxhi[0] = boxlo[0] + aprime[0];
  boxhi[1] = boxlo[1] + bprime[1];
  boxhi[2] = boxlo[2] + cprime[2];

  xy = bprime[0];
  xz = cprime[0];
  yz = cprime[1];
}

/* ----------------------------------------------------------------------
   compute rotation matrix to transform from general to restricted triclinic
   ABC = 3 general triclinic edge vectors
   rotmat = rotation matrix
   A`B`C` = 3 restricited triclinic edge vectors
------------------------------------------------------------------------- */

void Domain::general_to_restricted_rotation(double *a, double *b, double *c,
                                            double rotmat[3][3],
                                            double *aprime, double *bprime, double *cprime)
{
  // error checks
  // A,B,C cannot be co-planar
  // A x B must point in C direction (right-handed)

  double abcross[3];
  MathExtra::cross3(a,b,abcross);
  double dot = MathExtra::dot3(abcross,c);
  if (dot == 0.0)
    error->all(FLERR,"General triclinic edge vectors are co-planar");
  if (dot < 0.0)
    error->all(FLERR,"General triclinic edge vectors must be right-handed");

  // quat1 = convert A into A' along +x-axis
  // rot1 = unit vector to rotate A around
  // theta1 = angle of rotation calculated from
  //   A dot xunit = Ax = |A| cos(theta1)

  double rot1[3],quat1[4];
  double xaxis[3] = {1.0, 0.0, 0.0};

  double alen = MathExtra::len3(a);
  MathExtra::cross3(a,xaxis,rot1);
  MathExtra::norm3(rot1);
  double theta1 = acos(a[0]/alen);
  MathExtra::axisangle_to_quat(rot1,theta1,quat1);

  // rotmat1 = rotation matrix associated with quat1

  double rotmat1[3][3];
  MathExtra::quat_to_mat(quat1,rotmat1);

  // B1 = rotation of B by quat1 rotation matrix

  double b1[3];
  MathExtra::matvec(rotmat1,b,b1);

  // quat2 = rotation to convert B1 into B' in xy plane
  // Byz1 = projection of B1 into yz plane
  // +xaxis = unit vector to rotate B1 around
  // theta2 = angle of rotation calculated from
  //   Byz1 dot yunit = B1y = |Byz1| cos(theta2)
  // theta2 via acos() is positive (0 to PI)
  //   positive is valid if B1z < 0.0 else flip sign of theta2

  double byzvec1[3],quat2[4];
  MathExtra::copy3(b1,byzvec1);
  byzvec1[0] = 0.0;
  double byzvec1_len = MathExtra::len3(byzvec1);
  double theta2 = acos(b1[1]/byzvec1_len);
  if (b1[2] > 0.0) theta2 = -theta2;
  MathExtra::axisangle_to_quat(xaxis,theta2,quat2);

  // quat_single = rotation via single quat = quat2 * quat1
  // quat_r2g = rotation from restricted to general
  // rotmat = general to restricted rotation matrix

  double quat_single[4];
  MathExtra::quatquat(quat2,quat1,quat_single);
  MathExtra::quat_to_mat(quat_single,rotmat);

  // rotate general ABC to restricted triclinic A'B'C'

  MathExtra::matvec(rotmat,a,aprime);
  MathExtra::matvec(rotmat,b,bprime);
  MathExtra::matvec(rotmat,c,cprime);
}

/* ----------------------------------------------------------------------
   transform atom coords from general triclinic to restricted triclinic
------------------------------------------------------------------------- */

void Domain::general_to_restricted_coords(double *x)
{
  double xshift[3],xnew[3];

  xshift[0] = x[0] - boxlo[0];
  xshift[1] = x[1] - boxlo[1];
  xshift[2] = x[2] - boxlo[2];
  MathExtra::matvec(rotate_g2r,xshift,xnew);
  x[0] = xnew[0] + boxlo[0];
  x[1] = xnew[1] + boxlo[1];
  x[2] = xnew[2] + boxlo[2];
}

/* ----------------------------------------------------------------------
   transform atom coords from restricted triclinic to general triclinic
------------------------------------------------------------------------- */

void Domain::restricted_to_general_coords(double *x)
{
  double xshift[3],xnew[3];

  xshift[0] = x[0] - boxlo[0];
  xshift[1] = x[1] - boxlo[1];
  xshift[2] = x[2] - boxlo[2];
  MathExtra::matvec(rotate_r2g,xshift,xnew);
  x[0] = xnew[0] + boxlo[0];
  x[1] = xnew[1] + boxlo[1];
  x[2] = xnew[2] + boxlo[2];
}

void Domain::restricted_to_general_coords(double *x, double *xnew)
{
  double xshift[3];

  xshift[0] = x[0] - boxlo[0];
  xshift[1] = x[1] - boxlo[1];
  xshift[2] = x[2] - boxlo[2];
  MathExtra::matvec(rotate_r2g,xshift,xnew);
  xnew[0] += boxlo[0];
  xnew[1] += boxlo[1];
  xnew[2] += boxlo[2];
}

/* ----------------------------------------------------------------------
   transform atom vector from general triclinic to restricted triclinic
------------------------------------------------------------------------- */

void Domain::general_to_restricted_vector(double *v)
{
  double vnew[3];

  MathExtra::matvec(rotate_g2r,v,vnew);
  v[0] = vnew[0];
  v[1] = vnew[1];
  v[2] = vnew[2];
}

/* ----------------------------------------------------------------------
   transform atom vector from restricted triclinic to general triclinic
------------------------------------------------------------------------- */

void Domain::restricted_to_general_vector(double *v)
{
  double vnew[3];

  MathExtra::matvec(rotate_r2g,v,vnew);
  v[0] = vnew[0];
  v[1] = vnew[1];
  v[2] = vnew[2];
}

void Domain::restricted_to_general_vector(double *v, double *vnew)
{
  MathExtra::matvec(rotate_r2g,v,vnew);
}

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   if fix deform, remap velocity of fix group atoms by box edge velocities
   for triclinic, atoms must be in lamda coords (0-1) before pbc is called
   image = 10 or 20 bits for each dimension depending on sizeof(imageint)
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::pbc()
{
  int nlocal = atom->nlocal;
  if (!nlocal) return;
  int i;
  imageint idim,otherdims;
  double *lo,*hi,*period;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;

  // verify owned atoms have valid numerical coords
  // may not if computed pairwise force between 2 atoms at same location

  double *coord;
  int n3 = 3*nlocal;
  coord = &x[0][0];
  int flag = 0;
  for (i = 0; i < n3; i++)
    if (!std::isfinite(*coord++)) flag = 1;
  if (flag) error->one(FLERR,"Non-numeric atom coords - simulation unstable");

  // setup for PBC checks

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }

  // apply PBC to each owned atom

  for (i = 0; i < nlocal; i++) {
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
        x[i][0] += period[0];
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] += h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
        x[i][0] -= period[0];
        x[i][0] = MAX(x[i][0],lo[0]);
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] -= h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
    }

    if (yperiodic) {
      if (x[i][1] < lo[1]) {
        x[i][1] += period[1];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[5];
          v[i][1] += h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      if (x[i][1] >= hi[1]) {
        x[i][1] -= period[1];
        x[i][1] = MAX(x[i][1],lo[1]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[5];
          v[i][1] -= h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
    }

    if (zperiodic) {
      if (x[i][2] < lo[2]) {
        x[i][2] += period[2];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[4];
          v[i][1] += h_rate[3];
          v[i][2] += h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      if (x[i][2] >= hi[2]) {
        x[i][2] -= period[2];
        x[i][2] = MAX(x[i][2],lo[2]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[4];
          v[i][1] -= h_rate[3];
          v[i][2] -= h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   check that point is inside box boundaries, in [lo,hi) sense
   return 1 if true, 0 if false
------------------------------------------------------------------------- */

int Domain::inside(double* x)
{
  double *lo,*hi;
  double lamda[3];

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;

    if (x[0] < lo[0] || x[0] >= hi[0] ||
        x[1] < lo[1] || x[1] >= hi[1] ||
        x[2] < lo[2] || x[2] >= hi[2]) return 0;
    else return 1;

  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;

    x2lamda(x,lamda);

    if (lamda[0] < lo[0] || lamda[0] >= hi[0] ||
        lamda[1] < lo[1] || lamda[1] >= hi[1] ||
        lamda[2] < lo[2] || lamda[2] >= hi[2]) return 0;
    else return 1;
  }
}

/* ----------------------------------------------------------------------
   check that point is inside nonperiodic boundaries, in [lo,hi) sense
   return 1 if true, 0 if false
------------------------------------------------------------------------- */

int Domain::inside_nonperiodic(double* x)
{
  double *lo,*hi;
  double lamda[3];

  if (xperiodic && yperiodic && zperiodic) return 1;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;

    if (!xperiodic && (x[0] < lo[0] || x[0] >= hi[0])) return 0;
    if (!yperiodic && (x[1] < lo[1] || x[1] >= hi[1])) return 0;
    if (!zperiodic && (x[2] < lo[2] || x[2] >= hi[2])) return 0;
    return 1;

  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;

    x2lamda(x,lamda);

    if (!xperiodic && (lamda[0] < lo[0] || lamda[0] >= hi[0])) return 0;
    if (!yperiodic && (lamda[1] < lo[1] || lamda[1] >= hi[1])) return 0;
    if (!zperiodic && (lamda[2] < lo[2] || lamda[2] >= hi[2])) return 0;
    return 1;
  }
}

/* ----------------------------------------------------------------------
   warn if image flags of any bonded atoms are inconsistent
   could be a problem when using replicate or fix rigid
------------------------------------------------------------------------- */

void Domain::image_check()
{
  int i,j,k,n,imol,iatom;
  tagint tagprev;

  // only need to check if system is molecular and some dimension is periodic
  // if running verlet/split, don't check on KSpace partition since
  //    it has no ghost atoms and thus bond partners won't exist

  if (atom->molecular == Atom::ATOMIC) return;
  if (!xperiodic && !yperiodic && (dimension == 2 || !zperiodic)) return;
  if (strncmp(update->integrate_style,"verlet/split",12) == 0 &&
      universe->iworld != 0) return;

  // communicate unwrapped position of owned atoms to ghost atoms

  double **unwrap;
  memory->create(unwrap,atom->nmax,3,"domain:unwrap");

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    unmap(x[i],image[i],unwrap[i]);

  comm->forward_comm_array(3,unwrap);

  // compute unwrapped extent of each bond
  // flag if any bond component is longer than 1/2 of periodic box length
  // flag if any bond component is longer than non-periodic box length
  //   which means image flags in that dimension were different

  int molecular = atom->molecular;

  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  double delx,dely,delz;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;

  int flag = 0;
  for (i = 0; i < nlocal; i++) {
    if (molecular == Atom::MOLECULAR) n = num_bond[i];
    else {
      if (molindex[i] < 0) continue;
      imol = molindex[i];
      iatom = molatom[i];
      n = onemols[imol]->num_bond[iatom];
    }

    for (j = 0; j < n; j++) {
      if (molecular == Atom::MOLECULAR) {
        if (bond_type[i][j] <= 0) continue;
        k = atom->map(bond_atom[i][j]);
      } else {
        if (onemols[imol]->bond_type[iatom][j] < 0) continue;
        tagprev = tag[i] - iatom - 1;
        k = atom->map(onemols[imol]->bond_atom[iatom][j]+tagprev);
      }

      if (k == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR)
          error->one(FLERR,"Bond atom missing in image check");
        continue;
      }

      delx = fabs(unwrap[i][0] - unwrap[k][0]);
      dely = fabs(unwrap[i][1] - unwrap[k][1]);
      delz = fabs(unwrap[i][2] - unwrap[k][2]);

      if (xperiodic && delx > xprd_half) flag = 1;
      if (yperiodic && dely > yprd_half) flag = 1;
      if (dimension == 3 && zperiodic && delz > zprd_half) flag = 1;
      if (!xperiodic && delx > xprd) flag = 1;
      if (!yperiodic && dely > yprd) flag = 1;
      if (dimension == 3 && !zperiodic && delz > zprd) flag = 1;
    }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"Inconsistent image flags");

  if (lostbond == Thermo::WARN) {
    int all;
    MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
    if (all && comm->me == 0)
      error->warning(FLERR,"Bond atom missing in image check");
  }

  memory->destroy(unwrap);
}

/* ----------------------------------------------------------------------
   warn if end atoms in any bonded interaction
     are further apart than half a periodic box length
   could cause problems when bonded neighbor list is built since
     closest_image() could return wrong image
------------------------------------------------------------------------- */

void Domain::box_too_small_check()
{
  int i,j,k,n,imol,iatom;
  tagint tagprev;

  // only need to check if system is molecular and some dimension is periodic
  // if running verlet/split, don't check on KSpace partition since
  //    it has no ghost atoms and thus bond partners won't exist

  if (atom->molecular == Atom::ATOMIC) return;
  if (!xperiodic && !yperiodic && (dimension == 2 || !zperiodic)) return;
  if (strncmp(update->integrate_style,"verlet/split",12) == 0 &&
      universe->iworld != 0) return;

  // maxbondall = longest current bond length
  // if periodic box dim is tiny (less than 2 * bond-length),
  //   minimum_image() itself may compute bad bond lengths
  // in this case, image_check() should warn,
  //   assuming 2 atoms have consistent image flags

  int molecular = atom->molecular;

  double **x = atom->x;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  int nlocal = atom->nlocal;

  double delx,dely,delz,rsq;
  double maxbondme = 0.0;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;

  for (i = 0; i < nlocal; i++) {
    if (molecular == Atom::MOLECULAR) n = num_bond[i];
    else {
      if (molindex[i] < 0) continue;
      imol = molindex[i];
      iatom = molatom[i];
      n = onemols[imol]->num_bond[iatom];
    }

    for (j = 0; j < n; j++) {
      if (molecular == Atom::MOLECULAR) {
        if (bond_type[i][j] <= 0) continue;
        k = atom->map(bond_atom[i][j]);
      } else {
        if (onemols[imol]->bond_type[iatom][j] < 0) continue;
        tagprev = tag[i] - iatom - 1;
        k = atom->map(onemols[imol]->bond_atom[iatom][j]+tagprev);
      }

      if (k == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR)
          error->one(FLERR,"Bond atom missing in box size check");
        continue;
      }

      delx = x[i][0] - x[k][0];
      dely = x[i][1] - x[k][1];
      delz = x[i][2] - x[k][2];
      minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      maxbondme = MAX(maxbondme,rsq);
    }
  }

  if (lostbond == Thermo::WARN) {
    int all;
    MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
    if (all && comm->me == 0)
      error->warning(FLERR,"Bond atom missing in box size check");
  }

  double maxbondall;
  MPI_Allreduce(&maxbondme,&maxbondall,1,MPI_DOUBLE,MPI_MAX,world);
  maxbondall = sqrt(maxbondall);

  // maxdelta = furthest apart 2 atoms in a bonded interaction can be
  // include BONDSTRETCH factor to account for dynamics

  double maxdelta = maxbondall * BONDSTRETCH;
  if (atom->nangles) maxdelta = 2.0 * maxbondall * BONDSTRETCH;
  if (atom->ndihedrals) maxdelta = 3.0 * maxbondall * BONDSTRETCH;

  // warn if maxdelta > than half any periodic box length
  // since atoms in the interaction could rotate into that dimension

  int flag = 0;
  if (xperiodic && maxdelta > xprd_half) flag = 1;
  if (yperiodic && maxdelta > yprd_half) flag = 1;
  if (dimension == 3 && zperiodic && maxdelta > zprd_half) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,
                   "Bond/angle/dihedral extent > half of periodic box length");
}

/* ----------------------------------------------------------------------
  check warn if any proc's subbox is smaller than thresh
    since may lead to lost atoms in comm->exchange()
  current callers set thresh = neighbor skin
------------------------------------------------------------------------- */

void Domain::subbox_too_small_check(double thresh)
{
  int flag = 0;
  if (!triclinic) {
    if (subhi[0]-sublo[0] < thresh || subhi[1]-sublo[1] < thresh) flag = 1;
    if (dimension == 3 && subhi[2]-sublo[2] < thresh) flag = 1;
  } else {
    double delta = subhi_lamda[0] - sublo_lamda[0];
    if (delta*prd[0] < thresh) flag = 1;
    delta = subhi_lamda[1] - sublo_lamda[1];
    if (delta*prd[1] < thresh) flag = 1;
    if (dimension == 3) {
      delta = subhi_lamda[2] - sublo_lamda[2];
      if (delta*prd[2] < thresh) flag = 1;
    }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"Proc sub-domain size < neighbor skin, "
                   "could lead to lost atoms");
}

/* ----------------------------------------------------------------------
   minimum image convention in periodic dimensions
   use 1/2 of box size as test
   for triclinic, also add/subtract tilt factors in other dims as needed
   changed "if" to "while" to enable distance to
     far-away ghost atom returned by atom->map() to be wrapped back into box
     could be problem for looking up atom IDs when cutoff > boxsize
   should be used for most cases where the difference in the image count
     is small (usually 0 or 1)
   use minimum_image_big() when a large difference between image counts is expected
------------------------------------------------------------------------- */

static constexpr double MAXIMGCOUNT = 16;

void Domain::minimum_image(double &dx, double &dy, double &dz) const
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(dx) > (MAXIMGCOUNT * xprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dx);
      while (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(dy) > (MAXIMGCOUNT * yprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dy);
      while (fabs(dy) > yprd_half) {
        if (dy < 0.0) dy += yprd;
        else dy -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(dz) > (MAXIMGCOUNT * zprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dz);
      while (fabs(dz) > zprd_half) {
        if (dz < 0.0) dz += zprd;
        else dz -= zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (fabs(dz) > (MAXIMGCOUNT * zprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dz);
      while (fabs(dz) > zprd_half) {
        if (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        } else {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      }
    }
    if (yperiodic) {
      if (fabs(dy) > (MAXIMGCOUNT * yprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dy);
      while (fabs(dy) > yprd_half) {
        if (dy < 0.0) {
          dy += yprd;
          dx += xy;
        } else {
          dy -= yprd;
          dx -= xy;
        }
      }
    }
    if (xperiodic) {
      if (fabs(dx) > (MAXIMGCOUNT * xprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dx);
      while (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   minimum image convention in periodic dimensions
   use 1/2 of box size as test
   for triclinic, also add/subtract tilt factors in other dims as needed
   allow multiple box lengths to enable distance to
     far-away ghost atom returned by atom->map() to be wrapped back into box
     could be problem for looking up atom IDs when cutoff > boxsize
   this should be used when there is a large image count difference possible
     this applies for example to fix rigid/small
------------------------------------------------------------------------- */

void Domain::minimum_image_big(double &dx, double &dy, double &dz) const
{
  if (triclinic == 0) {
    if (xperiodic) {
      double dfactor = dx/xprd + 0.5;
      if (dx < 0) dfactor -= 1.0;
      if (dfactor > MAXSMALLINT)
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dx);
      dx -= xprd * static_cast<int>(dfactor);
    }
    if (yperiodic) {
      double dfactor = dy/yprd + 0.5;
      if (dy < 0) dfactor -= 1.0;
      if (dfactor > MAXSMALLINT)
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dy);
      dy -= yprd * static_cast<int>(dfactor);
    }
    if (zperiodic) {
      double dfactor = dz/zprd + 0.5;
      if (dz < 0) dfactor -= 1.0;
      if (dfactor > MAXSMALLINT)
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dz);
      dz -= zprd * static_cast<int>(dfactor);
    }

  } else {
    if (zperiodic) {
      double dfactor = dz/zprd + 0.5;
      if (dz < 0) dfactor -= 1.0;
      if (dfactor > MAXSMALLINT)
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dz);
      int factor = static_cast<int>(dfactor);
      dz -= zprd * factor;
      dy -= yz * factor;
      dx -= xz * factor;
    }
    if (yperiodic) {
      double dfactor = dy/yprd + 0.5;
      if (dy < 0) dfactor -= 1.0;
      if (dfactor > MAXSMALLINT)
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dy);
      int factor = static_cast<int>(dfactor);
      dy -= yprd * factor;
      dx -= xy * factor;
    }
    if (xperiodic) {
      double dfactor = dx/xprd + 0.5;
      if (dx < 0) dfactor -= 1.0;
      if (dfactor > MAXSMALLINT)
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dx);
      dx -= xprd * static_cast<int>(dfactor);
    }
  }
}

/* ----------------------------------------------------------------------
   return local index of atom J or any of its images that is closest to atom I
   if J is not a valid index like -1, just return it
------------------------------------------------------------------------- */

int Domain::closest_image(int i, int j)
{
  if (j < 0) return j;

  int *sametag = atom->sametag;
  double **x = atom->x;
  double *xi = x[i];

  int closest = j;
  double delx = xi[0] - x[j][0];
  double dely = xi[1] - x[j][1];
  double delz = xi[2] - x[j][2];
  double rsqmin = delx*delx + dely*dely + delz*delz;
  double rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = xi[0] - x[j][0];
    dely = xi[1] - x[j][1];
    delz = xi[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }

  return closest;
}

/* ----------------------------------------------------------------------
   return local index of atom J or any of its images that is closest to pos
   if J is not a valid index like -1, just return it
------------------------------------------------------------------------- */

int Domain::closest_image(const double * const pos, int j)
{
  if (j < 0) return j;

  const int * const sametag = atom->sametag;
  const double * const * const x = atom->x;

  int closest = j;
  double delx = pos[0] - x[j][0];
  double dely = pos[1] - x[j][1];
  double delz = pos[2] - x[j][2];
  double rsqmin = delx*delx + dely*dely + delz*delz;
  double rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = pos[0] - x[j][0];
    dely = pos[1] - x[j][1];
    delz = pos[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }

  return closest;
}

/* ----------------------------------------------------------------------
   find and return Xj image = periodic image of Xj that is closest to Xi
   for triclinic, add/subtract tilt factors in other dims as needed
   called by ServerMD class and LammpsInterface in lib/atc
------------------------------------------------------------------------- */

void Domain::closest_image(const double * const xi, const double * const xj, double * const xjimage)
{
  double dx = xj[0] - xi[0];
  double dy = xj[1] - xi[1];
  double dz = xj[2] - xi[2];

  if (triclinic == 0) {
    if (xperiodic) {
      if (dx < 0.0) {
        while (dx < 0.0) dx += xprd;
        if (dx > xprd_half) dx -= xprd;
      } else {
        while (dx > 0.0) dx -= xprd;
        if (dx < -xprd_half) dx += xprd;
      }
    }
    if (yperiodic) {
      if (dy < 0.0) {
        while (dy < 0.0) dy += yprd;
        if (dy > yprd_half) dy -= yprd;
      } else {
        while (dy > 0.0) dy -= yprd;
        if (dy < -yprd_half) dy += yprd;
      }
    }
    if (zperiodic) {
      if (dz < 0.0) {
        while (dz < 0.0) dz += zprd;
        if (dz > zprd_half) dz -= zprd;
      } else {
        while (dz > 0.0) dz -= zprd;
        if (dz < -zprd_half) dz += zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (dz < 0.0) {
        while (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        }
        if (dz > zprd_half) {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      } else {
        while (dz > 0.0) {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
        if (dz < -zprd_half) {
          dz += zprd;
          dy += yz;
          dx += xz;
        }
      }
    }
    if (yperiodic) {
      if (dy < 0.0) {
        while (dy < 0.0) {
          dy += yprd;
          dx += xy;
        }
        if (dy > yprd_half) {
          dy -= yprd;
          dx -= xy;
        }
      } else {
        while (dy > 0.0) {
          dy -= yprd;
          dx -= xy;
        }
        if (dy < -yprd_half) {
          dy += yprd;
          dx += xy;
        }
      }
    }
    if (xperiodic) {
      if (dx < 0.0) {
        while (dx < 0.0) dx += xprd;
        if (dx > xprd_half) dx -= xprd;
      } else {
        while (dx > 0.0) dx -= xprd;
        if (dx < -xprd_half) dx += xprd;
      }
    }
  }

  xjimage[0] = xi[0] + dx;
  xjimage[1] = xi[1] + dy;
  xjimage[2] = xi[2] + dz;
}

/* ----------------------------------------------------------------------
   remap the point into the periodic box no matter how far away
   adjust 3 image flags encoded in image accordingly
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before doing remap
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::remap(double *x, imageint &image)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];
  imageint idim,otherdims;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }

  if (xperiodic) {
    while (coord[0] < lo[0]) {
      coord[0] += period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim--;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    while (coord[0] >= hi[0]) {
      coord[0] -= period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim++;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) {
      coord[1] += period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    while (coord[1] >= hi[1]) {
      coord[1] -= period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) {
      coord[2] += period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    while (coord[2] >= hi[2]) {
      coord[2] -= period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    coord[2] = MAX(coord[2],lo[2]);
  }

  if (triclinic) lamda2x(coord,x);
}

/* ----------------------------------------------------------------------
   remap the point into the periodic box no matter how far away
   no image flag calculation
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before remap
------------------------------------------------------------------------- */

void Domain::remap(double *x)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }

  if (xperiodic) {
    while (coord[0] < lo[0]) coord[0] += period[0];
    while (coord[0] >= hi[0]) coord[0] -= period[0];
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) coord[1] += period[1];
    while (coord[1] >= hi[1]) coord[1] -= period[1];
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) coord[2] += period[2];
    while (coord[2] >= hi[2]) coord[2] -= period[2];
    coord[2] = MAX(coord[2],lo[2]);
  }

  if (triclinic) lamda2x(coord,x);
}

/* ----------------------------------------------------------------------
   remap all points into the periodic box no matter how far away
   adjust 3 image flags encoded in image accordingly
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before doing remap
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::remap_all()
{
  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) remap(x[i],image[i]);
}

/* ----------------------------------------------------------------------
   remap xnew to be within half box length of xold
   do it directly, not iteratively, in case is far away
   for triclinic, both points are converted to lamda coords (0-1) before remap
------------------------------------------------------------------------- */

void Domain::remap_near(double *xnew, double *xold)
{
  int n;
  double *coordnew,*coordold,*period,*half;
  double lamdanew[3],lamdaold[3];

  if (triclinic == 0) {
    period = prd;
    half = prd_half;
    coordnew = xnew;
    coordold = xold;
  } else {
    period = prd_lamda;
    half = prd_half_lamda;
    x2lamda(xnew,lamdanew);
    coordnew = lamdanew;
    x2lamda(xold,lamdaold);
    coordold = lamdaold;
  }

  // iterative form
  // if (xperiodic) {
  //   while (coordnew[0]-coordold[0] > half[0]) coordnew[0] -= period[0];
  //   while (coordold[0]-coordnew[0] > half[0]) coordnew[0] += period[0];
  // }

  if (xperiodic) {
    if (coordnew[0]-coordold[0] > period[0]) {
      n = static_cast<int> ((coordnew[0]-coordold[0])/period[0]);
      coordnew[0] -= n*period[0];
    }
    while (coordnew[0]-coordold[0] > half[0]) coordnew[0] -= period[0];
    if (coordold[0]-coordnew[0] > period[0]) {
      n = static_cast<int> ((coordold[0]-coordnew[0])/period[0]);
      coordnew[0] += n*period[0];
    }
    while (coordold[0]-coordnew[0] > half[0]) coordnew[0] += period[0];
  }

  if (yperiodic) {
    if (coordnew[1]-coordold[1] > period[1]) {
      n = static_cast<int> ((coordnew[1]-coordold[1])/period[1]);
      coordnew[1] -= n*period[1];
    }
    while (coordnew[1]-coordold[1] > half[1]) coordnew[1] -= period[1];
    if (coordold[1]-coordnew[1] > period[1]) {
      n = static_cast<int> ((coordold[1]-coordnew[1])/period[1]);
      coordnew[1] += n*period[1];
    }
    while (coordold[1]-coordnew[1] > half[1]) coordnew[1] += period[1];
  }

  if (zperiodic) {
    if (coordnew[2]-coordold[2] > period[2]) {
      n = static_cast<int> ((coordnew[2]-coordold[2])/period[2]);
      coordnew[2] -= n*period[2];
    }
    while (coordnew[2]-coordold[2] > half[2]) coordnew[2] -= period[2];
    if (coordold[2]-coordnew[2] > period[2]) {
      n = static_cast<int> ((coordold[2]-coordnew[2])/period[2]);
      coordnew[2] += n*period[2];
    }
    while (coordold[2]-coordnew[2] > half[2]) coordnew[2] += period[2];
  }

  if (triclinic) lamda2x(coordnew,xnew);
}

/* ----------------------------------------------------------------------
   remap the point to specific image flags
   x overwritten with result, reset image flag
   for triclinic, use h[] to add in tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::unmap_inv(double *x, imageint image)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  if (triclinic == 0) {
    x[0] -= xbox*xprd;
    x[1] -= ybox*yprd;
    x[2] -= zbox*zprd;
  } else {
    x[0] -= h[0]*xbox + h[5]*ybox + h[4]*zbox;
    x[1] -= h[1]*ybox + h[3]*zbox;
    x[2] -= h[2]*zbox;
  }
}

/* ----------------------------------------------------------------------
   unmap the point via image flags
   x overwritten with result, don't reset image flag
   for triclinic, use h[] to add in tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::unmap(double *x, imageint image)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  if (triclinic == 0) {
    x[0] += xbox*xprd;
    x[1] += ybox*yprd;
    x[2] += zbox*zprd;
  } else {
    x[0] += h[0]*xbox + h[5]*ybox + h[4]*zbox;
    x[1] += h[1]*ybox + h[3]*zbox;
    x[2] += h[2]*zbox;
  }
}

/* ----------------------------------------------------------------------
   unmap the point via image flags
   result returned in y, don't reset image flag
   for triclinic, use h[] to add in tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::unmap(const double *x, imageint image, double *y)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  if (triclinic == 0) {
    y[0] = x[0] + xbox*xprd;
    y[1] = x[1] + ybox*yprd;
    y[2] = x[2] + zbox*zprd;
  } else {
    y[0] = x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    y[1] = x[1] + h[1]*ybox + h[3]*zbox;
    y[2] = x[2] + h[2]*zbox;
  }
}

/* ----------------------------------------------------------------------
   adjust image flags due to triclinic box flip
   flip operation is changing box vectors A,B,C to new A',B',C'
     A' = A              (A does not change)
     B' = B + mA         (B shifted by A)
     C' = C + pB + nA    (C shifted by B and/or A)
   this requires the image flags change from (a,b,c) to (a',b',c')
   so that x_unwrap for each atom is same before/after
     x_unwrap_before = xlocal + aA + bB + cC
     x_unwrap_after = xlocal + a'A' + b'B' + c'C'
   this requires:
     c' = c
     b' = b - cp
     a' = a - (b-cp)m - cn = a - b'm - cn
   in other words, for xy flip, change in x flag depends on current y flag
   this is b/c the xy flip dramatically changes which tiled image of
     simulation box an unwrapped point maps to
------------------------------------------------------------------------- */

void Domain::image_flip(int m, int n, int p)
{
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    int xbox = (image[i] & IMGMASK) - IMGMAX;
    int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    int zbox = (image[i] >> IMG2BITS) - IMGMAX;

    ybox -= p*zbox;
    xbox -= m*ybox + n*zbox;

    image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
      (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
      (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
  }
}

/* ----------------------------------------------------------------------
   return 1 if this proc owns atom with coords x, else return 0
   x is returned remapped into periodic box
   if image flag is passed, flag is updated via remap(x,image)
   if image = nullptr is passed, no update with remap(x)
   if shrinkexceed, atom can be outside shrinkwrap boundaries
   called from create_atoms() in library.cpp
------------------------------------------------------------------------- */

int Domain::ownatom(int /*id*/, double *x, imageint *image, int shrinkexceed)
{
  double lamda[3];
  double *coord,*blo,*bhi,*slo,*shi;

  if (image) remap(x,*image);
  else remap(x);

  // if triclinic, convert to lamda coords (0-1)
  // for periodic dims, resulting coord must satisfy 0.0 <= coord < 1.0

  if (triclinic) {
    x2lamda(x,lamda);
    if (xperiodic && (lamda[0] < 0.0 || lamda[0] >= 1.0)) lamda[0] = 0.0;
    if (yperiodic && (lamda[1] < 0.0 || lamda[1] >= 1.0)) lamda[1] = 0.0;
    if (zperiodic && (lamda[2] < 0.0 || lamda[2] >= 1.0)) lamda[2] = 0.0;
    coord = lamda;
  } else coord = x;

  // box and subbox bounds for orthogonal vs triclinic

  if (triclinic == 0) {
    blo = boxlo;
    bhi = boxhi;
    slo = sublo;
    shi = subhi;
  } else {
    blo = boxlo_lamda;
    bhi = boxhi_lamda;
    slo = sublo_lamda;
    shi = subhi_lamda;
  }

  if (coord[0] >= slo[0] && coord[0] < shi[0] &&
      coord[1] >= slo[1] && coord[1] < shi[1] &&
      coord[2] >= slo[2] && coord[2] < shi[2]) return 1;

  // check if atom did not return 1 only b/c it was
  //   outside a shrink-wrapped boundary

  if (shrinkexceed) {
    int outside = 0;
    if (coord[0] < blo[0] && boundary[0][0] > 1) outside = 1;
    if (coord[0] >= bhi[0] && boundary[0][1] > 1) outside = 1;
    if (coord[1] < blo[1] && boundary[1][0] > 1) outside = 1;
    if (coord[1] >= bhi[1] && boundary[1][1] > 1) outside = 1;
    if (coord[2] < blo[2] && boundary[2][0] > 1) outside = 1;
    if (coord[2] >= bhi[2] && boundary[2][1] > 1) outside = 1;
    if (!outside) return 0;

    // newcoord = coords pushed back to be on shrink-wrapped boundary
    // newcoord is a copy, so caller's x[] is not affected

    double newcoord[3];
    if (coord[0] < blo[0] && boundary[0][0] > 1) newcoord[0] = blo[0];
    else if (coord[0] >= bhi[0] && boundary[0][1] > 1) newcoord[0] = bhi[0];
    else newcoord[0] = coord[0];
    if (coord[1] < blo[1] && boundary[1][0] > 1) newcoord[1] = blo[1];
    else if (coord[1] >= bhi[1] && boundary[1][1] > 1) newcoord[1] = bhi[1];
    else newcoord[1] = coord[1];
    if (coord[2] < blo[2] && boundary[2][0] > 1) newcoord[2] = blo[2];
    else if (coord[2] >= bhi[2] && boundary[2][1] > 1) newcoord[2] = bhi[2];
    else newcoord[2] = coord[2];

    // re-test for newcoord inside my sub-domain
    // use <= test for upper-boundary since may have just put atom at boxhi

    if (newcoord[0] >= slo[0] && newcoord[0] <= shi[0] &&
        newcoord[1] >= slo[1] && newcoord[1] <= shi[1] &&
        newcoord[2] >= slo[2] && newcoord[2] <= shi[2]) return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   create a lattice
------------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char **arg)
{
  if (lattice) delete lattice;
  lattice = nullptr;
  lattice = new Lattice(lmp,narg,arg);
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "region", error);

  if (strcmp(arg[1],"delete") == 0) {
    delete_region(arg[0]);
    return;
  }

  if (strcmp(arg[1],"none") == 0)
    error->all(FLERR,"Unrecognized region style 'none'");

  if (get_region_by_id(arg[0])) error->all(FLERR,"Reuse of region ID {}", arg[0]);

  // create the Region
  Region *newregion = nullptr;

  if (lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      std::string estyle = std::string(arg[1]) + "/" + lmp->non_pair_suffix();
      if (region_map->find(estyle) != region_map->end()) {
        RegionCreator &region_creator = (*region_map)[estyle];
        newregion = region_creator(lmp, narg, arg);
      }
    }

    if (!newregion && lmp->suffix2) {
      std::string estyle = std::string(arg[1]) + "/" + lmp->suffix2;
      if (region_map->find(estyle) != region_map->end()) {
        RegionCreator &region_creator = (*region_map)[estyle];
        newregion = region_creator(lmp, narg, arg);
      }
    }
  }

  if (!newregion && (region_map->find(arg[1]) != region_map->end())) {
    RegionCreator &region_creator = (*region_map)[arg[1]];
    newregion = region_creator(lmp, narg, arg);
  }

  if (!newregion)
    error->all(FLERR,utils::check_packages_for_style("region",arg[1],lmp));

  // initialize any region variables via init()
  // in case region is used between runs, e.g. to print a variable

  newregion->init();
  regions.insert(newregion);
}

/* ----------------------------------------------------------------------
   delete a region
------------------------------------------------------------------------- */

void Domain::delete_region(Region *reg)
{
  if (!reg) return;

  regions.erase(reg);
  delete reg;
}

void Domain::delete_region(const std::string &id)
{
  auto reg = get_region_by_id(id);
  if (!reg) error->all(FLERR,"Delete region {} does not exist", id);
  delete_region(reg);
}

/* ----------------------------------------------------------------------
   return pointer to region name matches existing region ID
   return null if no match
------------------------------------------------------------------------- */

Region *Domain::get_region_by_id(const std::string &name) const
{
  for (auto &reg : regions)
    if (name == reg->id) return reg;
  return nullptr;
}

/* ----------------------------------------------------------------------
   look up pointers to regions by region style name
   return vector with matching pointers
------------------------------------------------------------------------- */

const std::vector<Region *> Domain::get_region_by_style(const std::string &name) const
{
  std::vector<Region *> matches;
  if (name.empty()) return matches;

  for (auto &reg : regions)
    if (name == reg->style)  matches.push_back(reg);

  return matches;
}

/* ----------------------------------------------------------------------
   return list of regions as vector
------------------------------------------------------------------------- */

const std::vector<Region *> Domain::get_region_list()
{
  return std::vector<Region *>(regions.begin(), regions.end());
}

/* ----------------------------------------------------------------------
   (re)set boundary settings
   flag = 0, called from the input script
   flag = 1, called from change box command
------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg, int flag)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command: expected 3 arguments but found {}", narg);

  char c;
  for (int idim = 0; idim < 3; idim++)
    for (int iside = 0; iside < 2; iside++) {
      if (iside == 0) c = arg[idim][0];
      else if (iside == 1 && strlen(arg[idim]) == 1) c = arg[idim][0];
      else c = arg[idim][1];

      if (c == 'p') boundary[idim][iside] = 0;
      else if (c == 'f') boundary[idim][iside] = 1;
      else if (c == 's') boundary[idim][iside] = 2;
      else if (c == 'm') boundary[idim][iside] = 3;
      else {
        if (flag == 0) error->all(FLERR,"Unknown boundary keyword: {}", c);
        if (flag == 1) error->all(FLERR,"Unknown change_box keyword: {}", c);
      }
    }

  for (int idim = 0; idim < 3; idim++)
    if ((boundary[idim][0] == 0 && boundary[idim][1]) ||
        (boundary[idim][0] && boundary[idim][1] == 0))
      error->all(FLERR,"Both sides of boundary must be periodic");

  if (boundary[0][0] == 0) xperiodic = 1;
  else xperiodic = 0;
  if (boundary[1][0] == 0) yperiodic = 1;
  else yperiodic = 0;
  if (boundary[2][0] == 0) zperiodic = 1;
  else zperiodic = 0;

  // record if we changed a periodic boundary to a non-periodic one

  int pflag=0;
  if ((periodicity[0] && !xperiodic)
      || (periodicity[1] && !yperiodic)
      || (periodicity[2] && !zperiodic)) pflag = 1;

  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
    nonperiodic = 1;
    if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
        boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
        boundary[2][0] >= 2 || boundary[2][1] >= 2) nonperiodic = 2;
  }

  // force non-zero image flags to zero for non-periodic dimensions
  // keep track if a change was made, so we can print a warning message

  if (pflag) {
    pflag = 0;
    for (int i=0; i < atom->nlocal; ++i) {
      int xbox = (atom->image[i] & IMGMASK) - IMGMAX;
      int ybox = (atom->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (atom->image[i] >> IMG2BITS) - IMGMAX;
      if ((!xperiodic) && (xbox != 0)) { xbox = 0; pflag = 1; }
      if ((!yperiodic) && (ybox != 0)) { ybox = 0; pflag = 1; }
      if ((!zperiodic) && (zbox != 0)) { zbox = 0; pflag = 1; }
      atom->image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }
    int flag_all;
    MPI_Allreduce(&pflag,&flag_all, 1, MPI_INT, MPI_SUM, world);
    if ((flag_all > 0) && (comm->me == 0))
      error->warning(FLERR,"Resetting image flags for non-periodic dimensions");
  }
}

/* ----------------------------------------------------------------------
   print box info, orthogonal or triclinic
------------------------------------------------------------------------- */

void Domain::print_box(const std::string &prefix)
{
  if (comm->me == 0) {
    std::string mesg = prefix;
    if (triclinic == 0) {
      mesg += fmt::format("orthogonal box = ({:.8} {:.8} {:.8}) to "
                          "({:.8} {:.8} {:.8})\n",boxlo[0],boxlo[1],
                          boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
    } else {
      mesg += fmt::format("triclinic box = ({:.8} {:.8} {:.8}) to "
                          "({:.8} {:.8} {:.8}) with tilt "
                          "({:.8} {:.8} {:.8})\n",boxlo[0],boxlo[1],
                          boxlo[2],boxhi[0],boxhi[1],boxhi[2],xy,xz,yz);
    }
    utils::logmesg(lmp,mesg);
  }
}

/* ----------------------------------------------------------------------
   format boundary string for output
   assume str is 9 chars or more in length
------------------------------------------------------------------------- */

void Domain::boundary_string(char *str)
{
  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
    for (int iside = 0; iside < 2; iside++) {
      if (boundary[idim][iside] == 0) str[m++] = 'p';
      else if (boundary[idim][iside] == 1) str[m++] = 'f';
      else if (boundary[idim][iside] == 2) str[m++] = 's';
      else if (boundary[idim][iside] == 3) str[m++] = 'm';
    }
    str[m++] = ' ';
  }
  str[8] = '\0';
}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for all N atoms
   x = H lamda + x0;
------------------------------------------------------------------------- */

void Domain::lamda2x(int n)
{
  double **x = atom->x;

  for (int i = 0; i < n; i++) {
    x[i][0] = h[0]*x[i][0] + h[5]*x[i][1] + h[4]*x[i][2] + boxlo[0];
    x[i][1] = h[1]*x[i][1] + h[3]*x[i][2] + boxlo[1];
    x[i][2] = h[2]*x[i][2] + boxlo[2];
  }
}

void Domain::lamda2x(int n, int groupbit)
{
  double **x = atom->x;
  int *mask = atom->mask;

  for (int i = 0; i < n; i++)
    if (mask[i] & groupbit) {
      x[i][0] = h[0]*x[i][0] + h[5]*x[i][1] + h[4]*x[i][2] + boxlo[0];
      x[i][1] = h[1]*x[i][1] + h[3]*x[i][2] + boxlo[1];
      x[i][2] = h[2]*x[i][2] + boxlo[2];
    }
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for all N atoms
   lamda = H^-1 (x - x0)
------------------------------------------------------------------------- */

void Domain::x2lamda(int n)
{
  double delta[3];
  double **x = atom->x;

  for (int i = 0; i < n; i++) {
    delta[0] = x[i][0] - boxlo[0];
    delta[1] = x[i][1] - boxlo[1];
    delta[2] = x[i][2] - boxlo[2];

    x[i][0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
    x[i][1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
    x[i][2] = h_inv[2]*delta[2];
  }
}

void Domain::x2lamda(int n, int groupbit)
{
  double delta[3];
  double **x = atom->x;
  int *mask = atom->mask;

  for (int i = 0; i < n; i++)
    if (mask[i] & groupbit) {
      delta[0] = x[i][0] - boxlo[0];
      delta[1] = x[i][1] - boxlo[1];
      delta[2] = x[i][2] - boxlo[2];

      x[i][0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
      x[i][1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
      x[i][2] = h_inv[2]*delta[2];
    }

}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for one atom
   x = H lamda + x0;
   lamda and x can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::lamda2x(double *lamda, double *x)
{
  x[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0];
  x[1] = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1];
  x[2] = h[2]*lamda[2] + boxlo[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for one atom
   lamda = H^-1 (x - x0)
   x and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::x2lamda(double *x, double *lamda)
{
  double delta[3];
  delta[0] = x[0] - boxlo[0];
  delta[1] = x[1] - boxlo[1];
  delta[2] = x[2] - boxlo[2];

  lamda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  lamda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  lamda[2] = h_inv[2]*delta[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for one atom
   use my_boxlo & my_h_inv stored by caller for previous state of box
   lamda = H^-1 (x - x0)
   x and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::x2lamda(double *x, double *lamda,
                     double *my_boxlo, double *my_h_inv)
{
  double delta[3];
  delta[0] = x[0] - my_boxlo[0];
  delta[1] = x[1] - my_boxlo[1];
  delta[2] = x[2] - my_boxlo[2];

  lamda[0] = my_h_inv[0]*delta[0] + my_h_inv[5]*delta[1] + my_h_inv[4]*delta[2];
  lamda[1] = my_h_inv[1]*delta[1] + my_h_inv[3]*delta[2];
  lamda[2] = my_h_inv[2]*delta[2];
}

/* ----------------------------------------------------------------------
   convert 8 lamda corner pts of lo/hi box to box coords
   return bboxlo/hi = bounding box around 8 corner pts in box coords
------------------------------------------------------------------------- */

void Domain::bbox(double *lo, double *hi, double *bboxlo, double *bboxhi)
{
  double x[3];

  bboxlo[0] = bboxlo[1] = bboxlo[2] = BIG;
  bboxhi[0] = bboxhi[1] = bboxhi[2] = -BIG;

  x[0] = lo[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
}

/* ----------------------------------------------------------------------
   compute 8 corner pts of my triclinic sub-box
   output is in corners, see ordering in lamda_box_corners
------------------------------------------------------------------------- */

void Domain::box_corners()
{
  lamda_box_corners(boxlo_lamda,boxhi_lamda);
}

/* ----------------------------------------------------------------------
   compute 8 corner pts of my triclinic sub-box
   output is in corners, see ordering in lamda_box_corners
------------------------------------------------------------------------- */

void Domain::subbox_corners()
{
  lamda_box_corners(sublo_lamda,subhi_lamda);
}

/* ----------------------------------------------------------------------
   compute 8 corner pts of any triclinic box with lo/hi in lamda coords
   8 output corners are ordered with x changing fastest, then y, finally z
   could be more efficient if just coded with xy,yz,xz explicitly
------------------------------------------------------------------------- */

void Domain::lamda_box_corners(double *lo, double *hi)
{
  corners[0][0] = lo[0]; corners[0][1] = lo[1]; corners[0][2] = lo[2];
  lamda2x(corners[0],corners[0]);
  corners[1][0] = hi[0]; corners[1][1] = lo[1]; corners[1][2] = lo[2];
  lamda2x(corners[1],corners[1]);
  corners[2][0] = lo[0]; corners[2][1] = hi[1]; corners[2][2] = lo[2];
  lamda2x(corners[2],corners[2]);
  corners[3][0] = hi[0]; corners[3][1] = hi[1]; corners[3][2] = lo[2];
  lamda2x(corners[3],corners[3]);
  corners[4][0] = lo[0]; corners[4][1] = lo[1]; corners[4][2] = hi[2];
  lamda2x(corners[4],corners[4]);
  corners[5][0] = hi[0]; corners[5][1] = lo[1]; corners[5][2] = hi[2];
  lamda2x(corners[5],corners[5]);
  corners[6][0] = lo[0]; corners[6][1] = hi[1]; corners[6][2] = hi[2];
  lamda2x(corners[6],corners[6]);
  corners[7][0] = hi[0]; corners[7][1] = hi[1]; corners[7][2] = hi[2];
  lamda2x(corners[7],corners[7]);
}
