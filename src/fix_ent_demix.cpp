
#include "fix_ent_demix.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define BIG MAXTAGINT

enum {ATOM,MOLECULE};


/* ---------------------------------------------------------------------- */

FixEntDemix::FixEntDemix(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR, "Illegal fix manyb command");

  alpha  = utils::numeric(FLERR, arg[3], false, lmp);

//   A_type  = utils::numeric(FLERR, arg[3], false, lmp);
//   B_type = utils::numeric(FLERR, arg[4], false, lmp);

  if (comm->me == 0) {
    fprintf(screen, "\n <<<<<<<< Demix settings >>>>>>>>>> \n");
    //fprintf(screen, "Demix atom types   : (%d, %d)\n", A_type, B_type);
  }

  //comm_forward = 7;  // For q_M 
  comm_reverse = 1;  // For coordination numbers and derivatives

  // Add these flags for proper NPT handling
  restart_global = 1;  // Support restart with changed box
  box_change     = 1;  // Tell LAMMPS this fix is aware of box changes
}

/* ---------------------------------------------------------------------- */

FixEntDemix::~FixEntDemix()
{

}

/* ---------------------------------------------------------------------- */

int FixEntDemix::setmask()
{
  int mask = 0;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::POST_FORCE;
  mask |= FixConst::MIN_PRE_FORCE;
  mask |= FixConst::POST_NEIGHBOR;
  mask |= FixConst::POST_INTEGRATE;
  mask |= FixConst::END_OF_STEP;
  mask |= FixConst::PRE_EXCHANGE;
  return mask;
}

void FixEntDemix::setup(int /*vflag*/)
{

}

void FixEntDemix::init()
{
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL)->set_cutoff((RC + DELTA_MANY_B)  + neighbor->skin); //  NeighConst::REQ_OCCASIONAL


}

/* ---------------------------------------------------------------------- */

void FixEntDemix::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixEntDemix::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEntDemix::post_neighbor()
{
  //fprintf(screen, "post_neighbor\n");
}

void FixEntDemix::post_integrate()
{

}

// 2. Enhanced pre_exchange to handle box deformations properly
void FixEntDemix::pre_exchange()
{

}


// 3. Enhanced end_of_step with proper box change handling
void FixEntDemix::end_of_step()
{

}

void FixEntDemix::pre_force(int vflag)
{

  neighbor->build_one(list);

  double *N = atom->N;
  //double *N_d = atom->N_d;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  for (int i = 0; i < nlocal + nghost; i++) { //nghost
    N[i] = 0.0;
    //N_d[i][0] = 0.0;
  }
  //comm->reverse_comm(this);
  //ncoord_derivative();

  count_bonds();
  comm_forward = 1;
  comm->forward_comm(this);

  ncoord_correction();
  comm_forward = 1;
  comm->forward_comm(this);
}

void FixEntDemix::post_force(int vflag)
{

}

/* ---------------------------------------------------------------------- */

void FixEntDemix::pre_force_respa(int vflag, int ilevel, int /*iloop*/)
{
}
void FixEntDemix::count_bonds()
{
  // Get number of owned atoms and ghost atoms.
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;

  double **x         = atom->x;
  int *type          = atom->type;
  double *N          = atom->N;
  //double *N_d        = atom->N_d;

  int newton_pair = force->newton_pair;
  double *special_coul = force->special_coul;
  double factor_coul;

  // Get the neighbor list data.
  int inum = list->inum; //  + list->gnum loop over owned and ghost atoms
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // Loop over all atoms in the full neighbor list.
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    double xi = x[i][0], yi = x[i][1], zi = x[i][2];

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      int type_j = type[j];
      double xj = x[j][0], yj = x[j][1], zj = x[j][2];
      double dx = xi - xj;
      double dy = yi - yj;
      double dz = zi - zj;
      // *** Apply minimum image correction ***
      //domain->minimum_image(dx, dy, dz);
      
      double r = sqrt(dx*dx + dy*dy + dz*dz);
      // fc(r) is the cutoff switching function used to count bonds.
      double fc_val = factor_coul * fc(r);
      //if (fc_val == 0.0) continue;
      double fcd_val = factor_coul * fc_d(r);

      // Compute derivative components for fc:
      double dfx = fcd_val * (dx / r);
      double dfy = fcd_val * (dy / r);
      double dfz = fcd_val * (dz / r);
      
      // In this model we only consider O–H interactions.
      // There are two cases: either i is oxygen and j is hydrogen, or vice-versa.
      if (type[i] == type[j]) {
        // Oxygen i: add fc_val to its coordination.
        N[i] += fc_val;
        // // For the derivatives
        // N_d[i][0] += dfx;
        // N_d[i][1] += dfy;
        // N_d[i][2] += dfz;
      }
    }
  }
}

double FixEntDemix::fc(double r)
{
  if (r < RC - DELTA_MANY_B) return 1.0;
  else if (r > RC + DELTA_MANY_B) return 0.0;
  else {
    return 0.5 * (1.0 + sin((M_PI/2.0)*((RC - r) / DELTA_MANY_B)));
  }
}

double FixEntDemix::fc_d(double r)
{
  if (r < RC - DELTA_MANY_B || r > RC + DELTA_MANY_B) return 0.0;
  else {
    return - (M_PI / (4.0 * DELTA_MANY_B)) * cos((M_PI/2.0)*((RC - r) / DELTA_MANY_B));
  }
}

void FixEntDemix::ncoord_correction()
{
  int nlocal      = atom->nlocal, nghost = atom->nghost;
  int *type       = atom->type;
  double *N       = atom->N;

  double N_coord_deviation  = 3;
  double N_coord_min =  0;//2 - deviation_oxygen;

  for (int i = 0; i < nlocal + nghost; i++) {
      N[i] = (N[i] < N_coord_min ? N_coord_min : N[i]);
      N[i] = (N[i] >= N_coord_deviation ? N_coord_deviation : N[i]);
      //N_d[i] *= 0;
  }
}


double FixEntDemix::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)sizeof(double)*3*atom->nghost;
  return bytes;
}

/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

int FixEntDemix::pack_forward_comm(int n, int *list, double *buf,
                                int /*pbc_flag*/, int * /*pbc*/)
{
  int i, idx, m;
  m = 0;
  for (i = 0; i < n; i++) {
    idx = list[i];
    //   buf[m++] = atom->N_d[idx][0];
    //   buf[m++] = atom->N_d[idx][1];
    //   buf[m++] = atom->N_d[idx][2];
      buf[m++] = atom->N[idx];
  }
  return m;
}

void FixEntDemix::unpack_forward_comm(int n, int first, double *buf)
{
  int i, last, m;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
      atom->N[i] = buf[m++];
    //   atom->ncoord[i][0] = buf[m++];
    //   atom->ncoord[i][1] = buf[m++];
    //   atom->ncoord[i][2] = buf[m++];
  }
}

int FixEntDemix::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  double *N = atom->N;
  //double **N_d = atom->N_d;
  // Only pack 3 ncoord and 3 ncoord_d, removing ncoord[i][3]
  for (int i = first; i < first + n; i++) {
    buf[m++] = N[i];
    // buf[m++] = N_d[i][0];
    // buf[m++] = N_d[i][1];
    // buf[m++] = N_d[i][2];
  }
  return m;
}

void FixEntDemix::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  double *N = atom->N;
  //double **N_d = atom->N_d;
  for (int i = 0; i < n; i++) {
    int idx = list[i];
    N[idx] += buf[m++];

    // N_d[idx][0] += buf[m++];
    // N_d[idx][1] += buf[m++];
    // N_d[idx][2] += buf[m++];
  }
}