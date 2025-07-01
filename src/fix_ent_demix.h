
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(manyb,FixEntDemix)

#else

#ifndef LMP_FIX_ENT_DEMIX_H
#define LMP_FIX_ENT_DEMIX_H

#include "fix.h"
#include "pair.h"

namespace LAMMPS_NS {

class FixEntDemix : public Fix {
 friend class Neighbor; 
 public:
  FixEntDemix(class LAMMPS *, int, char **);
  ~FixEntDemix() override;
  int setmask() override;
  void setup(int ) override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void pre_force(int) override;
  void pre_exchange() override;
  void post_force(int) override;
  void pre_force_respa(int, int, int) override;
  void min_pre_force(int) override;
  void post_neighbor() override;
  void post_integrate() override;
  void end_of_step() override;
  double memory_usage() override;


  double alpha;

  int comm_flag;

 protected:
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  NeighList *list;
  double RC, DELTA_MANY_B;
  double CHI_MANYB, LAMBDA_MANYB;
  double fc(double r);
  double fc_d(double r);
  void count_bonds();
  void charges_shift();
  void ncoord_correction();

   // Add this function in your class or as a helper function
   inline int sbmask(int j) {
   return (j >> SBBITS) & 3;
   }

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix indent does not exist

Self-explanatory.

E: Variable for fix indent is invalid style

Only equal-style variables can be used.

E: Variable for fix indent is not equal style

Only equal-style variables can be used.

*/
