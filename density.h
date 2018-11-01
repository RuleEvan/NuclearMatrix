#ifndef DENSITY_H
#define DENSITY_H
#include "wave_function.h"

typedef struct sBasisCoeff 
{
  int64_t p;
  int64_t b;
  double *wave;
} BasisCoeff;

typedef struct wfnData
{
  int n_proton_i, n_proton_f, n_neutron_i, n_neutron_f;
  int n_states_i, n_states_f, n_shells, n_orbits, n_data;
  int n_eig_i, n_eig_f;
  BasisCoeff *bc_i, *bc_f;
  int *n_shell, *l_shell, *j_shell, *jz_shell, *tz_shell;
  int *n_orb, *l_orb;
  double *j_orb;
  double jz_i, jz_f;
  double *e_nuc_i, *e_nuc_f, *j_nuc_i, *j_nuc_f, *t_nuc_i, *t_nuc_f;

} wfnData;

void one_body_density(int j_op, int t_op);
void two_body_density(int j_op, int t_op);
void exp_from_wfn();
wfnData* read_wfn_data();

int s_compare(const void * a, const void * b);

#endif
