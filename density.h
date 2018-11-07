#ifndef DENSITY_H
#define DENSITY_H
#include "wave_function.h"

typedef struct sBasisCoeff 
{
  double *wave;
} BasisCoeff;

typedef struct sd_list
{
  int pi, pn;
  int phase;
  struct sd_list *next;
} sd_list;

typedef struct wf_list
{
  unsigned int p;
  unsigned int index;
  struct wf_list *next;
} wf_list;

sd_list* create_sd_node(int pi, int pf, int phase, sd_list* next);
sd_list* sd_append(sd_list* head, int pi, int pf, int phase);
wf_list* create_wf_node(unsigned int b, unsigned int index, wf_list* next);
wf_list* wf_append(wf_list* head, unsigned int b, unsigned int index);

typedef struct wfnData
{
  int n_proton_i, n_proton_f, n_neutron_i, n_neutron_f;
  int n_states_i, n_states_f, n_shells, n_orbits, n_data;
  int n_eig_i, n_eig_f;
  int n_sds_n_i, n_sds_n_f, n_sds_p_i, n_sds_p_f;
  wf_list **p_hash_i, **p_hash_f, **n_hash_i, **n_hash_f;
  BasisCoeff *bc_i, *bc_f;
  int *n_shell, *l_shell, *j_shell, *jz_shell, *tz_shell;
  int *n_orb, *l_orb;
  double *j_orb;
  double jz_i, jz_f;
  double *e_nuc_i, *e_nuc_f, *j_nuc_i, *j_nuc_f, *t_nuc_i, *t_nuc_f;
  unsigned int *p_sd_i, *p_sd_f, *n_sd_i, *n_sd_f;
  int same_basis;
} wfnData;

void one_body_density(int j_op, int t_op);
void two_body_density(int j_op, int t_op);
wfnData* read_wfn_data();


#endif
