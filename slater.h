#ifndef SLATER_H
#define SLATER_H
#include "lanczos.h"
int p_step(int n_s, int n_p, int *m_p);
int n_choose_k(int n, int k);
void orbitals_from_p(int p, int n_s, int n_p, int* orbitals);
void generate_single_particle_states();
int bin_from_p(int n_s, int n_p, int p);
int m_from_p(int p, int n_s, int n_p, int* m_shell);
int a_op(int n_s, int n_p, int p, int n_op, int* phase);
int a_dag_op(int n_s, int n_p, int p, int n_op, int* phase);
int bin_phase_from_p(int n_s, int n_p, int p, int n_op, int* phase);
void initialize_orbitals(int* n_shell, int* l_shell, int* j_shell, int* m_shell);
unsigned int count_set_bits(unsigned int n);
int p_from_binary(int n_s, int n_p, int b);
void orbitals_from_binary(int n_s, int n_p, int b, int *orbitals);
#endif
