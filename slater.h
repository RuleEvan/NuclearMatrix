#ifndef SLATER_H
#define SLATER_H
#include "lanczos.h"
int64_t p_step(int n_s, int n_p, int *m_p);
int64_t n_choose_k(int n, int k);
void orbitals_from_p(int64_t p, int n_s, int n_p, int* orbitals);
void generate_single_particle_states();
int64_t bin_from_p(int n_s, int n_p, int64_t p);
int m_from_p(int64_t p, int n_s, int n_p, int* m_shell);
int64_t a_op(int n_s, int n_p, int64_t p, int n_op, int* phase);
int64_t a_dag_op(int n_s, int n_p, int64_t p, int n_op, int* phase);
int64_t a_dag_a_op(int n_s, int n_p, int64_t p, int n_a, int n_b, int* phase);
int64_t bin_phase_from_p(int n_s, int n_p, int64_t p, int n_op, int* phase);
int64_t bin_from_orbitals(int n_s, int n_p, int *orbitals);
void initialize_orbitals(int* n_shell, int* l_shell, int* j_shell, int* m_shell);
unsigned int count_set_bits(unsigned int n);
int64_t p_from_binary(int n_s, int n_p, int64_t b);
void orbitals_from_binary(int n_s, int n_p, int64_t b, int *orbitals);
void generate_binomial_file();
#endif
