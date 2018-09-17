#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "stdint.h"

#define ALPHA_FS 0.007297352
#define R_NUC 1.1 // [fm]
#define A_NUC 76 // Atomic Mass
#define M_NUC 2.19 // nuclear matrix element
#define M_ELECTRON 0.511 // [MeV]
#define M_NEUTRON 939.57 // [MeV]
#define Z_ATOM 52 // Atomic Number
#define E_BETA 3.038 // [MeV]
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))
#define XI_MIX_0 0.0016 // Mean mixing angle
#define G_AXIAL 1.275 // Axial coupling
#define E_HADRON 1 // Hadronization scale [GeV] 
#define DELTA_SUPP 1.0/30.0 // Suppression factor
#define M_W1 80.0 // W1 gauge boson mass [GeV]
#define M_W2_MIN 0.715 // Minimum W2 gauge boson mass [TeV]

#define DENSITY_FILE "GE76_DENSITY.DAT"
#define NUM_SHELLS 99
#define POTENTIAL 1
#define A_FACTOR 9.155 // [MeV] Average nuclear excitation energy
#define B_OSC 2.927 // [fm] Oscillator parameter (times sqrt(2))
#define COR_FAC 1
#define PION_MASS 139.57 // [MeV] Charged pion mass
 
double Fermi(int Z, double T);
double RombergIntegrator(double (*f)(double), double a, double b, double tol); 
double Romberg3Vars(double (*f)(double, double, int), double a, double b, double p, int iv, double tol);
double talmi(double p, int iv);
double talmi_integrand(double q, double p, int iv);
double phase_integrand_1(double E);
double phase_integrand_2(double E);
double three_j(double j1, double j2, double j3, double m1, double m2, double m3);
double six_j(double j1, double j2, double j3, double j4, double j5, double j6);
double nine_j(double j11, double j12, double j13, double j21, double j22, double j23, double j31, double j32, double j33);
double triangle(double a, double b, double c);
double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m);
double g_factor(double a, double b, double c);
double h_factor(double a, double b, double c);
double brody_mosh_zero(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int l1, int l2); 
double brody_mosh(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int n1, int l1, int n2, int l2);
double compute_matrix_element();
double matrix_element_sigma_dot_sigma(double s, double ms, double t, double sp, double msp, double tp);
double b_coeff(double n, double l, double np, double lp, double p);
double compute_potential(double n, double np, double l, double lp, int iv);
double a_coeff(double n, double l, double k);
double mat_gt(double s, double sp);
double mat_f(double s, double sp);
double compute_matrix_element_M(int m_sw, int iv);
double compute_matrix_element_M_GT();
double compute_matrix_element_M_F();
double compute_matrix_element_TT(int iv);
double compute_matrix_element_M_0();
double v_light_limit(double q);
double v_pion_f1(double r);
double v_pion_f2(double r);
double v_light_limit_d(double r);
double v_pion_g1(double r);
double v_pion_g2(double r);
double v_pion_NN_g1(double r);
double v_pion_NN_g2(double r);
void square_matrix_times_vector(double *mat, double *vec, double *result, int dim);
void vector_times_scalar(double *v1, double c1, int dim);
void vector_difference(double *v1, double *v2, double *v3, int dim);
void lanczos(double *h, double *alpha, double *beta, int dim, int n_iter);
double vector_dot_product(double *v1, double *v2, int dim);
void matrix_elements(double *h, double *h_tri, int dim, int n_iter, double* alpha, double* beta);
void get_orthogonal_vector(double *h_tri, double *v, int dim, int n_done);
double characteristic_polynomial(double lambda, double *alpha, double *beta, int dim, int *num_sign);
void vector_max(double* v, double* max, int dim, int* n_max);
void vector_min(double* v, double* min, int dim, int* n_min);
void lanczos_eigenvalues(double*h, int dim, int n_iter);
int agree(double a, double b);
int p_step(int n_s, int n_p, int *m_p);
int n_choose_k(int n, int k);
void orbitals_from_p(int p, int n_s, int n_p, int* orbitals);
void generate_single_particle_states();
int bin_from_p(int n_s, int n_p, int p);
int a_op(int n_s, int n_p, int p, int n_op);
int a_dag_op(int n_s, int n_p, int p, int n_op);
int bin_phase_from_p(int n_s, int n_p, int p, int n_op, int* phase);
