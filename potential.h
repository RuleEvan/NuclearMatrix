#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "phase_space.h"
double b_coeff(double n, double l, double np, double lp, double p);
double compute_potential(double n, double np, double l, double lp, int iv);
double a_coeff(double n, double l, double k);
double v_light_limit(double q);
double v_pion_f1(double r);
double v_pion_f2(double r);
double v_light_limit_d(double r);
double v_pion_g1(double r);
double v_pion_g2(double r);
double v_pion_NN_g1(double r);
double v_pion_NN_g2(double r);
double talmi(double p, int iv);
double talmi_integrand(double q, double p, int iv);
void potential_spline_init(double (*f)(double, double), double r_min, double r_max, int r_steps);
double v_spline(double r);
double g_axial(double q_sq);
double g_vector(double q_sq);
double g_magnetic(double q_sq);
double g_pseudo(double q_sq);
double h_AA_GT_q(double q, double r);
double h_AA_T_q(double q, double r);
double h_AP_GT_q(double q, double r);
double h_AP_T_q(double q, double r);
double h_PP_GT_q(double q, double r);
double h_PP_T_q(double q, double r);
double h_MM_GT_q(double q, double r);
double h_MM_T_q(double q, double r);
double h_F_q(double q, double r);
double h_AA_GT_q_sd(double q, double r);
double h_AA_T_q_sd(double q, double r);
double h_AP_GT_q_sd(double q, double r);
double h_AP_T_q_sd(double q, double r);
double h_PP_GT_q_sd(double q, double r);
double h_PP_T_q_sd(double q, double r);
double h_MM_GT_q_sd(double q, double r);
double h_MM_T_q_sd(double q, double r);
double h_F_q_sd(double q, double r);
double h_AA_GT(double r);
double h_AA_T(double r);
double h_AP_GT(double r);
double h_AP_T(double r);
double h_PP_GT(double r);
double h_PP_T(double r);
double h_F(double r);
#endif
