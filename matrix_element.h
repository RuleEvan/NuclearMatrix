#ifndef MATRIX_ELEMENT_H
#define MATRIX_ELEMENT_H
#include "potential.h"
double compute_matrix_element_M(int m_sw, int iv);
double compute_matrix_element_M_GT();
double compute_matrix_element_M_F();
double compute_matrix_element_TT(int iv);
double compute_matrix_element_M_0();
double compute_matrix_element();
double matrix_element_sigma_dot_sigma(double s, double ms, double t, double sp, double msp, double tp);
double compute_matrix_element_av18(int n1p, int l1p, int n2p, int l2p, int lambdap, int mup, int sp, int msp, int tp, int mtp, int n1, int l1, int n2, int l2, int lambda, int mu, int s, int ms, int t, int mt);
double compute_matrix_element_c(int n1p, int l1p, int n2p, int l2p, int lambdap, int mup, int sp, int msp, int tp, int mtp, int n1, int l1, int n2, int l2, int lambda, int mu, int s, int ms, int t, int mt);

#endif
