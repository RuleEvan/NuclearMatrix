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
double compute_matrix_element_av18(int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_c(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_tau(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_sigma(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_sigma_tau(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_ls(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_ls_tau(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_l2(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_l2_tau(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_l2_sigma(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);
double compute_matrix_element_l2_sigma_tau(int iv, int n1p, double j1p, int n2p, double j2p, int j12p, int t12p, int n1, double j1, int n2, double j2, int j12, int t12);

#endif
