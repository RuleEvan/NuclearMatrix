#ifndef BRODY_H
#define BRODY_H
#include "phase_space.h"
double g_factor(double a, double b, double c);
double h_factor(double a, double b, double c);
double brody_mosh_zero(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int l1, int l2); 
double brody_mosh(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int n1, int l1, int n2, int l2);
#endif
