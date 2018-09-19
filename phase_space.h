#ifndef PHASE_SPACE_H
#define PHASE_SPACE_H
#include "romberg.h"
double phase_integrand_1(double E);
double phase_integrand_2(double E);
double phase_integrand_3(double E);
double phase_integrand_4(double E);
double phase_integrand_6(double E);
double phase_integrand_9(double E);
double G_0(int k);
double Fermi(int Z, double T);
#endif
