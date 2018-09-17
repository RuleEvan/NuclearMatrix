#ifndef PHASE_SPACE_H
#define PHASE_SPACE_H
#include "romberg.h"
double phase_integrand_1(double E);
double phase_integrand_2(double E);
double Fermi(int Z, double T);
#endif
