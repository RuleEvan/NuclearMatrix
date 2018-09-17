#ifndef ROMBERG_H
#define ROMBERG_H
#include "angular.h"
double RombergIntegrator(double (*f)(double), double a, double b, double tol); 
double Romberg3Vars(double (*f)(double, double, int), double a, double b, double p, int iv, double tol);
#endif
