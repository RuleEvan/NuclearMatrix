#ifndef DENSITY_H
#define DENSITY_H
#include "wave_function.h"

typedef struct sBasisCoeff 
{
  int p;
  double *wave;
} BasisCoeff;

void one_body_density();
void two_body_density(int j_op, int t_op);
void exp_from_wfn();

int s_compare(const void * a, const void * b);

#endif
