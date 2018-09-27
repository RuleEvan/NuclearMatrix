#include "romberg.h"
double RombergIntegrator(double (*f)(double), double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(a) + (*f)(b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}

 double Romberg2Vars(double (*f)(double, double), double a, double b, double r, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(a, r) + (*f)(b, r));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(a + (k + k - 1)*h, r);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}
 

double Romberg3Vars(double (*f)(double, double, int), double a, double b, double p, int iv, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(a, p, iv) + (*f)(b, p, iv));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(a + (k + k - 1)*h, p, iv);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}

