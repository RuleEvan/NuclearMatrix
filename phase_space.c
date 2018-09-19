#include "phase_space.h"
double Fermi(int Z, double T) {
  // Computes the electron Fermi function for a given atomic number Z and electron energy T
  double s = sqrt(1 - pow(ALPHA_FS*Z, 2.0));
  double p = sqrt(pow(T, 2.0) - pow(M_ELECTRON, 2));
  double eta = ALPHA_FS*Z*T/p;
  double rho = R_NUC*pow(A_NUC, 1.0/3.0);
  gsl_sf_result a, b;

  double f = 2*(1+s)/pow(gsl_sf_gamma(1+2*s), 2.0);
//  double f = 4.0/pow(gsl_sf_gamma(1+2*s), 2.0);
  f *= pow(2*p*rho*8.065544*pow(10, -4), 2.0*s - 2.0);
  f *= exp(M_PI*eta);

  gsl_sf_lngamma_complex_e(s, eta, &a, &b);

  f *= exp(2.0*a.val);


  return f;
}

double phase_integrand_1(double E_1) {
  // Integrand for the type I phase space integral
  double E_2 = E_BETA - E_1;
  double p_1 = sqrt(pow(E_1, 2.0) - pow(M_ELECTRON, 2.0));
  double p_2 = sqrt(pow(E_2, 2.0) - pow(M_ELECTRON, 2.0));
  double phase_int = Fermi(Z_ATOM + 2, E_1)*Fermi(Z_ATOM + 2, E_2);
  phase_int *=  2.0*p_1*p_2*E_1*E_2;

  return phase_int;
}

double phase_integrand_2(double E_1) {
  // Integrand for the type II phase space integral
  double E_2 = E_BETA - E_1;
  double p_1 = sqrt(pow(E_1, 2.0) - pow(M_ELECTRON, 2.0));
  double p_2 = sqrt(pow(E_2, 2.0) - pow(M_ELECTRON, 2.0));
  double phase_int = Fermi(Z_ATOM + 2, E_1)*Fermi(Z_ATOM + 2, E_2);
  phase_int *=  p_1*p_2*pow((E_1 - E_2)/M_ELECTRON, 2.0)*(E_1*E_2 - pow(M_ELECTRON, 2.0));

  return phase_int;
}
double phase_integrand_3(double E_1) {
  // Integrand for the type II phase space integral
  double E_2 = E_BETA - E_1;
  double p_1 = sqrt(pow(E_1, 2.0) - pow(M_ELECTRON, 2.0));
  double p_2 = sqrt(pow(E_2, 2.0) - pow(M_ELECTRON, 2.0));
  double phase_int = Fermi(Z_ATOM + 2, E_1)*Fermi(Z_ATOM + 2, E_2);
  phase_int *=  2.0*p_1*p_2*pow(E_1 - E_2, 2.0);

  return phase_int;
}
double phase_integrand_4(double E_1) {
  // Integrand for the type II phase space integral
  double E_2 = E_BETA - E_1;
  double p_1 = sqrt(pow(E_1, 2.0) - pow(M_ELECTRON, 2.0));
  double p_2 = sqrt(pow(E_2, 2.0) - pow(M_ELECTRON, 2.0));
  double phase_int = Fermi(Z_ATOM + 2, E_1)*Fermi(Z_ATOM + 2, E_2);
  phase_int *=  2.0*p_1*p_2*(E_1*E_2 - pow(M_ELECTRON, 2.0));

  return phase_int;
}
double phase_integrand_6(double E_1) {
  // Integrand for the type II phase space integral
  double E_2 = E_BETA - E_1;
  double p_1 = sqrt(pow(E_1, 2.0) - pow(M_ELECTRON, 2.0));
  double p_2 = sqrt(pow(E_2, 2.0) - pow(M_ELECTRON, 2.0));
  double phase_int = Fermi(Z_ATOM + 2, E_1)*Fermi(Z_ATOM + 2, E_2);
  phase_int *=  4.0*p_1*p_2*M_ELECTRON*(E_1 + E_2);

  return phase_int;
}
double phase_integrand_9(double E_1) {
  // Integrand for the type II phase space integral
  double E_2 = E_BETA - E_1;
  double p_1 = sqrt(pow(E_1, 2.0) - pow(M_ELECTRON, 2.0));
  double p_2 = sqrt(pow(E_2, 2.0) - pow(M_ELECTRON, 2.0));
  double phase_int = Fermi(Z_ATOM + 2, E_1)*Fermi(Z_ATOM + 2, E_2);
  phase_int *=  4.0*p_1*p_2*(E_1*E_2 + pow(M_ELECTRON, 2.0));

  return phase_int;
}

double G_0(int k) {
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double g = 2.0*6.641*pow(10, -16)/pow(r_a, 2.0);
  double ep = pow(10, -6);
  if (k == 1) {
    g *= RombergIntegrator(phase_integrand_1, M_ELECTRON + ep, E_BETA - M_ELECTRON - ep, ep);
  } else if (k == 2) {
    g *= RombergIntegrator(phase_integrand_2, M_ELECTRON + ep, E_BETA - M_ELECTRON - ep, ep);
  } else if( k == 3) {
    g *= RombergIntegrator(phase_integrand_3, M_ELECTRON + ep, E_BETA - M_ELECTRON - ep, ep);
  } else if (k == 4) {
    g *= RombergIntegrator(phase_integrand_4, M_ELECTRON + ep, E_BETA - M_ELECTRON - ep, ep);
  } else if (k == 6) {
    g *= RombergIntegrator(phase_integrand_6, M_ELECTRON + ep, E_BETA - M_ELECTRON - ep, ep);
  } else if (k == 9) {
    g *= RombergIntegrator(phase_integrand_9, M_ELECTRON + ep, E_BETA - M_ELECTRON - ep, ep);
  } else {
    printf("Invalid phase argument\n");
    return 0.0;
  }
  return g;
}
