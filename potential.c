#include "potential.h"
double a_coeff(double n, double l, double k) {
  double a = pow(-1.0, k)/gsl_sf_gamma(k + 1.0);
  a *= sqrt(2.0*gsl_sf_gamma(n + 1.0)/gsl_sf_gamma(n + l + 1.5))*gsl_sf_gamma(n + l + 1.5)/(gsl_sf_gamma(n - k + 1.0)*gsl_sf_gamma(k + l + 1.5));
  
  return a;
}

double b_coeff(double n, double l, double np, double lp, double p) {
  double b = 0.5*gsl_sf_gamma(p + 1.5);
  int k, kp;
  double b_sum = 0.0;
  for (k = 0; k <= n; k++) {
    for (kp = 0; kp <= np; kp++) {
      if (0.5*(2*k + 2*kp + l + lp) != p) {continue;}
      b_sum += a_coeff(n,l,k)*a_coeff(np, lp, kp);
    }
  }
  b *= b_sum;

  return b;
}

double compute_potential(double n, double np, double l, double lp, int iv) {
  // Sum the required Talmi integrals with the corresponding B coefficients
  int p;
  double v;
  for (p = (int)((l + lp)/2); p <= (int)((l + lp)/2 + n + np); p++) {
    v += b_coeff(n, l, np, lp, p)*talmi(p, iv);
  }
  return v;
}

double talmi(double p, int iv) {
  // Compute the order p Talmi integral
  // Set the limits of the integral and the error tolerance
  double r_min = 0.001;
  double r_max = 5.0;
  double tol = pow(10, -6);
  double I_p = Romberg3Vars(&talmi_integrand, r_min, r_max, p, iv, tol);
  I_p *= 2.0/gsl_sf_gamma(p + 1.5);
  
  return I_p;
}
    
double talmi_integrand(double q, double p, int iv) {
  // The integrand of the Talmi integral
  // Plug in the required potential here
  double v = pow(q, 2.0*p + 2.0)*exp(-q*q);
  if (COR_FAC == 1) {
    double beta = exp(-1.1*pow(B_OSC*q, 2))*(1.0 - 0.68*pow(B_OSC*q,2.0));
    v *= pow(1.0 - beta, 2.0);
  }
  if (iv == 1) {
    v *= v_light_limit(q);
  } else if (iv == 2) {
    v *= v_pion_f1(q);
  } else if (iv == 3) {
    v *= v_pion_f2(q);
  } else if (iv == 4) {
    v *= v_light_limit_d(q);
  } else if (iv == 5) {
    v *= v_av18_c(q);
  } else if (iv == 6) {
    v *= v_av18_tau(q);
  } else if (iv == 7) {
    v *= v_av18_sigma(q);
  } else if (iv == 8) {
    v *= v_av18_sigma_tau(q);
  } else if (iv == 9) {
    v *= v_av18_t(q);
  } else if (iv == 10) {
    v *= v_av18_t_tau(q);
  } else if (iv == 11) {
    v *= v_av18_ls(q);
  } else if (iv == 12) {
    v *= v_av18_ls_tau(q);
  } else if (iv == 13) {
    v *= v_av18_l2(q);
  } else if (iv == 14) {
    v *= v_av18_l2_tau(q);
  } else if (iv == 15) {
    v *= v_av18_l2_sigma(q);
  } else if (iv == 16) {
    v *= v_av18_l2_sigma_tau(q);
  } else if (iv == 17) {
    v *= v_av18_ls2(q);
  } else if (iv == 18) {
    v *= v_av18_ls2_tau(q);
  }
  return v;
}

double v_light_limit(double r) {
  // Nuclear potential in the case of light neutrinos 
  // Lepton kinematics are ignored (see Haxton review sec. 3.4.2)
  double lb = 0.0050677*B_OSC*A_FACTOR;
  double v = 2.0/M_PI/B_OSC;   
  v *= 1/r;
  v *= gsl_sf_Ci(r*lb)*sin(r*lb) + (M_PI/2.0 - gsl_sf_Si(r*lb))*cos(r*lb);

  return v;
}

double v_light_limit_d(double r) {
  // Nuclear potential in the case of light neutrinos 
  // Lepton kinematics are ignored (see Haxton review sec. 3.4.2)
  double lb = 0.0050677*B_OSC*A_FACTOR;
  double v = (sin(lb*r)*gsl_sf_Ci(lb*r) + cos(lb*r)*(M_PI/2.0 - gsl_sf_Si(lb*r)))/(M_PI*r*B_OSC) - lb/B_OSC*(cos(lb*r)*gsl_sf_Ci(lb*r)-sin(lb*r)*(M_PI/2.0 - gsl_sf_Si(lb*r)))/M_PI;   
  v *= 2.0;

  return v;
}

double v_pion_f1(double r) {
  r *= B_OSC;
  double v = 1/r;
  double x = r*PION_MASS*0.0050677;
  v *= (x - 2.0)*exp(-x);
  
  return v;  
}

double v_pion_f2(double r) {
  r *= B_OSC;
  double v = 1/r;
  double x = r*PION_MASS*0.0050677;
  v *= (x + 1.0)*exp(-x);
  
  return v;  
}

double v_pion_g1(double r) {
  r *= B_OSC;
  double v = 1/pow(r, 3.0);
  double x = r*PION_MASS*0.0050677;
  v *= -1.0*pow(x, 2.0)/3.0*(4.0 - x)*exp(-x);
  
  return v;
}

double v_pion_g2(double r) {
  r *= B_OSC;
  double v = 1/pow(r, 3.0);
  double x = r*PION_MASS*0.0050677;
  v *= -1.0*(2.0 + 2.0*x + 1.0/3.0*(x*x - x*x*x))*exp(-x);

  return v;
}

double v_pion_NN_g1(double r) {
  r *= B_OSC;
  double v = 1/pow(r, 3.0);
  double x = r*PION_MASS*0.0050677;
  v *= -1.0/3.0*x*x*exp(-x);

  return v;
}

double v_pion_NN_g2(double r) {
  r *= B_OSC;
  double v = 1/pow(r, 3.0);
  double x = r*PION_MASS*0.0050677;
  v *= -1.0*(1.0 + x + 1.0/3.0*x*x)*exp(-x);

  return v;
}

