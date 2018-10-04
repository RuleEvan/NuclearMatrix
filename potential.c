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
  if (iv < 0) {v *= v_spline(q);}
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

double g_vector(double q_sq) {
  double g = pow(1.0 + q_sq/pow(LAMBDA_V, 2.0), -2.0);

  return g;
}

double g_axial(double q_sq) {
  double g = pow(G_AXIAL, 2.0)*pow(1.0 + q_sq/pow(LAMBDA_A, 2.0), -2.0);

  return g;
}

double g_magnetic(double q_sq) {
  double g = (1.0 + KAPPA_1)*g_vector(q_sq);
 
  return g;
}

double g_pseudo(double q_sq) {
  double g = -2.0*M_NEUTRON*g_axial(q_sq)/(q_sq + pow(PION_MASS, 2.0));

  return g;
}

double h_AA_GT_q(double q, double r) {
  double h = pow(g_axial(q*q), 2.0)/pow(G_AXIAL, 2.0)*gsl_sf_bessel_j0(q*r*0.00507);

  return h;
}

double h_AA_T_q(double q, double r) {
  double h = pow(g_axial(q*q), 2.0)/pow(G_AXIAL, 2.0)*gsl_sf_bessel_j2(q*r*0.00507);

  return h;
}

double h_AP_GT_q(double q, double r) {
  double h = g_pseudo(q*q)/(3*M_NEUTRON*pow(G_AXIAL, 2.0))*g_axial(q*q)*q*q*gsl_sf_bessel_j0(q*r*0.00507);

  return h;
}

double h_AP_T_q(double q, double r) {
  double h = -g_pseudo(q*q)/(3*M_NEUTRON*pow(G_AXIAL, 2.0))*g_axial(q*q)*q*q*gsl_sf_bessel_j2(q*r*0.00507);

  return h;
}
 
double h_PP_GT_q(double q, double r) {
  double h = pow(g_pseudo(q*q), 2.0)*pow(q, 4.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*12)*gsl_sf_bessel_j0(q*r*0.00507);

  return h;
}

double h_PP_T_q(double q, double r) {
  double h = -pow(g_pseudo(q*q), 2.0)*pow(q, 4.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*12)*gsl_sf_bessel_j2(q*r*0.00507);

  return h;
}

double h_MM_GT_q(double q, double r) {
  double h = pow(g_magnetic(q*q), 2.0)*pow(q, 2.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*6)*gsl_sf_bessel_j0(q*r*0.00507);

  return h;
}

double h_MM_T_q(double q, double r) {
  double h = pow(g_magnetic(q*q), 2.0)*pow(q, 2.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*12)*gsl_sf_bessel_j2(q*r*0.00507);

  return h;
}

double h_F_q(double q, double r) {
  double h = g_vector(q*q)*gsl_sf_bessel_j0(q*r*0.00507);

  return h;
}

double h_AA_GT_q_sd(double q, double r) {
  double h = pow(g_axial(q*q), 2.0)/pow(G_AXIAL, 2.0)*gsl_sf_bessel_j0(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}

double h_AA_T_q_sd(double q, double r) {
  double h = pow(g_axial(q*q), 2.0)/pow(G_AXIAL, 2.0)*gsl_sf_bessel_j2(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}

double h_AP_GT_q_sd(double q, double r) {
  double h = g_pseudo(q*q)/(3*M_NEUTRON*pow(G_AXIAL, 2.0))*g_axial(q*q)*q*q*gsl_sf_bessel_j0(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}

double h_AP_T_q_sd(double q, double r) {
  double h = -g_pseudo(q*q)/(3*M_NEUTRON*pow(G_AXIAL, 2.0))*g_axial(q*q)*q*q*gsl_sf_bessel_j2(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}
 
double h_PP_GT_q_sd(double q, double r) {
  double h = pow(g_pseudo(q*q), 2.0)*pow(q, 4.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*12)*gsl_sf_bessel_j0(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}

double h_PP_T_q_sd(double q, double r) {
  double h = -pow(g_pseudo(q*q), 2.0)*pow(q, 4.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*12)*gsl_sf_bessel_j2(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}

double h_MM_GT_q_sd(double q, double r) {
  double h = pow(g_magnetic(q*q), 2.0)*pow(q, 2.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*6)*gsl_sf_bessel_j0(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}

double h_MM_T_q_sd(double q, double r) {
  double h = pow(g_magnetic(q*q), 2.0)*pow(q, 2.0)/(pow(G_AXIAL*M_NEUTRON, 2.0)*12)*gsl_sf_bessel_j2(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}

double h_F_q_sd(double q, double r) {
  double h = g_vector(q*q)*gsl_sf_bessel_j0(q*r*0.00507)*pow(q/PION_MASS, 2.0);

  return h;
}


double h_AA_GT(double r) {
  r *= B_OSC;
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double v = r_a/r*pow(G_AXIAL, 2.0);

  return v;
}

double h_AA_T(double r) {
  r *= B_OSC;
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double v = r_a/(2*r)*pow(G_AXIAL, 2.0);
 
  return v;
}

double h_AP_GT(double r) {
  r *= B_OSC;
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double x = r*PION_MASS*0.0050677;
  double v = -2*r_a/(3*r)*pow(G_AXIAL, 2.0)*exp(-x);

  return v;
}

double h_AP_T(double r) {
  r *= B_OSC;
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double x = r*PION_MASS*0.0050677;
  double v = r_a/(r*x*x)*pow(G_AXIAL, 2.0)*(2.0/3.0)*exp(-x);
  v *= (-3.0 + 3*exp(x) - 3*x - x*x);

  return v;
}

double h_PP_GT(double r) {
  r *= B_OSC;
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double x = r*PION_MASS*0.0050677;
  double v = -r_a/(6*r)*pow(G_AXIAL, 2.0)*exp(-x);
  v *= (x-2);

  return v;
}

double h_PP_T(double r) {
  r *= B_OSC;
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double x = r*PION_MASS*0.0050677;
  double v = -r_a/(6*r)*pow(G_AXIAL, 2.0)*exp(-x);
  v *= (1+x);

  return v;
}

double h_F(double r) {
  r *= B_OSC;
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);
  double v = r_a/r*pow(G_AXIAL, 2.0);

  return v;
}

gsl_interp_accel *potential_acc;
gsl_spline *potential_spline;

void potential_spline_init(double (*f)(double, double), double r_min, double r_max, int r_steps) {
  potential_acc = gsl_interp_accel_alloc();
  potential_spline = gsl_spline_alloc(gsl_interp_cspline, r_steps);
  printf("Initializing spline\n");
  double* r_arr = (double*)malloc(sizeof(double)*r_steps);
  double* v_arr = (double*)malloc(sizeof(double)*r_steps);
  double dr = (r_max - r_min)/r_steps; 
  double r_a = R_NUC*pow(A_NUC, 1.0/3.0);

  for (int i = 0; i < r_steps; i++) {
    double r = r_min + i*dr;
    r_arr[i] = r;
    v_arr[i] = 0.00507*2.0*r_a/M_PI*Romberg2Vars(f, 0.001, 10000.0, r, pow(10, -4));
  }
  gsl_spline_init(potential_spline, r_arr, v_arr, r_steps);
  
  return;
}
   
double v_spline(double r) {
  r *= B_OSC;
  double v = gsl_spline_eval(potential_spline, r, potential_acc);
  return v;
} 
