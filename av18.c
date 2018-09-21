#include "av18.h"
double woods_saxon(double r) {
  double a = 0.2; // [fm]
  double r_0 = 0.5; // [fm]
  double w = 1.0 + exp((r - r_0)/a);
  w = 1.0/w;
  
  return w;
}

double t_mu(double r) {
  double x = MU*r;
  double c = 2.1; // [fm]^-2
  double t= (1.0 + 3.0/x + 3.0/pow(x, 2.0));
  t *= exp(-x)/x;
  t *= pow(1.0 - exp(-c*r*r), 2.0);
  
  return t;
}


double v_av18_c_01_pp(double r) {
  double I = -11.27028;
  double P = 3346.6874;
  double Q = 1859.5627;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_c_01_np(double r) {
  double I = -10.66788;
  double P = 3126.5542;
  double Q = 1746.4298;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_c_01_nn(double r) {
  double I = -11.27028;
  double P = 3342.7664;
  double Q = 1857.4367;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_l2_01(double r) {
  double I = 0.12472;
  double P = 16.7780;
  double Q = 9.0972;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_c_00(double r) {
  double I = -2.09971;
  double P = 1204.4301;
  double Q = 511.9380;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_l2_00(double r) {
  double I = -0.31452;
  double P = 217.4559;
  double Q = 117.9063;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_c_11_pp(double r) {
  double I = -7.62701;
  double P = 1815.4920;
  double Q = 969.3863;
  double R = 1847.8059;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_c_11_np(double r) {
  double I = -7.62701;
  double P = 1813.5315;
  double Q = 966.2483;
  double R = 1847.8059;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_c_11_nn(double r) {
  double I = -7.62701;
  double P = 1811.5710;
  double Q = 967.2603;
  double R = 1847.8059;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_l2_11(double r) {
  double I = 0.06709;
  double P = 342.0669;
  double Q = 185.4713;
  double R = -615.2339;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_t_11(double r) {
  double I = 1.07985;
  double P = 0.0;
  double Q = -190.0949;
  double R = -811.2040;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_ls_11(double r) {
  double I = -0.62697;
  double P = -570.5571;
  double Q = -309.3605;
  double R = 819.1222;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_ls2_11(double r) {
  double I = 0.74129;
  double P = 9.3418;
  double Q = 5.0652;
  double R = -376.4384;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_c_10(double r) {
  double I = -8.62770;
  double P = 2605.2682;
  double Q = 1459.6345;
  double R = 411.9733;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_l2_10(double r) {
  double I = -0.13201;
  double P = 253.4350;
  double Q = 137.4144;
  double R = -1.0076;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_t_10(double r) {
  double I = 1.485601;
  double P = 0.0;
  double Q = -1126.8359;
  double R = 370.1324;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_ls_10(double r) {
  double I = 0.10180;
  double P = 86.0658;
  double Q = 46.6655;
  double R = -356.5175;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_ls2_10(double r) {
  double I = 0.07357;
  double P = -217.5791;
  double Q = -117.9731;
  double R = 18.3935;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_av18_ci_11(double r) {
  double v = 1.0/3.0*(v_av18_c_11_pp(r) + v_av18_c_11_nn(r) + v_av18_c_11_np(r));
  
  return v;
}

double v_av18_ci_01(double r) {
  double v = 1.0/3.0*(v_av18_c_01_pp(r) + v_av18_c_01_nn(r) + v_av18_c_01_np(r));
  
  return v;
}
 
double v_av18_c(double r) {
  double v = 9.0*v_av18_ci_11(r) + 3.0*v_av18_c_10(r) + 3.0*v_av18_ci_01(r) + v_av18_c_00(r);
  v *= 1.0/16.0;
  
  return v;
}

double v_av18_tau(double r) {
  double v = 3.0*v_av18_ci_11(r) - 3.0*v_av18_c_10(r) + v_av18_ci_01(r) - v_av18_c_00(r);
  v *= 1.0/16.0;

  return v;
}

double v_av18_sigma(double r) {
   double v = 3.0*v_av18_ci_11(r) + v_av18_c_10(r) -3.0*v_av18_ci_01(r) - v_av18_c_00(r);
  v *= 1.0/16.0;

  return v;
}

double v_av18_sigma_tau(double r) {
   double v = v_av18_ci_11(r) - v_av18_c_10(r) - v_av18_ci_01(r) + v_av18_c_00(r);
  v *= 1.0/16.0;

  return v;
}

double v_av18_t(double r) {
  double v = 3.0*v_av18_t_11(r) + v_av18_t_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_ls(double r) {
  double v = 3.0*v_av18_ls_11(r) + v_av18_ls_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_ls2(double r) {
  double v = 3.0*v_av18_ls2_11(r) + v_av18_ls2_10(r);
  v *= 1.0/4.0;
  
  return v;
}

double v_av18_t_tau(double r) {
  double v = v_av18_t_11(r) - v_av18_t_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_l2(double r) {
  double v = 9.0*v_av18_l2_11(r) + 3.0*v_av18_l2_10(r) + 3.0*v_av18_l2_01(r) + v_av18_l2_00(r);
  v *= 1.0/16.0;

  return v;
}

double v_av18_l2_tau(double r) {
  double v = 3.0*v_av18_l2_11(r) - 3.0*v_av18_l2_10(r) + v_av18_l2_01(r) - v_av18_l2_00(r);
  v *= 1.0/16.0;

  return v;
}

double v_av18_l2_sigma(double r) {
  double v = 3.0*v_av18_l2_11(r) + v_av18_l2_10(r) - 3.0*v_av18_l2_01(r) - v_av18_l2_00(r);
  v *= 1.0/16.0;

  return v;
}

double v_av18_l2_sigma_tau(double r) {
  double v = v_av18_l2_11(r) - v_av18_l2_10(r) - v_av18_l2_01(r) + v_av18_l2_00(r);
  v *= 1.0/16.0;

  return v;
}

double v_av18_ls_tau(double r) {
  double v = v_av18_ls_11(r) - v_av18_ls_10(r);
  v *= 1.0/4.0;

return v;
}

double v_av18_ls2_tau(double r) {
  double v = v_av18_ls2_11(r) - v_av18_ls2_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_cd_11(double r) {
  double v = 0.5*(v_av18_c_11_pp(r) + v_av18_c_11_nn(r)) - v_av18_c_11_np(r);
  v *= 1.0/6.0;

  return v;
}

double v_av18_cd_01(double r) {
  double v = 0.5*(v_av18_c_01_pp(r) + v_av18_c_01_nn(r)) - v_av18_c_01_np(r);
  v *= 1.0/6.0;

  return v;
}

double v_av18_cap_t(double r) {
  double v = 3.0*v_av18_cd_11(r) + v_av18_cd_01(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_sigma_cap_t(double r) {
  double v = v_av18_cd_11(r) - v_av18_cd_01(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_t_cap_t(double r) {
  double v = 0.5*(v_av18_t_11(r) + v_av18_t_11(r)) - v_av18_t_11(r);
  v *= 1.0/6.0;

  return v;
}

double v_av18_ca_11(double r) {
  double v = v_av18_c_11_pp(r) - v_av18_c_11_nn(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_ca_01(double r) {
  double v = v_av18_c_01_pp(r) - v_av18_c_01_nn(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_tau_z(double r) {
  double v = 3.0*v_av18_ca_11(r) + v_av18_ca_01(r);
  v *= 1.0/4.0;

  return v;
}

double v_av18_sigma_tau_z(double r) {
  double v = v_av18_ca_11(r) - v_av18_ca_01(r);
  v *= 1.0/4.0;

  return v;
}

double compute_matrix_element_c(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear matrix element of the operator unity with arbitrary radial function specified by iv
  int l1, l2, l1p, l2p;
  if (j12 != j12p) {return 0.0;}
  if (t12 != t12p) {return 0.0;}
  l1 = j1 + 0.5;
  l2 = j2 + 0.5;
  l1p = j1p + 0.5;
  l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

   
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double mat = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1.0);
      double m1 = sqrt(2*t12 + 1);
      m1 *= fact;
      // Now perform Brody-Moshinsky transformation
      m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
      mat += m1;
    }
  }
                        
  return mat;
}

double compute_matrix_element_tau(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear matrix element of the operator tau_i dot tau_j with arbitrary radial function specified by iv
  int l1, l2, l1p, l2p;
  if (j12 != j12p) {return 0.0;}
  if (t12 != t12p) {return 0.0;}
  l1 = j1 + 0.5;
  l2 = j2 + 0.5;
  l1p = j1p + 0.5;
  l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

   
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double mat = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1.0);
      double m1 = 0.5*(t12*(t12 + 1.0) - 1.5);
      m1 *= fact;
      m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
      mat += m1;
    }
  }
                        
  return mat;
}

double compute_matrix_element_sigma(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear matrix element of the operator sigma_i dot sigma_j with arbitrary radial function specified by iv
  int l1, l2, l1p, l2p;
  if (j12 != j12p) {return 0.0;}
  if (t12 != t12p) {return 0.0;}
  l1 = j1 + 0.5;
  l2 = j2 + 0.5;
  l1p = j1p + 0.5;
  l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

   
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double mat = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*t12 + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1.0);
      double m1 = 0.5*(s*(s + 1) - 1.5);
      m1 *= fact;
      m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
      mat += m1;
    }
  }
                        
  return mat;
}

double compute_matrix_element_sigma_tau(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear matrix element of the operator unity with arbitrary radial function specified by iv
  int l1, l2, l1p, l2p;
  if (j12 != j12p) {return 0.0;}
  if (t12 != t12p) {return 0.0;}
  l1 = j1 + 0.5;
  l2 = j2 + 0.5;
  l1p = j1p + 0.5;
  l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

   
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double mat = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1.0);
      double m1 = 0.5*(s*(s + 1) - 1.5)*0.5*(t12*(t12 + 1) - 1.5);
      m1 *= fact;
      m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
      mat += m1;
    }
  }
                        
  return mat;
}

double compute_matrix_element_tensor(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  
  double mat = 0.0;
  if ((j12 != j12p ) || (t12 != t12p)) {return 0.0;}
    
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, abs(lambda - 2));
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      double fact = sqrt((2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*lambda + 1)*sqrt(2*lambdap + 1);
      fact *= pow(-1.0, 1 + j12)*sqrt(15)*(-2.0)*3.0;
      fact *= six_j(j12, 1, lambdap, 2, lambda, 1);
      fact *= sqrt(2*t12 + 1)*sqrt(2*j12 + 1);
      fact *= compute_radial_matrix_element_y2(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, 1, t12);
      mat += fact;
    }
  }
                        
  return mat;
}

double compute_matrix_element_tensor_tau(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  
  double mat = 0.0;
  if ((j12 != j12p) || (t12 != t12p)) {return 0.0;}
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, abs(lambda - 2));
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      double fact = sqrt((2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= (2*lambda + 1)*(2*lambdap + 1);
      fact *= pow(-1.0, j12)*sqrt(30)*(-2.0)*3.0;
      fact *= six_j(j12, 1, lambdap, 2, lambda, 1);
      fact *= 0.5*(t12*(t12 + 1) - 1.5)*sqrt(2*j12 + 1);
      fact *= compute_radial_matrix_element_y2(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, 1, t12);
      mat += fact;
    }
  }
                        
  return mat;
}


double compute_matrix_element_ls(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear reduced matrix element L dot S
  
  double mat = 0.0;
  if ((j12 != j12p) || (t12 != t12p)) {return mat;}
    
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, abs(lambda - 2));
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      for (int s = 0; s <= 1; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambdap + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        fact *= pow(-1.0, s + j12 + lambda);
        fact *= six_j(j12, s, lambdap, 1, s, lambda);
        fact *= sqrt(s*(s + 1)*(2*s + 1));
        fact *= sqrt(2*j12 + 1);
        fact *= sqrt(2*t12 + 1);
        fact *= compute_radial_matrix_element_spin_orbit(iv, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, s, t12);
        mat += fact;
      }
    }
  }
                        
  return mat;
}

double compute_matrix_element_ls_tau(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear reduced matrix element L dot S
  
  double mat = 0.0;
  if ((j12 != j12p) || (t12 != t12p)) {return mat;}
    
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, abs(lambda - 2));
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
      for (int s = 0; s <= 1; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambdap + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        fact *= pow(-1.0, s + j12 + lambda);
        fact *= six_j(j12, s, lambdap, 1, s, lambdap);
        fact *= sqrt(s*(s + 1)*(2*s + 1));
        fact *= sqrt(2*j12 + 1);
        fact *= 0.5*(t12*(t12 + 1) - 1.5);
        fact *= compute_radial_matrix_element_spin_orbit(iv, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, s, t12);
        mat += fact;
      }
    }
  }
                        
  return mat;
}

double compute_matrix_element_l2(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear reduced matrix element L^2
  
  double mat = 0.0;
  if ((j12 != j12p) || (t12 != t12p)) {return mat;}
    
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    for (int s = 0; s <= 1; s++) {
    double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1)*sqrt(2*t12 + 1);
      fact *= compute_radial_matrix_element_l2(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
      mat += fact;
    }
  }
                        
  return mat;
}

double compute_matrix_element_l2_tau(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear reduced matrix element L^2
  
  double mat = 0.0;
  if ((j12 != j12p) || (t12 != t12p)) {return mat;}
    
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    for (int s = 0; s <= 1; s++) {
    double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1);
      fact *= 0.5*(t12*(t12 + 1) - 1.5);
      fact *= compute_radial_matrix_element_l2(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
      mat += fact;
    }
  }
                        
  return mat;
}

double compute_matrix_element_l2_sigma(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear reduced matrix element L^2
  
  double mat = 0.0;
  if ((j12 != j12p) || (t12 != t12p)) {return mat;}
    
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    for (int s = 0; s <= 1; s++) {
    double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
    fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
    fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
    fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
    if (fact == 0.0) {continue;}
    fact *= sqrt(2*j12 + 1);
    fact *= sqrt(2*t12 + 1)*2.0*(s*(s + 1) - 1.5);
    fact *= compute_radial_matrix_element_l2(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
    mat += fact;
    }
  }
                        
  return mat;
}

double compute_matrix_element_l2_sigma_tau(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear reduced matrix element L^2
  
  double mat = 0.0;
  if ((j12 != j12p) || (t12 != t12p)) {return mat;}
    
  int l1 = j1 + 0.5;
  int l2 = j2 + 0.5;
  int l1p = j1p + 0.5;
  int l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double m4 = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    for (int s = 0; s <= 1; s++) {
    double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
    fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
    fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
    fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
    if (fact == 0.0) {continue;}
    fact *= sqrt(2*j12 + 1);
    fact *= 0.5*(t12*(t12 + 1) - 1.5)*2.0*(s*(s + 1) - 1.5);
    fact *= compute_radial_matrix_element_l2(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
    mat += fact;
    }
  }
                        
  return mat;
}


double compute_matrix_element_ls2(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear matrix element of the operator unity with arbitrary radial function specified by iv

  int lambda, s;
  int l1, l2, l1p, l2p;
  if (j12 != j12p) {return 0.0;}
  if (t12 != t12p) {return 0.0;}
  l1 = j1 + 0.5;
  l2 = j2 + 0.5;
  l1p = j1p + 0.5;
  l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

   
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double mat = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
  int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, abs(lambda - 2));
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
 
      double fact = sqrt((2*lambda + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambdap + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      for (int j = 0 ; j <= 2; j++) {
        fact *= 3.0*sqrt(2*t12 + 1);
        fact *= sqrt((2*j12 + 1.0)*(2*j + 1)*(2*lambda + 1)*(2*lambdap + 1));
        fact *= pow(-1.0, j12 + 1);
        fact *= six_j(1, 1, j, 1, 1, 1);
        fact *= six_j(j12, 1, lambdap, j, lambda, 1);
        // Now perform Brody-Moshinsky transformation
        int n_rel, l_rel, n_cm, l_cm;
        int max = 2*n1 + 2*n2 + l1 + l2;
        int maxp = 2*n1p + 2*n2p + l1p + l2p;
        double radial_mat = 0.0;
        for (l_cm = 0; l_cm <= max; l_cm++) {
          for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
            if (pow(-1.0, l_cm + l_rel) != pow(-1.0, l1 + l2)) {continue;}
            double sym = 1.0 + pow(-1.0, s + l_rel);
            if (sym == 0.0) {continue;}
            for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
              n_rel = (max - l_rel - l_cm)/2 - n_cm;
              int n_relp = n_rel + (maxp - max)/2;
              if (n_relp < 0) {continue;}
              double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
              rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambda, n1p, l1p, n2p, l2p);
              rm *= sym;
              rm *= pow(-1.0, l_rel + l_cm)*(l_rel*(l_rel + 1)*(2*l_rel + 1));
              rm *= six_j(1, 1, j, l_rel, l_rel, l_rel);
              rm *= six_j(l_rel, lambdap, l_cm, lambda, l_rel, j);
              if (rm == 0.0) {continue;}
              rm *= compute_potential(n_rel, n_relp, l_rel, l_rel, iv);
              radial_mat += rm;
            }
          }
        }
        fact *= radial_mat;
        mat += fact;
      }
    }
  }
                        
  return mat;
}

double compute_matrix_element_ls2_tau(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear matrix element of the operator unity with arbitrary radial function specified by iv

  int lambda, s;
  int l1, l2, l1p, l2p;
  if (j12 != j12p) {return 0.0;}
  if (t12 != t12p) {return 0.0;}
  l1 = j1 + 0.5;
  l2 = j2 + 0.5;
  l1p = j1p + 0.5;
  l2p = j2p + 0.5;
  
  if (pow(-1.0, j1 + 0.5 - in1) == -1.0) {l1 -= 1;}
  if (pow(-1.0, j2 + 0.5 - in2) == -1.0) {l2 -= 1;}
  if (pow(-1.0, j1p + 0.5 - in1p) == -1.0) {l1p -= 1;}
  if (pow(-1.0, j2p + 0.5 - in2p) == -1.0) {l2p -= 1;}

   
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
 
  double mat = 0.0;
  // Convert from JJ to LS coupling (L is lambda)
 int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
  int lambda_max = MIN(l1 + l2, j12 + 1);
  int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
  int lambdap_max = MIN(l1p + l2p, j12p + 1);
  for (int lambda = lambda_min; lambda <= lambda_max; lambda++) {
    lambdap_min = MAX(lambdap_min, abs(lambda - 2));
    lambdap_max = MIN(lambdap_max, lambda + 2);
    if (lambdap_min > lambdap_max) {continue;}
    for (int lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
 
      double fact = sqrt((2*lambda + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambdap + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      for (int j = 0 ; j <= 2; j++) {
        fact *= 3.0*0.5*(t12*(t12 + 1) - 1.5);
        fact *= sqrt((2*j12 + 1.0)*(2*j + 1)*(2*lambda + 1)*(2*lambdap + 1));
        fact *= pow(-1.0, j12 + 1);
        fact *= six_j(1, 1, j, 1, 1, 1);
        fact *= six_j(j12, 1, lambdap, j, lambda, 1);
        // Now perform Brody-Moshinsky transformation
        int n_rel, l_rel, n_cm, l_cm;
        int max = 2*n1 + 2*n2 + l1 + l2;
        int maxp = 2*n1p + 2*n2p + l1p + l2p;
        double radial_mat = 0.0;
        for (l_cm = 0; l_cm <= max; l_cm++) {
          for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
            if (pow(-1.0, l_cm + l_rel) != pow(-1.0, l1 + l2)) {continue;}
            double sym = 1.0 + pow(-1.0, s + l_rel);
            if (sym == 0.0) {continue;}
            for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
              n_rel = (max - l_rel - l_cm)/2 - n_cm;
              int n_relp = n_rel + (maxp - max)/2;
              if (n_relp < 0) {continue;}
              double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
              rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambda, n1p, l1p, n2p, l2p);
              rm *= sym;
              rm *= pow(-1.0, l_rel + l_cm)*(l_rel*(l_rel + 1)*(2*l_rel + 1));
              rm *= six_j(1, 1, j, l_rel, l_rel, l_rel);
              rm *= six_j(l_rel, lambdap, l_cm, lambda, l_rel, j);
              if (rm == 0.0) {continue;}
              rm *= compute_potential(n_rel, n_relp, l_rel, l_rel, iv);
              radial_mat += rm;
            }
          }
        }
        fact *= radial_mat;
        mat += fact;
      }
    }
  }
                        
  return mat;
}

double compute_matrix_element_av18(int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  double v = 0.0;
  double v1 = compute_matrix_element_c(5, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v2 = compute_matrix_element_tau(6, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v3 = compute_matrix_element_sigma(7, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v4 = compute_matrix_element_sigma_tau(8, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v5 = compute_matrix_element_tensor(9, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v6 = compute_matrix_element_tensor_tau(10, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v7 = compute_matrix_element_ls(11, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v8 = compute_matrix_element_ls_tau(12, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v9 = compute_matrix_element_l2(13, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v10 = compute_matrix_element_l2_tau(14, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v11 = compute_matrix_element_l2_sigma(15, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v12 = compute_matrix_element_l2_sigma_tau(16, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v13 = compute_matrix_element_ls2(17, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v14 = compute_matrix_element_ls2_tau(18, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);


printf("%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14);
  v = v1 + v2 + v3 + v4;
  v += v5 + v6 + v7 + v8;
  v += v9 + v10 + v11 + v12;  
  v += v13 + v14;
  return v;
}

