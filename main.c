#include "nuc_mat.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main() {
  
//  printf("%g\n", (compute_matrix_element_TT(4)));
  printf("%g\n", compute_matrix_element_M_0());

  return 0;
}

double compute_matrix_element_M_0() {
  double mat1 = 2.0*R_NUC*pow(A_NUC, 1.0/3.0)*compute_matrix_element_M(2, 2);
  double mat2 = 2.0*R_NUC*pow(A_NUC, 1.0/3.0)*compute_matrix_element_TT(3);

  printf("S: %g T: %g\n", mat1, mat2);
  double mat = mat1 + mat2;

  return mat;
}

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
  }
  if (iv == 2) {
    v *= v_pion_f1(q);
  }
  if (iv == 3) {
    v *= v_pion_f2(q);
  }
  if (iv == 4) {
    v *= v_light_limit_d(q);
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
  

double compute_matrix_element_TT(int iv) {
  // Computes the total nuclear matrix element for the given operator
  // m_sw = 0 computes M2
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(DENSITY_FILE, "r");
  
  double mat = 0.0;
  int i;
  for (i = 0; i < NUM_SHELLS; i++) {
    // Each line of the file corresponds to a nuclear shell
    float density;
    fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density);

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
   
    int lambda, lambdap;
    int l1, l2, l1p, l2p;
    
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
 
    double m4 = 0.0;
    // Convert from JJ to LS coupling (L is lambda)
    int lambda_min = MAX(abs(l1 - l2), abs(j12 - 1));
    int lambda_max = MIN(l1 + l2, j12 + 1);
    int lambdap_min = MAX(abs(l1p - l2p), abs(j12p - 1));
    int lambdap_max = MIN(l1p + l2p, j12p + 1);
    for (lambda = lambda_min; lambda <= lambda_max; lambda++) {
      lambdap_min = MAX(lambdap_min, abs(lambda - 2));
      lambdap_max = MIN(lambdap_max, lambda + 2);
      if (lambdap_min > lambdap_max) {continue;}
      for (lambdap = lambdap_min; lambdap <= lambdap_max; lambdap++) {
        double fact = sqrt((2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, 1, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, 1, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        fact *= (2*lambda + 1)*(2*lambdap + 1);
        fact *= pow(-1.0, j12)*sqrt(30)*(-2.0)*3.0*3.0*3.0*3.0*3.0*3.0*3.0*3.0;
        fact *= six_j(j12, 1, lambdap, 2, lambda, 1);
        fact *= sqrt(5)*sqrt(2*j12 + 1);
        if ((in1 == in2) && (j1 == j2)) {fact *= 1/sqrt(2);}
        if ((in1p == in2p) && (j1p == j2p)) {fact *= 1/sqrt(2);}
        // Now perform Brody-Moshinsky transformation
        int n_rel, l_rel, n_cm, l_cm, l_relp;
        int max = 2*n1 + 2*n2 + l1 + l2;
        int maxp = 2*n1p + 2*n2p + l1p + l2p;
        double radial_mat = 0.0;
        for (l_cm = 0; l_cm <= max; l_cm++) {
          for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
            double sym = 1.0 + pow(-1.0, 1 + l_rel);
            if (sym == 0.0) {continue;}
            if (pow(-1.0, l_rel + l_cm) != pow(-1.0, l1 + l2)) {continue;}
            for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
              n_rel = (max - l_rel - l_cm)/2 - n_cm;
              for (l_relp = l_rel - 2; l_relp <= l_rel + 2; l_relp += 2) {
                if (l_relp < 0) {continue;}
                if (pow(-1.0, l_relp + l_cm) != pow(-1.0, l1p + l2p)) {continue;}
                int n_relp = n_rel + (maxp - max)/2 + (l_rel - l_relp)/2;
                if (n_relp < 0) {continue;}
                double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
                rm *= brody_mosh(n_relp, l_relp, n_cm, l_cm, lambdap, n1p, l1p, n2p, l2p);
                rm *= sym;
                rm *= sqrt(2*l_rel + 1)*sqrt(2*l_relp + 1)*three_j(l_relp, 2, l_rel, 0, 0, 0) *six_j(l_rel, l_relp, 2, lambdap, lambda, l_cm)*pow(-1.0, l_cm);
                if (rm == 0.0) {continue;}
                rm *= compute_potential(n_rel, n_relp, l_rel, l_relp, iv);
                radial_mat += rm;
              }
            }
          }
        }
        fact *= radial_mat;
        m4 += fact;
      }
    }
    mat += m4*density;
  }
  fclose(in_file);         
                        
  return mat;
}


double compute_matrix_element_M(int m_sw, int iv) {
  // Computes the total nuclear matrix element for the given operator
  // m_sw = 0 computes M2
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(DENSITY_FILE, "r");
  
  double mat = 0.0;
  int i;
  for (i = 0; i < NUM_SHELLS; i++) {
    // Each line of the file corresponds to a nuclear shell
    float density;
    fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density);

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
   
    int lambda, s;
    int l1, l2, l1p, l2p;
    
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
 
    double m4 = 0.0;
    // Convert from JJ to LS coupling (L is lambda)
    for (lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
      int s_max = MIN(lambda + j12, 1);
      int s_min = abs(lambda - j12);
      for (s = s_min; s <= s_max; s++) {
        double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        fact *= sqrt(2*j12 + 1.0);
        double m1 = 0.0;
        if (m_sw == 2) {
          m1 += pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
        }
        if (m_sw == 1) {
          m1 += 1.0;
        }
        m1 *= sqrt(5);
        m1 *= fact;
        if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m1 *= 1/sqrt(2);}
        if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m1 *= 1/sqrt(2);}
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
              if (rm == 0.0) {continue;}
              rm *= compute_potential(n_rel, n_relp, l_rel, l_rel, iv);
              radial_mat += rm;
            }
          }
        }
        m1 *= radial_mat;
        m4 += m1;
      }
    }
    mat += m4*density;
  }
  fclose(in_file);         
                        
  return mat;
}

double compute_matrix_element_M_GT() {
  // Computes the total nuclear matrix element for the given operator
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(DENSITY_FILE, "r");
  
  double mat = 0.0;
  int i;
  for (i = 0; i < NUM_SHELLS; i++) {
    // Each line of the file corresponds to a nuclear shell
    float density;
    fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density);

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
   
    int lambda, s;
    int l1, l2, l1p, l2p;
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
 
    double m4 = 0.0;
    // Convert from JJ to LS coupling (L is lambda)
    for (lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
      int s_max = MIN(lambda + j12, 1);
      int s_min = abs(lambda - j12);
      for (s = s_min; s <= s_max; s++) {
        double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        fact *= sqrt(2*j12 + 1.0);
        double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
        m1 *= fact;
        m1 *= sqrt(5);
        double anti_symm = 0.0;
        if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
        if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1p + l2p + lambda + s + 1);}
        m1 *= anti_symm;
        if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m1 *= 1/sqrt(2);}
        if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m1 *= 1/sqrt(2);}
        m4 += m1;
      }
    }
    mat += m4*density;
  }
  fclose(in_file);     
                        
  return mat;
}

double compute_matrix_element_M_F() {
  // Computes the total nuclear matrix element for the given operator
  // Uses the density matrix method to decompose into two-body matrix elements
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;

  // Open the file containing density matrix coefficients
  FILE *in_file;
  in_file = fopen(DENSITY_FILE, "r");
  
  double mat = 0.0;
  int i;
  for (i = 0; i < NUM_SHELLS; i++) {
    // Each line of the file corresponds to a nuclear shell
    float density;
    fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density);
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    double m4 = 0.0;
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    m4 *= sqrt(2.0*j12p + 1.0);
    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    mat += m4*density*sqrt(5.0);
  }
  fclose(in_file);
                        
  return mat;
}


double mat_gt(double s, double sp) {
  double mat = 0.0;
  if (s != sp) {return mat;}
  mat = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
  mat *= sqrt(5);
  return mat;
}
   
double mat_f(double s, double sp) {
  double mat = 0.0;
  if (s != sp) {return mat;}
  mat = 1.0;
  mat *= sqrt(5);
  return mat;
}

double triangle(double a, double b, double c) {
  // Computes the triangle coefficients necessary for sixj calculations
  // See Edmonds pg. 99
  double tri = 0.0;
  tri = gsl_sf_gamma(a + b - c + 1.0)*gsl_sf_gamma(a - b + c + 1.0)*gsl_sf_gamma(-a + b + c + 1.0);
  tri *= 1.0/gsl_sf_gamma(a + b + c + 2.0);
  tri = sqrt(tri);

  return tri;
}

double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m) {
  // Computes the Clebsch-Gordan coefficients between the uncoupled basis (j1, m1, j2, m2) and
  // the coupled basis (j1,j2; j,m)
  double cg = 0.0;
  if (((m1 + m2) != m) || (j > (j1 + j2)) || (j < abs(j1-j2))) {return cg;}
  int s;
  double f = 0.0;
  for (s = 0; s < 100; s++) {
    double d1 = j1 - m1 -s;
    double d2 = j - m - s;
    double d3 = j2 - j + m1 + s;
    if ((d1 < 0) || (d2 < 0) || (d3 < 0)){continue;}
    f += pow(-1.0, s + j1 - m1)*gsl_sf_gamma(j1 + m1 + s +1.0)*gsl_sf_gamma(j2 + j - m1 - s + 1.0)/(gsl_sf_gamma(s + 1.0)*gsl_sf_gamma(d1 + 1.0)*gsl_sf_gamma(d2 + 1.0)*gsl_sf_gamma(d3 + 1.0));
  }
  double w = (2.0*j + 1.0)*gsl_sf_gamma(j1 + j2 - j + 1.0)*gsl_sf_gamma(j1 - m1 + 1.0)*gsl_sf_gamma(j2 - m2 + 1.0)*gsl_sf_gamma(j + m + 1.0)*gsl_sf_gamma(j - m + 1.0)/(gsl_sf_gamma(j1 + j2 + j + 2.0)*gsl_sf_gamma(j1 - j2 + j + 1.0)*gsl_sf_gamma(-j1 + j2 + j + 1.0)*gsl_sf_gamma(j1 + m1 + 1.0)*gsl_sf_gamma(j2 + m2 + 1.0));
  cg = sqrt(w)*f;
  
  return cg;
} 

double three_j(double j1, double j2, double j3, double m1, double m2, double m3) {
  // Computes the Wigner 3J symbol from the corresponding Clebsch-Gordan coefficient
  double three_j = pow(-1.0, j1 - j2 - m3)/sqrt(2.0*j3 + 1.0)*clebsch_gordan(j1, j2, j3, m1, m2, -m3);

  return three_j;
}  

double six_j(double j1, double j2, double j3, double j4, double j5, double j6) {
  // Computes the Wigner six-j symbol using the Racah formual (Edmonds pg. 99)
  double six_j = 0;
  if ((j1 < 0) || (j2 < 0) || (j3 < 0) || (j4 < 0) || (j5 < 0) || (j6 < 0)) {
    printf("Unallowed quantum numbers: %g, %g, %g, %g, %g, %g\n", j1, j2, j3, j4, j5, j6);
    return six_j;
  }
  if ((j1 < abs(j2 - j3)) || (j1 > (j2 + j3)) || (j1 < abs(j5 - j6)) || (j1 > (j5 + j6)) || (j4 < abs(j2 - j6)) || (j4 > (j2 + j6)) || (j4 < abs(j5 - j3)) || (j4 > (j5 + j3))) {return six_j;}
  double t1 = triangle(j1, j2, j3);
  double t2 = triangle(j1, j5, j6);
  double t3 = triangle(j4, j2, j6);
  double t4 = triangle(j4, j5, j3);

  double w = 0.0;
  
  int z;
  for (z = 0; z < 100; z++) {
    double d1 = z - j1 - j2 - j3;
    double d2 = z - j1 - j5 - j6;
    double d3 = z - j4 - j2 - j6;
    double d4 = z - j4 - j5 - j3;
    double d5 = j1 + j2 + j4 + j5 - z;
    double d6 = j2 + j3 + j5 + j6 - z;
    double d7 = j3 + j1 + j6 + j4 - z;
    if ((d1 < 0) || (d2 < 0) || (d3 < 0) || (d4 < 0) || (d5 <0) || (d6 < 0) || (d7 < 0)){continue;}
    double f = gsl_sf_gamma(d1 + 1.0)*gsl_sf_gamma(d2 + 1.0)*gsl_sf_gamma(d3 + 1.0)*gsl_sf_gamma(d4 + 1.0)*gsl_sf_gamma(d5 + 1.0)*gsl_sf_gamma(d6 + 1.0)*gsl_sf_gamma(d7 + 1.0);
    w += pow(-1.0, z)*gsl_sf_gamma(z + 2.0)/f;
  }
  six_j = t1*t2*t3*t4*w;
 
  return six_j;
}

double nine_j(double j11, double j12, double j13, double j21, double j22, double j23, double j31, double j32, double j33) {
  // Computes the Wigner 9J-symbol from the necessary 6J-symbols
  double nine_j = 0.0;
  int i;
 
  for (i = 0; i < 100; i++) {
    double k = i/2.0;
    nine_j += pow(-1.0, i)*(i + 1.0)*six_j(j11, j21, j31, j32, j33, k)*six_j(j12, j22, j32, j21, k, j23)*six_j(j13, j23, j33, k, j11, j12);
  }
  
  return nine_j;
}

double g_factor(double a, double b, double c) {
  // Computes the G(l',l'',l) coefficient required to compute the BM coefficients
  // See Moshinsky (1959) Eq. 48
  double g = pow(-1.0, b)*sqrt(4.0*M_PI)*sqrt(gsl_sf_gamma(2.0*c + 2.0)/(gsl_sf_gamma(2.0*a + 2.0)*gsl_sf_gamma(2.0*b + 2.0)));
  return g;
}

double h_factor(double a, double b, double c) {
  // Computes the H(l',l'',l) coefficient required to compute the BM coefficients
  // See Moshinsky (1959) Eq. 54
  double h = sqrt(((2.0*a + 1.0)*(2.0*b + 1.0))/(4.0*M_PI*(2.0*c + 1.0)))*clebsch_gordan(a,b,c,0,0,0);
  return h;
}

double brody_mosh_zero(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int l1, int l2) {
  // Computes the n1 = n2 = 0 Brody-Moshinksky brackets
  // See Moshinsky (1959)
  double bm = (2.0*l_rel + 1.0)*(2.0*l_cm + 1.0)*(2.0*l1 + 1.0)*(2.0*l2 + 1.0)/(gsl_sf_gamma(l1 + 1.5)*gsl_sf_gamma(l2 + 1.5)*gsl_sf_gamma(n_rel + 1.0)*gsl_sf_gamma(n_rel + l_rel + 1.5)*gsl_sf_gamma(n_cm + 1.0)*gsl_sf_gamma(n_cm + l_cm + 1.5));
  bm = sqrt(bm);
  bm *= pow(2.0, -0.5*(l1 + l2));
  double f = 0.0;
  int l_rel_1, l_rel_2, l_cm_1, l_cm_2;
  for (l_rel_1 = 0; l_rel_1 <= l1; l_rel_1++) {
    for (l_rel_2 = 0; l_rel_2 <= l2; l_rel_2++) {
      if ((l_rel_1 + l_rel_2) != 2*n_rel + l_rel) {continue;}
      for (l_cm_1 = 0; l_cm_1 <= l1; l_cm_1++) {
        if ((l_rel_1 + l_cm_1) != l1) {continue;}
        for (l_cm_2 = 0; l_cm_2 <= l2; l_cm_2++) {
          if ((l_rel_2 + l_cm_2) != l2) {continue;}
          if ((l_cm_1 + l_cm_2) != 2*n_cm + l_cm) {continue;}
            f += pow(-1.0, l_rel_1)*g_factor(l_cm_1, l_rel_1, l1)*g_factor(l_cm_2, l_rel_2, l2)*h_factor(l_rel_1, l_rel_2, l_rel)*h_factor(l_cm_1, l_cm_2, l_cm)*nine_j(l_rel_1, l_rel_2, l_rel, l_cm_1, l_cm_2, l_cm, l1, l2, l_tot)*pow(-1.0, n_rel)*gsl_sf_gamma(0.5*(l_rel_1 + l_rel_2 - l_rel) + 1.0)*gsl_sf_gamma(0.5*(l_rel_1 + l_rel_2 + l_rel + 3.0))*pow(-1.0, n_cm)*gsl_sf_gamma(0.5*(l_cm_1 + l_cm_2 - l_cm) + 1.0)*gsl_sf_gamma(0.5*(l_cm_1 + l_cm_2 + l_cm + 3.0));
        }
      }
    }
  }
  bm *= f;

  return bm;
} 

double brody_mosh(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int n1, int l1, int n2, int l2) {
  // Uses the r^2 recursion relation of Moshinsky to compute the Brody-Moshinsky coefficients for  // arbitrary n1, n2
  // See Moshinsky (1959)
  double bm = 0.0;
  if ((n1 == 0) && (n2 == 0)) {
    bm = brody_mosh_zero(n_rel, l_rel, n_cm, l_cm, l_tot, l1, l2);
    return bm;
  } else if (n2 == 0) {
    if (n_rel > 0) {
      bm += 0.5*sqrt(n_rel*(n_rel + l_rel + 0.5))*brody_mosh(n_rel - 1, l_rel, n_cm, l_cm, l_tot, n1 - 1, l1, n2, l2);
      if (n_cm > 0) {
        bm += sqrt(n_rel*n_cm*(l_rel + 1.0)*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm - 1, l_cm + 1, l_tot, n1 - 1, l1, n2, l2);
      }
      if (l_cm > 0) {
        bm += sqrt(n_rel*(n_cm + l_cm + 0.5)*(l_rel + 1.0)*l_cm)*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm, l_cm - 1, l_tot, n1 - 1, l1, n2, l2);
      }
    }
    if (n_cm > 0) {
      bm += 0.5*sqrt(n_cm*(n_cm + l_cm + 0.5))*brody_mosh(n_rel, l_rel, n_cm - 1, l_cm, l_tot, n1 - 1, l1, n2, l2);
      if (l_rel > 0) {
        bm += sqrt((n_rel + l_rel + 0.5)*n_cm*l_rel*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel - 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm - 1, l_cm + 1, l_tot, n1 - 1, l1, n2, l2);
      }
    }
    if ((l_rel > 0) && (l_cm > 0)) {
      bm += sqrt((n_rel + l_rel + 0.5)*(n_cm + l_cm + 0.5)*l_rel*l_cm)*pow(-1.0, l_tot + l_rel + l_cm)*six_j(l_rel, l_rel - 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm, l_cm - 1, l_tot, n1 - 1, l1, n2, l2);
    }
    bm *= 1.0/sqrt((n1)*(n1 + l1 + 0.5));
  } else { 
   if (n_rel > 0) {
      bm += 0.5*sqrt(n_rel*(n_rel + l_rel + 0.5))*brody_mosh(n_rel - 1, l_rel, n_cm, l_cm, l_tot, n1, l1, n2 - 1, l2);
      if (n_cm > 0) {
        bm -= sqrt(n_rel*n_cm*(l_rel + 1.0)*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm - 1, l_cm + 1, l_tot, n1, l1, n2 - 1, l2);
      }
      if (l_cm > 0) {
        bm -= sqrt(n_rel*(n_cm + l_cm + 0.5)*(l_rel + 1.0)*l_cm)*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm, l_cm - 1, l_tot, n1, l1, n2 - 1, l2);
      }
    }
    if (n_cm > 0) {
      bm += 0.5*sqrt(n_cm*(n_cm + l_cm + 0.5))*brody_mosh(n_rel, l_rel, n_cm - 1, l_cm, l_tot, n1, l1, n2 - 1, l2);
      if (l_rel > 0) {
        bm -= sqrt((n_rel + l_rel + 0.5)*n_cm*l_rel*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel - 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm - 1, l_cm + 1, l_tot, n1, l1, n2 - 1, l2);
      }
    }
    if ((l_rel > 0) && (l_cm > 0)) {
      bm -= sqrt((n_rel + l_rel + 0.5)*(n_cm + l_cm + 0.5)*l_rel*l_cm)*pow(-1.0, l_tot + l_rel + l_cm)*six_j(l_rel, l_rel - 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm, l_cm - 1, l_tot, n1, l1, n2 - 1, l2);
    }
    bm *= 1.0/sqrt((n2)*(n2 + l2 + 0.5));
  }
  return bm;
}


double Fermi(int Z, double T) {
  // Computes the electron Fermi function for a given atomic number Z and electron energy T
  double s = sqrt(1 - pow(ALPHA_FS*Z, 2.0));
  double p = sqrt(pow(T, 2.0) - pow(M_ELECTRON, 2));
  double eta = ALPHA_FS*Z*T/p;
  double rho = R_NUC*pow(A_NUC, 1.0/3.0);
  gsl_sf_result a, b;

  double f = 2*(1+s)/pow(gsl_sf_gamma(1+2*s), 2.0);
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
  phase_int *=  p_1*p_2*E_1*E_2;

  return phase_int;
}

double phase_integrand_2(double E_1) {
  // Integrand for the type II phase space integral
  double E_2 = E_BETA - E_1;
  double p_1 = sqrt(pow(E_1, 2.0) - pow(M_ELECTRON, 2.0));
  double p_2 = sqrt(pow(E_2, 2.0) - pow(M_ELECTRON, 2.0));
  double phase_int = Fermi(Z_ATOM + 2, E_1)*Fermi(Z_ATOM + 2, E_2);
  phase_int *=  p_1*p_2*pow(M_ELECTRON, 2.0);

  return phase_int;
}
  

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

