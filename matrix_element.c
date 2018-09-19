#include "matrix_element.h"

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
      int n_rel, l_rel, n_cm, l_cm;
      int max = 2*n1 + 2*n2 + l1 + l2;
      int maxp = 2*n1p + 2*n2p + l1p + l2p;
      double radial_mat = 0.0;
      for (l_cm = 0; l_cm <= max; l_cm++) {
        for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
          if (pow(-1.0, l_cm + l_rel) != pow(-1.0, l1 + l2)) {continue;}
          for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
            n_rel = (max - l_rel - l_cm)/2 - n_cm;
            int n_relp = n_rel + (maxp - max)/2;
            if (n_relp < 0) {continue;}
            double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
            rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambda, n1p, l1p, n2p, l2p);
            if (rm == 0.0) {continue;}
            rm *= compute_potential(n_rel, n_relp, l_rel, l_rel, iv);
            radial_mat += rm;
          }
        }
      }
      m1 *= radial_mat;
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
      // Now perform Brody-Moshinsky transformation
      int n_rel, l_rel, n_cm, l_cm;
      int max = 2*n1 + 2*n2 + l1 + l2;
      int maxp = 2*n1p + 2*n2p + l1p + l2p;
      double radial_mat = 0.0;
      for (l_cm = 0; l_cm <= max; l_cm++) {
        for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
          if (pow(-1.0, l_cm + l_rel) != pow(-1.0, l1 + l2)) {continue;}
          for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
            n_rel = (max - l_rel - l_cm)/2 - n_cm;
            int n_relp = n_rel + (maxp - max)/2;
            if (n_relp < 0) {continue;}
            double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
            rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambda, n1p, l1p, n2p, l2p);
            if (rm == 0.0) {continue;}
            rm *= compute_potential(n_rel, n_relp, l_rel, l_rel, iv);
            radial_mat += rm;
          }
        }
      }
      m1 *= radial_mat;
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
      // Now perform Brody-Moshinsky transformation
      int n_rel, l_rel, n_cm, l_cm;
      int max = 2*n1 + 2*n2 + l1 + l2;
      int maxp = 2*n1p + 2*n2p + l1p + l2p;
      double radial_mat = 0.0;
      for (l_cm = 0; l_cm <= max; l_cm++) {
        for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
          if (pow(-1.0, l_cm + l_rel) != pow(-1.0, l1 + l2)) {continue;}
          for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
            n_rel = (max - l_rel - l_cm)/2 - n_cm;
            int n_relp = n_rel + (maxp - max)/2;
            if (n_relp < 0) {continue;}
            double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
            rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambda, n1p, l1p, n2p, l2p);
            if (rm == 0.0) {continue;}
            rm *= compute_potential(n_rel, n_relp, l_rel, l_rel, iv);
            radial_mat += rm;
          }
        }
      }
      m1 *= radial_mat;
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
      mat += m1;
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
        fact *= (2*lambda + 1)*(2*lambdap + 1);
        fact *= pow(-1.0, s + j12 + 1 + lambda + lambdap);
        fact *= six_j(s, s, 1, lambda, lambdap, j12);
        fact *= sqrt(s*(s + 1)*(2*s + 1));
        fact *= sqrt(2*j12 + 1);
        fact *= sqrt(2*t12 + 1);
        // Now perform Brody-Moshinsky transformation
        int n_rel, l_rel, n_cm, l_cm;
        int max = 2*n1 + 2*n2 + l1 + l2;
        int maxp = 2*n1p + 2*n2p + l1p + l2p;
        double radial_mat = 0.0;
        for (l_cm = 0; l_cm <= max; l_cm++) {
          for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
            if (pow(-1.0, l_rel + l_cm) != pow(-1.0, l1 + l2)) {continue;}
            for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
              n_rel = (max - l_rel - l_cm)/2 - n_cm;
              int n_relp = n_rel + (maxp - max)/2;
              if (n_relp < 0) {continue;}
              double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
              rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambdap, n1p, l1p, n2p, l2p);
              rm *= sqrt(l_rel*(l_rel + 1)*(2*l_rel + 1))*six_j(l_rel, l_rel, 1, lambdap, lambda, l_cm)*pow(-1.0, l_cm - l_rel);
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
        fact *= (2*lambda + 1)*(2*lambdap + 1);
        fact *= pow(-1.0, s + j12 + 1 + lambda + lambdap);
        fact *= six_j(s, s, 1, lambda, lambdap, j12);
        fact *= sqrt(s*(s + 1)*(2*s + 1));
        fact *= sqrt(2*j12 + 1);
        fact *= 0.5*(t12*(t12 + 1) - 1.5);
        // Now perform Brody-Moshinsky transformation
        int n_rel, l_rel, n_cm, l_cm;
        int max = 2*n1 + 2*n2 + l1 + l2;
        int maxp = 2*n1p + 2*n2p + l1p + l2p;
        double radial_mat = 0.0;
        for (l_cm = 0; l_cm <= max; l_cm++) {
          for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
            if (pow(-1.0, l_rel + l_cm) != pow(-1.0, l1 + l2)) {continue;}
            for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
              n_rel = (max - l_rel - l_cm)/2 - n_cm;
              int n_relp = n_rel + (maxp - max)/2;
              if (n_relp < 0) {continue;}
              double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
              rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambdap, n1p, l1p, n2p, l2p);
              rm *= sqrt(l_rel*(l_rel + 1)*(2*l_rel + 1))*six_j(l_rel, l_rel, 1, lambdap, lambda, l_cm)*pow(-1.0, l_cm - l_rel);
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
        fact *= (2*lambda + 1)*(2*lambdap + 1);
        fact *= sqrt(2*j12 + 1)*sqrt(2*j12p + 1);
        fact *= sqrt(2*t12 + 1)*sqrt(2*s + 1);
        fact *= nine_j(lambdap, lambda, 0, s, s, 0, j12p, j12, 0);
        // Now perform Brody-Moshinsky transformation
        int n_rel, l_rel, n_cm, l_cm;
        int max = 2*n1 + 2*n2 + l1 + l2;
        int maxp = 2*n1p + 2*n2p + l1p + l2p;
        double radial_mat = 0.0;
        for (l_cm = 0; l_cm <= max; l_cm++) {
          for (l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
            if (pow(-1.0, l_rel + l_cm) != pow(-1.0, l1 + l2)) {continue;}
            for (n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
              n_rel = (max - l_rel - l_cm)/2 - n_cm;
              int n_relp = n_rel + (maxp - max)/2;
              if (n_relp < 0) {continue;}
              double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
              rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambdap, n1p, l1p, n2p, l2p);
              rm *= l_rel*(l_rel + 1)*six_j(l_rel, l_rel, 1, lambdap, lambda, l_cm)*pow(-1.0, l_cm - l_rel);
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


double compute_matrix_element_(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
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
  for (lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (s = s_min; s <= s_max; s++) {
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      fact *= sqrt(2*j12 + 1.0)/sqrt(2*lambda + 1.0);
      double m1 = 1.0;
      m1 *= fact;
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
      mat += m1;
    }
  }
                        
  return mat;
}


double compute_matrix_element_tau_plus(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
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
 
  double mat = 0.0;
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
      double m1 = sqrt(5);
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
      mat += m1;
    }
  }
                        
  return mat;
}


double compute_matrix_element_sigma_tau_plus(int iv, int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv

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
 
  double mat = 0.0;
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
      mat += m1;
    }
  }
                        
  return mat;
}


double compute_matrix_element_M_GT() {
  // no radial part
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
  // No radial part
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

double compute_matrix_element_av18(int in1p, double j1p, int in2p, double j2p, int j12p, int t12p, int in1, double j1, int in2, double j2, int j12, int t12) {
  double v = 0.0;
  double v1 = compute_matrix_element_c(5, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v2 = compute_matrix_element_tau(6, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v3 = compute_matrix_element_sigma(7, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v4 = compute_matrix_element_sigma_tau(8, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
  double v7 = compute_matrix_element_ls(11, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
 double v8 = compute_matrix_element_ls_tau(12, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);
 double v9 = compute_matrix_element_l2(13, in1p, j1p, in2p, j2p, j12p, t12p, in1, j1, in2, j2, j12, t12);


printf("%g, %g, %g, %g, %g, %g, %g\n", v1, v2, v3, v4, v7, v8, v9);
  v = v1 + v2 + v3 + v4;
  v += v7 + v8 + v9;  
  return v;
}
