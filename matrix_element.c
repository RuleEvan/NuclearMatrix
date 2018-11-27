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
  if (iv == -1) {potential_spline_init(h_AA_T_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -2) {potential_spline_init(h_AA_T_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -3) {potential_spline_init(h_AP_T_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -4) {potential_spline_init(h_AP_T_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);}
  if (iv == -5) {potential_spline_init(h_PP_T_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);}
  if (iv == -6) {potential_spline_init(h_PP_T_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -7) {potential_spline_init(h_MM_T_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -8) {potential_spline_init(h_MM_T_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
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
        fact *= sqrt(2*lambda + 1)*sqrt(2*lambdap + 1);
        fact *= pow(-1.0, j12 + lambda + 1)*sqrt(15)*sqrt(32*M_PI/5)*3.0;
        fact *= six_j(j12, 1, lambdap, 2, lambda, 1);
        fact *= sqrt(5)*sqrt(2*j12 + 1);
        if ((in1 == in2) && (j1 == j2)) {fact *= 1/sqrt(2);}
        if ((in1p == in2p) && (j1p == j2p)) {fact *= 1/sqrt(2);}  
        fact *= compute_radial_matrix_element_y2(iv, n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, 1, t12);
        m4 += fact;
      }
    }
    mat += m4*density;
  }
  fclose(in_file);         
                        
  return mat;
}


double compute_matrix_element_tau_plus(int iv) {
  double mat = 0.0;
  FILE *in_file;
  in_file = fopen(DENSITY_FILE, "r");
  if (iv == -1) {potential_spline_init(h_F_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -2) {potential_spline_init(h_F_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 

  for (int i = 0; i < NUM_SHELLS; i++) {
    int in1, in2, ij1, ij2, ij12, it12;
    int in1p, in2p, ij1p, ij2p, ij12p, it12p;

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
        double m1 = sqrt(5);
        m1 *= fact;
        if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m1 *= 1/sqrt(2);}
        if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m1 *= 1/sqrt(2);}
        m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
        m4 += m1;
      }
    }
    mat += density*m4;
  }
                        
  return mat;
}


double compute_matrix_element_sigma_tau_plus(int iv) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(DENSITY_FILE, "r");
  double mat = 0.0;
  if (iv == -1) {potential_spline_init(h_AA_GT_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -2) {potential_spline_init(h_AA_GT_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -3) {potential_spline_init(h_AP_GT_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -4) {potential_spline_init(h_AP_GT_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);}
  if (iv == -5) {potential_spline_init(h_PP_GT_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);}
  if (iv == -6) {potential_spline_init(h_PP_GT_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -7) {potential_spline_init(h_MM_GT_q, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
  if (iv == -8) {potential_spline_init(h_MM_GT_q_sd, R_SPL_MIN, R_SPL_MAX, R_SPL_STEPS);} 
 
  for (int i = 0; i < NUM_SHELLS; i++) {
    int in1, in2, ij1, ij2, ij12, it12;
    int in1p, in2p, ij1p, ij2p, ij12p, it12p;

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
        double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
        m1 *= sqrt(5);
        m1 *= fact;
        if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m1 *= 1/sqrt(2);}
        if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m1 *= 1/sqrt(2);}
        m1 *= compute_radial_matrix_element_scalar(iv, n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12);
        m4 += m1;
      }
    }
    mat += m4*density;
  }
                        
  return mat;
}


double compute_matrix_element_M_GT() {
  // no radial part
  // Computes the total nuclear matrix element for the given operator
  // Uses the density matrix method to decompose into two-body matrix elements
  int jop = 0;
  int top = 2;
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
    //    if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m1 *= 1/sqrt(2);}
    //    if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m1 *= 1/sqrt(2);}
        m4 += m1;
      }
    }
    mat += 0.5*m4*density;
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
    if (t12 != 1 || t12p != 1) {continue;}
    if (j12 != j12p) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    m4 *= sqrt(2.0*j12p + 1.0);
    
//    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
//    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    mat += m4*density*sqrt(5.0)/2.0;
  }
  fclose(in_file);
                        
  return mat;
}

double compute_matrix_element_M_JF() {
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
    if (t12 != t12p) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0,  j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    m4 *= sqrt(j1*(2*j1 + 1)*(j1 + 1))*sqrt(j2*(2*j2 + 1)*(j2 + 1))*nine_j(j1p, j1, 1, j2p, j2, 1, j12p, j12, 2);
    m4 *= sqrt(2.0*j12 + 1.0)*sqrt(2*j12p + 1.0)*sqrt(5.0);
    
//    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
//    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    mat += m4*density*sqrt(2*t12 + 1)/2.0;
  }
  fclose(in_file);
                        
  return mat;
}

double compute_matrix_element_M_JF2() {
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
    if (t12 != t12p) {continue;}
    if (j12 != j12p) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (j2 == j2p) && (in2 == in2p)) {
      m4 = 1.0; 
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (j2 == j1p) && (in2 == in1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    m4 *= sqrt(j12*(j12 + 1)*(2*j12 + 3)*(2*j12 - 1));
    
//    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
//    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    mat += m4*density*sqrt(2*t12 + 1)/2.0*sqrt(2)*2.0/3.0;
    printf("Adding %g\n", mat);
  }
  fclose(in_file);
                        
  return mat;
}


// One-body matrix elements
double cme_1_sigma() {
  FILE* in_file;
  in_file = fopen("ne20_density_1", "r");
  double mat = 0.0;

  for (int i = 0; i < 7; i++) {
    int n1p, l1p, ij1p;
    int n1, l1, ij1;
    double density;
    fscanf(in_file, "%d %d %d %d %d %d %lf\n", &n1p, &l1p, &ij1p, &n1, &l1, &ij1, &density); 
    double j1p = ij1p/2.0;
    double j1 = ij1/2.0;
    if (l1 != l1p) {continue;}
    if (n1 != n1p) {continue;}
   
    double m1 = sqrt(2)*sqrt(6.0)*sqrt((2*j1 + 1)*(2*j1p + 1))*pow(-1.0, l1 + 0.5 + j1p + 1)*six_j(0.5, j1p, l1, j1, 0.5, 1.0);
    printf("%d %d %g %d %d %g %g\n", n1p, l1p, j1p, n1, l1, j1, m1);
    mat += m1*density;
  } 
  return mat;
} 

/*******************************************/
//  Matrix elements for Master Equation    //
/*******************************************/


  
