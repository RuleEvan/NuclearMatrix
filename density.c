#include "density.h"

void one_body_density(int j_op, int t_op) {
  // Reads in BIGSTICK basis/wavefunction (.trwfn) files along with
  // orbit definition (.sp) files and constructs the one-body density matrices
  // for each initial and final state eigenfunction
  // Initial and final wave functions must share the same orbit file
  // J_op and T_op are the total spin and isospin of the operator
  wfnData *wd;
  wd = read_wfn_data();
  FILE* out_file;
  out_file = fopen("ge76_density_1", "w");
  printf("j_op: %d t_op: %d\n", j_op, t_op);
  // Loop over initial eigenstates
  for (int psi_i = 0; psi_i < wd->n_eig_i; psi_i++) {
    double ji = wd->j_nuc_i[psi_i];
    double ti = wd->t_nuc_i[psi_i];
    double mti = 0.5*(wd->n_proton_i - wd->n_neutron_i);
    // Loop over final eigenstates
    for (int psi_f = 0; psi_f < wd->n_eig_f; psi_f++) {
      double cg_j = 0.0;
      double cg_t = 0.0;
      double jf = wd->j_nuc_f[psi_f];
      double tf = wd->t_nuc_f[psi_f];
      double mtf = 0.5*(wd->n_proton_f - wd->n_neutron_f);
      cg_j = clebsch_gordan(j_op, ji, jf, 0, 0, 0);
      if (cg_j == 0.0) {continue;}
      cg_t = clebsch_gordan(t_op, ti, tf, 0, mti, mtf);
      if (cg_t == 0.0) {continue;}
      cg_j *= pow(-1.0, j_op + ji + jf)*sqrt(2*j_op + 1)/sqrt(2*jf + 1);
      cg_t *= pow(-1.0, t_op + ti + tf)*sqrt(2*t_op + 1)/sqrt(2*tf + 1);
      printf("Initial state: # %d J: %g T: %g Final state: # %d J: %g T: %g\n", psi_i + 1, ji, ti, psi_f + 1, jf, tf);
      // Loop over final state orbits
      for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
        double j1 = wd->j_orb[i_orb1];
        // Loop over initial state orbits
        for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
          double j2 = wd->j_orb[i_orb2];
          double total = 0.0;
          if ((j_op > j1 + j2) || (j_op < abs(j1 - j2))) {continue;}
          printf("Orbs: %d, %d\n", i_orb1, i_orb2);
          // Loop over initial state SDs
          for (int b = 0; b < wd->n_shells; b++) {
            if (wd->l_shell[b] != wd->l_orb[i_orb2]) {continue;}
            if (wd->n_shell[b] != wd->n_orb[i_orb2]) {continue;}
            if (wd->j_shell[b]/2.0 != j2) {continue;}
            double mj2 = wd->jz_shell[b]/2.0;
            double mt2 = wd->tz_shell[b]/2.0;
            int64_t b_b = pow(2, wd->n_shells - (b + 1));
            // Loop over initial state shells
            for (int a = 0; a < wd->n_shells; a++) {
              if (wd->l_shell[a] != wd->l_orb[i_orb1]) {continue;}
              if (wd->n_shell[a] != wd->n_orb[i_orb1]) {continue;}
              if (wd->j_shell[a]/2.0 != j1) {continue;}
              double mj1 = wd->jz_shell[a]/2.0;
              double mt1 = wd->tz_shell[a]/2.0;
              double d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
              if (d2 == 0.0) {continue;}

              for (int j = 0; j < wd->n_states_i; j++) {
                int64_t bi = wd->bc_i[j].b;
                if (!(b_b & bi)) {continue;} // Check if c_b|p> vanishes
                int phase = 1;
                if (a != b) {
                  bi = a2_op_b(wd->n_shells, bi, a + 1, b + 1, &phase);
                  if (bi == 0) {continue;}
                }
                
                int64_t i_min = 0;
                int64_t i_max = wd->n_states_f - 1;
                while (i_max >= i_min) { 
                  int64_t i = floor((i_max + i_min)/2.0);
                  int64_t bf = wd->bc_f[i].b;
                  if (bf == bi) {
                    total += wd->bc_f[i].wave[psi_f]*wd->bc_i[j].wave[psi_i]*phase*d2;
                    break;
                  } else if (bf > bi) {
                    i_max = i - 1;
                  } else {
                    i_min = i + 1;
                  }
                }
              }
            }
          }
          total /= cg_j*cg_t;
          if (total != 0.0) {
            //printf("%d %d %d %d %d %d %g\n", wd->n_shell[i_orb1], wd->l_shell[i_orb1], wd->j_shell[i_orb1], wd->n_shell[i_orb2], wd->l_shell[i_orb2], wd->j_shell[i_orb2], total);
            printf("%d %d %g\n", i_orb1 + 1, i_orb2 + 1, total);
          }
        }
      }
    }
  }          
  return;
  fclose(out_file);
}

void two_body_density(int j_op, int t_op) {
  // Read in data  
  wfnData *wd = read_wfn_data();
  double* j_store = (double*) malloc(sizeof(double)*4);
  // Loop over initial many-body wave functions
  for (int psi_i = 0; psi_i < wd->n_eig_i; psi_i++) {
    double ji = wd->j_nuc_i[psi_i];
    double ti = wd->t_nuc_i[psi_i];
    // Many-body states do not have mt = 0
    double mti = 0.5*(wd->n_proton_i - wd->n_neutron_i);
    // Loop over final many-body wave functions
    for (int psi_f = 0; psi_f < wd->n_eig_f; psi_f++) {
      double cg_j = 0.0;
      double cg_t = 0.0;
      double jf = wd->j_nuc_f[psi_f];
      double tf = wd->t_nuc_f[psi_f];
      // Many-body states do not have mt = 0
      double mtf = 0.5*(wd->n_proton_f - wd->n_neutron_f);
      double mt_op = mtf - mti;
      // CG factors for coupling initial and final states to operator
      cg_j = clebsch_gordan(j_op, ji, jf, 0, 0, 0);
      if (cg_j == 0.0) {continue;}
      cg_t = clebsch_gordan(t_op, ti, tf, mt_op, mti, mtf);
      if (cg_t == 0.0) {continue;}
      cg_j *= pow(-1.0, j_op + ji + jf)*sqrt(2*j_op + 1)/sqrt(2*jf + 1);
      cg_t *= pow(-1.0, t_op + ti + tf)*sqrt(2*t_op + 1)/sqrt(2*tf + 1);
      printf("Initial state: %d Final State: %d \n", psi_i + 1, psi_f + 1);
      double mat_test = 0.0;
      // Loop over orbital a
      for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
        double j1 = wd->j_orb[i_orb1];
        // Loop over orbital b
        for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
          double j2 = wd->j_orb[i_orb2];
          // Loop over orbital c
          for (int i_orb3 = 0; i_orb3 < wd->n_orbits; i_orb3++) {
            double j4 = wd->j_orb[i_orb3];
            // Loop over orbital d
            for (int i_orb4 = 0; i_orb4 < wd->n_orbits; i_orb4++) {
              double j3 = wd->j_orb[i_orb4];
              // Allocate storage for each j12 and j34
              int j_min_12 = abs(j1 - j2);
              int j_max_12 = j1 + j2;
              int j_dim_12 = j_min_12 + j_max_12 + 1;
              int j_min_34 = abs(j3 - j4);
              int j_max_34 = j3 + j4;
              int j_dim_34 = j_min_34 + j_max_34 + 1;
              int j_dim = j_dim_12*j_dim_34;
              j_store = realloc(j_store, sizeof(double)*4*j_dim);
              for (int k = 0; k < 4*j_dim; k++) {
                j_store[k] = 0.0;
              }
              // Loop over shells for orbit a
              for (int a = 0; a < wd->n_shells; a++) {
                if (wd->l_shell[a] != wd->l_orb[i_orb1]) {continue;}
                if (wd->n_shell[a] != wd->n_orb[i_orb1]) {continue;}
                if (wd->j_shell[a]/2.0 != j1) {continue;}
                double mj1 = wd->jz_shell[a]/2.0;
                double mt1 = wd->tz_shell[a]/2.0;
                 // Loop over shells for orbit b
                for (int b = 0; b < wd->n_shells; b++) {
                  if (a == b) {continue;}
                  if (wd->l_shell[b] != wd->l_orb[i_orb2]) {continue;}
                  if (wd->n_shell[b] != wd->n_orb[i_orb2]) {continue;}
                  if (wd->j_shell[b]/2.0 != j2) {continue;}
                  double mj2 = wd->jz_shell[b]/2.0;
                  double mt2 = wd->tz_shell[b]/2.0;
                  // Loop over coupled angular momentum j12
                  for (int j12 = (int) abs(j1 - j2); j12 <= (int) (j1 + j2); j12++) {
                    if ((mj1 + mj2 > j12) || (mj1 + mj2 < -j12)) {continue;}
                    // CG factor for coupling j1, j2 to j12
                    double cg_j12 = clebsch_gordan(j1, j2, j12, mj1, mj2, mj1 + mj2);
                    if (cg_j12 == 0.0) {continue;}
                    for (int t12 = 0; t12 <= 1; t12++) {
                      if ((mt1 + mt2 > t12) || (mt1 + mt2 < -t12)) {continue;}
                      double cg_t12 = clebsch_gordan(0.5, 0.5, t12, mt1, mt2, mt1 + mt2);
                      // Loop over shells for orbit c
                      for (int c = 0; c < wd->n_shells; c++) {
                        if (wd->l_shell[c] != wd->l_orb[i_orb3]) {continue;}
                        if (wd->n_shell[c] != wd->n_orb[i_orb3]) {continue;}
                        if (wd->j_shell[c]/2.0 != j4) {continue;}
                        double mj4 = wd->jz_shell[c]/2.0;
                        double mt4 = wd->tz_shell[c]/2.0;
                        // Loop over shells for orbit d
                        for (int d = 0; d < wd->n_shells; d++) {
                          if (c == d) {continue;}
                          if (wd->l_shell[d] != wd->l_orb[i_orb4]) {continue;}
                          if (wd->n_shell[d] != wd->n_orb[i_orb4]) {continue;}
                          if (wd->j_shell[d]/2.0 != j3) {continue;}
                          double mj3 = wd->jz_shell[d]/2.0;
                          double mt3 = wd->tz_shell[d]/2.0;
                          // Loop over coupled angular momentum j34
                          for (int j34 = abs(j3 - j4); j34 <= j3 + j4; j34++) {
                            if ((mj3 + mj4 > j34) || (mj3 + mj4 < -j34)){continue;}
                            double cg_j34 = clebsch_gordan(j3, j4, j34, -mj3, -mj4, -mj3 - mj4);
                            if (cg_j34 == 0.0) {continue;}
                            if (mj1 + mj2 -mj3 - mj4 != 0.0) {continue;}
                            double cg_jop = clebsch_gordan(j12, j34, j_op, mj1 + mj2, -mj3 - mj4, 0);
                            if (cg_jop == 0.0) {continue;}
                            for (int t34 = 0; t34 <= 1; t34++) {
                              if ((mt3 + mt4 > t34) || (mt3 + mt4 < -t34)){continue;}
                              double cg_t34 = clebsch_gordan(0.5, 0.5, t34, -mt3, -mt4, -mt3 - mt4);
                              if (cg_t34 == 0.0) {continue;}
                              double cg_top = clebsch_gordan(t12, t34, t_op, mt1 + mt2, -mt3 - mt4, mt_op);
                              if (cg_top == 0.0) {continue;}
                              double d2 = pow(-1.0, j3 + j4 + mj3 + mj4)*cg_j12*cg_j34*cg_jop;
                              d2 *= pow(-1.0, 1 + mt3 + mt4)*cg_t12*cg_t34*cg_top;
                              double d1 = 0.0;
                              // Loop over initial wave function basis
                              for (int j = 0; j < wd->n_states_i; j++) {
                                int64_t bi = wd->bc_i[j].b;
                                int64_t b_c = pow(2, wd->n_shells - (c + 1));
                                if (!(b_c & bi)) {continue;}
                                int64_t b_d = pow(2, wd->n_shells - (d + 1));
                                if (!(b_d & bi)) {continue;}
                                if ((a != c) && (a != d)) {
                                  int64_t b_a = pow(2, wd->n_shells - (a + 1));
                                  if (b_a & bi) {continue;}
                                }
                                if ((b != c) && (b != d)) {
                                  int64_t b_b = pow(2, wd->n_shells - (b + 1));
                                  if (b_b & bi) {continue;}
                                }
                                int phase;
                                bi = a4_op_b(wd->n_shells, bi, a + 1, b + 1, c + 1, d + 1, &phase);
                                if (bi == 0) {continue;}

                                /*pt = a_op(wd->n_shells, wd->n_data, pt, c + 1, &phase1);
                                if (pt == 0) {printf("is zero\n");}
                                pt = a_op(wd->n_shells, wd->n_data - 1, pt, d + 1, &phase2);
                                if (pt == 0) {printf("is zero\n");}
                                pt = a_dag_op(wd->n_shells, wd->n_data - 2, pt, b + 1, &phase3);
                                if (pt == 0) {printf("is zero\n");}
                                pt = a_dag_op(wd->n_shells, wd->n_data - 1, pt, a + 1, &phase4);
                                if (pt == 0) {printf("is zero\n");}
                                if (phase != phase1*phase2*phase3*phase4) {printf("%d %d %d %d %ld %ld %d %d\n", a, b, c, d, pi, pt, phase, phase1*phase2*phase3*phase4);} */

                                int i_min = 0;
                                int i_max = wd->n_states_f - 1;
                                while (i_max >= i_min) {
                                  int i = floor((i_max + i_min)/2.0);
                                  int64_t bf = wd->bc_f[i].b;
                                  if (bf == bi) {
                                    d1 += wd->bc_f[i].wave[psi_f]*wd->bc_i[j].wave[psi_i]*phase;
                                    break;
                                  } else if (bf > bi) {
                                    i_max = i-1;
                                  } else {
                                    i_min = i+1;
                                  }
                                }
                              }
                              d2 *= d1; 
                              //printf("About to write to jstore: %d, %d\n", j_dim, t12 + 2*(t34 + 2*(j12*h_dim + j34)));
                              j_store[t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34)))] += d2/(cg_t*cg_j);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              for (int t12 = 0; t12 <= 1; t12++) {
                for (int t34 = 0; t34 <= 1; t34++) {
                  for (int ij12 = 0; ij12 < j_dim_12; ij12++) {
                    for (int ij34 = 0; ij34 < j_dim_34; ij34++) {
                      //if (fabs(j_store[t12 + 2*(t34 + 2*(j12*h_dim + j34))]) < pow(10, -6)) {j_store[t12 + 2*(t34 + 2*(j12*h_dim + j34))] = 0.0; continue;}
                      if (fabs(j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]) > pow(10, -8)) {printf("%d %g %d %g %d %d %d %g %d %g %d %d %g \n", wd->n_orb[i_orb1], wd->j_orb[i_orb1], wd->n_orb[i_orb2], wd->j_orb[i_orb2], ij12 + j_min_12, t12, wd->n_orb[i_orb3], wd->j_orb[i_orb3], wd->n_orb[i_orb4], wd->j_orb[i_orb4], ij34 + j_min_34, t34, j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]);}
                      if ((i_orb1 == i_orb3) && (i_orb2 == i_orb4)) {
                        mat_test += pow(-1.0, j3 + j4 + j_min_12 + ij12 + t12)*j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]*sqrt(2*(j_min_12 + ij12) + 1)*sqrt(2*t12 + 1);
                      }
                      if ((i_orb1 == i_orb4) && (i_orb2 == i_orb3)) {
                        mat_test += j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]*sqrt(2*(j_min_12 + ij12) + 1)*sqrt(2*t12 + 1);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      printf("%g\n", 0.5*mat_test);
    }
  } 
  free(j_store); 
  return;
}

wfnData* read_wfn_data() {
  wfnData *wd = malloc(sizeof(*wd));
  FILE* in_file;
  // Read in initial wavefunction data
  in_file = fopen(WFN_FILE_INITIAL, "r");
  fscanf(in_file, "%d\n", &wd->n_proton_i);
  fscanf(in_file, "%d\n", &wd->n_neutron_i);
  char buffer[100];
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &wd->n_shells);
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &wd->n_states_i);
  printf("Initial state contains %d protons and %d neutrons\n", wd->n_proton_i, wd->n_neutron_i);
  printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_i);

  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%lf", &wd->jz_i);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &wd->n_eig_i);
  printf("Initial state file contains %d eigenstates\n", wd->n_eig_i);
  fgets(buffer, 100, in_file);
  wd->e_nuc_i =  (double*) malloc(sizeof(double)*wd->n_eig_i);
  wd->j_nuc_i =  (double*) malloc(sizeof(double)*wd->n_eig_i);
  wd->t_nuc_i =  (double*) malloc(sizeof(double)*wd->n_eig_i);
  for (int i = 0; i < wd->n_eig_i; i++) {
    fscanf(in_file, "%lf", &wd->e_nuc_i[i]);
    fscanf(in_file, "%lf", &wd->j_nuc_i[i]);
    fscanf(in_file, "%lf", &wd->t_nuc_i[i]);
    wd->j_nuc_i[i] = round(fabs(wd->j_nuc_i[i]));
    wd->t_nuc_i[i] = round(fabs(wd->t_nuc_i[i]));
    if (wd->n_data % 2) {
      wd->j_nuc_i[i] -= 0.5;
      wd->t_nuc_i[i] -= 0.5;
    }

    fgets(buffer, 100, in_file);
  }
  wd->n_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->j_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->l_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->jz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->tz_shell = (int*) malloc(sizeof(int)*wd->n_shells);

  // Read in shell quantum numbers
  for (int i = 0; i < wd->n_shells; i++) {
    fscanf(in_file, "%*d %d %d %d %d %d\n", &wd->n_shell[i], &wd->l_shell[i], &wd->j_shell[i], &wd->jz_shell[i], &wd->tz_shell[i]);
  }
  wd->n_data = wd->n_proton_i + wd->n_neutron_i;
  int* orbitals = (int*) malloc(sizeof(int)*wd->n_data);
  wd->bc_i = malloc(sizeof(BasisCoeff)*wd->n_states_i);
  printf("Reading in initial state wavefunction coefficients\n");
  double norm = 0.0;
  for (int i = 0; i < wd->n_states_i; i++) {
    int64_t b = 0;
    for (int j = 0; j < wd->n_data; j++) {
      fscanf(in_file, "%d", &orbitals[j]);
      b += pow(2, wd->n_shells - orbitals[j]);
    }
    int64_t p = p_step(wd->n_shells, wd->n_data, orbitals);
    wd->bc_i[i].p = p;
    wd->bc_i[i].b = b;
    wd->bc_i[i].wave = (double*) malloc(sizeof(double)*wd->n_eig_i);
    fgets(buffer, 100, in_file);
    for (int j = 0; j < wd->n_eig_i; j++) {
      fscanf(in_file, "%lf\n", &wd->bc_i[i].wave[j]);
    }
    norm += pow(wd->bc_i[i].wave[0], 2.0);
  }
  printf("Norm: %g\n", norm);
  qsort(wd->bc_i, wd->n_states_i, sizeof(BasisCoeff), s_compare);  

  fclose(in_file);
  
    

  if (strcmp(WFN_FILE_INITIAL, WFN_FILE_FINAL) == 0) {
    printf("Initial and final states are identical\n");
    wd->n_proton_f = wd->n_proton_i;
    wd->n_neutron_f = wd->n_neutron_i;
    wd->n_states_f = wd->n_states_i;
    wd->jz_f = wd->jz_i;
    wd->n_eig_f = wd->n_eig_i;    
    wd->e_nuc_f =  (double*) malloc(sizeof(double)*wd->n_eig_f);
    wd->j_nuc_f =  (double*) malloc(sizeof(double)*wd->n_eig_f);
    wd->t_nuc_f =  (double*) malloc(sizeof(double)*wd->n_eig_f);
    wd->e_nuc_f = wd->e_nuc_i;
    wd->j_nuc_f = wd->j_nuc_i;
    wd->t_nuc_f = wd->t_nuc_i;
    wd->bc_f = malloc(sizeof(BasisCoeff)*wd->n_states_f);
    wd->bc_f = wd->bc_i;
  } else {
    in_file = fopen(WFN_FILE_FINAL, "r");
    int n_shells_test;
    fscanf(in_file, "%d\n", &wd->n_proton_f);
    fscanf(in_file, "%d\n", &wd->n_neutron_f);
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%d", &n_shells_test);
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%d", &wd->n_states_f);
    printf("Final state contains %d protons and %d neutrons\n", wd->n_proton_f, wd->n_neutron_f);
    if (wd->n_shells != n_shells_test) {printf("Error: number of shells does not agree between initial and final state model spaces\n"); exit(0);}
    printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_f);
    if (wd->n_proton_f + wd->n_neutron_f != wd->n_proton_i + wd->n_neutron_i) {printf("Error: total number of nucleons is not constant\n"); exit(0);}
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%lf", &wd->jz_f);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%d", &wd->n_eig_f);
    printf("Final state file contains %d eigenvalues\n", wd->n_eig_f);
    fgets(buffer, 100, in_file);
    wd->e_nuc_f =  (double*) malloc(sizeof(double)*wd->n_eig_f);
    wd->j_nuc_f =  (double*) malloc(sizeof(double)*wd->n_eig_f);
    wd->t_nuc_f =  (double*) malloc(sizeof(double)*wd->n_eig_f);

    for (int i = 0; i < wd->n_eig_f; i++) {
      fscanf(in_file, "%lf", &wd->e_nuc_f[i]);
      fscanf(in_file, "%lf", &wd->j_nuc_f[i]);
      fscanf(in_file, "%lf", &wd->t_nuc_f[i]);
      wd->j_nuc_f[i] = round(fabs(wd->j_nuc_f[i]));
      wd->t_nuc_f[i] = round(fabs(wd->t_nuc_f[i]));
      if (wd->n_data % 2) {
        wd->j_nuc_f[i] -= 0.5;
        wd->t_nuc_f[i] -= 0.5;
      }
      fgets(buffer, 100, in_file);
    }
    for (int i = 0; i < wd->n_shells; i++) {
      int n_test, l_test, j_test, jz_test, tz_test;
      fscanf(in_file, "%*d %d %d %d %d %d\n", &n_test, &l_test, &j_test, &jz_test, &tz_test);
      if ((n_test != wd->n_shell[i]) || (l_test != wd->l_shell[i]) || (j_test != wd->j_shell[i]) || (jz_test != wd->jz_shell[i]) || (tz_test != wd->tz_shell[i])) {
        printf("Error shell structure is not identical between initial and final state model spaces\n");
        exit(0);
      }
    }
    wd->bc_f = malloc(sizeof(BasisCoeff)*wd->n_states_f);
    for (int i = 0; i < wd->n_states_f; i++) {
      int64_t b = 0;
      for (int j = 0; j < wd->n_data; j++) {
        fscanf(in_file, "%d", &orbitals[j]);
        b += pow(2, wd->n_shells - orbitals[j]);
      }
      int64_t p = p_step(wd->n_shells, wd->n_data, orbitals);
      wd->bc_f[i].b = b;
      wd->bc_f[i].p = p;
      wd->bc_f[i].wave = (double*) malloc(sizeof(double)*wd->n_eig_f);
      fgets(buffer, 100, in_file);
      for (int j = 0; j < wd->n_eig_f; j++) {
        fscanf(in_file, "%lf\n", &wd->bc_f[i].wave[j]);
      }
    }
    qsort(wd->bc_f, wd->n_states_f, sizeof(BasisCoeff), s_compare);  
  
    fclose(in_file);
  }

  in_file = fopen(ORBIT_FILE, "r");
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d\n", &wd->n_orbits);
  printf("The model space contains %d orbitals\n", wd->n_orbits);
  wd->n_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->l_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->j_orb = (double*) malloc(sizeof(double)*wd->n_orbits);

  for (int i = 0; i < wd->n_orbits; i++) {
    double n_orb_f, l_orb_f;
    fscanf(in_file, "%lf %lf %lf %*d", &n_orb_f, &l_orb_f, &wd->j_orb[i]);
    wd->n_orb[i] = (int) n_orb_f;
    wd->l_orb[i] = (int) l_orb_f;
    printf("%d, %d, %g\n", wd->n_orb[i], wd->l_orb[i], wd->j_orb[i]);
  }
  fclose(in_file);
  return wd;
}

int s_compare(const void * a, const void * b) {
  BasisCoeff *basisa = (BasisCoeff *)a;
  BasisCoeff *basisb = (BasisCoeff *)b;
  int64_t pa = basisa->b;
  int64_t pb = basisb->b;
  int s = 0;
  if (pa - pb > 0) {s = 1;}
  if (pa - pb < 0) {s = -1;}

  return s;
}

void exp_from_wfn() {

  FILE* in_file;
  in_file = fopen("ge76_basis.trwfn", "r");
  int n_proton, n_neutron, n_states, n_shells;
  // Read in number of protons and neutrons
  fscanf(in_file, "%d\n", &n_proton);
  fscanf(in_file, "%d\n", &n_neutron);
  // Define buffer for skipping through input files
  char buffer[100];
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file); 
  fgets(buffer, 100, in_file);
  // Read in number of shells
  fscanf(in_file, "%d", &n_shells);
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  // Read in number of Slater determinants in the basis
  fscanf(in_file, "%d", &n_states);
  printf("protons: %d neutrons: %d shells: %d states: %d\n", n_proton, n_neutron, n_shells, n_states);

  int *n_shell = (int*) malloc(sizeof(int)*n_shells);
  int *j_shell = (int*) malloc(sizeof(int)*n_shells);
  int *l_shell = (int*) malloc(sizeof(int)*n_shells);
  int *jz_shell = (int*) malloc(sizeof(int)*n_shells);
  int *tz_shell = (int*) malloc(sizeof(int)*n_shells);

  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  // Read magnetic quantum number of many-body wave function
  int jz;
  fscanf(in_file, "%d", &jz);
  fgets(buffer, 100, in_file);
  // Read number of many-body eigenstates
  int n_eig;
  fscanf(in_file, "%d", &n_eig);
  fgets(buffer, 100, in_file);
  // Read in the energy, spin and isopin of each energy eigenstate
  double *e_nuc = (double*) malloc(sizeof(double)*n_eig);
  int *j_nuc = (int*) malloc(sizeof(int)*n_eig);
  int *t_nuc = (int*) malloc(sizeof(int)*n_eig);
  for (int i = 0; i < n_eig; i++) {
    double j_float, t_float;
    fscanf(in_file, "%lf", &e_nuc[i]);
    fscanf(in_file, "%lf", &j_float);
    fscanf(in_file, "%lf", &t_float);
    j_nuc[i] = (int) j_float;
    t_nuc[i] = (int) t_float;
    fgets(buffer, 100, in_file);
  }
  // Read in the quantum numbers of each shell
  for (int i = 0; i < n_shells; i++) {
    fscanf(in_file, "%*d %d %d %d %d %d\n", &n_shell[i], &l_shell[i], &j_shell[i], &jz_shell[i], &tz_shell[i]);
  }
  int n_data = n_proton + n_neutron;
  
  // Read in the basis states and corresponding coefficients, store states according to their p-value (Whitehead)
  int *orbitals = (int*) malloc(sizeof(int)*n_data);
 
  BasisCoeff* p_wave = malloc(sizeof(BasisCoeff)*n_states);
  for (int i = 0; i < n_states; i++) {
    for (int j = 0; j < n_data; j++) {
      fscanf(in_file, "%d", &orbitals[j]);
    }
    int p = p_step(n_shells, n_data, orbitals);
    p_wave[i].p = p;
    p_wave[i].wave = (double*) malloc(sizeof(double)*n_eig);
    fgets(buffer, 100, in_file);
    for (int j = 0; j < n_eig; j++) {
      fscanf(in_file, "%lf\n", &p_wave[i].wave[j]);
    }
  }
  fclose(in_file);
  qsort(p_wave, n_states, sizeof(BasisCoeff), s_compare);  
  int ji = 0;
  int jf = 1;
  double cg1 = pow(-1.0, 1.0 + jf + ji)*sqrt(2*jf + 1)/clebsch_gordan(1, ji, jf, 0, 0, 0);
  int *orb_i = (int*) malloc(sizeof(int)*n_shells);
  double mat = 0.0;
  for (int j = 0; j < n_states; j++) {
    int pi = p_wave[j].p;
    orbitals_from_p(pi, n_shells, n_data, orb_i);
    for (int k = 0; k < n_data; k++) {
      double jk = j_shell[orb_i[k] - 1]/2.0;
      int lk = l_shell[orb_i[k] - 1];
      double mjk = jz_shell[orb_i[k] - 1]/2.0;
      int nk = n_shell[orb_i[k] - 1];
      int mtk = tz_shell[orb_i[k] - 1];
      // Case where jp = jk
      double cg2 = 0.0;
      for (int ml = -lk; ml <= lk; ml++) {
        for (int ms = -1; ms <= 1; ms += 2) {
          cg2 += ms*pow(clebsch_gordan(lk, 0.5, jk, ml, ms/2.0, mjk), 2.0);   
        }
      }
      mat += cg2*p_wave[j].wave[0]*p_wave[j].wave[24];
      // Now case where jp = jk +/- 1
      double jp = 0;
      if (jk == lk + 0.5) {
        jp = jk - 1.0;
      } else if (jk == lk - 0.5) {
        jp = jk + 1.0;
      } else {
        printf("Error: %g %d\n", jk, lk);
      }
      if (mjk > jp || mjk < -jp) {continue;}
      cg2 = 0.0;
      for (int ml = -lk; ml <= lk; ml++) {
        for (int ms = -1; ms <= 1; ms += 2) {
          cg2 += ms*clebsch_gordan(lk, 0.5, jp, ml, ms/2.0, mjk)*clebsch_gordan(lk, 0.5, jk, ml, ms/2.0, mjk);
        }
      }
      if (cg2 == 0.0) {continue;}
      int n_op = 0;
      for (int n = 0; n < n_shells; n++) {
        if ((n_shell[n] == nk) && (l_shell[n] == lk) && (j_shell[n]/2.0 == jp) && (tz_shell[n] == mtk) && (jz_shell[n]/2.0 == mjk)) {
          n_op = n + 1;
          break;
        }
      }
      if (n_op == 0) {continue;}
      int phase1, phase2;
      int pn = a_op(n_shells, n_data, pi, orb_i[k], &phase1);
      if (pn == 0) {printf("Something wrong\n");}
      pn = a_dag_op(n_shells, n_data - 1, pn, n_op, &phase2);
      if (pn == 0) {continue;}
      for (int i = 0; i < n_states; i++) {
        if (p_wave[i].p == pn) {
          mat += cg2*p_wave[j].wave[0]*p_wave[i].wave[24]*phase1*phase2;
          break;
        }
      }
    }
  }
  mat *= cg1;
  printf("%g\n", mat*mat);
  return;
}


