#include "density.h"

void one_body_density() {
  // Reads in BIGSTICK basis/wavefunction (.trwfn) files along with
  // orbit definition (.sp) files and constructs the one-body density matrices
  // for each initial and final state eigenfunction
  // Initial and final wave functions must share the same orbit file
  
  FILE* in_file;
  in_file = fopen("ne20_basis.trwfn", "r");
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
  in_file = fopen("sd.sps", "r");
  fgets(buffer, 100, in_file);
  qsort(p_wave, n_states, sizeof(BasisCoeff), s_compare);  
  // Read in the n, l, j values of each orbit
  int n_orbits;
  fscanf(in_file, "%d\n", &n_orbits);
  int* n_orb = (int*) malloc(sizeof(int)*n_orbits);
  int* l_orb = (int*) malloc(sizeof(int)*n_orbits);
  double* j_orb = (double*) malloc(sizeof(double)*n_orbits);
  for (int i = 0; i < n_orbits; i++) {
    double n_orb_f, l_orb_f;
    fscanf(in_file, "%lf %lf %lf %*d", &n_orb_f, &l_orb_f, &j_orb[i]);
    n_orb[i] = (int) n_orb_f;
    l_orb[i] = (int) l_orb_f;
  }
  fclose(in_file);
  FILE* out_file;
  out_file = fopen("ne20_density_1", "w");
  // Loop over initial eigenstates
  for (int psi_i = 0; psi_i < n_eig; psi_i++) {
    int ji = j_nuc[psi_i];
    int ti = t_nuc[psi_i];
    // Loop over final eigenstates
    for (int psi_f = 0; psi_f < n_eig; psi_f++) {
      double cg_j = 0.0;
      double cg_t = 0.0;
      int jf = j_nuc[psi_f];
      int tf = t_nuc[psi_f];
      fprintf("Initial state: %d Final state: %d\n", psi_i, psi_f);
      for (int j_op = abs(ji - jf); j_op <= ji + jf; j_op++) {
        for (int t_op = abs(ti - tf); t_op <= ti + tf; t_op++) {
      cg_j = clebsch_gordan(j_op, ji, jf, 0, 0, 0);
      if (cg_j == 0.0) {continue;}
      cg_t = clebsch_gordan(t_op, ti, tf, 0, 0, 0);
      if (cg_t == 0.0) {continue;}
      fprintf(out_file, "j_op: %d t_op: %d\n", j_op, t_op);
      cg_j *= pow(-1.0, j_op + ji + jf)*sqrt(2*j_op + 1)/sqrt(2*jf + 1);
      cg_t *= pow(-1.0, t_op + ti + tf)*sqrt(2*t_op + 1)/sqrt(2*tf + 1);
      // Loop over final state orbits
      for (int i_orb1 = 0; i_orb1 < n_orbits; i_orb1++) {
        // Loop over initial state orbits
        for (int i_orb2 = 0; i_orb2 < n_orbits; i_orb2++) {
          double total = 0.0;
          // Loop over final state shells
          for (int a = 0; a < n_shells; a++) {
            int ij1 = j_shell[a];
            double j1 = ij1/2.0;
            if (l_shell[a] != l_orb[i_orb1]) {continue;}
            if (n_shell[a] != n_orb[i_orb1]) {continue;}
            if (j1 != j_orb[i_orb1]) {continue;}
            double mj1 = jz_shell[a]/2.0;
            double mt1 = tz_shell[a]/2.0;
            // Loop over initial state shells
            for (int b = 0; b < n_shells; b++) {
              int ij2 = j_shell[b];
              double j2 = ij2/2.0;
              if (l_shell[b] != l_orb[i_orb2]) {continue;}
              if (n_shell[b] != n_orb[i_orb2]) {continue;}
              if (j2 != j_orb[i_orb2]) {continue;}
              double mj2 = jz_shell[b]/2.0;
              double mt2 = tz_shell[b]/2.0;
              if ((j_op > j1 + j2) || (j_op < abs(j1 - j2))) {continue;}
              double density = 0.0;
              double d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2);
              d2 *= pow(-1.0, 0.5 - mt2);
              if (d2 == 0.0) {continue;}
              double d1 = 0.0;
              // Loop over initial state Slater determinants
              for (int j = 0; j < n_states; j++) {
                int pi = p_wave[j].p;
                int phase1, phase2;
                pi = a_op(n_shells, n_data, pi, b + 1, &phase1);
                if (pi == 0) {continue;}
                pi = a_dag_op(n_shells, n_data - 1, pi, a + 1, &phase2);
                if (pi == 0) {continue;}
                // Loop over final state Slater determinants
                int i_min = 0;
                int i_max = n_states - 1;
                while (1 == 1) {
                  int i = floor((i_max + i_min)/2.0);
                  int pf = p_wave[i].p;
                  //printf("%d %d %d\n", i, i_min, i_max);
                  if (pf == pi) {
                    d1 += p_wave[i].wave[psi_f]*p_wave[j].wave[psi_i]*phase1*phase2;
                    break;
                  } else if (pf > pi) {
                    i_max = i - 1;
                  } else {
                    i_min = i + 1;
                  }
                }
              }
              d2 *= d1;
              density += d2;
              total += density;
            }
          }
          total /= cg_j*cg_t;
          if (total != 0.0) {
            fprintf(out_file, "%d %d %d %d %d %d %g\n", n_shell[i_orb1], l_shell[i_orb1], j_shell[i_orb1], n_shell[i_orb2], l_shell[i_orb2], j_shell[i_orb2], total);
          }
        }
      }
    }
    }
    }
  }          
  return;
  fclose(out_file);
}

void two_body_density(int j_op, int t_op) {

  FILE* in_file;
  in_file = fopen("ne20_basis.trwfn", "r");
  int n_proton, n_neutron, n_states, n_shells;
  fscanf(in_file, "%d\n", &n_proton);
  fscanf(in_file, "%d\n", &n_neutron);
  char buffer[100];
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &n_shells);
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &n_states);
  printf("%d, %d, %d, %d\n", n_proton, n_neutron, n_states, n_shells);

  int *n_shell = (int*) malloc(sizeof(int)*n_shells);
  int *j_shell = (int*) malloc(sizeof(int)*n_shells);
  int *l_shell = (int*) malloc(sizeof(int)*n_shells);
  int *jz_shell = (int*) malloc(sizeof(int)*n_shells);
  int *tz_shell = (int*) malloc(sizeof(int)*n_shells);

  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  int jz;
  fscanf(in_file, "%d", &jz);
  fgets(buffer, 100, in_file);
  int n_eig;
  fscanf(in_file, "%d", &n_eig);
  fgets(buffer, 100, in_file);
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
  for (int i = 0; i < n_shells; i++) {
    fscanf(in_file, "%*d %d %d %d %d %d\n", &n_shell[i], &l_shell[i], &j_shell[i], &jz_shell[i], &tz_shell[i]);
  }
  int n_data = n_proton + n_neutron;
  int* orbitals = (int*) malloc(sizeof(int)*n_data);
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
  qsort(p_wave, n_states, sizeof(BasisCoeff), s_compare);  

  fclose(in_file);
  in_file = fopen("sd.sps", "r");
  fgets(buffer, 100, in_file);
  int n_orbits;
  fscanf(in_file, "%d\n", &n_orbits);
  int* n_orb = (int*) malloc(sizeof(int)*n_orbits);
  int* l_orb = (int*) malloc(sizeof(int)*n_orbits);
  double* j_orb = (double*) malloc(sizeof(double)*n_orbits);

  for (int i = 0; i < n_orbits; i++) {
    double n_orb_f, l_orb_f;
    fscanf(in_file, "%lf %lf %lf %*d", &n_orb_f, &l_orb_f, &j_orb[i]);
    n_orb[i] = (int) n_orb_f;
    l_orb[i] = (int) l_orb_f;
  }
  fclose(in_file);
  for (int psi_i = 0; psi_i < n_eig; psi_i++) {
    // Loop over initial many-body wave functions
    int ji = j_nuc[psi_i];
    int ti = t_nuc[psi_i];
    for (int psi_f = 0; psi_f < n_eig; psi_f++) {
      // Loop over final many-body wave functions
      double cg_j = 0.0;
      double cg_t = 0.0;
      int jf = j_nuc[psi_f];
      int tf = t_nuc[psi_f];
      cg_j = clebsch_gordan(j_op, ji, jf, 0, 0, 0);
      if (cg_j == 0.0) {continue;}
      cg_t = clebsch_gordan(t_op, ti, tf, 0, 0, 0);
      if (cg_t == 0.0) {continue;}
      cg_j *= pow(-1.0, j_op + ji + jf)*sqrt(2*j_op + 1)/sqrt(2*jf + 1);
      cg_t *= pow(-1.0, t_op + ti + tf)*sqrt(2*t_op + 1)/sqrt(2*tf + 1);
      printf("Initial state: %d Final State: %d \n", psi_i + 1, psi_f + 1);
      // Loop over orbital a
      for (int i_orb1 = 0; i_orb1 < n_orbits; i_orb1++) {
        double j1 = j_orb[i_orb1];
        // Loop over orbital b
        for (int i_orb2 = 0; i_orb2 <= i_orb1; i_orb2++) {
          double j2 = j_orb[i_orb2];
          // Loop over orbital c
          for (int i_orb3 = 0; i_orb3 <= i_orb2; i_orb3++) {
            double j3 = j_orb[i_orb3];
            // Loop over orbital d
            for (int i_orb4 = 0; i_orb4 < i_orb3; i_orb4++) {
              double j4 = j_orb[i_orb4];
              // Allocate storage for each j12 and j34
              int j_dim = j1 + j2 - abs(j1 - j2);
              int h_dim = j3 + j4 - abs(j3 - j4);
              j_dim *= h_dim;
              double* j_store = (double*) malloc(sizeof(double)*4*j_dim);
              for (int k = 0; k < 4*j_dim; k++) {
                j_store[k] = 0.0;
              }
              // Loop over shells for orbit a
              for (int a = 0; a < n_shells; a++) {
                if (l_shell[a] != l_orb[i_orb1]) {continue;}
                if (n_shell[a] != n_orb[i_orb1]) {continue;}
                if (j_shell[a]/2.0 != j1) {continue;}
                double mj1 = jz_shell[a]/2.0;
                double mt1 = tz_shell[a]/2.0;
                 // Loop over shells for orbit b
                for (int b = 0; b < n_shells; b++) {
                  if (l_shell[b] != l_orb[i_orb2]) {continue;}
                  if (n_shell[b] != n_orb[i_orb2]) {continue;}
                  if (j_shell[b]/2.0 != j2) {continue;}
                  double mj2 = jz_shell[b]/2.0;
                  double mt2 = tz_shell[b]/2.0;
                  // Loop over coupled angular momentum j12
                  for (int j12 = (int) abs(j1 - j2); j12 <= (int) (j1 + j2); j12++) {
                    if ((mj1 + mj2 > j12) || (mj1 + mj2 < -j12)) {continue;}
                    double cg_j12 = clebsch_gordan(j1, j2, j12, mj1, mj2, mj1 + mj2);
                    if (cg_j12 == 0.0) {continue;}
                    for (int t12 = 0; t12 <= 1; t12++) {
                      if ((mt1 + mt2 > t12) || (mt1 + mt2 < -t12)) {continue;}
                      double cg_t12 = clebsch_gordan(0.5, 0.5, t12, mt1, mt2, mt1 + mt2);
                      // Loop over shells for orbit c
                      for (int c = 0; c < n_shells; c++) {
                        if (l_shell[c] != l_orb[i_orb3]) {continue;}
                        if (n_shell[c] != n_orb[i_orb3]) {continue;}
                        if (j_shell[c]/2.0 != j3) {continue;}
                        double mj3 = jz_shell[c]/2.0;
                        double mt3 = tz_shell[c]/2.0;
                        // Loop over shells for orbit d
                        for (int d = 0; d < n_shells; d++) {
                          if (l_shell[d] != l_orb[i_orb4]) {continue;}
                          if (n_shell[d] != n_orb[i_orb4]) {continue;}
                          if (j_shell[d]/2.0 != j4) {continue;}
                          double mj4 = jz_shell[d]/2.0;
                          double mt4 = tz_shell[d]/2.0;
                          // Loop over coupled angular momentum j34
                          for (int j34 = abs(j3 - j4); j34 <= j3 + j4; j34++) {
                            if ((mj3 + mj4 > j34) || (mj3 + mj4 < -j34)){continue;}
                            double cg_j34 = clebsch_gordan(j3, j4, j34, mj3, mj4, mj3 + mj4);
                            if (cg_j34 == 0.0) {continue;}
                            double cg_jop = clebsch_gordan(j12, j34, j_op, mj1 + mj2, mj3 + mj4, 0);
                            if (cg_jop == 0.0) {continue;}
                            for (int t34 = 0; t34 <= 1; t34++) {
                              if ((mt3 + mt4 > t34) || (mt3 + mt4 < -t34)){continue;}
                              double cg_t34 = clebsch_gordan(0.5, 0.5, t34, mt3, mt4, mt3 + mt4);
                              double d2 = pow(-1.0, j3 + j4 + mj3 + mj4)*cg_j12*cg_j34*cg_jop*cg_t12*cg_t34;
                              double d1 = 0.0;
                              // Loop over initial wave function basis
                              for (int j = 0; j < n_states; j++) {
                                int pi = p_wave[j].p;
                                int phase1, phase2, phase3, phase4;
                                pi = a_op(n_shells, n_data, pi, c + 1, &phase1);
                                if (pi == 0) {continue;}
                                pi = a_op(n_shells, n_data - 1, pi, d + 1, &phase2);
                                if (pi == 0) {continue;}
                                pi = a_dag_op(n_shells, n_data - 2, pi, b + 1, &phase3);
                                if (pi == 0) {continue;}
                                pi = a_dag_op(n_shells, n_data - 1, pi, a + 1, &phase4);
                                if (pi == 0) {continue;}
                                int i_min = 0;
                                int i_max = n_states - 1;
                                while (1 == 1) {
                                  int i = ceil((i_max + i_min)/2.0);
                                  int pf = p_wave[i].p;
                                //  printf("%d, %d, %d, %d\n", pi, i, i_min, i_max);
                                  if (pf == pi) {
                                    d1 += p_wave[i].wave[psi_f]*p_wave[j].wave[psi_i]*phase1*phase2*phase3*phase4;
                                    break;
                                  } else if ((p_wave[i-1].p < pi) && (p_wave[i+1].p > pi)) {
                                    break;
                                  } else if (pf > pi) {
                                    i_max = i-1;
                                  } else {
                                    i_min = i+1;
                                  }
                                }
                              }
                              d2 *= d1;
                              j_store[t12 + 2*(t34 + 2*(j12*h_dim + j34))] += d2/(cg_t*cg_j);
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
                  for (int j12 = abs(j1 - j2); j12 <= j1 + j2; j12++) {
                    for (int j34 = abs(j3 - j4); j34 <= j3 + j4; j34++) {
                      if (fabs(j_store[t12 + 2*(t34 + 2*(j12*h_dim + j34))]) < pow(10, -6)) {j_store[t12 + 2*(t34 + 2*(j12*h_dim + j34))] = 0.0; continue;}
                      printf("%d %g %d %g %d %d %d %g %d %g %d %d %g \n", n_orb[i_orb1], j_orb[i_orb1], n_orb[i_orb2], j_orb[i_orb2], j12, t12, n_orb[i_orb3], j_orb[i_orb3], n_orb[i_orb4], j_orb[i_orb4], j34, t34, j_store[t12 + 2*(t34 + 2*(j12*h_dim + j34))]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  
return;
}

int s_compare(const void * a, const void * b) {
  BasisCoeff *basisa = (BasisCoeff *)a;
  BasisCoeff *basisb = (BasisCoeff *)b;

  return (basisa->p - basisb->p);
}

void exp_from_wfn() {

  FILE* in_file;
  in_file = fopen("ne20_basis.trwfn", "r");
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
  int j_op = 1;
  double cg1 = pow(-1.0, 1.0 + jf + ji)*sqrt(2*jf + 1)/clebsch_gordan(1, ji, jf, 0, 0, 0);
  int *orb_i = (int*) malloc(sizeof(int)*n_shells);
  int *orb_f = (int*) malloc(sizeof(int)*n_shells);
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

void exp_from_wfn_2body() {

  FILE* in_file;
  in_file = fopen("ne20_basis.trwfn", "r");
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
  int j_op = 1;
  double cg1 = pow(-1.0, 1.0 + jf + ji)*sqrt(2*jf + 1)/clebsch_gordan(1, ji, jf, 0, 0, 0);
  int *orb_i = (int*) malloc(sizeof(int)*n_shells);
  int *orb_f = (int*) malloc(sizeof(int)*n_shells);
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
