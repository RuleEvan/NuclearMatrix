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
  int* max_state = (int*) malloc(sizeof(int)*(wd->n_proton_i - 1));
  for (int i = 0; i < wd->n_proton_i - 1; i++) {
    max_state[i] = wd->n_shells - (wd->n_proton_i - i - 2);
  }
  int n_sds_p_int = p_step(wd->n_shells, wd->n_proton_i - 1, max_state);
  printf("Num intermediate proton sds: %d\n", n_sds_p_int);
  max_state = realloc(max_state, sizeof(int)*(wd->n_neutron_i - 1));
  for (int i = 0; i < wd->n_neutron_i - 1; i++) {
    max_state[i] = wd->n_shells - (wd->n_neutron_i - i - 2);
  }
  int n_sds_n_int = p_step(wd->n_shells, wd->n_neutron_i - 1, max_state);
  printf("Num intermediate neutron sds: %d\n", n_sds_n_int);
  free(max_state);

  sd_list **p_list_i = (sd_list**) malloc(sizeof(sd_list*)*wd->n_shells);
  sd_list **n_list_i = (sd_list**) malloc(sizeof(sd_list*)*wd->n_shells);
  printf("Building initial state proton jumps...\n");
  // For each shell, loop over all proton SDs
  // If c_a | SD > is non-zero, store the resulting SD
  for (int a = 0; a < wd->n_shells; a++) {
    p_list_i[a] = create_sd_node(0,0,1,NULL);
    for (int j = 1; j <= wd->n_sds_p_i; j++) {
      int phase;
      unsigned int pn = a_op(wd->n_shells, wd->n_proton_i, j, a + 1, &phase);
      if (pn == 0) {continue;}
      sd_append(p_list_i[a], j, pn, phase);
    }
  }
  printf("Done.\n");
  printf("Building final state proton jumps...\n");
  int *p_list_f = (int*) calloc(wd->n_shells*n_sds_p_int, sizeof(int));
  for (int a = 0; a < wd->n_shells; a++) {
    for (int j = 1; j <= wd->n_sds_p_f; j++) {
      int phase;
      unsigned int pn = a_op(wd->n_shells, wd->n_proton_f, j, a + 1, &phase);
      if (pn == 0) {continue;}
      p_list_f[(pn - 1) + a*n_sds_p_int] = phase*j;
    }
  }
  printf("Done.\n");

  printf("Building initial state neutron jumps...\n");
  
  for (int a = 0; a < wd->n_shells; a++) {
    n_list_i[a] = create_sd_node(0,0,1,NULL);;
    for (int j = 1; j <= wd->n_sds_n_i; j++) {
      int phase;
      int pf = a_op(wd->n_shells, wd->n_neutron_i, j, a + 1, &phase);
      if (pf == 0) {continue;}
      sd_append(n_list_i[a], j, pf, phase);
    }
  }
  printf("Done.\n");
  int *n_list_f = (int*) calloc(wd->n_shells*n_sds_n_int, sizeof(int));
  printf("Building final state neutron jumps...\n");
  
  for (int a = 0; a < wd->n_shells; a++) {
    for (int j = 1; j <= wd->n_sds_n_f; j++) {
      int phase;
      unsigned int pn = a_op(wd->n_shells, wd->n_neutron_f, j, a + 1, &phase);
      if (pn == 0) {continue;}
      n_list_f[(pn - 1) + a*n_sds_n_int] = phase*j;
    }
  }
  printf("Done.\n");
 
  // Loop over initial eigenstates
  for (int psi_i = 0; psi_i < 9; psi_i++) {
    double ji = wd->j_nuc_i[psi_i];
    double ti = wd->t_nuc_i[psi_i];
    double mti = 0.5*(wd->n_proton_i - wd->n_neutron_i);
    // Loop over final eigenstates
    for (int psi_f = 0; psi_f < 9; psi_f++) {
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
          // Loop over initial state SDs
          for (int b = 0; b < wd->n_shells; b++) {
            if (wd->l_shell[b] != wd->l_orb[i_orb2]) {continue;}
            if (wd->n_shell[b] != wd->n_orb[i_orb2]) {continue;}
            if (wd->j_shell[b]/2.0 != j2) {continue;}
            double mj2 = wd->jz_shell[b]/2.0;
            // Loop over initial state shells
            for (int a = 0; a < wd->n_shells; a++) {
              if (wd->l_shell[a] != wd->l_orb[i_orb1]) {continue;}
              if (wd->n_shell[a] != wd->n_orb[i_orb1]) {continue;}
              if (wd->j_shell[a]/2.0 != j1) {continue;}
              double mj1 = wd->jz_shell[a]/2.0;
              double mt1 = 0.5;
              double mt2 = 0.5;
              double d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
                double d1 = 0.0;
 
              if (d2 != 0.0) {
                sd_list *node = p_list_i[b]->next;
                // Loop over SDs resulting from c_a| SD > 
                while (node != NULL) {
                  // Get the resulting SD
                  int ppn = node->pn;
                  // Check if the resulting SD appears in the 
                  // Corresponding list of SDs resulting from C_b | SD >
                  int ppf = p_list_f[(ppn - 1) + a*n_sds_p_int];
                  if (ppf == 0) {node = node->next; continue;}
                  int ppi = node->pi;
                  int phase1 = node->phase;
                  int phase2 = 1;
                  // Negative SDs indicate negative phase
                  if (ppf < 0) {
                    ppf *= -1;
                    phase2 = -1;
                  }
                  node = node->next;
                  wf_list *node2 = wd->p_hash_i[ppi - 1];
                  wf_list *node3 = wd->p_hash_f[ppf - 1];
                  unsigned int index_i;
                  unsigned int index_f;
                  while (node2 != NULL) {
                    index_i = node2->index;
                    node3 = wd->p_hash_f[ppf - 1];
                    while (node3 != NULL) {
                      index_f = node3->index;
                      if (node2->p == node3->p) {
                        d1 += wd->bc_i[index_i].wave[psi_i]*wd->bc_f[index_f].wave[psi_f]*phase1*phase2;
                      }
                      node3 = node3->next;
                    }
                    node2 = node2->next;
                  }
         
                }
                total += d2*d1;
              }
              mt1 = -0.5;
              mt2 = -0.5;
              d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
              if (d2 != 0.0) {

                d1 = 0.0;
                sd_list* node = n_list_i[b]->next;
                while (node != NULL) {
                  int pnn = node->pn;
                  int pnf = n_list_f[(pnn - 1) + a*n_sds_n_int];
                  if (pnf == 0) {node = node->next; continue;}
                  int pni = node->pi;
                  int phase1 = node->phase;
                  int phase2 = 1;
                  if (pnf < 0) {
                    pnf *= -1;
                    phase2 = -1;
                  }
                  node = node->next;
                  wf_list *node2 = wd->n_hash_i[pni - 1];
                  wf_list *node3 = wd->n_hash_f[pnf - 1];
                  unsigned int index_i;
                  unsigned int index_f;
                  while (node2 != NULL) {
                    index_i = node2->index;
                    node3 = wd->n_hash_f[pnf - 1];
                    while (node3 != NULL) {
                      index_f = node3->index;
                      if (node2->p == node3->p) {
                        d1 += wd->bc_i[index_i].wave[psi_i]*wd->bc_f[index_f].wave[psi_f]*phase1*phase2;
                      }
                      node3 = node3->next;
                    }
                    node2 = node2->next;
                  }
         
                }
                total += d2*d1;
              }
              mt1 = 0.5;
              mt2 = -0.5;
              d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
              if (d2 != 0.0) {

                d1 = 0.0;
                sd_list* node = n_list_i[b]->next;
                while (node != NULL) {
                  int ppf = node->pn;
                  int ppi = node->pi; 
                  int phase1 = node->phase;
                  // Get list of n_f associated to p_f
                  wf_list *node2 = wd->n_hash_f[ppf - 1];
                  unsigned int index_i;
                  unsigned int index_f;
                  // Loop over n_f 
                  while (node2 != NULL) {
                    index_i = node2->index;
                    unsigned int pnf = node2->p;
                    wf_list *node3 = wd->n_hash_i[ppi - 1];
                    while (node3 != NULL) {
                      unsigned int pni = node3->p;
                      int phase2 = 1;
                      if (n_list_f[(pnf - 1) + a*n_sds_n_int] == pni) {
                        index_f = node3->index;
                        d1 += wd->bc_i[index_i].wave[psi_i]*wd->bc_f[index_f].wave[psi_f]*phase1*phase2;
                      } else if (n_list_f[(pnf - 1) + a*n_sds_n_int] == -pni) {
                        phase2 = -1;
                        index_f = node3->index;
                        d1 += wd->bc_i[index_i].wave[psi_i]*wd->bc_f[index_f].wave[psi_f]*phase1*phase2;
                      }

                      node3 = node3->next;
                    }
                    node2 = node2->next;
                  }
                }  
              }
              total += d1*d2;
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
  int* max_state = (int*) malloc(sizeof(int)*(wd->n_proton_i - 1));
  for (int i = 0; i < wd->n_proton_i - 1; i++) {
    max_state[i] = wd->n_shells - (wd->n_proton_i - i - 2);
  }
  int n_sds_p_int = p_step(wd->n_shells, wd->n_proton_i - 1, max_state);
  printf("Num intermediate proton sds: %d\n", n_sds_p_int);
  max_state = realloc(max_state, sizeof(int)*(wd->n_neutron_i - 1));
  for (int i = 0; i < wd->n_neutron_i - 1; i++) {
    max_state[i] = wd->n_shells - (wd->n_neutron_i - i - 2);
  }
  int n_sds_n_int = p_step(wd->n_shells, wd->n_neutron_i - 1, max_state);
  printf("Num intermediate neutron sds: %d\n", n_sds_n_int);
  free(max_state);

  int ns = wd->n_shells;
  sd_list **p_list_i = (sd_list**) malloc(sizeof(sd_list*)*ns*(ns + 1)/2);
  sd_list **n_list_i = (sd_list**) malloc(sizeof(sd_list*)*ns*(ns + 1)/2);
  printf("Building initial state proton jumps...\n");
  for (int c = 0; c < ns; c++) {
    for (int d = 0; d < c; d++) {
      p_list_i[c + d*ns] = create_sd_node(0,0,1,NULL);
      for (int j = 1; j <= wd->n_sds_p_i; j++) {
        int phase1, phase2;
        int pf = a_op(ns, wd->n_proton_i, j, c + 1, &phase1);
        if (pf == 0) {continue;}
        pf = a_op(ns, wd->n_proton_i - 1, pf, d + 1, &phase2);
        if (pf == 0) {continue;}
        sd_append(p_list_i[c + d*ns], j, pf, phase1*phase2);
      }
    }
  }
  printf("Done.\n");
  printf("Building final state proton jumps...\n");
  int *p_list_f = (int*) calloc(ns*(ns + 1)/2*n_sds_p_int, sizeof(int));
  for (int a = 0; a < ns; a++) {
    for (int b = 0; b < a; b++) {
      for (int j = 1; j <= wd->n_sds_p_f; j++) {
        int phase1, phase2;
        int pf = a_op(ns, wd->n_proton_f, j, a + 1, &phase1);
        if (pf == 0) {continue;}
        pf = a_op(ns, wd->n_proton_f - 1, pf, b + 1, &phase2);
        if (pf == 0) {continue;}
        p_list_f[(pf - 1) + n_sds_p_int*(a + b*ns)] = phase1*phase2*j;
      }
    }
  }
  printf("Done.\n");

  printf("Building initial state neutron jumps...\n");
  
  for (int c = 0; c < ns; c++) {
    for (int d = 0; d < c; d++) {
      n_list_i[c + d*ns] = create_sd_node(0,0,1,NULL);;
      for (int j = 1; j <= wd->n_sds_n_i; j++) {
        int phase1, phase2;
        int pf = a_op(ns, wd->n_neutron_i, j, c + 1, &phase1);
        if (pf == 0) {continue;}
        pf = a_op(ns, wd->n_neutron_i - 1, pf, d + 1, &phase2);
        if (pf == 0) {continue;}
        sd_append(n_list_i[c + d*ns], j, pf, phase1*phase2);
      }
    }
  }
  printf("Done.\n");
  int *n_list_f = (int*) calloc(ns*(ns + 1)*n_sds_n_int, sizeof(int));
  printf("Building final state neutron jumps...\n");
  
  for (int a = 0; a < ns; a++) {
    for (int b = 0; b < a; b++) {
      for (int j = 1; j <= wd->n_sds_n_f; j++) {
        int phase1, phase2;
        int pf = a_op(ns, wd->n_neutron_f, j, a + 1, &phase1);
        if (pf == 0) {continue;}
        pf = a_op(ns, wd->n_neutron_f - 1, j, b + 1, &phase2);
        if (pf == 0) {continue;}
        n_list_f[(pf - 1) + n_sds_n_int*(a + b*ns)] = phase1*phase2*j;
      }
    }
  }
  printf("Done.\n");

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
                  int swp1 = 1;
                  if (b > a) {
                    int f = a;
                    a = b;
                    b = f;
                    swp1 = -1;
                  }
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
                          int swp2 = 1;
                          if (d > c) {
                            int f = c;
                            c = d;
                            d = f;
                            swp2 = -1;
                          }
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
                              //printf("About to write to jstore: %d, %d\n", j_dim, t12 + 2*(t34 + 2*(j12*h_dim + j34)));
                              j_store[t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34)))] += d1*d2/(cg_t*cg_j);
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
  wd->n_data = wd->n_proton_i + wd->n_neutron_i;
  char buffer[100];
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &wd->n_shells);
  wd->n_shells /= 2; // Mod out isospin
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
  //wd->tz_shell = (int*) malloc(sizeof(int)*wd->n_shells);

  // Read in shell quantum numbers
  for (int i = 0; i < wd->n_shells; i++) {
    fscanf(in_file, "%*d %d %d %d %d %*d\n", &wd->n_shell[i], &wd->l_shell[i], &wd->j_shell[i], &wd->jz_shell[i]);
  }
  // Skip shells with opposite isospin
  for (int i = 0; i < wd->n_shells; i++) {
    fscanf(in_file, "%*d %*d %*d %*d %*d %*d\n");
  }

  wd->bc_i = malloc(sizeof(BasisCoeff)*wd->n_states_i);
  printf("Reading in initial state wavefunction coefficients\n");
  double norm = 0.0;
  int* max_state = (int*) malloc(sizeof(int)*wd->n_proton_i);
  for (int i = 0; i < wd->n_proton_i; i++) {
    max_state[i] = wd->n_shells - (wd->n_proton_i - i - 1);
  }
  wd->n_sds_p_i = p_step(wd->n_shells, wd->n_proton_i, max_state);
  printf("Num initial proton SDs: %d\n", wd->n_sds_p_i);
  max_state = realloc(max_state, sizeof(int)*wd->n_neutron_i);
  for (int i = 0; i < wd->n_neutron_i; i++) {
    max_state[i] = wd->n_shells - (wd->n_neutron_i - i - 1);
  }
  wd->n_sds_n_i = p_step(wd->n_shells, wd->n_neutron_i, max_state);

  printf("Num initial neutron SDs: %d\n", wd->n_sds_n_i);
  free(max_state);
  wd->p_hash_i = (wf_list**) calloc(wd->n_sds_p_i, sizeof(wf_list*));
  wd->n_hash_i = (wf_list**) calloc(wd->n_sds_n_i, sizeof(wf_list*));
  int *p_orbitals = (int*) malloc(sizeof(int)*wd->n_proton_i);
  int *n_orbitals = (int*) malloc(sizeof(int)*wd->n_neutron_i);
  for (int i = 0; i < wd->n_states_i; i++) {
    int in = 0;
    int ip = 0;
    for (int j = 0; j < wd->n_data; j++) {
      int i_orb;
      fscanf(in_file, "%d", &i_orb);
      if (i_orb <= wd->n_shells) {
        p_orbitals[ip] = i_orb;
        ip++;
      } else {
        n_orbitals[in] = i_orb - wd->n_shells;
        in++;
      }
    }
    unsigned int pp = p_step(wd->n_shells, wd->n_proton_i, p_orbitals);
    unsigned int pn = p_step(wd->n_shells, wd->n_neutron_i, n_orbitals);
    if (wd->p_hash_i[pp - 1] == NULL) {
      wd->p_hash_i[pp - 1] = create_wf_node(pn, i, NULL);
    } else {
      wf_append(wd->p_hash_i[pp - 1], pn, i);
    }
    if (wd->n_hash_i[pn - 1] == NULL) {
      wd->n_hash_i[pn - 1] = create_wf_node(pp, i, NULL);
    } else {
      wf_append(wd->n_hash_i[pn - 1], pp, i);
    }

    wd->bc_i[i].wave = (double*) malloc(sizeof(double)*wd->n_eig_i);
    fgets(buffer, 100, in_file);
    for (int j = 0; j < wd->n_eig_i; j++) {
      fscanf(in_file, "%lf\n", &wd->bc_i[i].wave[j]);
    }
    norm += pow(wd->bc_i[i].wave[0], 2.0);
  }
  printf("Norm: %g\n", norm);
  free(p_orbitals);
  free(n_orbitals);
  fclose(in_file);
  
    

  if (strcmp(WFN_FILE_INITIAL, WFN_FILE_FINAL) == 0) {
    printf("Initial and final states are identical\n");
    wd->same_basis = 1;
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
//    wd->bc_f = malloc(sizeof(BasisCoeff)*wd->n_states_f);
    wd->bc_f = wd->bc_i;
    wd->p_hash_f = wd->p_hash_i;
    wd->n_hash_f = wd->n_hash_i;
    wd->n_sds_n_f = wd->n_sds_n_i;
    wd->n_sds_p_f = wd->n_sds_p_i;
    wd->n_sd_f = wd->n_sd_i;
    wd->p_sd_f = wd->p_sd_i;
  } else {
    wd->same_basis = 0;
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
    int* max_state = (int*) malloc(sizeof(int)*wd->n_proton_f);
    for (int i = 0; i < wd->n_proton_f; i++) {
      max_state[i] = wd->n_shells - (wd->n_proton_f - i - 1);
    }
    wd->n_sds_p_f = p_step(wd->n_shells, wd->n_proton_f, max_state);
    printf("Num final proton SDs: %d\n", wd->n_sds_p_f);
    max_state = realloc(max_state, sizeof(int)*wd->n_neutron_f);
    for (int i = 0; i < wd->n_neutron_f; i++) {
      max_state[i] = wd->n_shells - (wd->n_neutron_f - i - 1);
    }
    wd->n_sds_n_f = p_step(wd->n_shells, wd->n_neutron_f, max_state);
    printf("Num initial neutron SDs: %d\n", wd->n_sds_n_f);
    free(max_state);
    norm = 0.0;
    for (int i = 0; i < wd->n_states_f; i++) {
      unsigned int bn = 0;
      unsigned int bp = 0;
      for (int j = 0; j < wd->n_data; j++) {
        int i_orb;
        fscanf(in_file, "%d", &i_orb);
        if (i_orb <= wd->n_shells) {
          bp += pow(2, wd->n_shells - i_orb);
        } else {
          bn += pow(2, 2*wd->n_shells - i_orb);
        }
      }
      unsigned int bh = bp ^ bn;
      if (wd->p_hash_f[bh] == NULL) {
        wd->p_hash_f[bh] = create_wf_node(bp, i, NULL);
      } else {
        wf_append(wd->p_hash_f[bh], bn, i);
      }
      wd->bc_f[i].wave = (double*) malloc(sizeof(double)*wd->n_eig_f);
      fgets(buffer, 100, in_file);
      for (int j = 0; j < wd->n_eig_f; j++) {
        fscanf(in_file, "%lf\n", &wd->bc_f[i].wave[j]);
      }
      norm += pow(wd->bc_f[i].wave[0], 2.0);
    }
    printf("Final state norm: %g\n", norm);

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


sd_list* create_sd_node(int pi, int pn, int phase, sd_list* next) {
  sd_list* new_node = (sd_list*)malloc(sizeof(sd_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->pi = pi;
  new_node->pn = pn;
  new_node->phase = phase;
  new_node->next = next;

  return new_node;
}

sd_list* sd_append(sd_list* head, int pi, int pn, int phase) {
  sd_list* current = head;
  while (current->next != NULL) {
    current = current->next;
  }

  sd_list* new_node = create_sd_node(pi, pn, phase, NULL);
  current->next = new_node;

  return head;
}

wf_list* create_wf_node(unsigned int p, unsigned int index, wf_list* next) {
  wf_list* new_node = (wf_list*)malloc(sizeof(wf_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->p = p;
  new_node->index = index;
  new_node->next = next;
  
  return new_node;
}

wf_list* wf_append(wf_list* head, unsigned int p, unsigned int index) {
  wf_list* current = head;
  while (current->next != NULL) {
    current = current->next;
  }

  wf_list* new_node = create_wf_node(p, index, NULL);
  current->next = new_node;

  return head;
}
