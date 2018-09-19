#include "slater.h"
int p_step(int n_s, int n_p, int *m_p) {
  int p = gsl_sf_choose(n_s, n_p);
  for (int i = 0; i < n_p; i++) {
    p -= n_choose_k(n_s-m_p[i], n_p - i);
  }
  return p;
}
void orbitals_from_p(int p, int n_s, int n_p, int* orbitals) {
  int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        orbitals[k] = j;
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }
  return;
}

unsigned int count_set_bits(unsigned int n) {
  unsigned int count = 0;
  while (n) {
    n & (n-1);
    count++;
  }
  return count;
}

void initialize_orbitals(int* n_shell, int* l_shell, int* j_shell, int* m_shell) {
  int n_s = 0;
  for (int q = 0; q <= N_OSC_QUANTA; q++) {
    int n_max;
    if ((q % 2) == 0) {n_max = q/2;}
    else {n_max = (q - 1)/2;}
    for (int n = 0; n <= n_max; n++) {
      int l = q - 2*n;
      for (int jj = abs(2*l-1); jj <= 2*l + 1; jj += 2) {
        for (int mjj = -jj; mjj <= jj; mjj += 2) {
          n_shell[n_s] = n;
          l_shell[n_s] = l;
          j_shell[n_s] = jj;
          m_shell[n_s] = mjj; 
          n_s++;
        }
      }
    }
  }
  
  return;
}

void generate_single_particle_states() {
  unsigned state;
  state = 0;
  int n_s = (N_OSC_QUANTA + 1)*(N_OSC_QUANTA + 2)*(N_OSC_QUANTA + 3)/3; 
  printf("Number of single particle states: %d\n", n_s);
  int* n_shell = (int*)malloc(sizeof(int)*n_s);
  int* l_shell = (int*)malloc(sizeof(int)*n_s);
  int* j_shell = (int*)malloc(sizeof(int)*n_s);
  int* m_shell = (int*)malloc(sizeof(int)*n_s);
  initialize_orbitals(n_shell, l_shell, j_shell, m_shell);
  int* p_state = (int*)malloc(sizeof(int)*N_PROTON);
  for (int i = 0; i < N_PROTON; i++) {
    p_state[i] = n_s - (N_PROTON - i - 1);
  }
  int n_sds_p = p_step(n_s, N_PROTON, p_state);
  printf("proton SDs: %d\n", n_sds_p);
  int* n_state = (int*)malloc(sizeof(int)*N_NEUTRON);
  for (int i = 0; i < N_NEUTRON; i++) {
    n_state[i] = n_s - (N_NEUTRON - i - 1);
  }
  int n_sds_n = p_step(n_s, N_NEUTRON, n_state);
  printf("neutron SDs: %d\n", n_sds_n);
  int match = 0;
  for (int p1 = 0; p1 < n_sds_p; p1++) {
    int m_tot_p1 = m_from_p(p1, n_s, N_PROTON, m_shell);
    for (int n1 = 0; n1 < n_sds_n; n1++) {
      int m_tot_n1 = m_from_p(n1, n_s, N_NEUTRON, m_shell);
      if (m_tot_p1 + m_tot_n1 != M_BASIS) {continue;}
      match++;
    }
  }
  printf("%d\n", match);
  double *hamiltonian = (double*)malloc(sizeof(double)*match*match);
  for (int i = 0; i < match*match; i++) {hamiltonian[i] = 0.0;}
  int k,l;
  k = 0;
  for (int p1 = 0; p1 < n_sds_p; p1++) {
    int m_tot_p1 = m_from_p(p1, n_s, N_PROTON, m_shell);
    for (int n1 = 0; n1 < n_sds_n; n1++) {
      int m_tot_n1 = m_from_p(n1, n_s, N_NEUTRON, m_shell);
      if ((m_tot_p1 + m_tot_n1) != M_BASIS) {continue;}
      l = 0;
      for (int p2 = 0; p2 < n_sds_p; p2++) {
        int m_tot_p2 = m_from_p(p2, n_s, N_PROTON, m_shell);
        for (int n2 = 0; n2 < n_sds_n; n2++) {
          int m_tot_n2 = m_from_p(n2, n_s, N_NEUTRON, m_shell);
          if ((m_tot_p2 + m_tot_n2) != M_BASIS) {continue;}
          hamiltonian[k*match + l] = k*l;        
          l++;
        }
      }
      k++;
    }
  }
  lanczos_eigenvalues(hamiltonian, match, N_LANCZOS);  

  return;
}

int m_from_p(int p, int n_s, int n_p, int* m_shell) {
  int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  int m_tot = 0;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        m_tot += m_shell[j];
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }
  return m_tot;
}


int bin_from_p(int n_s, int n_p, int p) {
  int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  int bin = 0;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        bin += pow(2, n_s - j);
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }
  return bin;
}

int bin_phase_from_p(int n_s, int n_p, int p, int n_op, int* phase) {
  int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  int bin = 0;
  *phase = 1;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        bin += pow(2, n_s - j);
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        if (j < n_op) {*phase *= 1;}
        break;
      }
    }
  }
  return bin;
}

int a_op(int n_s, int n_p, int p, int n_op) {
  int phase;
  phase = 1;
  unsigned int bin = bin_phase_from_p(n_s, n_p, p, n_op, &phase);
  unsigned int check = pow(2, n_s - n_op);
  if (!(bin & check)) {return 0;}
  bin -= (bin & check);
  printf("bin: %d phase: %d\n", bin, phase);
  return bin;
}

int a_dag_op(int n_s, int n_p, int p, int n_op) {
  unsigned int bin = bin_from_p(n_s, n_p, p);
  unsigned int check = pow(2, n_s - n_op);
  if (bin & check) {return 0;}
  bin += check;
  return bin;
}

int n_choose_k(int n, int k) {
  int c = 0;
  if (k > n) {return c;}
  c = gsl_sf_choose(n, k);
  return c;
}
