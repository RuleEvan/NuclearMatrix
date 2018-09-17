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

void generate_single_particle_states() {
  unsigned state;
  state = 0;
  int n_s = 12;
  int n_p = 5;
  int n_n = 6;
  int* m_shell = (int*)malloc(sizeof(int)*n_s);
  m_shell[0]=5;
  m_shell[1]=3;
  m_shell[2]=3;
  m_shell[3]=1;
  m_shell[4]=1;
  m_shell[5]=1;
  m_shell[6]=-1;
  m_shell[7]=-1;
  m_shell[8]=-1;
  m_shell[9]=-3;
  m_shell[10]=-3;
  m_shell[11]=-5; 
  int* last_state_p = (int*)malloc(sizeof(int)*n_p);
  for (int i = 0; i < n_p; i++) {
    last_state_p[i] = n_s - (n_p - i - 1);
  }
  int n_sds_p = p_step(n_s, n_p, last_state_p);
  printf("proton SDs: %d\n", n_sds_p);
  int* last_state_n = (int*)malloc(sizeof(int)*n_n);
  for (int i = 0; i < n_n; i++) {
    last_state_n[i] = n_s - (n_n - i - 1);
  }
  int n_sds_n = p_step(n_s, n_n, last_state_n);
  printf("neutron SDs: %d\n", n_sds_n);
  printf("Bin: %d\n", a_op(n_s, n_p, 1, 2));
  return;
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
