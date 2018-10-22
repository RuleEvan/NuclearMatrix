#include "wave_function.h"

void generate_wave_function() {
  double hbarw;
  int n1, n2, n3, n4;
  int j1, j2, j3, j4;
  int jc, tc;
  double two_body_me, p_me, p_me_cm, p_me_rel;
  FILE *in_file;
  int num_shells;
  in_file = fopen("o16.dat", "r");
  fscanf(in_file, "%d", &num_shells);
  printf("%d\n", num_shells);
  int *n_shell = (int*) malloc(sizeof(int)*num_shells);
  int *j_shell = (int*) malloc(sizeof(int)*num_shells);
  for (int i = 0; i < num_shells; i++) {
    fscanf(in_file, "%d %d", &n_shell[i], &j_shell[i]);
  }
  int np, nn, mval, ipar, n_pre;
  double e_shift, c_jsq, c_tsq; 
  fscanf(in_file, "%d %d %d %d %d %lf %lf %lf\n", &np, &nn, &mval, &ipar, &n_pre, &e_shift, &c_jsq, &c_tsq);
  printf("%g\n", e_shift); 
  double xxcm;
  fscanf(in_file, "%lf\n", &xxcm);
  int *inph = (int*) malloc(sizeof(int)*3);
  for (int i = 0; i < 3; i++) {
    fscanf(in_file, "%d ", &inph[i]);
  }
  int itmax, nkeep, isaveb, istop9, isave4, iexp, itimer, ibasis, icoul, itcon;
  fscanf(in_file, "%d %d %d %d %d %d %d %d %d %d\n", &itmax, &nkeep, &isaveb, &istop9, &isave4, &iexp, &itimer, &ibasis, &icoul, &itcon);
  
  fclose(in_file);
  
}


