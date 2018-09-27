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
  fclose(in_file);
  
}
