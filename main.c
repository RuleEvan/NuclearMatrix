#include "master_eq.h"

int main(int argc, char *argv[]) {
 // one_body_density();
//  printf("%g\n", pow(cme_1_sigma(), 2.0));
  two_body_density(2, 0);
  printf("%g\n", compute_matrix_element_M_JF());
  return 0;
}
