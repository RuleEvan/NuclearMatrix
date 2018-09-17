#include "nuc_mat.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main() {
  
//  printf("%g\n", (compute_matrix_element_TT(4)));
  printf("%g\n", compute_matrix_element_M_0());

  return 0;
}
