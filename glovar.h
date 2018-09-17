#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "stdint.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define ALPHA_FS 0.007297352
#define R_NUC 1.1 // [fm]
#define A_NUC 76 // Atomic Mass
#define M_NUC 2.19 // nuclear matrix element
#define M_ELECTRON 0.511 // [MeV]
#define M_NEUTRON 939.57 // [MeV]
#define Z_ATOM 52 // Atomic Number
#define E_BETA 3.038 // [MeV]
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))
#define XI_MIX_0 0.0016 // Mean mixing angle
#define G_AXIAL 1.275 // Axial coupling
#define E_HADRON 1 // Hadronization scale [GeV] 
#define DELTA_SUPP 1.0/30.0 // Suppression factor
#define M_W1 80.0 // W1 gauge boson mass [GeV]
#define M_W2_MIN 0.715 // Minimum W2 gauge boson mass [TeV]
#define MU 0.6995 // [fm]^-1

#define DENSITY_FILE "GE76_DENSITY.DAT"
#define NUM_SHELLS 99
#define POTENTIAL 1
#define A_FACTOR 9.155 // [MeV] Average nuclear excitation energy
#define B_OSC 2.927 // [fm] Oscillator parameter (times sqrt(2))
#define COR_FAC 1
#define PION_MASS 139.57 // [MeV] Charged pion mass
