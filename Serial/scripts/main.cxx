/**
 * main.cxx
 *
 *  Description:
 *    File to show implementation of Metropolis.cxx while varying parameter beta
 *    for the diatomic hydrogen molecule.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <random>
#include <type_traits>
#include <typeinfo>

#include "H2/wavefunction.h"
#include "Metropolis.h"
#include "Vector3D.h"
#include "newton.h"

/**
 * Function to compute the ground state energy of the hydrogen
 * molecule given a struct of input parameters.
 */
void ComputeE(arg_t args) {

  FILE *fptr;
  fptr = fopen("beta_results.txt", "a");

  Metropolis alg;
  alg.SetParameters(args);
  alg.StartNewWalk();
  for (int j = 0; j < 1'000'000; j++) {
    alg.Step();
    alg.SamplePDF();
    alg.ApplyRejectionCriteria();
    alg.SumEnergy();
  }
  fprintf(fptr, "%f %f %f %.15f %.15f %f %f\n", alg.GetS(), alg.GetBeta(),
          alg.GetDelta(), alg.GetE0_Final(), alg.GetE0_Final_STDDEV(),
          alg.GetNE(), alg.GetN_A());

  fclose(fptr);
}

int main() {

  FILE *fptr = fopen("beta_results.txt", "w");
  fprintf(fptr, "#s beta delta Efinal stddev NumE N_A\n");
  fclose(fptr);

  std::vector<arg_t> params;

  // Start Loop to Vary beta and delta
  double b_min = 0.2, b_max = 2., b_step = 0.001;
  int N_b = 1 + (b_max - b_min) / b_step;
  double d_min = 1.0, d_step = 0.1;
  double P_A, d_tmp;

  double s = 0.74 / 0.529, alpha = 2.0, delta, beta;

  params.reserve(N_b);

  for (int j = 0; j < N_b; j++) {
    beta = j * b_step + b_min;
    delta = 2.0;
    d_tmp = 3.;
    while (d_tmp > d_min) {
      P_A = ComputeAcceptanceProbability(d_tmp, beta, alpha, s);
      if (P_A >= .45 && P_A <= .55) {
        delta = d_tmp;
        break;
      }
      d_tmp -= d_step;
    }
    params.emplace_back(beta, delta, s, alpha);
  }

#pragma omp parallel for
  for (unsigned i = 0; i < params.size(); ++i) {
    ComputeE(params[i]);
  }

  return 0;
}