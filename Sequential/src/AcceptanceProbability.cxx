/**
 * AcceptanceProbability.cxx
 *
 *  Description:
 *      Function to compute the acceptance probability given a set of input
 *      parameters for the Metropolis algorithim.
 */

#include "H2/wavefunction.h"
#include "Metropolis.h"
#include "Vector3D.h"
#include "newton.h"
#include <cmath>
#include <cstdlib>
#include <random>

double ComputeAcceptanceProbability(double delta, double beta, double alpha,
                                    double s) {
  using PAR::psi;

  int Iterations = 300;

  double a = newton(PAR::func, 1., s);

  // This pseudo-rng should be thread safe
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> d(-0.5, 0.5); // to use: d(rng)
  std::uniform_real_distribution<double> P(0, 1);      // to use: d(rng)

  Vector3D r1(d(rng), d(rng), d(rng));
  Vector3D r2(d(rng), d(rng), d(rng));

  double PDF_k, PDF_t;
  PDF_k = pow(psi(r1, r2, a, alpha, beta, s), 2.0);

  int N_A = 0;

  for (int i = 0; i < Iterations; i++) {
    Vector3D r1t = r1 + delta * Vector3D(d(rng), d(rng), d(rng));
    Vector3D r2t = r2 + delta * Vector3D(d(rng), d(rng), d(rng));

    PDF_t = pow(psi(r1t, r2t, a, alpha, beta, s), 2.0);

    if ((PDF_t / PDF_k) >= 1) {
      r1 = r1t;
      r2 = r2t;
      PDF_k = PDF_t;
      N_A++;
    } else {
      double P0 = P(rng);
      if ((PDF_t / PDF_k) >= P0) {
        r1 = r1t;
        r2 = r2t;
        PDF_k = PDF_t;
        N_A++;
      }
    }
  }
  return ((1.0 * N_A) / Iterations);
}