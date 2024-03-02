#include "Metropolis.h"
#include "H2/wavefunction.h"
#include "Vector3D.h"

// Initialize variables to start new random walk
void Metropolis::StartNewWalk() {
  d = std::uniform_real_distribution<double>(-0.5, 0.5);
  P = std::uniform_real_distribution<double>(0, 1);
  PDF_t = 0;
  Ratio_PDF = 0;
  sumE = 0;
  sumE2 = 0;
  ACCEPT_TRIAL = 0;
  N_A = 0;
  NumE = 0;

  r1 = Vector3D(d(rng), d(rng), d(rng));
  r2 = Vector3D(d(rng), d(rng), d(rng));
  PDF_k = pow(H2::psi(r1, r2), 2.0);
}

void Metropolis::Step() {
  r1t = r1 + _delta * Vector3D(d(rng), d(rng), d(rng));
  r2t = r2 + _delta * Vector3D(d(rng), d(rng), d(rng));
}

void Metropolis::SamplePDF() {
  using H2::psi;
  // Computing psi is computationally expensive so this should be faster
  PDF_t = pow(psi(r1t, r2t), 2.0);
}

void Metropolis::ApplyRejectionCriteria() {
  Ratio_PDF = PDF_t / PDF_k;
  if (Ratio_PDF >= 1) {
    N_A++;
    ACCEPT_TRIAL = 1;
  } else {
    double P0 = P(rng);
    if (Ratio_PDF >= P0) {
      N_A++;
      ACCEPT_TRIAL = 1;
    } else
      ACCEPT_TRIAL = 0;
  }
}

void Metropolis::SumEnergy() {
  using H2::Energy;
  if (ACCEPT_TRIAL) {
    int istep = N_A % NUM_SKIP;
    if (istep == 0) {
      r1 = r1t;
      r2 = r2t;
      PDF_k = PDF_t;
      double E_local = Energy(r1, r2);
      sumE += E_local;
      sumE2 += E_local * E_local;
      NumE++;
    }
  } else {
    double E_local = Energy(r1, r2);
    sumE += E_local;
    sumE2 += E_local * E_local;
    NumE++;
  }
  ACCEPT_TRIAL = 0;
}