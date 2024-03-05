#ifndef _METROPOLIS_H
#define _METROPOLIS_H

#include <iostream>
#include <random>

#include "Vector3D.h"
#include "constants.h"
#include "wavefunction.h"

class Metropolis {
  // Uniformly distributed random number
  // generator only used to seed the Mersenne
  // Twister engine
  std::random_device rd;
  // Mersenne Twister engine is a fast pseudo-random
  // number engine used to feed the uniform_real_distribution
  std::mt19937 rng;
  // Uniform distributions for the random walk (d) and
  // Von Neumann Rejection (P)
  std::uniform_real_distribution<double> d; // to use: d(rng)
  std::uniform_real_distribution<double> P;
  // Electron position vectors and trial vectors
  Vector3D r1, r2, r1t, r2t;
  // Values of Omega for the electron positions and the trial positions
  double PDF_k, PDF_t, Ratio_PDF;
  // Variables to store energy functional
  double sumE, sumE2;
  // Boolean flag for whether the trial was accepted or not
  int ACCEPT_TRIAL;
  // Step Size
  double _delta;
  // Number of accepted trials and number of times energy functional was added
  int N_A, NumE;
  // Number of N_A to skip on
  int NUM_SKIP = 11;

public:
  Metropolis() : rng(rd()){};

  Metropolis(double seed) : rng(seed){};

  // Getters
  Vector3D GetR1() const { return r1; }
  Vector3D GetR2() const { return r2; }
  double GetPDF() const { return PDF_k; }
  double GetE() const { return sumE; }
  double GetE2() const { return sumE2; }
  double GetNE() const { return NumE; }
  double GetE0_Final() const { return (sumE / NumE + 1. / H2::s) * 2 * 13.1; }
  double GetE0_Final_STDDEV() const {
    return sqrt(((sumE2 / NumE) - (sumE / NumE) * (sumE / NumE)) / NumE) * 2 *
           13.1;
  }

  // Define the step size for the random walk
  void SetStepSize(double delta) { _delta = delta; }

  // Use this to check for correlations
  void PrintRNG() { std::cout << d(rng) << " " << P(rng) << '\n'; }

  // Initialize variables to start new random walk
  void StartNewWalk();

  // Take random step for new trial position
  void Step();

  // Sample the PDF for the current and trial positions
  void SamplePDF();

  // Apply the rejection criteria on the trial position
  void ApplyRejectionCriteria();

  // Sum up the energies to perform MC Integration
  void SumEnergy();

  void PrintState() const {
    std::cout << "Number of Iterations = " << NumE << '\n';
    std::cout << "beta = " << H2::beta << "\tdelta = " << H2::delta
              << "\ta = " << H2::a << "\ts = " << H2::s << '\n';
    std::cout << "\tEgs = " << (sumE / NumE + 1 / H2::s) * 2 * 13.1 << " +/- "
              << sqrt(((sumE2 / NumE) - (sumE / NumE) * (sumE / NumE)) / NumE) *
                     2 * 13.1
              << '\n';
  }
};

#endif