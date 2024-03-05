#include "H2/wavefunction.h"

namespace H2 {

using H2::_alpha;
using H2::a;
using H2::beta;
using H2::s;

// Function to get paramater a from
double func(double x) { return (1. / (1 + exp(-s / x)) - x); }

double phi(Vector3D r) { return (exp(-Mag(r.L()) / a) + exp(-Mag(r.R()) / a)); }

double f12(Vector3D r12) {
  return (exp(Mag(r12) / (_alpha * (1. + beta * Mag(r12)))));
}

double psi(Vector3D r1, Vector3D r2) {
  return (phi(r1) * phi(r2) * f12(r1 - r2));
}

// Function to compute the energy functional
double Energy(Vector3D r1, Vector3D r2) {
  // Define the separation vectors
  Vector3D r1L = r1.L();
  Vector3D r1R = r1.R();
  Vector3D r2L = r2.L();
  Vector3D r2R = r2.R();
  Vector3D r12 = r1 - r2;

  double eta = 1. / (1. + beta * Mag(r12));

  double KE = -1. / (a * a) - pow(eta, 3.) * (0.25 * eta + 1. / Mag(r12));

  KE += ((exp(-Mag(r1L) / a) / Mag(r1L) + exp(-Mag(r1R) / a) / Mag(r1R)) +
         0.5 * eta * eta *
             (exp(-Mag(r1L) / a) * Dot(r1L, r12) / (Mag(r1L) * Mag(r12)) +
              exp(-Mag(r1R) / a) * Dot(r1R, r12) / (Mag(r1R) * Mag(r12)))) /
        (a * phi(r1));

  KE += ((exp(-Mag(r2L) / a) / Mag(r2L) + exp(-Mag(r2R) / a) / Mag(r2R)) -
         0.5 * eta * eta *
             (exp(-Mag(r2L) / a) * Dot(r2L, r12) / (Mag(r2L) * Mag(r12)) +
              exp(-Mag(r2R) / a) * Dot(r2R, r12) / (Mag(r2R) * Mag(r12)))) /
        (a * phi(r2));

  double PE = -(1. / Mag(r1L)) - (1. / Mag(r1R)) - (1. / Mag(r2L)) -
              (1. / Mag(r2R)) + (1. / Mag(r12));

  return (KE + PE);
}

} // namespace H2