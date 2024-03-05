#ifndef _WAVEFUNCTION_H
#define _WAVEFUNCTION_H

#include "Vector3D.h"
#include "constants.h"
#include <cmath>

namespace H2 {

// Function to get paramater a from
double func(double x);

double phi(Vector3D r);

double f12(Vector3D r12);

double psi(Vector3D r1, Vector3D r2);

// Function to compute the energy functional
double Energy(Vector3D r1, Vector3D r2);

} // namespace H2

#endif
