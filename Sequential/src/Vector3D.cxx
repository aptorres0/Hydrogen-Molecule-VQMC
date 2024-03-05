#include "Vector3D.h"

#include <cmath>
#include <iostream>

#include "H2/constants.h"

Vector3D Vector3D::operator+(Vector3D const &rhs) const {
  Vector3D result;
  result.fX = fX + rhs.fX;
  result.fY = fY + rhs.fY;
  result.fZ = fZ + rhs.fZ;
  return result;
}

Vector3D Vector3D::operator-(Vector3D const &rhs) const {
  Vector3D result;
  result.fX = fX - rhs.fX;
  result.fY = fY - rhs.fY;
  result.fZ = fZ - rhs.fZ;
  return result;
}

Vector3D Vector3D::operator*(const double &rhs) const {
  return Vector3D(fX * rhs, fY * rhs, fZ * rhs);
}

Vector3D &Vector3D::operator=(const Vector3D &rhs) {
  // check for self-assignment
  if (this == &rhs)
    return *this;
  fX = rhs.fX;
  fY = rhs.fY;
  fZ = rhs.fZ;
  return *this;
}

void Vector3D::Print() const {
  std::cout << fX << ", " << fY << ", " << fZ << '\n';
}

// Get Separation Vector from Left proton
Vector3D Vector3D::L() const {
  using H2::s;
  Vector3D result;
  result.fX = fX + s / 2.;
  result.fY = fY;
  result.fZ = fZ;
  return result;
}

// Get Separation Vector from Right proton
Vector3D Vector3D::R() const {
  using H2::s;
  Vector3D result;
  result.fX = fX - s / 2.;
  result.fY = fY;
  result.fZ = fZ;
  return result;
}

double Mag(Vector3D r) {
  // return sqrtl(powl(r.GetX(), 2.) + powl(r.GetY(), 2.) + powl(r.GetZ(), 2.));
  return sqrt(r.GetX() * r.GetX() + r.GetY() * r.GetY() + r.GetZ() * r.GetZ());
}

double Dot(Vector3D r1, Vector3D r2) {
  return (r1.GetX() * r2.GetX() + r1.GetY() * r2.GetY() +
          r1.GetZ() * r2.GetZ());
}
