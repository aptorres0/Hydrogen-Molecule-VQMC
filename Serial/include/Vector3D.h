#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include <iostream>

class Vector3D {
private:
  double fX;
  double fY;
  double fZ;

public:
  // Constructors
  Vector3D() : fX(0), fY(0), fZ(0){};
  Vector3D(double x, double y, double z) : fX(x), fY(y), fZ(z){};
  // Copy Constructor
  Vector3D(const Vector3D &v) : fX(v.fX), fY(v.fY), fZ(v.fZ){};

  // assignment operation
  Vector3D &operator=(const Vector3D &rhs);

  // Binary Operators
  Vector3D operator+(Vector3D const &rhs) const;
  Vector3D operator-(Vector3D const &rhs) const;
  Vector3D operator*(const double &rhs) const;
  friend Vector3D operator*(const double &lhs, const Vector3D &rhs) {
    return rhs * lhs;
  }

  // Function to print out vector
  void Print() const;

  // Getters
  double GetX() const { return fX; }
  double GetY() const { return fY; }
  double GetZ() const { return fZ; }

  // Get Separation Vector from Left proton
  Vector3D L() const;
  // Get Separation Vector from Right proton
  Vector3D R() const;
};

// Function to return magnitude of a vector
double Mag(Vector3D r);
// Function to compute the dot product of two vectors
double Dot(Vector3D r1, Vector3D r2);

#endif