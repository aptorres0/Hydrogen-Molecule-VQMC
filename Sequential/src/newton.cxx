#include "newton.h"
#include <iostream>

#include <cmath>

// function to perform Newton's method to find zero of function
double newton(double (*func)(double), double xstart = 1.) {
  double tol = 0.001;
  int MaxIter = 1000;
  int i = 0;
  double xi = 0., x0 = xstart;
  double df, f = func(x0);
  double h = 0.01;

  // std::cout << "Initial condition: " << x0 << '\n';

  while (fabs(f) > tol && i < MaxIter) {
    f = func(x0);

    // Using forward difference scheme:
    df = (func(x0 + h) - func(x0)) / h;

    xi = x0 - f / df;
    x0 = xi;
    i++;

    // std::cout << "i = " << i << " : " << xi << '\n';

    // Backtracking
    //    If we're not in the linear portion of the function
    //    approximate at the mean between the two values
    while ((func(xi) * func(xi)) > (func(x0) * func(x0))) {
      xi = (x0 + xi) / 2.;
    }
  }

  return xi;
}