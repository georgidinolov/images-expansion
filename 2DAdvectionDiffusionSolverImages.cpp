// This is the implementation of the Images solver.
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include "2DAdvectionDiffusionSolverImages.hpp"

using namespace std;

namespace{ 
  inline double square(double x) {
    return x * x;
  }
  
  inline double pi() {
    return std::atan(1)*4;
  }
}

TwoDAdvectionDiffusionSolverImages::TwoDAdvectionDiffusionSolverImages()
  : solver_x_(OneDAdvectionDiffusionSolverImages()),
    solver_y_(OneDAdvectionDiffusionSolverImages())
{}

TwoDAdvectionDiffusionSolverImages::
TwoDAdvectionDiffusionSolverImages(double mu_x,
				   double mu_y,
				   double sigma_x,
				   double sigma_y,
				   int order,
				   double IC_x,
				   double IC_y,
				   double a, 
				   double b,
				   double c, 
				   double d)
  : solver_x_(OneDAdvectionDiffusionSolverImages(mu_x,
						 sigma_x,
						 order,
						 IC_x,
						 a,b)),
    solver_y_(OneDAdvectionDiffusionSolverImages(mu_y,
						 sigma_y,
						 order,
						 IC_y,
						 c,d))
{}

// The summation is done in the log-scale because the eigenvalues for
// the problem are getting very small.
double TwoDAdvectionDiffusionSolverImages::solve(double t, double x, double y)
{
  double out = solver_x_.solve(t,x) *
    solver_y_.solve(t,y);
  return out;
}

double TwoDAdvectionDiffusionSolverImages::solve(double t, 
						 double x, double y,
						 double a, double b,
						 double c, double d)
{
  double out_1 = solver_x_.solve(t,x,a,b);
  double out_2 = solver_y_.solve(t,y,c,d);
  double out = out_1 * out_2;
  return out;
}

double TwoDAdvectionDiffusionSolverImages::
numerical_likelihood(double t, double x, double y, double h) {

  double derivative = solver_x_.numerical_likelihood(t,x,h)*
    solver_y_.solve(t,y);

  return derivative;
}

double TwoDAdvectionDiffusionSolverImages::
likelihood(double t, double x, double y) {
  // double out = solver_x_.likelihood(t,x) *
  //   solver_y_.likelihood(t,y);

  double out = solver_x_.likelihood(t,x) *
    solver_y_.likelihood(t,y);

  return out;
}
