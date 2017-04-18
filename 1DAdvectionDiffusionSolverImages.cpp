// This is the implementation of the Images solver.
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include "1DAdvectionDiffusionSolverImages.hpp"

using namespace std;

namespace{ 
  inline double square(double x) {
    return x * x;
  }
  
  inline double pi() {
    return std::atan(1.0)*4.0;
  }
}

OneDAdvectionDiffusionSolverImages::OneDAdvectionDiffusionSolverImages()
  : mu_(0),
    sigma_(1),
    order_(1),
    IC_(0.5),
    L_(1)
{
  // TODO(gdinolov): This should throw an exception.
  if (order_<1) {
    cout << "order_ < 1!!!" << endl;
  }
}

OneDAdvectionDiffusionSolverImages::OneDAdvectionDiffusionSolverImages(
  double mu,
  double sigma,
  int order,
  double IC,
  double a,
  double b)
  : mu_(mu),
    sigma_(sigma),
    order_(order),
    IC_(IC),
    a_(a),
    b_(b),
    L_(b_-a_)
{}

// The summation is done in the log-scale because the eigenvalues for
// the problem are getting very small.
double OneDAdvectionDiffusionSolverImages::solve(double t, double x)
{
  double out = solve(t,x,a_,b_);
  return out;
}

double OneDAdvectionDiffusionSolverImages::solve(double t, 
						  double x,
						  double a,
						  double b)
{
  int nn = 2*order_ + 1;
  double out = 0.0;

  for (int i=0; i<nn; ++i) {
    int n = i-order_;
    double d_1_i = pow(x-IC_-2.0*n*(b-a), 2) / (2*pow((sigma_*sqrt(t)),2));
    double d_2_i = pow(x+IC_-2.0*a-2.0*n*(b-a), 2) / (2*pow((sigma_*sqrt(t)),2));

    out = out +
      exp(-d_1_i) - exp(-d_2_i);
  }
  out = out * 1.0 / (sqrt(2*pi()) * (sigma_*sqrt(t))) * 
    exp( -(pow((mu_*t),2) - 2*(mu_*t)*(x-IC_)) / (2*pow((sigma_*sqrt(t)),2)) );

  return out;
}

double OneDAdvectionDiffusionSolverImages::
numerical_likelihood(double t, double x, double h) {

  double h2_inv = square(1.0/h);

  double sol_1 = solve(t,
		       x,
		       a_-h,
		       b_+h);

  double sol_2 = solve(t,
		       x,
		       a_,
		       b_+h);

  double sol_3 = solve(t,
		       x,
		       a_-h,
		       b_);

  double sol_4 = solve(t,
		       x,
		       a_,
		       b_);

  double derivative = -1 * h2_inv *
    (-1.0*sol_1 + sol_2 + sol_3 - sol_4);

  // double derivative = -1.0 / h *
  //    (sol_4 - sol_3);
  
  return derivative;
}

double OneDAdvectionDiffusionSolverImages::
likelihood(double t, double x) {
  int nn = 2*order_ + 1;

  // std::cout << "solution before diff = "
  // 	    << solve(t,x,a_,b_) << std::endl;
  
  double out = 0.0;
  for (int i=0; i<nn; ++i) {
    int n = i-order_;
    double d_1_i = pow(x-IC_-2.0*n*(b_-a_), 2) / (2*pow((sigma_*sqrt(t)),2));
    double d_2_i = pow(x+IC_-2.0*a_-2.0*n*(b_-a_), 2) / (2*pow((sigma_*sqrt(t)),2));

    out = out +
      4*pow(n,2)*(2*d_1_i-1)*exp(-d_1_i) - 
      4*n*(n-1)*(2*d_2_i-1)*exp(-d_2_i);
  }
  out = out * 1.0 / (sqrt(2*pi()) * pow((sigma_*sqrt(t)),3)) * 
    exp( -(pow((mu_*t),2) - 2*(mu_*t)*(x-IC_)) / (2*pow((sigma_*sqrt(t)),2)) );

  return out;
}

std::vector<double> OneDAdvectionDiffusionSolverImages::
solve(double t, 
      std::vector<double>x) {
  std::vector<double> out;
  int x_size = x.size();
  for (int i=0; i<x_size; ++i) {
    out.push_back(solve(t,x[i]));
  }
  return out;
}

std::vector<std::vector<double>> OneDAdvectionDiffusionSolverImages::
solve(std::vector<double> t, 
      std::vector<double> x) 
{ 
  std::vector<std::vector<double>> out;
  int t_size = t.size();
  for (int i=0; i<t_size; ++i) {
    out.push_back(solve(t[i], x));
  }
  return(out);
}
