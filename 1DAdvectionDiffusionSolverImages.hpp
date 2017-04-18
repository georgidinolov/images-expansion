// This is the header file 1DAdvectionDiffusionSolverFourier.hpp.  It
// contains the interface for the method for solving the advection
// diffusion equation in 1D using the method of images.

#include <vector>

class OneDAdvectionDiffusionSolverImages
{
public:
  OneDAdvectionDiffusionSolverImages();
  OneDAdvectionDiffusionSolverImages(double mu,
				     double sigma,
				     int order,
				     double IC,
				     double a,
				     double b);
  double solve(double t, double x);
  double solve(double t, double x, double a, double b);
  double numerical_likelihood(double t, double x, double h);
  double likelihood(double t, double x);

  std::vector<double> solve(double t, std::vector<double> x);
  std::vector<std::vector<double>> solve(std::vector<double> t, 
					 std::vector<double> x);

private:
  double mu_;
  double sigma_;
  double order_;
  double IC_;
  double a_;
  double b_;
  double L_;
};
