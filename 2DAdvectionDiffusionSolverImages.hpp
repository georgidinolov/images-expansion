// This is the header file 2DAdvectionDiffusionSolverFourier.hpp.  It
// contains the interface for the method for solving the advection
// diffusion equation in 2D using the method of images. 
// THERE IS NO CORRELATION!

#include "1DAdvectionDiffusionSolverImages.hpp"

using namespace std;

class TwoDAdvectionDiffusionSolverImages

{
public:
  TwoDAdvectionDiffusionSolverImages();
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
				     double d);

  double solve(double t, double x, double y);
  double solve(double t,
	       double x, double y,
	       double a, double b,
	       double c, double d);
  double numerical_likelihood(double t, double x, double y, 
  			      double h);
  double likelihood(double t, double x, double y);

private:
  OneDAdvectionDiffusionSolverImages solver_x_;
  OneDAdvectionDiffusionSolverImages solver_y_;
};
