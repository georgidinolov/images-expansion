#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "1DAdvectionDiffusionSolverImages.hpp"

using namespace std;

int main ()
{
  OneDAdvectionDiffusionSolverImages *images_solver;

  int order = 1000;
  double sigma_x = 1.0;
  double a = -0.5;
  double b = 0.5;
  double x = 0.5;
  double t = 1;
  images_solver = 
    new OneDAdvectionDiffusionSolverImages(0,
					   sigma_x,
					   order,
					   0,
					   a,
					   b);
					   

  std::vector<double> xs = {-0.5, -0.25, 0 , 0.25, 0.5};
  for (unsigned i=0; i<xs.size(); ++i) {
    x = xs[i];
    double likelihood = images_solver->likelihood(t,x);
    std::cout << "likelihood = " << likelihood<< std::endl;
  }

  delete images_solver;
  return 0;
}
