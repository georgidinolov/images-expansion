#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "2DAdvectionDiffusionSolverImages.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 10) {
    printf("You must provide input\n");
    printf("The input is: \n sigma_x, \n sigma_y, \n t,\n a,\n x_T,\n b,\n c,\n y_T,\n d\n"); 
    exit(0);
  }
  
  TwoDAdvectionDiffusionSolverImages *images_solver;

  int order = 2000;
  double mu_x = 0;
  double mu_y = 0;

  double sigma_x = std::stod(argv[1]);
  double sigma_y = std::stod(argv[2]);
  double t = std::stod(argv[3]);
  double a = std::stod(argv[4]);
  double x = std::stod(argv[5]);
  double b = std::stod(argv[6]);
    
  double c = std::stod(argv[7]);
  double y = std::stod(argv[8]);
  double d = std::stod(argv[9]);

  images_solver = 
    new TwoDAdvectionDiffusionSolverImages(mu_x,
					   mu_y,
					   sigma_x,
					   sigma_y,
					   order,
					   0,
					   0,
					   a, 
					   b,
					   c, 
					   d);
  double h = 1.0/128.0;
  
  double solution = images_solver->solve(t,x,y);
  std::cout << "solution = " << solution << std::endl;
  std::cout << "solution_minus_h = " << images_solver->solve(t,x,y,a-h,b,c,d)
	    << std::endl;
  
  double likelihood = images_solver->likelihood(t,x,y);
  std::cout << "likelihood = " << likelihood << std::endl;

  
  likelihood = images_solver->likelihood(t,x,y);
  std::cout << "likelihood_numeric = "
	    << images_solver->numerical_likelihood(t,x,y, h) << std::endl;

  delete images_solver;
  return 0;
}
