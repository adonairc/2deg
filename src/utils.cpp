#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_expint.h>

int main (int argc, char* argv[])
{

  std::ofstream fout;
  // fout.open("data.dat");

  double N = 100;
  double xi = -10;
  double xf = 10;
  double dx = (xf-xi+1)/N;


  for(int i = 0; i < N; i++){
    double x = xi + dx*i;
    double si = gsl_sf_Si(x);
    //double ci = gsl_sf_Ci(x);
    // fout << x << "," << si << "," << ci << result << std::endl;
    std::cout << x << "," << si << std::endl;
  }

  // fout.close();
  return 0;
}