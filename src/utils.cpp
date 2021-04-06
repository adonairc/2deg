#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <string>
#include <vector>
#include <algorithm>
template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

int main (int argc, char* argv[])
{

  int N;
  std::vector<double> jx;
  std::vector<double> jy;

  std::ifstream infile(argv[1]);

  std::vector<double> x_coords;
  std::vector<double> y_coords;



  for( std::string line; getline( infile, line ); )
  {
    if (line[0] != '#'){
      std::stringstream ss(line);
      std::vector<double> v;
      while(ss.good()) {
        std::string substr;
        getline(ss, substr, ','); //get first string delimited by comma
        v.push_back(stod(substr));
      }
      x_coords.push_back(v[0]);
      y_coords.push_back(v[1]);
      jx.push_back(v[2]);
      jy.push_back(v[3]);
    }
  }

  N = (int)sqrt(x_coords.size());
  double max_x = *std::max_element(x_coords.begin(), x_coords.end());
  double min_x = *std::min_element(x_coords.begin(), x_coords.end());

  double max_y = *std::max_element(y_coords.begin(), y_coords.end());
  double min_y = *std::min_element(y_coords.begin(), y_coords.end());
  // double dx = abs(x_coords[1] - x_coords[0]);
  double dy = abs(y_coords[1] - y_coords[0]);
  double dx = dy;
  

 //
 // Current and derivatives
//
double m = 0.0;
for (int i = 0; i < N-1; i++){
  for (int j = 0; j < N-1; j++){
    m += 0.5*(x_coords[j + N*i]*jy[j + N*i]- y_coords[j + N*i]*jx[j + N*i])*dx*dy;
  }
}
std::cout << "Total magnetic moment: " << m << std::endl;
  return 0;
}