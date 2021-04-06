#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <string>
#include <vector>
#include<algorithm>
using namespace std;

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
  std::vector<std::vector<double> > input;

  if (argc < 2) {
    printf("Error. Missing current file.");
    exit(1);
  }

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
        // std::cout << substr << std::endl;
        v.push_back(stod(substr));
      }
      input.push_back(v);
      x_coords.push_back(v[0]);
      y_coords.push_back(v[1]);
    }
  }

  N = (int)sqrt(input.size());
  // std::vector<double> jx(N*N,0.0);
  // std::vector<double> jy(N*N,0.0);
  // std::vector<double> m(N*N,0.0);

  complex<double> ii(0.0,1.0);
  vector<complex<double>> jx_k(N*N);
  vector<complex<double>> jy_k(N*N);
  
  vector<complex<double>> mag(N*N);
  vector<complex<double>> mag_k(N*N);
  
  double max_x = *std::max_element(x_coords.begin(), x_coords.end());
  double min_x = *std::min_element(x_coords.begin(), x_coords.end());

  double max_y = *std::max_element(y_coords.begin(), y_coords.end());
  double min_y = *std::min_element(y_coords.begin(), y_coords.end());

  // double max_z = *std::max_element(z_coords.begin(), z_coords.end());
  // double min_z = *std::min_element(z_coords.begin(), z_coords.end());


  double Lx = abs(max_x - min_x);
  double Ly = abs(max_y - min_y);

  for (int ki = 0; ki < N; ki++){
    for (int kj = 0; kj < N; kj++){
      for (int i = 0; i < N; i++ ){
        for (int j = 0; j< N; j++){
          jx_k[kj +  N*ki] += complex<double>(input[j+N*i][2],0.0)*exp(-2.0*ii*M_PI*(double)(ki*i+kj*j)/(double)N);
          jy_k[kj +  N*ki] += complex<double>(input[j+N*i][3],0.0)*exp(-2.0*ii*M_PI*(double)(ki*i+kj*j)/(double)N);
        }
      } 
    }
  }
  for (int ki = 0; ki < N; ki++){
    for (int kj = 0; kj < N; kj++){
      if (ki == 0 && kj ==0){
        mag_k[kj +  N*ki] = complex<double>(0.0,0.0);
      } else{
        mag_k[kj +  N*ki] =ii*((double)kj*jx_k[kj + N*ki] - (double)ki*jy_k[kj + N*ki])/(double)(ki*ki + kj*kj);
      }
    }
  }

  
  for (int i = 0; i < N; i++ ){
    for (int j = 0; j< N; j++){
      for (int ki = 0; ki < N; ki++){
        for (int kj = 0; kj < N; kj++){
          mag[j +  N*i] += (1.0/(double)(N*N))*mag_k[kj+N*ki]*exp(2.0*ii*M_PI*(double)(ki*i+kj*j)/(double)N);
        }
      } 
    }
  }
  for (int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++){
      std::cout << x_coords[j + N*i] << "," <<  y_coords[j + N*i] << ","  << real(mag[j + N*i]) << std::endl;
    }
  }

  
  return 0;
}