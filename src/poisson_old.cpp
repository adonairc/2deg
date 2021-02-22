#include <iostream>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <string>
#include <vector>
#include <fftw3.h>
#define N 100

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
  double a,b,c;
  double xi = -10;
  double xf = 10;
  std::vector<double> xs = linspace(xi,xf,N);

  double yi = -10;
  double yf =10;
  std::vector<double> ys = linspace(yi,yf,N);

  double zi =-10;
  double zf = 10;
  std::vector<double> zs = linspace(zi,zf,N);


  double L = xf - xi;
  int l;


  std::cout << "L = " << L << std::endl;

  a = 1.0;
  b = 0.0;
  c = 0.5;

  fftw_plan p1,p2;

  std::string filename = "field.dat";
  std::ofstream fmag;
  fmag.open(filename);

  std::string filename2 = "charge.dat";
  std::ofstream fcharge;
  fcharge.open(filename2);
 
  std::vector<double> charge_real(N*N*N,0.0);
  std::vector<double> charge_k(N*N*N,0.0);

  std::vector<double> field_real(N*N*N,0.0);
  std::vector<double> field_k(N*N*N,0.0);


  p1 = fftw_plan_r2r_3d(N, N, N, charge_real.data(), charge_k.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
  p2 = fftw_plan_r2r_3d(N, N, N, field_k.data(), field_real.data(),FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

  // fftw_complex *charge_real;
  // fftw_complex *charge_k;
  // fftw_complex *field_k;
  // fftw_complex *field_real;


  // charge_real = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));
  // charge_k = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));

  // field_k = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));
  // field_real = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));


  l = 0; //k + N* (j + N * i)
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        charge_real[l] = a*exp(-(pow(xs[i],2) + pow(ys[j],2) + pow(zs[k],2))/(2*pow(c,2)));
        fcharge <<  xs[i] << "," << ys[j] << "," << zs[k] << "," << charge_real[l] << std::endl;
        l++;
      }
    }
  }

  fftw_execute(p1);

  
  //
  // Solve Poisson equation
  // 
  l = 0;
  for (int p = 1; p < N; p++) {
    for (int q = 1; q < N; q++) {
      for (int r = 1; r < N; r++) {
        //std::cout << charge_k[l] << std::endl;
        field_k[l]= charge_k[l]/(pow(2.0*M_PI*p/L,2)+pow(2.0*M_PI*q/L,2)+pow(2.0*M_PI*r/L,2));
        l++;
      }
    }
  }

  // Inverse Fourier transform of magnetization
  fftw_execute(p2);

  // Normalization
  l = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        field_real[l] *= 1.0/pow(N,3);
        l++;
      }
    }
  }

// Output
  l = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        fmag << xs[i] << "," << ys[j] << "," << zs[k] << "," << field_real[l] << std::endl;
        l++;
      }
    }
  }

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);

  fmag.close();
  fcharge.close();
  // fftw_free(charge_real);
  // fftw_free(charge_k);
  // fftw_free(field_real);
  // fftw_free(field_k);
 
  return 0;
}