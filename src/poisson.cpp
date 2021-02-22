#include <iostream>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <string>
#include <vector>
#include <fftw3.h>
#define N 200

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
  double xi = -100;
  double xf = 100;
  std::vector<double> xs = linspace(xi,xf,N);

  double yi = -100;
  double yf =100;
  std::vector<double> ys = linspace(yi,yf,N);

  double zi =-100;
  double zf = 100;
  std::vector<double> zs = linspace(zi,zf,N);

  double L = xf - xi;
  int l;


  std::cout << "L = " << L << std::endl;

  a = 1.0;
  b = 0.0;
  c = 1.0;
  double R = 10.0;

  fftw_plan p1,p2;

  std::string filename = "field.dat";
  std::ofstream fmag;
  fmag.open(filename);

  std::string filename2 = "charge.dat";
  std::ofstream fcharge;
  fcharge.open(filename2);
 
  // std::vector<double> charge_real(N*N*N,0.0);
  // std::vector<double> charge_k(N*N*N,0.0);

  // std::vector<double> field_real(N*N*N,0.0);
  // std::vector<double> field_k(N*N*N,0.0);


  fftw_complex *charge_real;
  fftw_complex *charge_k;
  fftw_complex *field_k;
  fftw_complex *field_real;


  charge_real = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));
  charge_k = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));

  field_k = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));
  field_real = (fftw_complex*) fftw_malloc(N*N*N*sizeof(fftw_complex));

  p1 = fftw_plan_dft_3d(N, N, N, charge_real, charge_k, FFTW_FORWARD, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_3d(N, N, N, field_k, field_real, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
          // charge_real[l][1] = 0.0;
        if (sqrt(pow(xs[i],2) + pow(ys[j],2) + pow(zs[k],2)) < R){

          charge_real[ k + N* (j + N * i)][0] = a*exp(-(pow(xs[i],2) + pow(ys[j],2) + pow(zs[k],2))/(2*pow(c,2)));
          
        }
        else{
          charge_real[ k + N* (j + N * i)][0] = 0.0;
        }
        charge_real[ k + N* (j + N * i)][1] = 0.0;
        
        fcharge <<  xs[i] << "," << ys[j] << "," << zs[k] << "," << charge_real[ k + N* (j + N * i)][0] << std::endl;
      }
    }
  }
  fftw_execute(p1);
  
  //
  // Solve Poisson equation
  // 
  l = 0;
  for (int p = 0; p < N; p++) {
    for (int q = 0; q < N; q++) {
      for (int r = 0; r < N; r++) {
        if (p == 0 && q == 0 && r == 0){
          field_k[r + N* (q + N * p)][0] = 0.0;
          field_k[r + N* (q + N * p)][1] = 0.0;
        } else {
              // field_k[r + N* (q + N * p)][0]= -(2.0*M_PI*r/L)*charge_k[r + N* (q + N * p)][1]/(pow(2*M_PI*p/L,2) + pow(2*M_PI*q/L,2) + pow(2*M_PI*r/L,2));
              // field_k[r + N* (q + N * p)][1]= (2.0*M_PI*r/L)*charge_k[r + N* (q + N * p)][0]/(pow(2*M_PI*p/L,2) + pow(2*M_PI*q/L,2) + pow(2*M_PI*r/L,2));

          field_k[r + N* (q + N * p)][0]= charge_k[r + N* (q + N * p)][0]/(pow(sin(2.0*M_PI*p/L),2) + pow(sin(2.0*M_PI*q/L),2) + pow(sin(2.0*M_PI*r/L),2));
          field_k[r + N* (q + N * p)][1]= charge_k[r + N* (q + N * p)][1]/(pow(sin(2.0*M_PI*p/L),2) + pow(sin(2.0*M_PI*q/L),2) + pow(sin(2.0*M_PI*r/L),2));

        }

        //     field_k[r + N* (q + N * p)][1]= charge_k[r + N* (q + N * p)][1]/fact;
        // l = r + N* (q + N * p);
        // double fact = 0.0;
        // field_k[r + N* (q + N * p)][0] = 0.0;
        // field_k[r + N* (q + N * p)][1] = 0.0;

        // if (2*p < N){
        //   fact = ((double)p*p);
        // } else {
        //   fact = ((double)(N-p)*(N-p));
        // }

        // if (2*q < N){
        //   fact *= ((double)q*q);
        // } else {
        //   fact *= ((double)(N-q)*(N-q));
        // }

        // if (2*r < N){
        //   fact *= ((double)r*r);
        // } else {
        //   fact *= ((double)(N-r)*(N-r));
        // }

        // if(fact!=0) {
        //     field_k[r + N* (q + N * p)][0]= charge_k[r + N* (q + N * p)][0]/fact;
        //     field_k[r + N* (q + N * p)][1]= charge_k[r + N* (q + N * p)][1]/fact;
        // } else {
        //     field_k[r + N* (q + N * p)][0]= 0;
        //     field_k[r + N* (q + N * p)][1]= 0;

        // }
        // std::cout << charge_k[l] << std::endl;
        // field_k[l][0] = (2.0*M_PI*r/L)*charge_k[l][1]/(pow(2.0*M_PI*p/L,2)+pow(2.0*M_PI*q/L,2)+pow(2.0*M_PI*r/L,2));;
        // field_k[l][1] = -(2.0*M_PI*r/L)*charge_k[l][0]/(pow(2.0*M_PI*p/L,2)+pow(2.0*M_PI*q/L,2)+pow(2.0*M_PI*r/L,2));
        // field_k[l][0] = charge_k[l][0]/(pow(2.0*M_PI*p/L,2)+pow(2.0*M_PI*q/L,2)+pow(2.0*M_PI*r/L,2));
        // field_k[l][1] = charge_k[l][1]/(pow(2.0*M_PI*p/L,2)+pow(2.0*M_PI*q/L,2)+pow(2.0*M_PI*r/L,2));
        

        // l++;
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
        l = k + N* (j + N * i);

        field_real[l][0] *= 1.0/pow(N,3);
        fmag << xs[i] << "," << ys[j] << "," << zs[k] << "," << field_real[l][0] << std::endl;

        // field_real[l][1] *= 1.0/pow(N,3);
        // l++;
      }
    }
  }

// Output
  // l = 0;
  // double dx = abs(xs[1] - xs[0]);
  // double dy = abs(ys[1] - ys[0]);
  // double dz = abs(zs[1] - zs[0]);
  // for (int i = 1; i < N-1; i++) {
  //   for (int j = 1; j < N-1; j++) {
  //     for (int k = 1; k < N-1; k++) {

  //       // Central finite differences
  //       double f_dx = field_real[k + N* (j + N * (i+1))];
  //       double f_mdx = field_real[k + N* (j + N * (i-1))];

  //       double f_dy = field_real[k + N* ((j+1) + N * i)];
  //       double f_mdy = field_real[k + N* ((j-1) + N * i)];

  //       double f_dz = field_real[(k+1) + N* (j + N * i)];
  //       double f_mdz = field_real[(k-1) + N* (j + N * i)];

  //       double bx = -(f_dx - f_mdx)/(2.0*dx);
  //       double by = -(f_dy - f_mdy)/(2.0*dy);
  //       double bz = -(f_dz - f_mdz)/(2.0*dz);
  //       fmag << xs[i] << "," << ys[j] << "," << zs[k] << "," << sqrt(pow(bx,2)+pow(by,2)+pow(bz,2)) << std::endl;
  //       // fmag << xs[i] << "," << ys[j] << "," << zs[k] << "," << field_real[l] << std::endl;
  //       l++;
  //     }
  //   }
  // }

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);

  fmag.close();
  fcharge.close();
  // fftw_free(charge_real);
  fftw_free(charge_k);
  // fftw_free(field_real);
  fftw_free(field_k);
 
  return 0;
}