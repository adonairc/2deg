#include <iostream>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <string>
#include <vector>
#include <fftw3.h>

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
  // const int Nx = 100;
  // const int Ny = 100;
  // const int Nz = 10;


  // double xi = -1000;
  // double xf = 1000;
  // std::vector<double> xs = linspace(xi,xf,Nx);

  // double yi = -1000;
  // double yf =1000;
  // std::vector<double> ys = linspace(yi,yf,Ny);

  // double zi =-20.0;
  // double zf = 20.0;
  // std::vector<double> zs = linspace(zi,zf,Nz);


  if (argc < 2) {
    printf("Error. Missing current file.");
    exit(1);
  }

  fftw_plan plan_x, plan_y,plan_m;

  std::string filename = "magmom_fft.dat";
  std::ofstream fmag;
  fmag.open(filename);
  
  std::ifstream infile(argv[1]);
  std::vector<double> x_coords;
  std::vector<double> y_coords;
  std::vector<double> jxs;
  std::vector<double> jys;

  fftw_complex *j_x;
  fftw_complex *j_y;
  fftw_complex *j_kx;
  fftw_complex *j_ky;

  fftw_complex *m;
  fftw_complex *m_k;



  std::cout << "[*] Reading file " << argv[1] << std::endl;

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
      jxs.push_back(v[2]*cos(v[0]) - v[3]*sin(v[1]));
      jys.push_back(v[2]*sin(v[0]) + v[3]*cos(v[1]));
      y_coords.push_back(v[0]*cos(v[1]));
      x_coords.push_back(v[0]*sin(v[1]));
    }
  }

  int N = sqrt(jxs.size());

  
  double max_x = *std::max_element(x_coords.begin(), x_coords.end());
  double min_x = *std::min_element(x_coords.begin(), x_coords.end());

  double max_y = *std::max_element(y_coords.begin(), y_coords.end());
  double min_y = *std::min_element(y_coords.begin(), y_coords.end());

  std::cout << "Total grid points from file: " << jxs.size() << std::endl;
  std::cout << "Field-of-view available from the current file:" << std::endl;
  std::cout << "x = [" << min_x  << "," << max_x << "]" << std::endl;
  std::cout << "y = [" << min_y  << "," << max_y << "]" << std::endl;
  std::cout << "N = " << N   << std::endl;

  j_x = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  j_y = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));

  j_kx = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  j_ky = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));

  m = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  m_k = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      j_x[j + N*(i)][0] = jxs[i + N*j];
      j_x[j + N*(i)][1] = 0.0;

      j_y[j + N*(i)][0] = jys[i + N*j];
      j_y[j + N*(i)][1] = 0.0;
      std::cout << x_coords[i + N*j] << "," << y_coords[i + N*(j)] << "," << jxs[i + N*(j)] << "," << jys[i + N*(j)] <<std::endl;

    }
  }

  plan_x = fftw_plan_dft_2d(N, N, j_x, j_kx, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan_x);

  plan_y = fftw_plan_dft_2d(N, N, j_y, j_ky, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan_y);
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     std::cout << "[" << j_kx[j + N*i][0]<< "," << j_kx[j + N*i][1] << "] [" << j_ky[j + N*i][0]<< "," << j_ky[j + N*i][1] << "]" << std::endl;
  //   }
  // }
  
  //
  // Solve for magnetization in k-space
  // 
  for (int i = 1; i < N; i++) {
    for (int j = 1; j < N; j++) {

      std::complex<double> jkx(j_kx[j + N*(i)][0],j_kx[j + N*(i)][1]);
      std::complex<double> jky(j_ky[j + N*(i)][0],j_ky[j + N*(i)][1]);
      std::complex<double> mk =   (std::complex<double>(0,2*M_PI*j/N)*jkx - std::complex<double>(0,2*M_PI*i/N)*jky)/(pow((2*M_PI*i/N),2)+ pow((2*M_PI*j/N),2));
      m_k[j + N*(i)][0] = std::real(mk);
      m_k[j + N*(i)][1] = std::imag(mk);
      // std::cout << mk << std::endl;
    }
  }
  // Inverse Fourier transform of magnetization
  plan_m = fftw_plan_dft_2d(N, N, m_k, m, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan_m);

  // Normalization
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      m[j + N*(i)][0] *= 1.0/N;
      m[j + N*(i)][1] *= 1.0/N;
      // std::cout << m[j + N*(i)][0] << std::endl;
    }
  }

// Output
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
     fmag << x_coords[i + N*j] << "," << y_coords[i + N*j] << "," << m[j + N*(i)][0] << std::endl;
    }
  }


  fmag.close();
  fftw_destroy_plan(plan_x);
  fftw_destroy_plan(plan_y);
  fftw_destroy_plan(plan_m);
  fftw_free(j_x);
  fftw_free(j_y);
  fftw_free(j_kx);
  fftw_free(j_ky);
 
  return 0;
}