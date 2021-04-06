#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <string>
#include <vector>
#include <fftw3.h>
#include<algorithm>
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

  fftw_complex *jx;
  fftw_complex *jx_k;


  fftw_complex *jy;
  fftw_complex *jy_k;

  fftw_complex *mag;
  fftw_complex *mag_k;

  jx = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  jx_k = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  jy = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  jy_k = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));

  mag = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  mag_k = (fftw_complex*) fftw_malloc(N*N*sizeof(fftw_complex));
  


  fftw_plan p1,p2, q;
  int i,j;
  p1 = fftw_plan_dft_2d(N,N, jx, jx_k , FFTW_FORWARD, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_2d(N,N, jy, jy_k , FFTW_FORWARD, FFTW_ESTIMATE);
  q = fftw_plan_dft_2d(N,N, mag_k, mag, FFTW_BACKWARD, FFTW_ESTIMATE);


  
  double max_x = *std::max_element(x_coords.begin(), x_coords.end());
  double min_x = *std::min_element(x_coords.begin(), x_coords.end());

  double max_y = *std::max_element(y_coords.begin(), y_coords.end());
  double min_y = *std::min_element(y_coords.begin(), y_coords.end());

  // double max_z = *std::max_element(z_coords.begin(), z_coords.end());
  // double min_z = *std::min_element(z_coords.begin(), z_coords.end());


  double L1 = abs(max_x - min_x);
  double L2 = abs(max_y - min_y);

  double invL1s=1.0/(L1*L1);
  double invL2s=1.0/(L2*L2);



  for (int i = 0; i < N; i++ ){
    for (int j = 0; j< N; j++){
          jx[j +  N*i][0] = input[j+N*i][2];
          jx[j +  N*i][1] = 0.0;
          jy[j +  N*i][0] = input[j+N*i][3];
          jy[j +  N*i][1] = 0.0;
    }
  }


  fftw_execute(p1);
  fftw_execute(p2);

  int l=-1;
  for (int i = 0; i < N-1; i++){   // f = g / ( kx² + ky² )  
    for(int j = 0; j < N-1; j++){
          l=l+1;
          double fact=0;
          double prefact_jx = 0.0;
          double prefact_jy = 0.0;

          // field_k[l][0] =0;
          // field_k[l][1] =0;

          // if(2*i<N){
              fact=((double)i*i)*invL1s;
              prefact_jy = (double)i;
          // }else{
          //     fact=((double)(N-i)*(N-i))*invL1s;
          //     prefact_jy = (double)(N-i);
          // }

          // if(2*j<N){
              fact+=((double)j*j)*invL2s;
              prefact_jx = (double)j;
          // }else{
          //     fact+=((double)(N-j)*(N-j))*invL2s;
          //     prefact_jx = (double)(N-j);

          // }

          if(fact!=0){
              mag_k[l][0] = -(prefact_jx*jx_k[l][1]/fact) + (prefact_jy*jy_k[l][1]/fact);
              mag_k[l][1] = (prefact_jx*jx_k[l][0]/fact) - (prefact_jy*jy_k[l][0]/fact);
          }
          else{
              mag_k[l][0] = 0.0;
              mag_k[l][1] = 0.0;
          }
    }
  }

  fftw_execute(q);
  l=-1;

  for (int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++){
      l=l+1;
      std::cout << x_coords[l] << "," <<  y_coords[l] << ","  << mag[l][0]  << std::endl;
    }
  }

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(q);
  fftw_cleanup();




  
  return 0;
}