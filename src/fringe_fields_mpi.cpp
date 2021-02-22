#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <string>
#include <vector>
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;


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
  mpi::environment env;
  mpi::communicator world;
  
  const int Nx = 100;
  const int Ny = 100;
  const int Nz = 40;


  double xi = -100;
  double xf = 100;
  std::vector<double> xs = linspace(xi,xf,Nx);

  double yi = -100;
  double yf =100;
  std::vector<double> ys = linspace(yi,yf,Ny);

  double zi =-20.0;
  double zf = 20.0;
  std::vector<double> zs = linspace(zi,zf,Nz);

  std::vector<std::vector<double> > ms;


std::vector<std::vector<std::vector<double>>> divs(Nx-1, std::vector<std::vector<double>>(Ny-1, std::vector<double>(Nz-1)));


  // for (int i = 0; i < Nx - 1; i ++)
  //   for (int j = 0; i < Ny - 1; j ++)
  //     for (int k = 0; i < Nz - 1; k ++)
  //       divs[i][j][k] = 0.0;


  if (argc < 2) {
    printf("Error. Missing orbital magnetic moment file.");
    exit(1);
  }

  std::ifstream infile(argv[1]);
  std::vector<double> x_coords;
  std::vector<double> y_coords;
  std::vector<double> z_coords;

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
      ms.push_back(v);
      x_coords.push_back(v[0]);
      y_coords.push_back(v[1]);
      z_coords.push_back(v[2]);
    }
  }
  
  double max_x = *std::max_element(x_coords.begin(), x_coords.end());
  double min_x = *std::min_element(x_coords.begin(), x_coords.end());

  double max_y = *std::max_element(y_coords.begin(), y_coords.end());
  double min_y = *std::min_element(y_coords.begin(), y_coords.end());

  double max_z = *std::max_element(z_coords.begin(), z_coords.end());
  double min_z = *std::min_element(z_coords.begin(), z_coords.end());
  std::cout << "Total grid points from file: " << ms.size() << std::endl;
  std::cout << "Field-of-view available from the magnetization file:" << std::endl;
  std::cout << "x = [" << min_x  << "," << max_x << "]" << std::endl;
  std::cout << "y = [" << min_y  << "," << max_y << "]" << std::endl;
  std::cout << "z = [" << min_z  << "," << max_z << "]" << std::endl;

  double size = x_coords.size();

  // Calculate the magnetic scalar potential function
  double dx = std::abs(xs[1] - xs[0]);
  double dy = std::abs(ys[1] - ys[0]);
  double dz = std::abs(zs[1] - zs[0]);

  double dV = dx*dy*dz;


  // Perpendicular view
  // std::string filename = "bfield_yz.dat";
  // std::ofstream fmag;
  // fmag.open(filename);
  // std::cout << "Calculating fringe-field perpendicular view in Y-Z plane ..." << std::endl;
  // double x = 0.0;
  // for (double y: ys){
  //   for (double z: zs) {
  //     // x = - y;
  //     double phi = 0.0;
  //     double phi_diff_y = 0.0;
  //     double phi_diff_z = 0.0;
      
  //     int sign_y = 1;
  //     int sign_z = 1;


  //     // Finite difference points
  //     // if (x != xf ) {
  //     //   x_diff = x + dx;
  //     // } else {
  //     //   x_diff = x - dx;
  //     //   sign_x = -1;
  //     // }

  //     if (y != yf ) {
  //       y_diff = y + dy;
  //     } else {
  //       y_diff = y - dy;
  //       sign_y = -1;
  //     }

  //     if (z != zf ) {
  //       z_diff = z + dz;
  //     } else {
  //       z_diff = z - dz;
  //       sign_z = -1;
  //     }
  //     for (auto m: ms) {
  //       double r = sqrt(pow(x-m[0],2)+pow(y-m[1],2)+pow(z,2));
  //       double r_diff_y = sqrt(pow(x-m[0],2)+pow(y_diff-m[1],2)+pow(z,2));
  //       double r_diff_z = sqrt(pow(x-m[0],2)+pow(y-m[1],2)+pow(z_diff,2));
  //       phi += (1.0/4.0*M_PI)*(m[2]*da/r);
  //       phi_diff_y += (1.0/4.0*M_PI)*(m[2]*da/r_diff_y);
  //       phi_diff_z += (1.0/4.0*M_PI)*(m[2]*da/r_diff_z);
  //     }

  //     double by = -sign_y*(phi_diff_y -  phi)/2*dy;
  //     double bz = -sign_z*(phi_diff_z -  phi)/2*dz;
  //     fmag  << y << "," << z << "," << by << "," << bz << std::endl;
  //   }
  // }

  // Top-view
  // std::string filename = "bfield_xy.dat";
  // std::ofstream fmag;
  // fmag.open(filename);
  // double z = 20;
  // std::cout << "Calculating fringe-field magnitude in X-Y plane at z = " << z << std::endl;
  // for (double x: xs){
  //   for (double y: ys) {
  //     double phi = 0.0;
  //     double phi_diff_x = 0.0;
  //     double phi_diff_y = 0.0;
  //     double phi_diff_z = 0.0;
      
  //     int sign_x = 1;
  //     int sign_y = 1;
  //     int sign_z = 1;


  //     // Finite difference points
  //     if (x != xf ) {
  //       x_diff = x + dx;
  //     } else {
  //       x_diff = x - dx;
  //       sign_x = -1;
  //     }

  //     if (y != yf ) {
  //       y_diff = y + dy;
  //     } else {
  //       y_diff = y - dy;
  //       sign_y = -1;
  //     }

  //     if (z != zf ) {
  //       z_diff = z + dz;
  //     } else {
  //       z_diff = z - dz;
  //       sign_z = -1;
  //     }

  //     for (auto m: ms) {
  //       double r = sqrt(pow(x-m[0],2)+pow(y-m[1],2)+pow(z,2));
  //       double r_diff_x = sqrt(pow(x_diff-m[0],2)+pow(y-m[1],2)+pow(z,2));
  //       double r_diff_y = sqrt(pow(x-m[0],2)+pow(y_diff-m[1],2)+pow(z,2));
  //       double r_diff_z = sqrt(pow(x-m[0],2)+pow(y-m[1],2)+pow(z_diff,2));

  //       phi += (1.0/4.0*M_PI)*(m[2]/r);
  //       phi_diff_x += (1.0/4.0*M_PI)*(m[2]/r_diff_x);
  //       phi_diff_y += (1.0/4.0*M_PI)*(m[2]/r_diff_y);
  //       phi_diff_z += (1.0/4.0*M_PI)*(m[2]/r_diff_z);
  //     }

  //     double bx = -(1.0/27.724)*sign_x*(phi_diff_x -  phi)/2*dx;
  //     double by = -(1.0/27.724)*sign_y*(phi_diff_y -  phi)/2*dy;
  //     double bz = -(1.0/27.724)*sign_z*(phi_diff_z -  phi)/2*dz;
  //     double b = sqrt(pow(bx,2)+pow(by,2)+pow(bz,2));
  //     fmag << x << "," << y << "," << b << std::endl;
  //   }
  // }

  // VOLUME
    // Perpendicular view
  std::string filename = "bfield_vol_yz.dat";
  std::ofstream fmag;
  fmag.open(filename);
  std::cout << "Calculating fringe-field perpendicular view in Y-Z plane ..." << std::endl;

 
 std::cout << "[*] Calculating divergence of magnetization..." << std::endl;
 for (int i = 0; i < Nx-1; i++){
    for (int j = 0; j < Ny-1; j++){
      for (int k = 0; k < Nz-1; k++) {
        double m_x = 0.0;
        double m_y = 0.0;
        double m_z = 0.0;

        // Integrals evaluated at the advanced point
        double m_dx = 0.0;
        double m_dy = 0.0;
        double m_dz = 0.0;
        for (auto m: ms) {
          // double r = sqrt(pow(xs[i]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k]-m[2],2));
          // double r_dx = sqrt(pow(xs[i+1]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k]-m[2],2));
          // double r_dy = sqrt(pow(xs[i]-m[0],2)+pow(ys[j+1]-m[1],2)+pow(zs[k]-m[2],2));
          // double r_dz = sqrt(pow(xs[i]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k+1]-m[2],2));

          m_x  +=  dV*m[3]/sqrt(pow(xs[i]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k]-m[2],2));
          m_dx  +=  dV*m[3]/sqrt(pow(xs[i+1]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k]-m[2],2));

          m_y  +=  dV*m[4]/sqrt(pow(xs[i]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k]-m[2],2));
          m_dy  +=  dV*m[4]/sqrt(pow(xs[i]-m[0],2)+pow(ys[j+1]-m[1],2)+pow(zs[k+1]-m[2],2));

          m_z  +=  dV*m[5]/sqrt(pow(xs[i]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k]-m[2],2));
          m_dz  +=  dV*m[5]/sqrt(pow(xs[i]-m[0],2)+pow(ys[j]-m[1],2)+pow(zs[k+1]-m[2],2));

        }
        divs[i][j][k] = (m_dx - m_x)/dx + (m_dy - m_y)/dy + (m_dz - m_z)/dz;
        // std::cout << "divs = " << divs[i][j][k] << std::endl;
      }
    }
  }


  std::cout << "[*] Calculating gradient of the potential..." << std::endl;

  int i = floor(Nx/2.0);
  
  for (int j = 0; j < Ny-2; j++){
    for (int k = 0; k < Nz-2; k++) {
      double div = divs[i][j][k];
      double div_dx = divs[i+1][j][k];
      double div_dy = divs[i][j+1][k];
      double div_dz = divs[i][j][k+1];
     
      fmag  << ys[j] << "," << zs[k] << "," << (div_dx - div)/dx  << "," << (div_dy - div)/dy  << "," << (div_dz - div)/dz  << std::endl;;
    }
  }


  std::cout << "Finished. Magnetic field written at " << filename << std::endl;

  fmag.close();
  return 0;
}