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

  int N,Nz;
  std::vector<std::vector<double> > js;
  std::vector<double> dxjy;
  std::vector<double> dyjx;
  std::vector<double> jx;
  std::vector<double> jy;

  std::ifstream infile(argv[1]);
  std::ofstream field;
  std::ofstream dxjy_file;
  std::ofstream dyjx_file;

  field.open("fringe_field_line.dat");
  dxjy_file.open("derivative_djy_dx.dat");
  dyjx_file.open("derivative_djx_dy.dat");
  std::vector<double> x_coords;
  std::vector<double> y_coords;

  std::vector<double> full_x_coords;
  std::vector<double> full_y_coords;

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
      js.push_back(v);
      full_x_coords.push_back(v[0]);
      full_y_coords.push_back(v[1]);
    }
  }

  N = (int)sqrt(js.size());
  double max_x = *std::max_element(full_x_coords.begin(), full_x_coords.end());
  double min_x = *std::min_element(full_x_coords.begin(), full_x_coords.end());

  double max_y = *std::max_element(full_y_coords.begin(), full_y_coords.end());
  double min_y = *std::min_element(full_y_coords.begin(), full_y_coords.end());
  // double dx = abs(x_coords[1] - x_coords[0]);
  double dy = abs(full_y_coords[1] - full_y_coords[0]);
  double dx =dy;
  
  // Slab
  Nz = 20;
  double d = 10.0;
  std::vector<double> z_coords = linspace(-d/2.0,d/2.0,Nz);

 //
 // Current and derivatives
 for (int i = 0; i < N-1; i++){
    for (int j = 0; j < N-1; j++){
      dxjy.push_back((js[j+N*(i+1)][3] - js[j+N*i][3])/dx);
      dyjx.push_back((js[(j+1)+N*i][2] - js[j+N*i][2])/dy);
      jx.push_back(js[j+N*i][2]);
      jy.push_back(js[j+N*i][3]);
      x_coords.push_back(js[j+N*i][0]);
      y_coords.push_back(js[j+N*i][1]);
      dxjy_file << js[j+N*i][0] << "," <<js[j+N*i][1] <<"," << (js[j+N*(i+1)][3] - js[j+N*i][3])/dx << std::endl;
      dyjx_file << js[j+N*i][0] << "," <<js[j+N*i][1] <<"," << (js[(j+1)+N*i][2] - js[j+N*i][2])/dy << std::endl;
    }
  }
  dxjy_file.close();
  dyjx_file.close();

  // double max_x = *std::max_element(x_coords.begin(), x_coords.end());
  // double min_x = *std::min_element(x_coords.begin(), x_coords.end());

  // double max_y = *std::max_element(y_coords.begin(), y_coords.end());
  // double min_y = *std::min_element(y_coords.begin(), y_coords.end());
  // // double dx = abs(x_coords[1] - x_coords[0]);
  // double dy = abs(y_coords[1] - y_coords[0]);
  // double dx =dy;
  std::cout << "dx = " << dx << std::endl;
  std::cout << "dy = " << dy << std::endl;



  // Field calculation cell
  int Ncell = 100;
  double min_x_cell = -50.0;
  double max_x_cell = 50.0;
  double min_y_cell = -50.0;
  double max_y_cell = 50.0;
  
  std::vector<double> xs_cell = linspace(min_x_cell,max_x_cell,Ncell);
  std::vector<double> ys_cell = linspace(max_y_cell,min_x_cell,Ncell);
  // std::vector<double> ys_cell = linspace(min_y_cell,max_x_cell,Ncell);

 

  
  dx = abs(xs_cell[1]-xs_cell[0]);
  dy = abs(ys_cell[1]-ys_cell[0]);

   // To calculate the z-component of the B-field 
  std::vector<double> zs_cell;
  // double dz = 1e-4;
  double dz = dx; // Same spacing as the data
  double height = 10.0 + d/2.0;  // 10 nm above the slab
  zs_cell.push_back(height);
  zs_cell.push_back(height+dz);


  // X-Y plane  
  std::cout << "Calculating B-field at height " << height << std::endl;
  double displacement = 4.0;
  double z = zs_cell[0];

  for (int i = 0; i < Ncell; i++) {
    // for (int j = 0; j < Ncell; j++) {

      double phi = 0.0;
      double phi_dx = 0.0;
      double phi_dy = 0.0;
      double phi_dz = 0.0;
     
      double x = xs_cell[i];
      double y = ys_cell[i];

      double xdx = xs_cell[i+1];
      double ydy = ys_cell[i+1];
      double zdz = zs_cell[1];

      // double bx_biot = 0.0;
      // double by_biot = 0.0;
      // double bz_biot = 0.0;

      // Primed variables
      double xp, yp, zp;

      // Slab integral
      for (int ip = 0; ip < N -1; ip++){
        for (int jp = 0; jp < N - 1; jp++){
          for (int kp = 0; kp < Nz; kp++){
            double z = zs_cell[0];

            xp = x_coords[jp+N*ip];
            yp = y_coords[jp+N*ip];
            zp = z_coords[kp];

            double div = -(1.0/d)*cos(M_PI*zp/d)*cos(M_PI*zp/d)*((zp-displacement) * dxjy[jp + N*ip]) + (1.0/d)*cos(M_PI*zp/d)*cos(M_PI*zp/d)*((zp-displacement) * dyjx[jp + N*ip]) - (2*M_PI/(d*d))*(xp*jy[jp + N*ip]-yp*jx[jp + N*ip])*sin(M_PI*zp/d)*cos(M_PI*zp/d);
            // double div = - (2*M_PI/(d*d))*(xp*jy[jp + N*ip]-yp*jx[jp + N*ip])*sin(M_PI*zp/d)*cos(M_PI*zp/d);
            phi += div/sqrt(pow(x-xp,2)+pow(y-yp,2)+pow(z-zp,2));
            // std::cout <<  dxjy[jp + N*ip]<< std::endl;

            phi_dx += div/sqrt(pow(xdx-xp,2)+pow(y-yp,2)+pow(z-zp,2));
            phi_dy += div/sqrt(pow(x-xp,2)+pow(ydy-yp,2)+pow(z-zp,2));
            phi_dz += div/sqrt(pow(x-xp,2)+pow(y-yp,2)+pow(zdz-zp,2));

            // // Biot-Savart Law
            // bx_biot += (height-d/2)*jy[jp + N*ip]/sqrt(pow(pow(x-xp,2)+pow(y-yp,2)+pow((height-d/2),2),3));
            // by_biot += -(height-d/2)*jx[jp + N*ip]/sqrt(pow(pow(x-xp,2)+pow(y-yp,2)+pow((height-d/2),2),3));
            // bz_biot += ((y-yp)*jx[jp + N*ip]-(x-xp)*jy[jp + N*ip])/sqrt(pow(pow(x-xp,2)+pow(y-yp,2)+pow((height-d/2),2),3));

          }
        }
      }
      double bx = -(4*M_PI*1e2)*(phi_dx - phi)/dx;
      double by = -(4*M_PI*1e2)*(phi_dy - phi)/dy;
      double bz = -(4*M_PI*1e2)*(phi_dz - phi)/dz;
      field << sqrt(2.0)*xs_cell[i] << "," << bx << "," << by << "," <<bz << std::endl;
      // field << xs_cell[i] << "," << ys_cell[j] << "," << sqrt(bx_biot*bx_biot + by_biot*by_biot + bz_biot*bz_biot) << std::endl;
      // std::cout << xs_cell[i] << "," << ys_cell[j] << "," << sqrt(bx*bx + by*by + bz*bz) << std::endl;

    // }
  }

  field.close();


  return 0;
}