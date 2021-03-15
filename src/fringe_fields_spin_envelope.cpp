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
  std::vector<std::vector<double> > s;
  std::vector<double> dxsx;
  std::vector<double> dysy;
  std::vector<double> sz;

  std::ifstream infile(argv[1]);
  std::ofstream field;
  std::ofstream fdxsx;
  std::ofstream fdysy;

  field.open("fringe_fields_spin_envelope.dat");
  fdxsx.open("spin_derivative_dsx_dx.dat");
  fdysy.open("spin_derivative_dsy_dy.dat");

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
      full_x_coords.push_back(v[0]);
      full_y_coords.push_back(v[1]);
      s.push_back(v);
    }
  }

   N = (int)sqrt(s.size());
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
      dxsx.push_back((s[j+N*(i+1)][2] - s[j+N*i][2])/dx);
      dysy.push_back((s[(j+1)+N*i][3] - s[j+N*i][3])/dy);
      
      sz.push_back(s[j+N*i][4]);
      x_coords.push_back(s[j+N*i][0]);
      y_coords.push_back(s[j+N*i][1]);
    }
  }

  // Field calculation cell
  int Ncell = 100;
  double min_x_cell = -50.0;
  double max_x_cell = 50.0;
  double min_y_cell = -50.0;
  double max_y_cell = 50.0;
  
  std::vector<double> xs_cell = linspace(min_x_cell,max_x_cell,Ncell);
  std::vector<double> ys_cell = linspace(min_y_cell,max_x_cell,Ncell);

  // To calculate the z-component of the B-field 
  std::vector<double> zs_cell;
  double dz = 1e-4;
  double height = 10.0 + d/2.0;  // 10 nm above the slab
  zs_cell.push_back(height);
  zs_cell.push_back(height+dz);

  
  dx = abs(xs_cell[1]-xs_cell[0]);
  dy = abs(ys_cell[1]-ys_cell[0]);


  // X-Y plane  
  std::cout << "Calculating B-field at height " << height << std::endl;
  double z = zs_cell[0];

  for (int i = 0; i < Ncell; i++) {
    for (int j = 0; j < Ncell; j++) {

      double phi = 0.0;
      double phi_dx = 0.0;
      double phi_dy = 0.0;
      double phi_dz = 0.0;
     
      double x = xs_cell[i];
      double y = ys_cell[j];
      double z = zs_cell[0];

      double xdx = xs_cell[i+1];
      double ydy = ys_cell[j+1];
      double zdz = zs_cell[1];

      // Primed variables
      double xp, yp, zp;

      // Slab integral
      for (int ip = 0; ip < N -1; ip++){
        for (int jp = 0; jp < N - 1; jp++){
          for (int kp = 0; kp < Nz; kp++){
            

            xp = x_coords[jp+N*ip];
            yp = y_coords[jp+N*ip];
            zp = z_coords[kp];

            double div = (1.0/d)*cos(M_PI*zp/d)*cos(M_PI*zp/d)*dxsx[jp + N*ip] + (1.0/d)*cos(M_PI*zp/d)*cos(M_PI*zp/d)*dysy[jp + N*ip] - (M_PI/(d*d))*sz[jp + N*ip]*sin(2*M_PI*zp/d);
            phi += div/sqrt(pow(x-xp,2)+pow(y-yp,2)+pow(z-zp,2));


            
            // std::cout <<  dxjy[jp + N*ip]<< std::endl;

            phi_dx += div/sqrt(pow(xdx-xp,2)+pow(y-yp,2)+pow(z-zp,2));
            phi_dy += div/sqrt(pow(x-xp,2)+pow(ydy-yp,2)+pow(z-zp,2));
            phi_dz += div/sqrt(pow(x-xp,2)+pow(y-yp,2)+pow(zdz-zp,2));

          }
        }
      }
      double bx = -(phi_dx - phi)/dx;
      double by = -(phi_dy - phi)/dy;
      double bz = -(phi_dz - phi)/dz;
      field << xs_cell[i] << "," << ys_cell[j] << "," << sqrt(bx*bx + by*by + bz*bz) << std::endl;
      // std::cout << xs_cell[i] << "," << ys_cell[j] << "," << sqrt(bx*bx + by*by + bz*bz) << std::endl;

    }
  }

  field.close();
  fdxsx.close();
  fdysy.close();


  return 0;
}