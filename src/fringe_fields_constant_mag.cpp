#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <string>
#include <vector>

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
  std::vector<std::vector<double> > ms;

  // if (argc < 2) {
  //   printf("Error. Missing orbital magnetic moment file.");
  //   exit(1);
  // }

  std::ifstream infile(argv[1]);
  std::ofstream field;

  field.open("fringe_fields_cons_mag.dat");

  std::vector<double> x_coords;
  std::vector<double> y_coords;

  // std::cout << "[*] Reading 2D magnetization file " << argv[1] << std::endl;
  // double c = 10.0;

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
      
      x_coords.push_back(v[0]);
      y_coords.push_back(v[1]);
      // mocking up
      // v[2] = 1.0;
      // v[2] = 1.0*exp(-(pow(v[0],2)+pow(v[1],2))/(2*pow(c,2)));
      ms.push_back(v);
      // std::cout << v[0] << ","  << v[1] << "," << v[2] << std::endl;
    }
  }

  N = (int)sqrt(ms.size());
  double max_x = *std::max_element(x_coords.begin(), x_coords.end());
  double min_x = *std::min_element(x_coords.begin(), x_coords.end());

  double max_y = *std::max_element(y_coords.begin(), y_coords.end());
  double min_y = *std::min_element(y_coords.begin(), y_coords.end());

  // Test
  // N = 200;

  // std::vector<double> xs2deg = linspace(-500,500,N);


  // std::vector<double> ys2deg = linspace(-500,500,N);

  // for (auto x: xs2deg){
  //   for (auto y: ys2deg){
  //     std::vector<double> v;
  //     v.push_back(x);
  //     v.push_back(y);
  //     v.push_back(1.0*exp(-(pow(v[0],2)+pow(v[1],2))/(2*pow(c,2))));
  //     std::cout << v[0] << ","  << v[1] << "," << v[2] << std::endl;
  //     ms.push_back(v);
  //   }
  // }

  
  

  


  // Super-cell z dimension
  int Ncell = 101;
  double min_x_cell = -5.0;
  double max_x_cell = 5.0;

  double d = 1.0;

  double min_y_cell = -5.0;
  double max_y_cell = 5.0;

  double min_z_cell = 0.0;
  double max_z_cell = 5.0;
  
  std::vector<double> xs_cell = linspace(min_x_cell,max_x_cell,Ncell);
  std::vector<double> ys_cell = linspace(min_y_cell,max_x_cell,Ncell);
  std::vector<double> zs_cell = linspace(min_z_cell,max_z_cell,Ncell);

  double dx = abs(xs_cell[1]-xs_cell[0]);
  double dy = abs(ys_cell[1]-ys_cell[0]);
  double dz = abs(zs_cell[1]-zs_cell[0]);

  // // X-Z plane
  // double y = 0.0;
  // for (int i = 0; i < Ncell-1; i++) {
  //   // for (int j = 0; j < Ncell-1; j++) {
  //     for (int k = 0; k < Ncell-1; k++) {

  //       double phi = 0.0;
  //       double phi_dx = 0.0;
  //       double phi_dz = 0.0;

  //       for (auto m: ms){
  //         //
  //         phi += -(1.0/(4.0*M_PI)*m[2]*dx*dy*(1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(y-m[1],2)+pow(zs_cell[k]+d/2.0,2)) - 1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(y-m[1],2)+pow(zs_cell[k]-d/2.0,2))));
  //         phi_dx += -(1.0/(4.0*M_PI)*m[2]*dx*dy*(1.0/sqrt(pow(xs_cell[i+1]-m[0],2)+pow(y-m[1],2)+pow(zs_cell[k]+d/2.0,2)) - 1.0/sqrt(pow(xs_cell[i+1]-m[0],2)+pow(y-m[1],2)+pow(zs_cell[k]-d/2.0,2))));
  //         phi_dz += -(1.0/(4.0*M_PI)*m[2]*dx*dy*(1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(y-m[1],2)+pow(zs_cell[k+1]+d/2.0,2)) - 1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(y-m[1],2)+pow(zs_cell[k+1]-d/2.0,2))));
  //       }

  //       double bx = -(phi_dx - phi)/dx;
  //       double bz = -(phi_dz - phi)/dz;
  //       field << xs_cell[i] << "," << zs_cell[k] << "," << bx << "," << bz << std::endl;
  //     }
  //   // }
  // }

  // X-Y plane
  double z = 2;
  int k = -1;
  while (zs_cell[k+1] < z){
    k++;
  }
  
  std::cout << "z = " << zs_cell[k] << std::endl;
  for (int i = 0; i < Ncell-1; i++) {
    for (int j = 0; j < Ncell-1; j++) {

        double phi = 0.0;
        double phi_dx = 0.0;
        double phi_dy = 0.0;
        double phi_dz = 0.0;

        for (auto m: ms){
          //
          phi += -(1.0/(4.0*M_PI)*m[2]*dx*dy*(1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(ys_cell[j]-m[1],2)+pow(zs_cell[k]+d/2.0,2)) - 1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(ys_cell[j]-m[1],2)+pow(zs_cell[k]-d/2.0,2))));
          phi_dx += -(1.0/(4.0*M_PI)*m[2]*dx*dy*(1.0/sqrt(pow(xs_cell[i+1]-m[0],2)+pow(ys_cell[j]-m[1],2)+pow(zs_cell[k]+d/2.0,2)) - 1.0/sqrt(pow(xs_cell[i+1]-m[0],2)+pow(ys_cell[j]-m[1],2)+pow(zs_cell[k]-d/2.0,2))));
          phi_dy += -(1.0/(4.0*M_PI)*m[2]*dx*dy*(1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(ys_cell[j+1]-m[1],2)+pow(zs_cell[k]+d/2.0,2)) - 1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(ys_cell[j+1]-m[1],2)+pow(zs_cell[k]-d/2.0,2))));
          phi_dz += -(1.0/(4.0*M_PI)*m[2]*dx*dy*(1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(ys_cell[j]-m[1],2)+pow(zs_cell[k+1]+d/2.0,2)) - 1.0/sqrt(pow(xs_cell[i]-m[0],2)+pow(ys_cell[j]-m[1],2)+pow(zs_cell[k+1]-d/2.0,2))));
        }

        double bx = -(phi_dx - phi)/dx;
        double by = -(phi_dy - phi)/dy;
        double bz = -(phi_dz - phi)/dz;
        field << xs_cell[i] << "," << ys_cell[j] << "," << sqrt(bx*bx + by*by + bz*bz) << std::endl;
      }
  }
  field.close();



  return 0;
}