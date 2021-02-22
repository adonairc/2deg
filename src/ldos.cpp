#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <string>
#include <vector>
#include <chrono> 


#define PI M_PI
double k0 = 1.0;

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

typedef struct {
  double r;
  double theta_r;
  double ke;
  double tau;
} params_t;

typedef struct{
  //G0
  double diag_real;
  double diag_imag;
  double up_dn_real;
  double up_dn_imag;
  double dn_up_real;
  double dn_up_imag;

} results_t;

// Integration parameters
gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
double error;
double control = 1e-7;

void error_handler (const char * reason, const char * file, int line, int gsl_errno) {
  std::cout << "Error (" << gsl_errno << ") : " << reason << ", file = " << file  << std::endl;
}


// G0
double radial_dn_up_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double ke = (params->ke);
  double tau = (params->tau);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(qplus*rho);
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(qminus*rho);

  double S = qplus*(sin(theta_q+tau)*R2plus - (PI/2.0)*cos(theta_q-tau)*sin(qplus*rho)) - qminus*(sin(theta_q+tau)*R2minus - (PI/2.0)*cos(theta_q-tau)*sin(qminus*rho));
  return(S/(Q*sqrt(f)));
}

double radial_dn_up_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double ke = (params->ke);
  double tau = (params->tau);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));


  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(qplus*rho);
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(qminus*rho);

  double S = qplus*(cos(theta_q-tau)*R2plus + (PI/2.0)*sin(theta_q+tau)*sin(qplus*rho)) - qminus*(cos(theta_q-tau)*R2minus+(PI/2.0)*sin(theta_q+tau)*sin(qminus*rho));
  return(S/(Q*sqrt(f)));
}


double radial_up_dn_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double ke = (params->ke);
  double tau = (params->tau);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(qplus*rho);
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(qminus*rho);

  double S = -qplus*(sin(theta_q+tau)*R2plus + (PI/2.0)*cos(theta_q-tau)*sin(qplus*rho)) + qminus*(sin(theta_q+tau)*R2minus + (PI/2.0)*cos(theta_q-tau)*sin(qminus*rho));
  return(S/(Q*sqrt(f)));
}

double radial_up_dn_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double ke = (params->ke);
  double tau = (params->tau);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(qplus*rho);
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(qminus*rho);

  double S = qplus*(cos(theta_q-tau)*R2plus - (PI/2.0)*sin(theta_q+tau)*sin(qplus*rho)) - qminus*(cos(theta_q-tau)*R2minus - (PI/2.0)*sin(theta_q+tau)*sin(qminus*rho));
  return(S/(Q*sqrt(f)));
}


double radial_diagonal_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double ke = (params->ke);
  double tau = (params->tau);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double Iplus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double Iminus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);
  return(-(qplus*Iplus+qminus*Iminus)/Q);
}

double radial_diagonal_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double ke = (params->ke);
  double tau = (params->tau);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  return(-(PI/2)*(qplus*cos(qplus*rho) + qminus*cos(qminus*rho))/Q);
}

results_t computeGreensFunction(params_t params){
  double r = (params.r);
  double theta_r = (params.theta_r);

  double result_diag_real;
  double result_diag_imag;
  double result_up_dn_real;
  double result_up_dn_imag;
  double result_dn_up_real;
  double result_dn_up_imag;

  // Diagonal
  gsl_function F_diag_real;
  F_diag_real.function = &radial_diagonal_real;
  F_diag_real.params = &params;
  gsl_integration_qags (&F_diag_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_diag_real, &error);

  gsl_function F_diag_imag;
  F_diag_imag.function = &radial_diagonal_imag;
  F_diag_imag.params = &params;
  gsl_integration_qags (&F_diag_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_diag_imag, &error);

  // Up - Down
  gsl_function F_up_dn_real;
  F_up_dn_real.function = &radial_up_dn_real;
  F_up_dn_real.params = &params;
  gsl_integration_qags (&F_up_dn_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_up_dn_real, &error);

  gsl_function F_up_dn_imag;
  F_up_dn_imag.function = &radial_up_dn_imag;
  F_up_dn_imag.params = &params;
  gsl_integration_qags (&F_up_dn_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_up_dn_imag, &error);

  // Down - Up
  gsl_function F_dn_up_real;
  F_dn_up_real.function = &radial_dn_up_real;
  F_dn_up_real.params = &params;
  gsl_integration_qags (&F_dn_up_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_dn_up_real, &error);

  gsl_function F_dn_up_imag;
  F_dn_up_imag.function = &radial_dn_up_imag;
  F_dn_up_imag.params = &params;
  gsl_integration_qags (&F_dn_up_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_dn_up_imag, &error);

  results_t res { result_diag_real, result_diag_imag, result_up_dn_real, result_up_dn_imag, result_dn_up_real, result_dn_up_imag};
  return res;
}


int main (int argc, char* argv[])
{
  
  int N,EnergyPoints;
  double tau,initial_coord,final_coord,kef;
  bool calc_line, calc_density,calc_total;

  calc_line = false;
  calc_density = false;
  calc_total = false;

  if (argc < 14) {
  // Command-line arguments
    std::cerr << "Usage: ldos CALCULATION OPTIONS" << std::endl;
    std::cerr << "CALCULATION:" << std::endl;
    std::cerr << "line - Energy integrated line plot" << std::endl;
    std::cerr << "density - Energy integrated LDOS density plot" << std::endl;
    std::cerr << "OPTIONS:" << std::endl;
    std::cerr << "-N <number of grid points>" << std::endl;
    std::cerr << "-En <number of energy points>" << std::endl;
    std::cerr << "-t <spin-orbit ratio comma-separated list> e.g. 0.4,0.5,0.6" << std::endl;
    std::cerr << "-i <initial coordinate> (in units of lambda_so)" << std::endl;
    std::cerr << "-f <final coordinate> (in units of lambda_so)" << std::endl;
    std::cerr << "-kf <fermi wave-vector> (in units of k0)" << std::endl;
    // std::cerr << "-tu <t up>" << std::endl;
    // std::cerr << "-td <t down>" << std::endl;
    return 1;
  }

  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-N") {
      if (i + 1 < argc) {
        N  = std::stoi(argv[i+1]); 
        i++;
      } else { 
        std::cerr << "-N option requires one argument." << std::endl;
        return 1;
      }  
    } else if (std::string(argv[i]) == "-t") {
      if (i + 1 < argc) {
        tau  = std::stof(argv[i+1]);
        i++;
      } else { 
        std::cerr << "-t option requires one argument." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "-Ne") {
      if (i + 1 < argc) {
        EnergyPoints  = std::stoi(argv[i+1]);
        i++;
      } else { 
        std::cerr << "-Ne option requires one argument." << std::endl;
        return 1;
      }  
    }
     else if (std::string(argv[i]) == "-i") {
      if (i + 1 < argc) {
        initial_coord  = std::stof(argv[i+1])*(2*M_PI); // Units of lambda_so
        i++; 
      } else { 
        std::cerr << "-i option requires one argument." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "-f") {
      if (i + 1 < argc) {
        final_coord  = std::stof(argv[i+1])*(2*M_PI); // Units of lambda_so
        i++; 
      } else { 
        std::cerr << "-f option requires one argument." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "-kf") {
      if (i + 1 < argc) {
        kef = std::stof(argv[i+1]); 
        i++;
      } else { 
        std::cerr << "-kf option requires one argument." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "line") {
        calc_line = true;
        i++;
    }
     else if (std::string(argv[i]) == "density") {
        calc_density = true;
        i++;
    }
    
  }

  gsl_set_error_handler_off();
  // gsl_set_error_handler (&error_handler);

  
  double xi = initial_coord;
  double xf = final_coord;
  std::vector<double> xs = linspace(xi,xf,N);



  double yi = initial_coord;
  double yf = final_coord;
  std::vector<double> ys = linspace(yi,yf,N);

  std::vector<double> rs = linspace(xi,xf,N);

  double ti = 0.0;
  double tf = 2*M_PI;
  std::vector<double> thetas = linspace(ti,tf,N);

  
  double kei = 0.001;
  std::vector<double> kes = linspace(kei,kef,EnergyPoints);

  double t_up = 0.0;
  double t_dn = -(1.0/(2*PI));

  const std::complex<double> ii(0, 1);

  
  std::cout << "\n---------------[Parameters]---------------" << std::endl;
  std::cout << "tau : " << tau << std::endl;
  std::cout << "kef : "<< kef << std::endl;
  std::cout << "grid size :  [" << N << "," << N << "] points"<< std::endl;
  std::cout << "energy integration points : "<< EnergyPoints << std::endl;
  std::cout << "x :  [" << xi << "," << xf << "] (in units of 1/k0)"<< std::endl;
  std::cout << "y :  [" << yi << "," << yf << "] (in units of 1/k0) "<< std::endl;
  std::cout << "------------------------------------------" << std::endl;
  auto start = std::chrono::steady_clock::now(); 

  if(calc_density) {
    std::string filename = "ldos_energy_integrated_ke_"+std::to_string(kef)+"_tau_"+std::to_string(tau)+".dat";
    std::ofstream fldos;
    fldos.open(filename);
    for( double x: xs){
      for (double y: ys){

      double r = sqrt(pow(x,2)+pow(y,2));
      double theta_r = atan2(y,x);
      double ldos = 0.0;
      

      for (double ke: kes){
        params_t params = {r,theta_r,ke,tau};
        results_t greens_functions = computeGreensFunction(params);
        std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
        std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
        std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);

        ldos +=  (-1.0/PI)*std::imag((t_up-t_dn)*(pow(G_diag,2) - G_dn_up*G_up_dn));
      }
      fldos << x << "," << y << "," << ldos << std::endl;
      }  
    }
    fldos.close();
  }

  gsl_integration_workspace_free (w);
  return 0;
}