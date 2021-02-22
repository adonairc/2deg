#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <string>

#define PI 3.14159265

// Model parameters
double k0 = 1.0;
double ke = 1.2;
double tau = 1.31;
double t_up = 0.0;
double t_dn = -(1.0/(2*PI));


typedef struct {
  double r;
  double theta_r;
} params_t;

typedef struct{
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
double control = 1e-5;

void error_handler (const char * reason, const char * file, int line, int gsl_errno) {
  std::cout << "Error (" << gsl_errno << ") : " << reason << ", file = " << file  << std::endl;
}

double radial_up_dn_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(qplus*rho);
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(qminus*rho);

  double S = qplus*(sin(theta_q+tau)*R2plus - (PI/2.0)*cos(theta_q-tau)*sin(qplus*rho)) - qminus*(sin(theta_q+tau)*R2minus-(PI/2.0)*cos(theta_q-tau)*sin(qminus*rho));
  return(S/(Q*sqrt(f)));
}

double radial_up_dn_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);

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


double radial_dn_up_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);

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

double radial_dn_up_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);

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

  //std::cout << "r = "<<r << ", theta_r = "<< theta_r << std::endl;
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

  results_t res { result_diag_real, result_diag_imag, result_up_dn_real, result_up_dn_imag, result_dn_up_real, result_dn_up_imag };
  return res;
}


int main (int argc, char* argv[])
{
  
  int N = 20;
  const std::complex<double> ii(0, 1);
  // gsl_set_error_handler_off();
  gsl_set_error_handler (&error_handler);

  double xi = -10;
  double xf = 10;
  double dx = (xf-xi+1)/N;

  double yi = -10;
  double yf = 10;
  double dy = (yf-yi+1)/N;

  //double taus[] = {0.785, 0.871, 0.971, 1.15, 1.31, -1.31};
  double taus[] = {PI/4, PI/2 };

  std::cout << "Calculating Green's functions" << std::endl;


  double r = 0.0;
  double theta_r = 0.0;
  params_t my_params = {r,theta_r};

  

  for(int t = 0; t < 2; t++){
    tau = taus[t];
    std::string filename_fint = "intensity_tau_"+std::to_string(tau)+".dat";
    std::string filename_fpol = "polarization_tau_"+std::to_string(tau)+".dat";
    std::string filename_fldos_imp = "ldos_imp_tau_"+std::to_string(tau)+".dat";

    std::string filename_flineldos_imp = "ldos_line_imp_tau_"+std::to_string(tau)+".dat";
    std::string filename_flinepol_imp_sx = "polarization_line_imp_sx_tau_"+std::to_string(tau)+".dat";
    std::string filename_flinepol_imp_sy = "polarization_line_imp_sy_tau_"+std::to_string(tau)+".dat";
    std::string filename_flinepol_imp_sz = "polarization_line_imp_sz_tau_"+std::to_string(tau)+".dat";

    std::string filename_fpol_imp_sz = "polarization_imp_sz_tau_"+std::to_string(tau)+".dat";
    std::string filename_fpol_imp_sy = "polarization_imp_sy_tau_"+std::to_string(tau)+".dat";
    std::string filename_fpol_imp_sx = "polarization_imp_sx_tau_"+std::to_string(tau)+".dat";

    std::ofstream fint;
    std::ofstream fpol;
    std::ofstream fldos;

    std::ofstream flineldos;

    std::ofstream flinepol_imp_sx;
    std::ofstream flinepol_imp_sy;
    std::ofstream flinepol_imp_sz;

    std::ofstream fpol_imp_sx;    
    std::ofstream fpol_imp_sy;
    std::ofstream fpol_imp_sz;

    fint.open(filename_fint);
    fpol.open(filename_fpol);
    fldos.open(filename_fldos_imp);
    flineldos.open(filename_flineldos_imp);

    flinepol_imp_sx.open(filename_flinepol_imp_sx);
    flinepol_imp_sy.open(filename_flinepol_imp_sy);
    flinepol_imp_sz.open(filename_flinepol_imp_sz);

    fpol_imp_sz.open(filename_fpol_imp_sz);
    fpol_imp_sy.open(filename_fpol_imp_sy);
    fpol_imp_sx.open(filename_fpol_imp_sx);

    results_t results = computeGreensFunction(my_params);
    std::complex<double> G_diag (results.diag_real , results.diag_imag);
    double homo_dos = (-2.0/PI)*std::imag(G_diag);
    

    std::cout << "##### PARAMETERS ######" << std::endl;
    std::cout << "ke = " << ke << std::endl;
    std::cout << "tau = " << tau <<  std::endl;
    std::cout << "t_up = " << t_up  << std::endl;
    std::cout << "t_dn = " << t_dn  << std::endl;
    std::cout << "Homogeneous density = " << homo_dos << std::endl;
    std::cout << "Calculating inhomogeneous Green's functions..." << std::endl;

    // Line LDOS plots
    // for (int i = 0; i < 1000 ; i++){
    //   double r = 0.01*i;
    //   double theta_r = 3.0*PI/4.0;
    //   params_t my_params = {r,theta_r};
    //   results_t greens_functions = computeGreensFunction(my_params);
    //   std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
    //   std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
    //   std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
    //   flineldos << r << "," << r*(-1.0/PI)*std::imag((t_up+t_dn)*(pow(G_diag,2) - G_dn_up*G_up_dn)) << std::endl;
    //   flinepol_imp_sx << r << "," << (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up)) << std::endl;
    //   flinepol_imp_sy << r  << "," << (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_dn-t_up)*(G_diag*G_up_dn + G_diag*G_dn_up)) << std::endl;
    //   flinepol_imp_sz << r << "," << (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_diag + G_dn_up*G_up_dn)) << std::endl;

    // }
    for(int i = 0; i < N; i++){
      for(int j = 0; j < N; j++){
        double x = xi + dx*i;
        double y = yi + dy*j;

        //std::cout << "x = " << x << ", y = " << y << std::endl;
        double r = sqrt(pow(x,2)+pow(y,2));
        double theta_r = atan2(y,x);
        params_t my_params = {r,theta_r};

        results_t greens_functions = computeGreensFunction(my_params);
        std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
        std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
        std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);

        double intensity = (2*pow(abs(G_diag),2) + pow(abs(G_up_dn),2) + pow(abs(G_dn_up),2));
        double polarization = (pow(abs(G_diag),2) - pow(abs(G_up_dn),2))/(pow(abs(G_diag),2) + pow(abs(G_up_dn),2));

        fint << x << "," << y << "," << r*intensity << std::endl;
        fpol << x << "," << y << "," << polarization << std::endl;
        fldos << x << "," << y << "," << (-1.0/PI)*std::imag((t_up-t_dn)*(pow(G_diag,2) - G_dn_up*G_up_dn)) << std::endl;
        fpol_imp_sz << x << "," << y << "," << (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_diag + G_dn_up*G_up_dn)) << std::endl;
        fpol_imp_sy << x << "," << y << "," << (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_dn-t_up)*(G_diag*G_up_dn + G_diag*G_dn_up)) << std::endl;
        fpol_imp_sx << x << "," << y << "," << (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up)) << std::endl;


         // Spin-orbit current

        
        std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
        std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

        std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
        std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

        double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
        double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 

        std::cout << "J_theta_so = " << J_theta_so << ", J_r_so = " << J_r_so << std::endl;

        //std::cout << "("<< x << "," << y << ") = " << (-1.0/PI)*std::imag((t_up-t_dn)*(pow(G_diag,2))) << ","<< (1.0/PI)*std::imag((t_up-t_dn)*G_dn_up*G_up_dn)<< std::endl;
        //fldos << x << "," << y << "," << (-1.0/PI)*std::imag(t_up*(pow(abs(G_diag),2) + G_dn_up*std::conj(G_up_dn)) + t_dn*(pow(abs(G_diag),2) + G_up_dn*std::conj(G_dn_up))) << std::endl;
      }
    }
    fint.close();
    fpol.close();
    fpol_imp_sz.close();
    fpol_imp_sy.close();
    fpol_imp_sx.close();
    fldos.close();
    flineldos.close();
    flinepol_imp_sx.close();
    flinepol_imp_sy.close();
    flinepol_imp_sz.close();


  }
  gsl_integration_workspace_free (w);
  return 0;
}