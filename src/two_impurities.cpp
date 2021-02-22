#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <string>
#include <vector>
#include <Eigen/Eigen>


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
double control = 1e-4;

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
  
  int N = 200;
  gsl_set_error_handler_off();
  //gsl_set_error_handler (&error_handler);


  double xi = -5.0;
  double xf = 5.0;
  std::vector<double> xs = linspace(xi,xf,N);


  double yi = -5.0;
  double yf = 5.0;
  std::vector<double> ys = linspace(yi,yf,N);


  double ri = 0.01;
  double rf = 5;

  double ti = 0.0;
  double tf = 2*M_PI;
  std::vector<double> thetas = linspace(ti,tf,500);


  Eigen::MatrixXcf Id(2,2);
  Id(0,0) = std::complex<double>(1.0,0);
  Id(1,0) = std::complex<double>(0,0);
  Id(0,1) = std::complex<double>(0,0);
  Id(1,1) = std::complex<double>(1.0,0);

  // std::cout << "thetas size: " << thetas.size() << std::endl;
  // for (double d : thetas)
  //   std::cout << d << " ";
  // std::cout << std::endl;


  // Model parameters
  double ke = 1.2;

  double a = 2.0; // distance between impurities in units of 1/k0

  // Impurity 0
  double t0_up = 0.0;
  double t0_dn = -(1.0/(2*M_PI));;

  Eigen::MatrixXcf t0(2,2);
  t0(0,0) = t0_up;
  t0(1,0) = std::complex<double>(0,0);
  t0(0,1) = std::complex<double>(0,0);
  t0(1,1) = t0_dn;

  double r0 = 1.0;
  double theta_0 = M_PI;
  
  // Impurity 1
  double t1_up = -(1.0/(2*PI));
  double t1_dn = 0.0;


  Eigen::MatrixXcf t1(2,2);
  t1(0,0) = t1_up;
  t1(1,0) = std::complex<double>(0,0);
  t1(0,1) = std::complex<double>(0,0);
  t1(1,1) = t1_dn;


  double r1 = 1.0;
  double theta_1 = 0.0;

  const std::complex<double> ii(0, 1);

  // double taus[] = {PI/4.0, 0.871, 0.971, 1.15, 1.31, -1.31, PI/2.0};
  double taus[] = {PI/4.0,1.31,PI/2.0};

  std::cout << "Calculating LDOS..." << std::endl;

  for(int t = 0; t < 1; t++){
    double tau = taus[t];

    Eigen::MatrixXcf G_r1_r0(2,2);
    Eigen::MatrixXcf G_r0_r1(2,2);
    Eigen::MatrixXcf T00(2,2),T10(2,2),T01(2,2),T11(2,2);


    std::string filename = "two_imp_ldos_tau_"+std::to_string(tau)+".dat";
    std::ofstream fldos;
    fldos.open(filename);

    // Intra-impurity scaterring
    double r01 = std::abs(r0 - r1);
    double theta_01 = theta_0 - theta_1;

    // std::cout << "r01 = " << r01 << std::endl;

    params_t params_01 = {r01,theta_01,ke,tau};
    results_t greens_functions_01 = computeGreensFunction(params_01);
    G_r0_r1(0,0) = std::complex<double>(greens_functions_01.diag_real , greens_functions_01.diag_imag);
    G_r0_r1(1,0) = std::complex<double>(greens_functions_01.up_dn_real , greens_functions_01.up_dn_imag);
    G_r0_r1(0,1) = std::complex<double>(greens_functions_01.dn_up_real , greens_functions_01.dn_up_imag);
    G_r0_r1(1,1) = std::complex<double>(greens_functions_01.diag_real , greens_functions_01.diag_imag);

    double r10 = std::abs(r1 - r0);
    double theta_10 = theta_1 - theta_0;

    // std::cout << "r10 = " << r10 << std::endl;

    params_t params_10 = {r10,theta_10,ke,tau};
    results_t greens_functions_10 = computeGreensFunction(params_10);
    G_r1_r0(0,0) = std::complex<double>(greens_functions_10.diag_real , greens_functions_10.diag_imag);
    G_r1_r0(1,0) = std::complex<double> (greens_functions_10.up_dn_real , greens_functions_10.up_dn_imag);
    G_r1_r0(0,1) = std::complex<double> (greens_functions_10.dn_up_real , greens_functions_10.dn_up_imag);
    G_r1_r0(1,1) = std::complex<double>(greens_functions_10.diag_real , greens_functions_10.diag_imag);

    T00 = t0*(Id - t0*G_r0_r1*t1*G_r1_r0).inverse();
    T10 = t1*G_r1_r0*t0*(Id - t0*G_r0_r1*t1*G_r1_r0).inverse();

    T01 = t0*G_r0_r1*t1*(Id - t1*G_r1_r0*t0*G_r0_r1).inverse();
    T11 = t1*(Id - t1*G_r1_r0*t0*G_r0_r1).inverse();

    for( double x: xs){
      for (double y: ys){
          double r = sqrt(pow(x,2)+pow(y,2));
          double theta = atan2(y,x);

          double r_minus_r0 = std::abs(r - r0);
          double theta_minus_theta_0 = theta - theta_0;
          // std::cout << "r_minus_r0 = " << r_minus_r0 << std::endl;

          Eigen::MatrixXcf G_r_r0(2,2);
          Eigen::MatrixXcf G_r0_r(2,2);
          Eigen::MatrixXcf G_r_r1(2,2);
          Eigen::MatrixXcf G_r1_r(2,2);
          Eigen::MatrixXcf Delta_G(2,2);

          params_t params_0 = {r_minus_r0,theta_minus_theta_0,ke,tau};

          results_t greens_functions_0 = computeGreensFunction(params_0);
          G_r_r0(0,0) = std::complex<double>(greens_functions_0.diag_real , greens_functions_0.diag_imag);
          G_r_r0(1,0) = std::complex<double>(greens_functions_0.up_dn_real , greens_functions_0.up_dn_imag);
          G_r_r0(0,1) = std::complex<double>(greens_functions_0.dn_up_real , greens_functions_0.dn_up_imag);
          G_r_r0(1,1) = std::complex<double>(greens_functions_0.diag_real , greens_functions_0.diag_imag);

          G_r0_r(0,0) = G_r_r0(0,0);
          G_r0_r(1,0) = -G_r_r0(1,0);
          G_r0_r(0,1) = -G_r_r0(0,1);
          G_r0_r(1,1) = G_r_r0(1,1);

          double r_minus_r1 = std::abs(r - r1);
          double theta_minus_theta_1 = theta - theta_1;

          // std::cout << "r_minus_r1 = " << r_minus_r1 << std::endl;

          params_t params_1 = {r_minus_r1,theta_minus_theta_1,ke,tau};

          results_t greens_functions_1 = computeGreensFunction(params_1);
          G_r_r1(0,0) = std::complex<double>(greens_functions_1.diag_real , greens_functions_1.diag_imag);
          G_r_r1(1,0) = std::complex<double>(greens_functions_1.up_dn_real , greens_functions_1.up_dn_imag);
          G_r_r1(0,1) = std::complex<double>(greens_functions_1.dn_up_real , greens_functions_1.dn_up_imag);
          G_r_r1(1,1) = std::complex<double>(greens_functions_1.diag_real , greens_functions_1.diag_imag);

          G_r1_r(0,0) = G_r_r1(0,0);
          G_r1_r(1,0) = -G_r_r1(1,0);
          G_r1_r(0,1) = -G_r_r1(0,1);
          G_r1_r(1,1) = G_r_r1(1,1);

          Delta_G = G_r_r0*T00*G_r0_r + G_r_r0*T01*G_r1_r + G_r_r1*T10*G_r0_r + G_r_r1*T11*G_r1_r;

          double ldos = -(1/M_PI)*std::imag(Delta_G(0,0) + Delta_G(1,1));


          fldos << x << "," << y << "," << ldos << std::endl;
          std::cout << x << "," << y << ", - G00 = " << Delta_G(0,0)  << ", G11 = " << Delta_G(1,1)<< std::endl;
         
        } 
      }
    fldos.close();
  }


  gsl_integration_workspace_free (w);
  return 0;
}