#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <boost/program_options.hpp>
#include <string>
#include <vector>
#include <chrono> 

#define PI M_PI
double k0 = 0.01241; // GaAs (beta ~ 10 meV.nm)
namespace po = boost::program_options;


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
  // Gradients
  double grad_r_diag_real;
  double grad_r_diag_imag;
  double grad_r_up_dn_real;
  double grad_r_up_dn_imag;
  double grad_r_dn_up_real;
  double grad_r_dn_up_imag;

  double grad_theta_diag_real;
  double grad_theta_diag_imag;
  double grad_theta_up_dn_real;
  double grad_theta_up_dn_imag;
  double grad_theta_dn_up_real;
  double grad_theta_dn_up_imag;

} results_t;



void error_handler (const char * reason, const char * file, int line, int gsl_errno) {
  if (gsl_errno == 22){
    std::cout << "Error (" << gsl_errno << ") : " << reason << ", file = " << file  << ", line = " << line << std::endl;
  }
}

// Radial gradients
double grad_r_diag_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double ke = (params->ke);
  double tau = (params->tau);

  double rho = abs(r*cos(theta_q-theta_r));
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f) ;
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(qplus*rho);
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(qminus*rho);


  double T1 = (qplus*qplus*R2plus + qminus*qminus*R2minus)*cos(theta_q-theta_r);
  double T2 = (qminus + qplus)/r;

  return((1.0/(2*M_PI*M_PI))*(T1+T2)/Q);
}

double grad_r_diag_imag(double theta_q, void * p) {
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


  double T1 = (PI/2.0)*(qplus*qplus*sin(qplus*rho) + qminus*qminus*sin(qminus*rho))*cos(theta_q-theta_r);

  return((1.0/(2*M_PI*M_PI))*T1/Q);
}


double grad_r_dn_up_real(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_r_dn_up_imag(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = (PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_r_up_dn_real(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = -sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_r_up_dn_imag(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = -(PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));


  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}


// Angular gradients

double grad_theta_diag_real(double theta_q, void * p) {
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


  double T1 = (qplus*qplus*R2plus + qminus*qminus*R2minus)*sin(theta_q-theta_r);
  double T2 = -(qminus + qplus)*tan(theta_q-theta_r)/r;

  return(-(1.0/(2*M_PI*M_PI))*(T1+T2)/Q);
}

double grad_theta_diag_imag(double theta_q, void * p) {
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

  double T1 = (PI/2.0)*(qplus*qplus*sin(qplus*rho) + qminus*qminus*sin(qplus*rho))*sin(theta_q-tau);
  return(-(1.0/(2*M_PI*M_PI))*T1/Q);
}


double grad_theta_dn_up_real(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_theta_dn_up_imag(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = (PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_theta_up_dn_real(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = -sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_theta_up_dn_imag(double theta_q, void * p) {
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

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(qplus*rho)-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(qminus*rho)-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = -(PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
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
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
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
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
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
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
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
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
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
  // std::cout << "radial_diagonal_real() Q = " << Q << std::endl;

  return(-(1.0/(2*M_PI*M_PI))*(qplus*Iplus+qminus*Iminus)/Q);

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
  // std::cout << "radial_diagonal_imag() Q = " << Q << std::endl;

  return(-(1.0/(2*M_PI*M_PI))*(PI/2)*(qplus*cos(qplus*rho) + qminus*cos(qminus*rho))/Q);

}

results_t computeGreensFunction(params_t params){
  // Integration parameters
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double error;
  double control = 1e-6;

  double r = (params.r);
  double theta_r = (params.theta_r);

  double result_diag_real;
  double result_diag_imag;
  double result_up_dn_real;
  double result_up_dn_imag;
  double result_dn_up_real;
  double result_dn_up_imag;

  double result_grad_r_diag_real;
  double result_grad_r_diag_imag;
  double result_grad_r_up_dn_real;
  double result_grad_r_up_dn_imag;
  double result_grad_r_dn_up_real;
  double result_grad_r_dn_up_imag;

  double result_grad_theta_diag_real;
  double result_grad_theta_diag_imag;
  double result_grad_theta_up_dn_real;
  double result_grad_theta_up_dn_imag;
  double result_grad_theta_dn_up_real;
  double result_grad_theta_dn_up_imag;

  // std::cout << "computeGreensFunction() r = "<<r << ", theta_r = "<< theta_r << ", ke = "<< params.ke << ", tau = "<< params.tau << std::endl;
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



  // Radial Gradients

  gsl_function Grad_r_diag_real;
  Grad_r_diag_real.function = &grad_r_diag_real;
  Grad_r_diag_real.params = &params;
  gsl_integration_qags (&Grad_r_diag_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_diag_real, &error);

  gsl_function Grad_r_diag_imag;
  Grad_r_diag_imag.function = &grad_r_diag_imag;
  Grad_r_diag_imag.params = &params;
  gsl_integration_qags (&Grad_r_diag_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_diag_imag, &error);

  gsl_function Grad_r_up_dn_real;
  Grad_r_up_dn_real.function = &grad_r_up_dn_real;
  Grad_r_up_dn_real.params = &params;
  gsl_integration_qags (&Grad_r_up_dn_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_up_dn_real, &error);

  gsl_function Grad_r_up_dn_imag;
  Grad_r_up_dn_imag.function = &grad_r_up_dn_imag;
  Grad_r_up_dn_imag.params = &params;
  gsl_integration_qags (&Grad_r_up_dn_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_up_dn_imag, &error);

  gsl_function Grad_r_dn_up_real;
  Grad_r_dn_up_real.function = &grad_r_dn_up_real;
  Grad_r_dn_up_real.params = &params;
  gsl_integration_qags (&Grad_r_dn_up_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_dn_up_real, &error);

  gsl_function Grad_r_dn_up_imag;
  Grad_r_dn_up_imag.function = &grad_r_dn_up_imag;
  Grad_r_dn_up_imag.params = &params;
  gsl_integration_qags (&Grad_r_dn_up_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_dn_up_imag, &error);
   

  // Angular Gradients

  gsl_function Grad_theta_diag_real;
  Grad_theta_diag_real.function = &grad_theta_diag_real;
  Grad_theta_diag_real.params = &params;
  gsl_integration_qags (&Grad_theta_diag_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_diag_real, &error);

  gsl_function Grad_theta_diag_imag;
  Grad_theta_diag_imag.function = &grad_theta_diag_imag;
  Grad_theta_diag_imag.params = &params;
  gsl_integration_qags (&Grad_theta_diag_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_diag_imag, &error);

  gsl_function Grad_theta_up_dn_real;
  Grad_theta_up_dn_real.function = &grad_theta_up_dn_real;
  Grad_theta_up_dn_real.params = &params;
  gsl_integration_qags (&Grad_theta_up_dn_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_up_dn_real, &error);

  gsl_function Grad_theta_up_dn_imag;
  Grad_theta_up_dn_imag.function = &grad_theta_up_dn_imag;
  Grad_theta_up_dn_imag.params = &params;
  gsl_integration_qags (&Grad_theta_up_dn_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_up_dn_imag, &error);

  gsl_function Grad_theta_dn_up_real;
  Grad_theta_dn_up_real.function = &grad_theta_dn_up_real;
  Grad_theta_dn_up_real.params = &params;
  gsl_integration_qags (&Grad_theta_dn_up_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_dn_up_real, &error);

  gsl_function Grad_theta_dn_up_imag;
  Grad_theta_dn_up_imag.function = &grad_theta_dn_up_imag;
  Grad_theta_dn_up_imag.params = &params;
  gsl_integration_qags (&Grad_theta_dn_up_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_dn_up_imag, &error);

  results_t res { result_diag_real, result_diag_imag, result_up_dn_real, result_up_dn_imag, result_dn_up_real, result_dn_up_imag, result_grad_r_diag_real,result_grad_r_diag_imag, result_grad_r_up_dn_real,result_grad_r_up_dn_imag,result_grad_r_dn_up_real,result_grad_r_dn_up_imag,result_grad_theta_diag_real,result_grad_theta_diag_imag,result_grad_theta_up_dn_real,result_grad_theta_up_dn_imag,result_grad_theta_dn_up_real,result_grad_theta_dn_up_imag};
  // std::cout << "abserr: "<<error <<std::endl;
  gsl_integration_workspace_free (w);
  return res;
}


results_t computeGreensFunctionEnergyIntegration(params_t params){
  // Integration parameters
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double error;
  double control = 1e-6;

  double r = (params.r);
  double theta_r = (params.theta_r);

  double result_diag_real;
  double result_diag_imag;
  double result_up_dn_real;
  double result_up_dn_imag;
  double result_dn_up_real;
  double result_dn_up_imag;

  double result_grad_r_diag_real;
  double result_grad_r_diag_imag;
  double result_grad_r_up_dn_real;
  double result_grad_r_up_dn_imag;
  double result_grad_r_dn_up_real;
  double result_grad_r_dn_up_imag;

  double result_grad_theta_diag_real;
  double result_grad_theta_diag_imag;
  double result_grad_theta_up_dn_real;
  double result_grad_theta_up_dn_imag;
  double result_grad_theta_dn_up_real;
  double result_grad_theta_dn_up_imag;

  // std::cout << "computeGreensFunction() r = "<<r << ", theta_r = "<< theta_r << ", ke = "<< params.ke << ", tau = "<< params.tau << std::endl;
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



  // Radial Gradients

  gsl_function Grad_r_diag_real;
  Grad_r_diag_real.function = &grad_r_diag_real;
  Grad_r_diag_real.params = &params;
  gsl_integration_qags (&Grad_r_diag_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_diag_real, &error);

  gsl_function Grad_r_diag_imag;
  Grad_r_diag_imag.function = &grad_r_diag_imag;
  Grad_r_diag_imag.params = &params;
  gsl_integration_qags (&Grad_r_diag_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_diag_imag, &error);

  gsl_function Grad_r_up_dn_real;
  Grad_r_up_dn_real.function = &grad_r_up_dn_real;
  Grad_r_up_dn_real.params = &params;
  gsl_integration_qags (&Grad_r_up_dn_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_up_dn_real, &error);

  gsl_function Grad_r_up_dn_imag;
  Grad_r_up_dn_imag.function = &grad_r_up_dn_imag;
  Grad_r_up_dn_imag.params = &params;
  gsl_integration_qags (&Grad_r_up_dn_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_up_dn_imag, &error);

  gsl_function Grad_r_dn_up_real;
  Grad_r_dn_up_real.function = &grad_r_dn_up_real;
  Grad_r_dn_up_real.params = &params;
  gsl_integration_qags (&Grad_r_dn_up_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_dn_up_real, &error);

  gsl_function Grad_r_dn_up_imag;
  Grad_r_dn_up_imag.function = &grad_r_dn_up_imag;
  Grad_r_dn_up_imag.params = &params;
  gsl_integration_qags (&Grad_r_dn_up_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_r_dn_up_imag, &error);
   

  // Angular Gradients

  gsl_function Grad_theta_diag_real;
  Grad_theta_diag_real.function = &grad_theta_diag_real;
  Grad_theta_diag_real.params = &params;
  gsl_integration_qags (&Grad_theta_diag_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_diag_real, &error);

  gsl_function Grad_theta_diag_imag;
  Grad_theta_diag_imag.function = &grad_theta_diag_imag;
  Grad_theta_diag_imag.params = &params;
  gsl_integration_qags (&Grad_theta_diag_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_diag_imag, &error);

  gsl_function Grad_theta_up_dn_real;
  Grad_theta_up_dn_real.function = &grad_theta_up_dn_real;
  Grad_theta_up_dn_real.params = &params;
  gsl_integration_qags (&Grad_theta_up_dn_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_up_dn_real, &error);

  gsl_function Grad_theta_up_dn_imag;
  Grad_theta_up_dn_imag.function = &grad_theta_up_dn_imag;
  Grad_theta_up_dn_imag.params = &params;
  gsl_integration_qags (&Grad_theta_up_dn_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_up_dn_imag, &error);

  gsl_function Grad_theta_dn_up_real;
  Grad_theta_dn_up_real.function = &grad_theta_dn_up_real;
  Grad_theta_dn_up_real.params = &params;
  gsl_integration_qags (&Grad_theta_dn_up_real, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_dn_up_real, &error);

  gsl_function Grad_theta_dn_up_imag;
  Grad_theta_dn_up_imag.function = &grad_theta_dn_up_imag;
  Grad_theta_dn_up_imag.params = &params;
  gsl_integration_qags (&Grad_theta_dn_up_imag, theta_r-PI/2.0, theta_r+PI/2.0, 0, control, 1000, w, &result_grad_theta_dn_up_imag, &error);

  results_t res { result_diag_real, result_diag_imag, result_up_dn_real, result_up_dn_imag, result_dn_up_real, result_dn_up_imag, result_grad_r_diag_real,result_grad_r_diag_imag, result_grad_r_up_dn_real,result_grad_r_up_dn_imag,result_grad_r_dn_up_real,result_grad_r_dn_up_imag,result_grad_theta_diag_real,result_grad_theta_diag_imag,result_grad_theta_up_dn_real,result_grad_theta_up_dn_imag,result_grad_theta_dn_up_real,result_grad_theta_dn_up_imag};
  // std::cout << "abserr: "<<error <<std::endl;
  gsl_integration_workspace_free (w);
  return res;
}

results_t computeGreensFunctionLDOS(params_t params){
  // Integration parameters
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double error;
  double control = 1e-6;

  double r = (params.r);
  double theta_r = (params.theta_r);

  double result_diag_real = 0.0;
  double result_diag_imag  = 0.0;
  double result_up_dn_real  = 0.0;
  double result_up_dn_imag  = 0.0;
  double result_dn_up_real = 0.0;
  double result_dn_up_imag = 0.0;

  double result_grad_r_diag_real = 0.0;
  double result_grad_r_diag_imag = 0.0;
  double result_grad_r_up_dn_real = 0.0;
  double result_grad_r_up_dn_imag = 0.0;
  double result_grad_r_dn_up_real = 0.0;
  double result_grad_r_dn_up_imag = 0.0;

  double result_grad_theta_diag_real = 0.0;
  double result_grad_theta_diag_imag = 0.0;
  double result_grad_theta_up_dn_real = 0.0;
  double result_grad_theta_up_dn_imag = 0.0;
  double result_grad_theta_dn_up_real = 0.0;
  double result_grad_theta_dn_up_imag = 0.0;

  // std::cout << "computeGreensFunction() r = "<<r << ", theta_r = "<< theta_r << ", ke = "<< params.ke << ", tau = "<< params.tau << std::endl;
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

  results_t res { result_diag_real, result_diag_imag, result_up_dn_real, result_up_dn_imag, result_dn_up_real, result_dn_up_imag, result_grad_r_diag_real,result_grad_r_diag_imag, result_grad_r_up_dn_real,result_grad_r_up_dn_imag,result_grad_r_dn_up_real,result_grad_r_dn_up_imag,result_grad_theta_diag_real,result_grad_theta_diag_imag,result_grad_theta_up_dn_real,result_grad_theta_up_dn_imag,result_grad_theta_dn_up_real,result_grad_theta_dn_up_imag};
  // std::cout << "abserr: "<<error <<std::endl;
  gsl_integration_workspace_free (w);
  return res;
}



// void output_parameters(std::ostream& output){

// }

int main (int argc, char* argv[])
{
  int Nx,Ny,Nz, Nt,EnergyPoints;
  double ratio,tau,initial_coord,final_coord,kef;
  double eff_mass,Ed,Eu,G,Ef;
  bool calc_line, calc_density,calc_total,calc_current,calc_current_int,calc_volume,calc_ldos;
  std::complex<double> t_up,t_dn;
  const std::complex<double> ii(0, 1);


  std::vector<double> x_coords;
  std::vector<double> y_coords;
  std::vector<double> z_coords;
  // GaAs
  eff_mass = 0.067;

  calc_line = false;
  calc_density = false;
  calc_total = false;
  calc_current = false;
  calc_current_int = false;
  calc_volume = false;
  calc_ldos = false;


  // Declare the supported options.
  po::options_description desc("Allowed Options");
  desc.add_options()
      ("help","Print program options")
      ("x-coords,x",po::value<std::vector<double>>(&x_coords)->multitoken()->required(),"X coordinates in nm (e.g -x -10 10)")
      ("y-coords,y",po::value<std::vector<double>>(&y_coords)->multitoken()->required(),"Y coordinates in nm (e.g -y -10 10)")
      ("z-coords,z",po::value<std::vector<double>>(&z_coords)->multitoken()->required(),"Z coordinates in nm (e.g -z -10 10)")
      ("Nx", po::value<int>(&Nx)->required(), "Number of grid points along x-axis")
      ("Ny", po::value<int>(&Ny)->required(), "Number of grid points along y-axis")
      ("Nz", po::value<int>(&Nz)->required(), "Number of grid points along z-axis")
      ("Nt", po::value<int>(&Nt), "Number of angular points (theta)")
      ("Ne", po::value<int>(&EnergyPoints)->required(), "Energy points")
      ("Ef", po::value<double>(&Ef)->required(), "Fermi energy (eV)")
      ("effective-mass,m", po::value<double>(&eff_mass)->required(), "Electron effective mass")
      ("ratio,r", po::value<double>(&ratio), "Spin-orbit strength ratio (alpha/beta)")
      ("tau,t", po::value<double>(&tau), "Spin-orbit strength phase (tau)")
      ("k0", po::value<double>(&k0)->required(), "Spin-orbit wave-vector (1/nm)")
      ("broadening,g", po::value<double>(&G)->required(), "Broadening energy (eV)")
      ("e_up", po::value<double>(&Eu)->required(), "Impurity resonance energy for spin up (eV)")
      ("e_dn", po::value<double>(&Ed)->required(), "Impurity resonance energy for spin down (eV)")
      ("line", po::bool_switch(&calc_line), "Calculate line plot of magnetization (angularly integrated)")
      ("density", po::bool_switch(&calc_density), "Calculate density plot of magnetization")
      ("current", po::bool_switch(&calc_current), "Calculate local current")
      ("current-int", po::bool_switch(&calc_current_int), "Calculate local current integrated in energy")
      ("volume", po::bool_switch(&calc_volume), "Calculate the volume magnetization")
      ("ldos", po::bool_switch(&calc_ldos), "Calculate the LDOS")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help") || argc == 1) {
    std::cout << desc << "\n";
    return 1;
  }

  if (!vm.count("tau")){
    tau = std::atan(ratio);
  }


  // Converts Fermi energy to wave-number in units of 1/nm
  kef = sqrt(2.0*eff_mass*13.097767*Ef); // m_0/\hbar^{2} = 13.097767 eV^{-1}nm^{-2}

  
  gsl_set_error_handler_off();
  // gsl_set_error_handler (&error_handler);

  double xi = x_coords[0];
  double xf = x_coords[1];
  std::vector<double> xs = linspace(xi,xf,Nx);

  double yi = y_coords[0];
  double yf = y_coords[1];
  std::vector<double> ys = linspace(yi,yf,Ny);

  double zi = z_coords[0];
  double zf = z_coords[1];
  std::vector<double> zs = linspace(zi,zf,Nz);

  // double ri = 0.01*(2*M_PI);
  // double rf = 20*(2*M_PI);
  std::vector<double> rs = linspace(xi,xf,Nx);

 
  // std::cout << "xs size: " << xs.size() << std::endl;
  // for (double d : thetas)
  //   std::cout << d << " ";
  // std::cout << std::endl;

  // Model parameters
  // int EnergyPoints = 50;
  double kei = 0.001;
  // double kef = 1.2; // in units of k0
  std::vector<double> kes = linspace(kei,kef,EnergyPoints);


  double dE = (kes[1]*kes[1]-kes[0]*kes[0])/(2.0*13.097767*eff_mass) ;
  std::vector<std::complex<double> > t_ups;
  std::vector<std::complex<double> > t_dns;

  for (int i = 0; i < EnergyPoints; i++){
    // double E = ke*ke/(2.0*13.097767*eff_mass);
    // Energy dependent T-matrix in \hbar^{2}/m^{*}

    t_ups.push_back(-ii*(1.0 + std::exp(std::complex<double>(0, 2)*atan(((kes[i]*kes[i]/(2.0*13.097767*eff_mass))-Eu)/G))));
    t_dns.push_back(-ii*(1.0 + std::exp(std::complex<double>(0, 2)*atan(((kes[i]*kes[i]/(2.0*13.097767*eff_mass))-Ed)/G))));
    // t_dns.push_back(std::complex<double>(0, 0));

  }
  
  // Down - delta = 0 -> t_dn = 0
  // Up  - delta pi/2 -> t_up = -2i
  t_up = std::complex<double>(0.0,-2.0);
  t_dn = std::complex<double>(0.0,0.0);

  // //double tau = 1.31*k0;
  // double t_up = 0.0;
  // double t_dn = -(1.0/(2*PI));



  // double taus[] = {0.785, 0.871, 0.971, 1.15, 1.31, -1.31, PI/2.0};
  // double taus[] = {PI/4.0,0.871, 0.971, 1.15, 1.31, -1.31, PI/2.0};
  // double taus[] = {0.0,M_PI/8.0,3.0*M_PI/8.0,M_PI/2.0};

  
  // for(int t = 0; t < 4; t++){
  //   double tau = taus[t];
  std::cout << "\n---------------[ 2DEG ]---------------" << std::endl;
  std::cout << "ratio : " << ratio << std::endl;
  std::cout << "tau : " << tau << std::endl;
  std::cout << "fermi energy : "<< Ef << " (eV)" << std::endl;
  std::cout << "kef : "<< kef << std::endl;
  std::cout << "grid size :  [" << Nx << "," << Ny << "," << Nz <<"] points"<< std::endl;
  std::cout << "energy integration points : "<< EnergyPoints << std::endl;
  std::cout << "[";
  for (auto ke: kes){
    std::cout << ke*ke/(2.0*13.097767*eff_mass) << " ";
  }
  std::cout << "]" << std::endl;
  std::cout << "dE = " << dE << " eV" << std::endl;
  std::cout << "x :  [" << xi << "," << xf << "] (nm)"<< std::endl;
  std::cout << "y :  [" << yi << "," << yf << "] (nm) "<< std::endl;
  std::cout << "z :  [" << zi << "," << zf << "] (nm) "<< std::endl;
  std::cout << "\n---------------[ Impurity ]---------------" << std::endl;
  std::cout << "e_up : "<< Eu << " (eV)" << std::endl;
  std::cout << "e_down : "<< Ed << " (eV)" << std::endl;
  std::cout << "broadening : "<< G << " (eV)" << std::endl;
  std::cout << "------------------------------------------" << std::endl;

  auto start = std::chrono::steady_clock::now();


   // Current energy integrated plot
  if (calc_current_int){
    std::cout << "Calculating current integrated in energy..." << std::endl;
    std::string filename = "current_energy_integrated_r_"+std::to_string(ratio)+".dat";
    std::ofstream fcur;
    fcur.open(filename);
  
          
    for( double x: xs){
      for (double y: ys){

        double r = sqrt(pow(x,2)+pow(y,2));
        double theta_r = atan2(y,x);
        
        // Current
        double J_r = 0.0;
        double J_theta = 0.0;


        for (int i = 0; i < EnergyPoints; i++){

          params_t params = {r,theta_r,kes[i],tau};

          //double E = ke*ke/(2.0*13.097767*eff_mass);

          // Energy dependent T-matrix in \hbar^{2}/m^{*}
          t_up = t_ups[i];
          t_dn = t_dns[i];
          results_t greens_functions = computeGreensFunction(params);
          std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
          std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
          std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
          std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
          std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
          std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
          std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
          std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
          std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


          std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
          std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

          std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
          std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

          double J_r_grad = 1.937025*(-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
          double J_theta_grad = 1.937025*(-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

          double J_r_so = 3.87405*(-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
          double J_theta_so = 3.87405*(-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


          J_r +=  (J_r_grad + J_r_so)*dE;
          J_theta +=  (J_theta_grad + J_theta_so)*dE;

          // double J_x = J_r*cos(theta_r) - J_theta*sin(theta_r);
          // double J_y = J_r*sin(theta_r) + J_theta*cos(theta_r);
                    
        }
        fcur << r << "," << theta_r << "," << J_r << "," << J_theta << std::endl;
      }
    }
    std::cout << "Energy-integrated current density was written in " << filename << std::endl;
    fcur.close();
    
  }

   // Current plot
  if (calc_current){
    std::cout << "Calculating current ..." << std::endl;
    std::string filename = "current_energy_"+std::to_string(Ef)+"_r_"+std::to_string(ratio)+".dat";
    std::ofstream fcur;
    fcur.open(filename);
    
    // double xi = initial_coord;
    // double xf = final_coord;
    // std::vector<double> xvec = linspace(xi,xf,25);


    // // double yi = -2*(2*M_PI);
    // // double yf = 2*(2*M_PI);
    // double yi = initial_coord;
    // double yf = final_coord;
    // std::vector<double> yvec = linspace(yi,yf,25);

    double ke = sqrt((2.0*eff_mass)*Ef/0.076348);
    // Energy dependent T-matrix in \hbar^{2}/m^{*}
    t_up = -(ii + std::exp(std::complex<double>(0, 2)*atan((Ef-Eu)/G)));
    t_dn = -(ii + std::exp(std::complex<double>(0, 2)*atan((Ef-Ed)/G)));

    for( double x: xs){
      for (double y: ys){

        double r = sqrt(pow(x,2)+pow(y,2));
        double theta_r = atan2(y,x);

        params_t params = {r,theta_r,ke,tau};

        results_t greens_functions = computeGreensFunction(params);
        std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
        std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
        std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
        std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
        std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
        std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
        std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
        std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
        std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


        std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
        std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

        std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
        std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

        double J_r_grad = 1.937025*(-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
        double J_theta_grad = 1.937025*(-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

        double J_r_so = 3.87405*(-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
        double J_theta_so = 3.87405*(-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


        double J_r =  J_r_grad + J_r_so;
        double J_theta =  J_theta_grad + J_theta_so;


        // double J_x = J_r*cos(theta_r) - J_theta*sin(theta_r);
        // double J_y = J_r*sin(theta_r) + J_theta*cos(theta_r);
        
        fcur << r << "," << theta_r << "," << J_r << "," << J_theta << std::endl;
        
      } 
    }
    std::cout << "Current written in " << filename << std::endl;
    fcur.close();
    
  }

  // Line plot
  if (calc_line){
    std::cout << "Nt = " << Nt << std::endl;
    std::cout << "Calculating magnetization line plot (angularly integrated)..." << std::endl;
    auto start = std::chrono::steady_clock::now(); 
    std::string filename_jint = "magmom_line_r_"+std::to_string(ratio)+".dat";
    std::ofstream fjint;
    fjint.open(filename_jint);

    double ti = 0.0;
    double tf = 2*M_PI;
    std::vector<double> thetas = linspace(ti,tf,Nt);
    
    for (double r: rs){
      double mz = 0.0;
      double msx = 0.0;
      double msy = 0.0;
      double msz = 0.0;
      
      for (double theta_r : thetas){
        double J_r = 0.0;
        double J_theta = 0.0;
        for (int i = 0; i < EnergyPoints; i++){          
          params_t params = {r,theta_r,kes[i],tau};

          // Energy dependent T-matrix in \hbar^{2}/m^{*}
          t_up = t_ups[i];
          t_dn = t_dns[i];

          results_t greens_functions = computeGreensFunction(params);
          std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
          std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
          std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
          std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
          std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
          std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
          std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
          std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
          std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


          std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
          std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

          std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
          std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

          double J_r_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
          double J_theta_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

          double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
          double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


          J_r +=  J_r_grad + J_r_so;
          J_theta +=  J_theta_grad + J_theta_so;

          msx += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_diag + G_up_dn*G_dn_up));
          msy += (-1.0/PI)*std::imag(-ii*(t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up));
          msz += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up));


        }
        double J_x = J_r*cos(theta_r) - J_theta*sin(theta_r);
        double J_y = J_r*sin(theta_r) + J_theta*cos(theta_r);
        mz += (r*cos(theta_r)*J_y - r*sin(theta_r)*J_x);
        // std::cout << mz  << std::endl; 
      }
      fjint << r << "," << mz  << "," << msx  << "," << msy  << "," << msz << std::endl; 
      std::cout << r << "," << mz  << "," << msx  << "," << msy  << "," << msz << std::endl; 
    }
    
    fjint.close();
    std::cout << "Calculations written in " << filename_jint << std::endl;
  }

  // Density plot
  if (calc_density){
    std::cout << "Calculating density plots for orbital and spin magnetization ..." << std::endl;
    std::string filename = "magmom_energy_integrated_r_"+std::to_string(ratio)+".dat";
    std::string filename_spin = "spinmag_energy_integrated_r_"+std::to_string(ratio)+".dat";
    std::ofstream fmag;
    std::ofstream fspin;
    fmag.open(filename);
    fspin.open(filename_spin);

    double sum_orb  = 0.0;
    double sum_spin  = 0.0;

    
    for( double x: xs){
      for (double y: ys){

        double r = sqrt(pow(x,2)+pow(y,2));
        double theta_r = atan2(y,x);

        // Current
        double J_r = 0.0;
        double J_theta = 0.0;

        // Spin magnetization
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;

       

        for (int i = 0; i < EnergyPoints; i++){
          params_t params = {r,theta_r,kes[i],tau};


          // Energy dependent T-matrix in \hbar^{2}/m^{*}
          

          results_t greens_functions = computeGreensFunction(params);
          std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
          std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
          std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
          std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
          std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
          std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
          std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
          std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
          std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


          std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
          std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

          std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
          std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

          double J_r_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
          double J_theta_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

          double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
          double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


          J_r +=  (J_r_grad + J_r_so)*dE;
          J_theta +=  (J_theta_grad + J_theta_so)*dE;

          sx += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_diag + G_up_dn*G_dn_up))*dE;
          sy += (-1.0/PI)*std::imag(-ii*(t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up))*dE;
          sz += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up))*dE;
        
        }
        double J_x = J_r*cos(theta_r) - J_theta*sin(theta_r);
        double J_y = J_r*sin(theta_r) + J_theta*cos(theta_r);

        sum_orb += (x*J_y - y*J_x);
        sum_spin += sqrt(sx*sx + sy*sy + sz*sz);
        
        fmag << x << "," << y << "," << (x*J_y - y*J_x) << std::endl; // In units of \mu_{B}.nm^{-2}
        fspin << x << "," << y << "," << sx << "," << sy << "," << sz << std::endl;
        
      } 
    }
    std::cout << "Orbital magnetization written in " << filename << std::endl;
    std::cout << "Spin magnetization written in " << filename_spin << std::endl;
    std::cout << "Total orbital magnetization = " <<sum_orb << std::endl;
    std::cout << "Total spin magnetization = " << sum_spin << std::endl;
    fmag.close();
    fspin.close();
  }

  // Volume density plot
  if (calc_volume){
    std::cout << "Calculating volume density for orbital and spin magnetization..." << std::endl;
    std::string filename = "magmom_vol_energy_integrated_r_"+std::to_string(ratio)+".dat";
    std::string filename_spin = "spinmag_vol_energy_integrated_r_"+std::to_string(ratio)+".dat";

    std::ofstream fmag;
    std::ofstream fspin;


    fmag.open(filename);
    fspin.open(filename_spin);

    fmag << "#---------------[Parameters]---------------" << std::endl;
    fmag << "# ratio : " << ratio << std::endl;
    fmag << "# tau : " << tau << std::endl;
    fmag << "# fermi energy : "<< Ef << " (eV)" << std::endl;
    fmag << "# kf : "<< kef << " (1/nm)" << std::endl;
    fmag << "# grid size :  [" << Nx << "," << Ny << "," << Nz <<"] points"<< std::endl;
    fmag << "# energy integration points : "<< EnergyPoints << std::endl;
    fmag << "# x :  [" << xi << "," << xf << "] (in units of nm)"<< std::endl;
    fmag << "# y :  [" << yi << "," << yf << "] (in units of nm) "<< std::endl;
    fmag << "# z :  [" << zi << "," << zf << "] (in units of nm) "<< std::endl;
    fmag << "#------------------------------------------" << std::endl;

    fspin << "#---------------[Parameters]---------------" << std::endl;
    fspin << "# ratio : " << ratio << std::endl;
    fspin << "# tau : " << tau << std::endl;
    fspin << "# fermi energy : "<< Ef << " (eV)" << std::endl;
    fspin << "# kf : "<< kef << " (1/nm)" << std::endl;
    fspin << "# grid size :  [" << Nx << "," << Ny << "," << Nz <<"] points"<< std::endl;
    fspin << "# energy integration points : "<< EnergyPoints << std::endl;
    fspin << "# x :  [" << xi << "," << xf << "] (in units of nm)"<< std::endl;
    fspin << "# y :  [" << yi << "," << yf << "] (in units of nm) "<< std::endl;
    fspin << "# z :  [" << zi << "," << zf << "] (in units of nm) "<< std::endl;
    fspin << "#------------------------------------------" << std::endl;

    double L = abs(z_coords[1] - z_coords[0]);

    std::cout << "Quantum-well height (z) : "<< L <<  " nm" << std::endl;
    double jx[Nx][Ny];
    double jy[Nx][Ny];
    

    double sx[Nx][Ny];
    double sy[Nx][Ny];
    double sz[Nx][Ny];


    for( int i = 0; i < Nx; i++){
      for ( int j = 0; j < Ny; j++){

        double r = sqrt(pow(xs[i],2)+pow(ys[j],2));
        double theta_r = atan2(ys[j],xs[i]);

        // Partial variables
        double J_r = 0.0;
        double J_theta = 0.0;

        double psx = 0.0;
        double psy = 0.0;
        double psz = 0.0;

    
        for (int e = 0; e < EnergyPoints; e++){
          params_t params = {r,theta_r,kes[e],tau};

          // Energy dependent T-matrix in \hbar^{2}/m^{*}
          t_up = t_ups[e];
          t_dn = t_dns[e];

          results_t greens_functions = computeGreensFunction(params);
          std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
          std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
          std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
          std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
          std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
          std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
          std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
          std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
          std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


          std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
          std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

          std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
          std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

          double J_r_grad = 1.937025*(-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
          double J_theta_grad = 1.937025*(-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

          double J_r_so = 3.87405*(-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
          double J_theta_so = 3.87405*(-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


          J_r +=  (J_r_grad + J_r_so)*dE;
          J_theta +=  (J_theta_grad + J_theta_so)*dE;

          psx += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_diag + G_up_dn*G_dn_up))*dE;
          psy += (-1.0/PI)*std::imag(-ii*(t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up))*dE;
          psz += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up))*dE;
        }
        jx[i][j] = J_r*cos(theta_r) - J_theta*sin(theta_r);
        jy[i][j] = J_r*sin(theta_r) + J_theta*cos(theta_r);
        sx[i][j] = psx;
        sy[i][j] = psy;
        sz[i][j] = psz;
        // double J_x = (J_r*cos(theta_r) - J_theta*sin(theta_r))*(2.0/L)*pow(cos(M_PI*z/L),2);
        // double J_y = (J_r*sin(theta_r) + J_theta*cos(theta_r))*(2.0/L)*pow(cos(M_PI*z/L),2);

        // sx = sx*(2.0/L)*pow(cos(M_PI*z/L),2);
        // sy = sy*(2.0/L)*pow(cos(M_PI*z/L),2);
        // sz = sz*(2.0/L)*pow(cos(M_PI*z/L),2);
        
      } 
    }
    for( int i = 0; i < Nx; i++){
      for ( int j = 0; j < Ny; j++){
        for ( int k = 0; k < Nz; k++){
          double mx = -(1.0/L)*zs[k]*jy[i][j]*pow(cos(M_PI*zs[k]/L),2)/2.206271475;
          double my = (1.0/L)*zs[k]*jx[i][j]*pow(cos(M_PI*zs[k]/L),2)/2.206271475;
          double mz = (1.0/L)*(xs[i]*jy[i][j] - ys[j]*jx[i][j])*pow(cos(M_PI*zs[k]/L),2)/2.206271475;

          sx[i][j] *=(1.0/L)*pow(cos(M_PI*zs[k]/L),2);
          sy[i][j] *=(1.0/L)*pow(cos(M_PI*zs[k]/L),2);
          sz[i][j] *= (1.0/L)*pow(cos(M_PI*zs[k]/L),2);
          fmag << xs[i] << "," << ys[j]  << "," << zs[k] << "," << mx << "," << my << "," << mz << std::endl; // In units of \mu_{B}.nm^{-3}          
          fspin << xs[i] << "," << ys[j] << "," << zs[k] << "," << sx[i][j] << "," << sy[i][j] << "," << sz[i][j] << std::endl;
        }
      }
    }

    std::cout << "Orbital magnetization written in " << filename << std::endl;
    std::cout << "Spin magnetization written in " << filename_spin << std::endl;
    fmag.close();
    fspin.close();
  }

  // Total orbital magnetization
  if (calc_total){
    std::cout << "Calculating total orbital magnetization..." << std::endl;
    std::string integrated_magmom_filename = "integrated_magmom_r_" + std::to_string(ratio) + ".dat";
    std::ofstream fint;
    fint.open(integrated_magmom_filename);

    double magmon = 0.0;   
    for( double x: xs){
      for (double y: ys){
        double r = sqrt(pow(x,2)+pow(y,2));
        double theta_r = atan2(y,x);

        double J_r = 0.0;
        double J_theta = 0.0;
        for (double ke: kes){
           double E = 0.076348*ke*ke/(2.0*eff_mass);

          // Energy dependent T-matrix in \hbar^{2}/m^{*}
          t_up = -(ii + std::exp(std::complex<double>(0, 2)*atan((E-Eu)/G)));
          t_dn = -(ii + std::exp(std::complex<double>(0, 2)*atan((E-Ed)/G)));

          params_t params = {r,theta_r,ke,tau};

          results_t greens_functions = computeGreensFunction(params);
          std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
          std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
          std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
          std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
          std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
          std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
          std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
          std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
          std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


          std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
          std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

          std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
          std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

          double J_r_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
          double J_theta_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

          double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
          double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


          J_r += J_r_grad + J_r_so;
          J_theta += J_theta_grad + J_theta_so;

        }

        double J_x = J_r*cos(theta_r) - J_theta*sin(theta_r);
        double J_y = J_r*sin(theta_r) + J_theta*cos(theta_r);
        magmon += x*J_y - y*J_x;
      }
    }
    fint << tau << "," << magmon << std::endl;
    std::cout << "Calculations written in " << integrated_magmom_filename << std::endl;
  }

  // LDOS
  if (calc_ldos){
    std::cout << "Calculating LDOS ..." << std::endl;
    std::string filename = "ldos_r_"+std::to_string(ratio)+".dat";
    
    std::ofstream fldos;
    
    fldos.open(filename);

    for( double x: xs){
      for (double y: ys){

        double r = sqrt(pow(x,2)+pow(y,2));
        double theta_r = atan2(y,x);

        double ldos = 0.0;

        for (int i = 0; i < EnergyPoints; i++){
          params_t params = {r,theta_r,kes[i],tau};


          // Energy dependent T-matrix in \hbar^{2}/m^{*}
          t_up = t_ups[i];
          t_dn = t_dns[i];

          results_t greens_functions = computeGreensFunctionLDOS(params);
          std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
          std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
          std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
          
          ldos +=  (-1.0/PI)*std::imag((t_up-t_dn)*(pow(G_diag,2) - G_dn_up*G_up_dn));

        }
       
       
        
      fldos << x << "," << y << "," << ldos << std::endl; // In units of \mu_{B}.nm^{-2}
        
      } 
    }
    std::cout << "LDOS written in " << filename << std::endl;
    fldos.close();
  }

  auto end = std::chrono::steady_clock::now(); 

  std::chrono::duration<double> elapsed_seconds = end-start;

  std::cout << "Total computing time: " << elapsed_seconds.count() << "s" << std::endl; 
  
  
  return 0;
}