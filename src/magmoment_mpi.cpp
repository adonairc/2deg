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
#include <boost/mpi.hpp>


#define PI M_PI
double k0 = 0.01241; // GaAs (beta ~ 10 meV.nm)
namespace po = boost::program_options;
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

typedef struct {
  double r;
  double theta_r;
  double e;
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
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(e+pow(k0,2)*f) ;
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(abs(qplus*rho));
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(abs(qminus*rho));


  double T1 = (qplus*qplus*R2plus + qminus*qminus*R2minus)*cos(theta_q-theta_r);
  double T2 = (qminus + qplus)/r;

  return((1.0/(2*M_PI*M_PI))*(T1+T2)/Q);
}

double grad_r_diag_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));


  double T1 = (PI/2.0)*(qplus*qplus*sin(qplus*rho) + qminus*qminus*sin(qminus*rho))*cos(theta_q-theta_r);

  return((1.0/(2*M_PI*M_PI))*T1/Q);
}


double grad_r_dn_up_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_r_dn_up_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = (PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_r_up_dn_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = -sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_r_up_dn_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = -(PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));


  return((1.0/(2*M_PI*M_PI))*(T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
}


// Angular gradients

double grad_theta_diag_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(abs(qplus*rho));
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(abs(qminus*rho));


  double T1 = (qplus*qplus*R2plus + qminus*qminus*R2minus)*sin(theta_q-theta_r);
  double T2 = -(qminus + qplus)*tan(theta_q-theta_r)/r;

  return(-(1.0/(2*M_PI*M_PI))*(T1+T2)/Q);
}

double grad_theta_diag_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double T1 = (PI/2.0)*(qplus*qplus*sin(qplus*rho) + qminus*qminus*sin(qplus*rho))*sin(theta_q-tau);
  return(-(1.0/(2*M_PI*M_PI))*T1/Q);
}


double grad_theta_dn_up_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_theta_dn_up_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = (PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_theta_up_dn_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = (PI/2.0)*cos(theta_q-tau)*(-qplus*qplus*cos(qplus*rho) + qminus*qminus*cos(qminus*rho));
  double T2 = -sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

double grad_theta_up_dn_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R1plus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double R1minus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);


  double T1 = cos(theta_q-tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);
  
  double T2 = -(PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((1.0/(2*M_PI*M_PI))*(T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

// G0
double radial_dn_up_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(abs(qplus*rho));
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(abs(qminus*rho));

  double S = qplus*(sin(theta_q+tau)*R2plus - (PI/2.0)*cos(theta_q-tau)*sin(qplus*rho)) - qminus*(sin(theta_q+tau)*R2minus - (PI/2.0)*cos(theta_q-tau)*sin(qminus*rho));
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
}

double radial_dn_up_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));


  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(abs(qplus*rho));
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(abs(qminus*rho));

  double S = qplus*(cos(theta_q-tau)*R2plus + (PI/2.0)*sin(theta_q+tau)*sin(qplus*rho)) - qminus*(cos(theta_q-tau)*R2minus+(PI/2.0)*sin(theta_q+tau)*sin(qminus*rho));
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
}


double radial_up_dn_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(abs(qplus*rho));
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(abs(qminus*rho));

  double S = -qplus*(sin(theta_q+tau)*R2plus + (PI/2.0)*cos(theta_q-tau)*sin(qplus*rho)) + qminus*(sin(theta_q+tau)*R2minus + (PI/2.0)*cos(theta_q-tau)*sin(qminus*rho));
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
}

double radial_up_dn_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(abs(qplus*rho));
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(abs(qminus*rho));

  double S = qplus*(cos(theta_q-tau)*R2plus - (PI/2.0)*sin(theta_q+tau)*sin(qplus*rho)) - qminus*(cos(theta_q-tau)*R2minus - (PI/2.0)*sin(theta_q+tau)*sin(qminus*rho));
  return((1.0/(2*M_PI*M_PI))*S/(Q*sqrt(f)));
}


double radial_diagonal_real(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));


  double Iplus = -cos(qplus*rho)*gsl_sf_Ci(abs(qplus*rho))-sin(qplus*rho)*gsl_sf_Si(qplus*rho);
  double Iminus = -cos(qminus*rho)*gsl_sf_Ci(abs(qminus*rho))-sin(qminus*rho)*gsl_sf_Si(qminus*rho);
  // std::cout << "radial_diagonal_real() Q = " << Q << std::endl;

  return(-(1.0/(2*M_PI*M_PI))*(qplus*Iplus+qminus*Iminus)/Q);

}

double radial_diagonal_imag(double theta_q, void * p) {
  params_t * params = (params_t *) p;
  double r = (params->r);
  double theta_r = (params->theta_r);
  double e = (params->e);
  double tau = (params->tau);

  double rho = r*cos(theta_q-theta_r);
  double f = 1 + sin(2*tau)*sin(2*theta_q);
  double Q = sqrt(abs(e+pow(k0,2)*f));
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));


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

//
// Fermi-Dirac distribution function
//
// k = 8.617333262145×10−5 eV.K^{-1}

// double fermi(double ke, double kf, double eff_mass, double beta){
//   // ke*ke/(2.0*eff_mass*13.097767) = E;
//   // kf*kf/(2.0*eff_mass*13.097767) = E;
//   return 1.0/(exp(beta*((ke*ke/(2.0*eff_mass*13.097767))-(kf*kf/(2.0*eff_mass*13.097767))))+1.0);
// }
int main (int argc, char* argv[])
{
  int Nx,Ny,Nz, Nt,EnergyPoints;
  double tau,ratio,alpha,beta,initial_coord,final_coord,kf,ef;
  double eff_mass,Ed,Eu,G,Ef;
  std::complex<double> t_up,t_dn;
  const std::complex<double> ii(0, 1);

  mpi::environment env;
  mpi::communicator world;

  std::vector<double> x_coords;
  std::vector<double> y_coords;
  std::vector<double> z_coords;
  // GaAs
  eff_mass = 0.067;




  // Declare the supported options.
  po::options_description desc("Allowed Options");
  desc.add_options()
      ("help","Print program options")
      ("x-coords,x",po::value<std::vector<double>>(&x_coords)->multitoken()->required(),"X coordinates in nm (e.g -x -10 10)")
      ("y-coords,y",po::value<std::vector<double>>(&y_coords)->multitoken()->required(),"Y coordinates in nm (e.g -y -10 10)")
      // ("z-coords,z",po::value<std::vector<double>>(&z_coords)->multitoken(),"Z coordinates in nm (e.g -z -10 10)")
      ("Nx", po::value<int>(&Nx)->required(), "Number of grid points along x-axis")
      ("Ny", po::value<int>(&Ny)->required(), "Number of grid points along y-axis")
      // ("Nz", po::value<int>(&Nz), "Number of grid points along z-axis")
      ("Nt", po::value<int>(&Nt), "Number of angular points (theta)")
      ("Ne", po::value<int>(&EnergyPoints)->required(), "Energy points")
      ("Ef", po::value<double>(&Ef), "Fermi energy (eV)")
      ("effective-mass,m", po::value<double>(&eff_mass)->required(), "Electron effective mass")
      ("ratio,r", po::value<double>(&ratio), "Spin-orbit strength ratio (alpha/beta)")
      ("tau,t", po::value<double>(&tau), "Rashba spin-orbit strength (alpha)")
      ("k0", po::value<double>(&k0), "Dresselhaus spin-orbit strength (beta)")
      ("alpha,a", po::value<double>(&alpha), "Rashba spin-orbit strength (alpha)")
      ("beta,b", po::value<double>(&beta), "Dresselhaus spin-orbit strength (beta)")
      ("kf", po::value<double>(&kf), "Fermi wave-vector (1/nm)")
      ("ef", po::value<double>(&ef), "Fermi energy (eV)")
      ("broadening,g", po::value<double>(&G), "Broadening energy (eV)")
      ("e_up", po::value<double>(&Eu), "Impurity resonance energy for spin up (eV)")
      ("e_dn", po::value<double>(&Ed), "Impurity resonance energy for spin down (eV)")
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

  
  gsl_set_error_handler_off();
  // gsl_set_error_handler (&error_handler);

  double xi = x_coords[0];
  double xf = x_coords[1];
  std::vector<double> xs = linspace(xi,xf,Nx);

  double yi = y_coords[0];
  double yf = y_coords[1];
  std::vector<double> ys = linspace(yi,yf,Ny);

  // double zi = z_coords[0];
  // double zf = z_coords[1];
  // std::vector<double> zs = linspace(zi,zf,Nz);

  // double ri = 0.01*(2*M_PI);
  // double rf = 20*(2*M_PI);
  std::vector<double> rs = linspace(xi,xf,Nx);

  std::vector< std::array<double,2> > points;
  std::vector< std::array<double,8> > partial_results;

  for (auto x : xs){
    for (auto y: ys){
      points.push_back(std::array<double,2>{x,y});
    }
  }

  int start = world.rank()*points.size()/world.size();
  int end;
  if (world.rank() == (world.size() - 1)){
    end = points.size();
  } else {
    end = start + points.size()/world.size();
  }

 
  if (vm.count("alpha") && vm.count("beta")){
    k0 = 13.097767*eff_mass*sqrt(alpha*alpha + beta*beta);
    tau = atan2(alpha,beta);
    ratio = alpha/beta;
  }
  else if (vm.count("tau") && vm.count("k0")){
    ratio = tan(tau);
    alpha = k0*sin(tau)/(13.097767*eff_mass);
    beta = k0*cos(tau)/(13.097767*eff_mass);
  }
  else {
    std::cout << "Please enter SO parameters" << std::endl;
    exit(1);
  }

  if (!vm.count("kf")){
      kf = sqrt(2.0*13.097767*eff_mass*ef);
  }
  

  t_up = std::complex<double>(1.0,0.0);
  t_dn = std::complex<double>(-1.0,0.0);

  if (world.rank() == 0) {

    std::cout << "\n---------------[ 2DEG ]---------------" << std::endl;
    std::cout << "alpha : " << alpha << std::endl;
    std::cout << "beta : " << beta << std::endl;
    std::cout << "tau: " << tau << std::endl;
    std::cout << "k0 : " << k0 << std::endl;
    std::cout << "kf : "<< kf << std::endl;
    std::cout << "grid size :  [" << Nx << "," << Ny << "] points"<< std::endl;
    std::cout << "energy integration points : "<< EnergyPoints << std::endl;
    std::cout << "x :  [" << xi << "," << xf << "] (nm)"<< std::endl;
    std::cout << "y :  [" << yi << "," << yf << "] (nm) "<< std::endl;

  }

  std::cout << "Process #" << world.rank() << " -> [" << start << "," << end -1 <<"] job size = " << (end-start) << " points" <<  std::endl;


  // Density plot
    // std::cout << "Calculating density plots for current, orbital and spin magnetization ..." << std::endl;
    
    double temperature = 150; // in Kelvin
    double boltzmann_factor = (8.617333262145e-5)*temperature; // in 1/eV

    // double e_initial = -2.0*k0*k0/(2.0*13.097767*eff_mass);
    // double e_final = kef*kef/(2.0*13.097767*eff_mass);
    // std::cout << "kf = " << kf << std::endl;
    double Einitial = -13.097767*eff_mass*pow(alpha+beta,2)/2.0+0.00001;
    // double Einitial = -2.0(*k0*k0);
    double Efinal = (kf*kf) ;
    std::vector<double> energies = linspace(Einitial,Efinal,EnergyPoints);
    // double dE = (kes[1]*kes[1]-kes[0]*kes[0])/(2.0*13.097767*eff_mass);
    double dE = abs(energies[1]/(2.0*13.097767*eff_mass) - energies[0]/(2.0*13.097767*eff_mass)); // eV (k in 1/nm)
    // std::cout << "dE = " << dE << std::endl;
    for (int i = start; i < end; i++){
        double x = points[i][0];
        double y = points[i][1];
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
          // double k = sqrt(energies[i]*(2.0*13.097767*eff_mass));
          // double energy = kes[i]*kes[i]/(2.0*13.097767*eff_mass);
          params_t params = {r,theta_r,energies[i],tau};
          // std::cout << "E = " << kes[i]*kes[i]/(2.0*13.097767*eff_mass) << std::endl;
          // std::cout << "k^2 = " << energies[i] << std::endl;

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


          // J_r +=  (J_r_grad + J_r_so)*(1.0/(exp(-beta*(energies[i]-(kf*kf))/(2.0*13.097767*eff_mass))))*dE;
          // J_theta +=  (J_theta_grad + J_theta_so)*(1.0/(exp(-beta*(energies[i]-(kf*kf))/(2.0*13.097767*eff_mass))))*dE;


          J_r +=  (J_r_grad + J_r_so)*dE;
          J_theta +=  (J_theta_grad + J_theta_so)*dE;

          // J_r +=  (J_r_grad)*dE;
          // J_theta +=  (J_theta_grad)*dE;

          sx += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_diag + G_up_dn*G_dn_up))*dE;
          sy += (-1.0/PI)*std::imag(-ii*(t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up))*dE;
          sz += (-1.0/PI)*std::imag((t_up-t_dn)*(G_diag*G_up_dn + G_diag*G_dn_up))*dE;
        
        }
        double J_x = J_r*cos(theta_r) - J_theta*sin(theta_r);
        double J_y = J_r*sin(theta_r) + J_theta*cos(theta_r);


        partial_results.push_back(std::array<double,8>{x,y,(x*J_y - y*J_x),sx,sy,sz,J_x,J_y});
        
        // fmag << x << "," << y << "," << (x*J_y - y*J_x) << std::endl; // In units of \mu_{B}.nm^{-2}
        // fspin << x << "," << y << "," << sx << "," << sy << "," << sz << std::endl;
        // fj << x << "," << y << "," << sqrt(J_x*J_x + J_y*J_y) << "," << sz << std::endl;
         
    }

    auto req = world.isend(0,0,partial_results);

    if (world.rank() == 0) {
      std::string filename = "magmom_alpha_"+std::to_string(alpha)+"_beta_"+std::to_string(beta)+".dat";
      std::string filename_spin = "spinmag_alpha_"+std::to_string(alpha)+"_beta_"+std::to_string(beta)+".dat";
      std::string filename_j = "current_alpha"+std::to_string(alpha)+"_beta_"+std::to_string(beta)+".dat";
      std::ofstream fmag;
      std::ofstream fspin;
      std::ofstream fj;
      fmag.open(filename);
      fspin.open(filename_spin);
      fj.open(filename_j);

      for (int p = 0; p < world.size(); p++){
        std::vector<std::array<double,8>>data_received;
        // std::cout << "Receiving from #" << p << "\n";
        world.recv(p,0,data_received);
        for (auto point: data_received){
          fmag << point[0] << "," <<  point[1] << "," <<  point[2] << std::endl;
          fspin <<  point[0] << "," <<  point[1] << "," << point[3] << "," <<  point[4] << "," <<  point[5] << std::endl;
          fj <<  point[0] << "," <<  point[1] << "," << point[6] << "," <<  point[7] <<  std::endl;
        }
      }

      std::cout << "Current written in " << filename_j << std::endl;
      std::cout << "Orbital magnetization written in " << filename << std::endl;
      std::cout << "Spin magnetization written in " << filename_spin << std::endl;
      fmag.close();
      fspin.close();
      fj.close();
    }
  
  
  return 0;
}