#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <string>

#define PI M_PI
double k0 = 1.0;

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

// Integration parameters
gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
double error;
double control = 1e-7;

void error_handler (const char * reason, const char * file, int line, int gsl_errno) {
  std::cout << "Error (" << gsl_errno << ") : " << reason << ", file = " << file  << std::endl;
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
  double Q = sqrt(pow(ke,2)+pow(k0,2)*f);
  double qplus = abs(Q+k0*sqrt(f));
  double qminus = abs(Q-k0*sqrt(f));

  double R2plus = cos(qplus*rho)*gsl_sf_Si(qplus*rho)-sin(qplus*rho)*gsl_sf_Ci(qplus*rho);
  double R2minus = cos(qminus*rho)*gsl_sf_Si(qminus*rho)-sin(qminus*rho)*gsl_sf_Ci(qminus*rho);


  double T1 = (qplus*qplus*R2plus + qminus*qminus*R2minus)*cos(theta_q-theta_r);
  double T2 = (qminus + qplus)/r;

  return((T1+T2)/Q);
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

  return(T1/Q);
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
  double T2 = sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
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
  
  double T2 = (PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
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
  double T2 = -sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
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
  
  double T2 = -(PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((T1+T2)*cos(theta_q-theta_r)/(Q*sqrt(f)));
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

  return(-(T1+T2)/Q);
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
  return(-T1/Q);
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
  double T2 = sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
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
  
  double T2 = (PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
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
  double T2 = -sin(theta_q+tau)*(qplus*qplus*R1plus - qminus*qminus*R1minus);

  return((T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
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
  
  double T2 = -(PI/2.0)*sin(theta_q+tau)*(qplus*qplus*cos(qplus*rho) - qminus*qminus*cos(qminus*rho));

  return((T1+T2)*sin(theta_q-theta_r)/(Q*sqrt(f)));
}

// G0
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

  double S = qplus*(sin(theta_q+tau)*R2plus - (PI/2.0)*cos(theta_q-tau)*sin(qplus*rho)) - qminus*(sin(theta_q+tau)*R2minus - (PI/2.0)*cos(theta_q-tau)*sin(qminus*rho));
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

  double S = qplus*(cos(theta_q-tau)*R2plus + (PI/2.0)*sin(theta_q+tau)*sin(qplus*rho)) - qminus*(cos(theta_q-tau)*R2minus+(PI/2.0)*sin(theta_q+tau)*sin(qminus*rho));
  return(S/(Q*sqrt(f)));
}


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

  double S = -qplus*(sin(theta_q+tau)*R2plus + (PI/2.0)*cos(theta_q-tau)*sin(qplus*rho)) + qminus*(sin(theta_q+tau)*R2minus + (PI/2.0)*cos(theta_q-tau)*sin(qminus*rho));
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

  // std::cout << "r = "<<r << ", theta_r = "<< theta_r << std::endl;
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
  return res;
}


int main (int argc, char* argv[])
{
  
  int N = 30;
  gsl_set_error_handler_off();
  //gsl_set_error_handler (&error_handler);

  double xi = -3.0;
  double xf = 3.0;
  double dx = (xf-xi+1)/N;

  double yi = -3.0;
  double yf = 3.0;
  double dy = (yf-yi+1)/N;


  double ri = 0.01;
  double rf = 5;
  double dr = (rf-ri+1)/N;


  double ti = 0.0;
  double tf = 2*PI;
  double dt = (tf-ti)/N;


  // Model parameters
  double ke = 1.2;
  //double tau = 1.31*k0;
  double t_up = 0.0;
  double t_dn = -(1.0/(2*PI));;

  const std::complex<double> ii(0, 1);

  double taus[] = {0.785, 0.871, 0.971, 1.15, 1.31, -1.31};
  // double taus[] = {PI/2.0};

  double rs[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

  std::cout << "Calculating charge current..." << std::endl;

  for(int t = 0; t < 6; t++){
    double tau = taus[t];

    std::string filename = "current_tau_"+std::to_string(tau)+".dat";
    // std::string filename_grad = "grad_tau_"+std::to_string(tau)+".dat";
    std::string filename_jint = "jintegral_tau_"+std::to_string(tau)+".dat";
  
    std::ofstream fcur;
    std::ofstream fldos;
    // std::ofstream fgrad;
    std::ofstream fjint;

    fcur.open(filename);
    fldos.open("ldos_test.dat");
    // fgrad.open(filename_grad);
    fjint.open(filename_jint);

    // for (int i = 0; i < N; i++){
    //   for (int j = 0; j < N; j++){
    //     double r = ri + dr*i;
    //     double theta_r = ti + dt*j;
    //     params_t params = {r,theta_r,ke,tau};

    //     results_t greens_functions = computeGreensFunction(params);
        
    //     std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
    //     std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
    //     std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);

    //     std::complex<double> Grad_r_diag (greens_functions.j_r_diag_real , greens_functions.j_r_diag_imag);
    //     std::complex<double> J_r_up_dn (greens_functions.j_r_up_dn_real , greens_functions.j_r_up_dn_imag);
    //     std::complex<double> J_r_dn_up (greens_functions.j_r_dn_up_real , greens_functions.j_r_dn_up_imag);

    //     std::complex<double> J_theta_diag (greens_functions.j_theta_diag_real , greens_functions.j_theta_diag_imag);
    //     std::complex<double> J_theta_up_dn (greens_functions.j_theta_up_dn_real , greens_functions.j_theta_up_dn_imag);
    //     std::complex<double> J_theta_dn_up (greens_functions.j_theta_dn_up_real , greens_functions.j_theta_dn_up_imag);

    //     std::complex<double> A_r_up_dn(cos(theta_r-tau),sin(theta_r+tau));
    //     std::complex<double> A_r_dn_up(cos(theta_r-tau),-sin(theta_r+tau));

    //     std::complex<double> A_t_up_dn(sin(theta_r-tau),cos(theta_r+tau));
    //     std::complex<double> A_t_dn_up(sin(theta_r-tau),-cos(theta_r+tau));

    //     double J_r = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(J_r_dn_up*G_up_dn - J_r_up_dn*G_dn_up) + k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
    //     double J_theta = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(J_theta_dn_up*G_up_dn - J_theta_up_dn*G_dn_up) + k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 



    //     double Grad_r = (-1.0/PI)*std::imag(J_r_diag);
    //     double Grad_theta = (-1.0/PI)*std::imag(J_theta_diag);

    //     //double Grad_r = pow(abs(J_r_diag),2);
    //     //double Grad_theta = pow(abs(J_theta_diag),2);

    //     fgrad << r << "," << theta_r << "," << Grad_r << "," << Grad_theta << std::endl;
    //     fcur << r << "," << theta_r << "," << J_r << "," << J_theta << std::endl;

    //   }
    // }

    // Circle integrations
    for (int i = 1; i < 50; i++){
      double r = float(i);

      double sum_jr = 0.0;
      double sum_jt = 0.0;
      for (int j = 0; j < N; j++){
        double theta_r = ti + dt*j;
        params_t params = {r,theta_r,ke,tau};

        results_t greens_functions = computeGreensFunction(params);
        
        std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
        // Inverted!!!
        std::complex<double> G_dn_up (greens_functions.up_dn_real , greens_functions.up_dn_imag);
        std::complex<double> G_up_dn (greens_functions.dn_up_real , greens_functions.dn_up_imag);
        // Inverted!!!
        std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
        std::complex<double> Grad_r_dn_up (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
        std::complex<double> Grad_r_up_dn (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
        // Inverted!!!
        std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
        std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
        std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


        std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
        std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

        std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
        std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);


        double J_r_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
        double J_theta_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

        double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
        double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


        double J_r = J_r_grad + J_r_so;
        double J_theta = J_theta_grad + J_theta_so;

        sum_jr += J_r;
        sum_jt += J_theta;

       }
       fjint << r << "," << sum_jr + sum_jt << std::endl;

    }
    fjint.close();

    // for(int i = 0; i < N; i++){
    //   for(int j = 0; j < N; j++){
    //     double x = xi + dx*i;
    //     double y = yi + dy*j;

 for(int i = 0; i < N; i++){
      for(int j = 0; j < N; j++){
        double x = xi + dx*i;
        double y = yi + dy*j;
        double r = sqrt(pow(x,2)+pow(y,2));
        double theta_r = atan2(y,x);

        // double theta_r = ti + dt*j
        // double theta_r = 0.0;
        // double r = 1.0;
        // double r = xi + dx*i;
        params_t params = {r,theta_r,ke,tau};

        results_t greens_functions = computeGreensFunction(params);
        std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
        // Inverted!!!
        std::complex<double> G_dn_up (greens_functions.up_dn_real , greens_functions.up_dn_imag);
        std::complex<double> G_up_dn (greens_functions.dn_up_real , greens_functions.dn_up_imag);
        // Inverted!!!
        std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
        std::complex<double> Grad_r_dn_up (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
        std::complex<double> Grad_r_up_dn (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
        // Inverted!!!
        std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
        std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
        std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


        std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
        std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

        std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
        std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);
        // std::complex<double> A_r_up_dn(cos(theta_r-tau),sin(theta_r+tau));
        // std::complex<double> A_r_dn_up(cos(theta_r-tau),-sin(theta_r+tau));

        // std::complex<double> A_t_up_dn(sin(theta_r-tau),-cos(theta_r+tau));
        // std::complex<double> A_t_dn_up(sin(theta_r-tau),-cos(theta_r+tau));

        double J_r_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
        double J_theta_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

        double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
        double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


        double J_r = J_r_grad + J_r_so;
        double J_theta = J_theta_grad + J_theta_so;

        // double J_r = J_r_grad ;
        // double J_theta = J_theta_grad;
        std::cout << "\nr = " << r << ", theta_r = " << theta_r << std::endl;
        std::cout <<"j_r_grad = (" << J_r_grad  << "," << J_theta_grad << ") , j_r_so = (" << J_r_so  <<  "," << J_theta_so << ")" <<  std::endl;
        // std::cout << "A_r_up_dn = " << A_r_up_dn << ", A_r_dn_up: " << A_r_dn_up << "A_t_up_dn = " << A_t_up_dn << ", A_t_dn_up: " << A_t_dn_up << std::endl;
        // std::cout << "DeltaG_up_dn = " << -(t_up-t_dn)*G_diag*G_up_dn << std::endl;
        // std::cout << "DeltaG_dn_up = " << (t_up-t_dn)*G_diag*G_dn_up << std::endl;

        // std::cout  << "G_up_dn: " << G_up_dn << std::endl;
        // std::cout  << "G_dn_up: " << G_dn_up << std::endl;

        //std::cout << x << "," << y << "," << Jx << "," << Jy << std::endl;
        //fcur << x << "," << y << "," << Jx << "," << Jy << std::endl;
        fcur << r << "," << theta_r << "," << J_r << "," << J_theta << std::endl;
        // fldos << x << "," << y << "," << (-1.0/PI)*std::imag((t_up-t_dn)*(pow(G_diag,2) - G_dn_up*G_up_dn)) << std::endl;

        // fgrad << x << "," << y << "," << Grad_x << "," << Grad_y << std::endl;
       
      } 
    }
    fcur.close();
    fldos.close();
    // fgrad.close();
  }


  gsl_integration_workspace_free (w);
  return 0;
}