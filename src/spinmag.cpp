#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <string>
#include <vector>

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


  results_t res { result_diag_real, result_diag_imag, result_up_dn_real, result_up_dn_imag, result_dn_up_real, result_dn_up_imag};
  return res;
}


int main (int argc, char* argv[])
{
  
  int N = 200;
  gsl_set_error_handler_off();
  //gsl_set_error_handler (&error_handler);

  double xi = -2*(2*M_PI);
  double xf = 2*(2*M_PI);
  std::vector<double> xs = linspace(xi,xf,N);


  double yi = -2*(2*M_PI);
  double yf = 2*(2*M_PI);
  std::vector<double> ys = linspace(yi,yf,N);


  double ri = 0.01*(2*M_PI);
  double rf = 15*(2*M_PI);
  std::vector<double> rs = linspace(ri,rf,N);

  double ti = 0.0;
  double tf = 2*M_PI;
  std::vector<double> thetas = linspace(ti,tf,N);

  // std::cout << "thetas size: " << thetas.size() << std::endl;
  // for (double d : thetas)
  //   std::cout << d << " ";
  // std::cout << std::endl;


  // Model parameters
  double kei = 0.001;
  double kef = 1.2; // in units of k0
  std::vector<double> kes = linspace(kei,kef,50);

  //double tau = 1.31*k0;
  double t_up = 0.0;
  double t_dn = -(1.0/(2*PI));

  const std::complex<double> ii(0, 1);

  // double taus[] = {0.785, 0.871, 0.971, 1.15, 1.31, -1.31, PI/2.0};
  // double taus[] = {PI/4.0,0.871, 0.971, 1.15, 1.31, -1.31, PI/2.0};
  double taus[] = {0.0,M_PI/8.0,3.0*M_PI/8.0,M_PI/2.0};
  std::cout << "Calculating spin magnetic density..." << std::endl;

  for(int t = 0; t < 4; t++){
    double tau = taus[t];

    std::cout << "\n------------[Parameters]------------" << std::endl;
    std::cout << "tau : " << tau << std::endl;
    std::cout << "grid :  [" << N << "," << N << "] points"<< std::endl;
    std::cout << "x :  [" << xi << "," << xf << "] (in units of 1/k0)"<< std::endl;
    std::cout << "y :  [" << yi << "," << yf << "] (in units of 1/k0) "<< std::endl;
    std::cout << "------------------------------------" << std::endl;
   

    // Circle integrations

    // std::string filename_jint = "magmom_lint_tau_"+std::to_string(tau)+".dat";
    // std::ofstream fjint;
    // fjint.open(filename_jint);
    // double sum_magmon = 0.0;
    
    // for (double r: rs){
    //   for(double theta_r: thetas){
    //     params_t params = {r,theta_r,ke,tau};

    //     results_t greens_functions = computeGreensFunction(params);
        
    //     std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
    //     std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
    //     std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
    //     std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
    //     std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
    //     std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
    //     std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
    //     std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
    //     std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


    //     std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
    //     std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

    //     std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
    //     std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);


    //     double J_r_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
    //     double J_theta_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

    //     double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
    //     double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


    //     double J_r = J_r_grad + J_r_so;
    //     double J_theta = J_theta_grad + J_theta_so;

    //     sum_magmon += r*J_theta - theta_r*J_r;

    //    }
    //    fjint << r << "," << sum_magmon << std::endl;

    // }
    // fjint.close();


   // Energy integration


    std::string filename_x = "spinmag_x_energy_integrated_ke_1.2_tau_"+std::to_string(tau)+".dat";
    std::ofstream fmag_x;
    fmag_x.open(filename_x);


    std::string filename_y = "spinmag_y_energy_integrated_ke_1.2_tau_"+std::to_string(tau)+".dat";
    std::ofstream fmag_y;
    fmag_y.open(filename_y);

    std::string filename_z = "spinmag_z_energy_integrated_ke_1.2_tau_"+std::to_string(tau)+".dat";
    std::ofstream fmag_z;
    fmag_z.open(filename_z);


    for( double x: xs){
      for (double y: ys){
          double r = sqrt(pow(x,2)+pow(y,2));
          double theta_r = atan2(y,x);

          double mag_x = 0.0;
          double mag_y = 0.0;
          double mag_z = 0.0;

          for (double ke: kes){
            params_t params = {r,theta_r,ke,tau};

            results_t greens_functions = computeGreensFunction(params);
            std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
            std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
            std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
           
            mag_x += (-1.0/PI)*std::imag((t_up-t_dn)*G_diag*(G_dn_up - G_up_dn));
            mag_y += (-1.0/PI)*std::imag(-ii*(t_up-t_dn)*G_diag*(G_dn_up + G_up_dn));
            mag_z += (-1.0/PI)*std::imag((t_up-t_dn)*(pow(G_diag,2) + G_dn_up*G_up_dn));

          }       
        fmag_x << x << "," << y << "," << mag_x << std::endl;
        fmag_y << x << "," << y << "," << mag_y << std::endl;
        fmag_z << x << "," << y << "," << mag_z << std::endl;

        } 
      }
    fmag_x.close();
    fmag_y.close();
    fmag_z.close();

  }



  //Surface integral of the orbital magnetization
  // std::vector<double> tau_vec = linspace(0.0,M_PI/2.0,20);
  // std::vector<double> tau_vec;
  // tau_vec.push_back(0.330694);
  // std::string filename = "integrated_magmom_tau_0.330694.dat";
  // std::ofstream fmag;
  // fmag.open(filename);

  // for( double tau: tau_vec){
  //   double magmon = 0.0;
  //   for( double x: xs){
  //     for (double y: ys){
  //         double r = sqrt(pow(x,2)+pow(y,2));
  //         double theta_r = atan2(y,x);

  //         params_t params = {r,theta_r,ke,tau};

  //         results_t greens_functions = computeGreensFunction(params);
  //         std::complex<double> G_diag (greens_functions.diag_real , greens_functions.diag_imag);
  //         std::complex<double> G_up_dn (greens_functions.up_dn_real , greens_functions.up_dn_imag);
  //         std::complex<double> G_dn_up (greens_functions.dn_up_real , greens_functions.dn_up_imag);
  //         std::complex<double> Grad_r_diag (greens_functions.grad_r_diag_real , greens_functions.grad_r_diag_imag);
  //         std::complex<double> Grad_r_up_dn (greens_functions.grad_r_up_dn_real , greens_functions.grad_r_up_dn_imag);
  //         std::complex<double> Grad_r_dn_up (greens_functions.grad_r_dn_up_real , greens_functions.grad_r_dn_up_imag);
  //         std::complex<double> Grad_theta_diag (greens_functions.grad_theta_diag_real , greens_functions.grad_theta_diag_imag);
  //         std::complex<double> Grad_theta_up_dn (greens_functions.grad_theta_up_dn_real , greens_functions.grad_theta_up_dn_imag);
  //         std::complex<double> Grad_theta_dn_up (greens_functions.grad_theta_dn_up_real , greens_functions.grad_theta_dn_up_imag);


  //         std::complex<double> A_r_up_dn = cos(theta_r)*std::exp(ii*tau) + ii*sin(theta_r)*std::exp(-ii*tau);
  //         std::complex<double> A_r_dn_up = cos(theta_r)*std::exp(-ii*tau) - ii*sin(theta_r)*std::exp(ii*tau);

  //         std::complex<double> A_t_up_dn = ii*cos(theta_r)*std::exp(-ii*tau)- sin(theta_r)*std::exp(ii*tau);
  //         std::complex<double> A_t_dn_up = -ii*cos(theta_r)*std::exp(ii*tau)- sin(theta_r)*std::exp(-ii*tau);

  //         double J_r_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_r_dn_up*G_up_dn - Grad_r_up_dn*G_dn_up) );
  //         double J_theta_grad = (-1.0/PI)*std::imag(std::complex<double>(0,1)*(t_up-t_dn)*(Grad_theta_dn_up*G_up_dn - Grad_theta_up_dn*G_dn_up) ); 

  //         double J_r_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_r_up_dn*G_dn_up - A_r_dn_up*G_up_dn) );
  //         double J_theta_so = (-1.0/PI)*std::imag(k0*(t_up-t_dn)*G_diag*(A_t_up_dn*G_dn_up - A_t_dn_up*G_up_dn) ); 


  //         double J_r = J_r_grad + J_r_so;
  //         double J_theta = J_theta_grad + J_theta_so;

  //         double J_x = J_r*cos(theta_r) - J_theta*sin(theta_r);
  //         double J_y = J_r*sin(theta_r) + J_theta*cos(theta_r);
  //         magmon += x*J_y - y*J_x;
         
  //       }
  //     }
  //    fmag << tau << "," << magmon << std::endl;
  //   }
  // fmag.close();
         


  gsl_integration_workspace_free (w);
  return 0;
}