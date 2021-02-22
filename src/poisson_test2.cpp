#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fftw3.h>
#include <iostream>
#include <vector>

using namespace std;
int main() {
int N1=200;
int N2=200;

double pi = 3.141592653589793;
double L1  = 20.0;
double dx = L1/(double)(N1-1);
double L2= 20.0;
double dy=L2/(double)(N2-1);

double invL1s=1.0/(L1*L1);
double invL2s=1.0/(L2*L2);

// std::vector<double> in1(N1*N2,0.0);
// std::vector<double> out2(N1*N2,0.0);

fftw_complex *in1;
fftw_complex *out2;

fftw_complex *in2;
fftw_complex *out1;

// std::vector<double> in2(N1*N2,0.0);
// std::vector<double> out1(N1*N2,0.0);

in1 = (fftw_complex*) fftw_malloc(N1*N2*sizeof(fftw_complex));
out2 = (fftw_complex*) fftw_malloc(N1*N2*sizeof(fftw_complex));

in2 = (fftw_complex*) fftw_malloc(N1*N2*sizeof(fftw_complex));
out1 = (fftw_complex*) fftw_malloc(N1*N2*sizeof(fftw_complex));


std::vector<double> X(N1,0.0);
std::vector<double> Y(N2,0.0);


fftw_plan p, q;
int i,j;
p = fftw_plan_dft_2d(N1,N2, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
q = fftw_plan_dft_2d(N1,N2, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);

double a =0.1;
double c = 1.0;


double R1 = -2.0;
double R2 = 5.0;
int l=-1;
for(i = 0;i <N1;i++){
    X[i] =-10.0+(double)i*dx ;           
    for(j = 0;j<N2;j++){
        l=l+1;
        Y[j] =-10.0+ (double)j*dy ;
        double r = sqrt(pow(X[i],2)+pow(Y[j],2));
        // printf("%f,%f - %f\n",X[i],Y[j],r);
        // if ( Y[j] > R1 &&  Y[j] < R2) {
        if (r < R2) {
            //in1[l][0] = pow(cos(pi*Y[j]/2.0),2);
            // in1[l][0] = (2.0*pi/4.0)*cos(pi*Y[j]/4.0)*sin(pi*Y[j]/4.0)*pow(cos(pi*X[i]/20.0),2);
            in1[l][0] =  a*exp(-(pow(X[i],2) + pow(Y[j],2))/(2*pow(c,2)));
            // in1[l][1] = 0.0; 
            in1[l][1] = 0.0;    
             // printf("%f,%f,1,1\n",X[i],Y[j]);
        }
        else {
            in1[l][0] = 0.0;
            in1[l][1] = 0.0;    
        }
        
        //in1[l]= sin(pi*X[i]) + sin(pi*Y[j]) ; // row major ordering
    }
}
// in1[int(N1/2) + N2*int(N2/2)][0] = 1.0;
fftw_execute(p);

l=-1;
for ( i = 0; i < N1; i++){   // f = g / ( kx² + ky² )  
    for( j = 0; j < N2; j++){
            l=l+1;
        double fact=0;

        // magnetization along x
        double num = 0;
        in2[l][0] =0;
        in2[l][1] =0;

        if(2*i<N1){
            fact=((double)i*i)*invL1s;
            // printf("num1 = %f\t i = %i\n",num,i);

        }else{
            fact=((double)(N1-i)*(N1-i))*invL1s;
            // printf("num2 = %f\n",num);
        }
        if(2*j<N2){
            fact+=((double)j*j)*invL2s;
            num = ((double)j)/L2;

        }else{
            fact+=((double)(N2-j)*(N2-j))*invL2s;
            num = ((double)(N2-j))/L1;

        }
        if(fact!=0){
            // printf("%f\n",num);
            in2[l][0] = out1[l][0]/fact;
            in2[l][1] = out1[l][1]/fact;
            // printf("in2 = (%f,%f)\n",);
        }else{
            in2[l][0] = 0.0;
            in2[l][1] = 0.0;
        }
    }
}

fftw_execute(q);
l=-1;

double erl1 = 0.;
for ( i = 0; i < N1-1; i++) {
    for( j = 0; j < N2-1; j++){
        l=l+1;

        //erl1 +=1.0/pi/pi*fabs( in1[l][0]-  0.25*out2[l][0]/((double)(N1-1))/((double)(N2-1))); 

        double bx = -0.25*((out2[j + N1*(i+1)][0]- out2[j + N1*i][0])/dx)/((double)(N1-1))/((double)(N2-1));
        double by = -0.25*((out2[j+1 + N1*i][0]- out2[j + N1*i][0])/dy)/((double)(N1-1))/((double)(N2-1));
        
        std::cout << -10.0+(double)i*dx << "," <<  -10.0+ (double)j*dy << "," <<  bx << "," <<  by << std::endl;
        // std::cout << -10.0+(double)i*dx << "," <<  -10.0+ (double)j*dy << "," <<  in1[j+N2*i][0] << "," <<  0.25*out2[j+N2*i][0]/((double)(N1-1))/((double)(N2-1)) << std::endl;

        // printf("%3d %10.5f %10.5f\n", l, in1[l],  0.25*out2[l]/((double)(N1-1))/((double)(N2-1)));

    }
}

// cout<<"error=" <<erl1 <<endl ;  
fftw_destroy_plan(p); fftw_destroy_plan(q); fftw_cleanup();

return 0;
}