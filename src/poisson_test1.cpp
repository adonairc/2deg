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
double L1  = 2.0;
double dx = L1/(double)(N1-1);
double L2= 2.0;
double dy=L2/(double)(N2-1);

double invL1s=1.0/(L1*L1);
double invL2s=1.0/(L2*L2);

std::vector<double> in1(N1*N2,0.0);
std::vector<double> in2(N1*N2,0.0);
std::vector<double> out1(N1*N2,0.0);
std::vector<double> out2(N1*N2,0.0);
std::vector<double> X(N1,0.0);
std::vector<double> Y(N2,0.0);


fftw_plan p, q;
int i,j;
p = fftw_plan_r2r_2d(N1,N2, in1.data(), out1.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
q = fftw_plan_r2r_2d(N1,N2, in2.data(), out2.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

int l=-1;
for(i = 0;i <N1;i++){
    X[i] =-1.0+(double)i*dx ;           
    for(j = 0;j<N2;j++){
        l=l+1;
        Y[j] =-1.0+ (double)j*dy ;
        in1[l] = 0.0;
        //in1[l]= sin(pi*X[i]) + sin(pi*Y[j]) ; // row major ordering
    }
}
in1[int(N1/2) + N2*int(N2/2)] = 1.0;
fftw_execute(p);

l=-1;
for ( i = 0; i < N1; i++){   // f = g / ( kx² + ky² )  
    for( j = 0; j < N2; j++){
            l=l+1;
        double fact=0;
        in2[l]=0;

        if(2*i<N1){
            fact=((double)i*i)*invL1s;;
        }else{
            fact=((double)(N1-i)*(N1-i))*invL1s;
        }
        if(2*j<N2){
            fact+=((double)j*j)*invL2s;
        }else{
            fact+=((double)(N2-j)*(N2-j))*invL2s;
        }
        if(fact!=0){
            in2[l] = out1[l]/fact;
        }else{
            in2[l] = 0.0;
        }
    }
}

fftw_execute(q);
l=-1;
double erl1 = 0.;
for ( i = 0; i < N1; i++) {
    for( j = 0; j < N2; j++){
        l=l+1;

        erl1 +=1.0/pi/pi*fabs( in1[l]-  0.25*out2[l]/((double)(N1-1))/((double)(N2-1))); 

        std::cout << -1.0+(double)i*dx << "," <<  -1.0+ (double)j*dy << "," <<  in1[l] << "," <<  0.25*out2[l]/((double)(N1-1))/((double)(N2-1)) << std::endl;

        // printf("%3d %10.5f %10.5f\n", l, in1[l],  0.25*out2[l]/((double)(N1-1))/((double)(N2-1)));

    }
}

// cout<<"error=" <<erl1 <<endl ;  
fftw_destroy_plan(p); fftw_destroy_plan(q); fftw_cleanup();

return 0;
}