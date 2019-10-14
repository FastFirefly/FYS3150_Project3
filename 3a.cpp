#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <omp.h>
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
#include "lib.h"

#define   ZERO       1.0E-10

using namespace std;
using namespace arma;

double func_cart(double x1, double y1, double z1, double x2, double y2, double z2){
  if  ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) != 0)
    return exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))) 
              / sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  else 
    return 0;
}

double int_function(double x1, double y1, double z1, double x2, double y2, double z2) {
    double alpha = 2;
    // evaluate the different terms of the exponential
    double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    //printf("%f \n", exp(exp1+exp2)/deno);
    return exp(exp1+exp2)/deno;
} // end of function to evaluate


double gauss_legendre(double a, double b, int N) {
    double *x = new double [N];
    double *w = new double [N];
    //  set up the mesh points and weights
    gauleg(a,b,x,w, N);

    //  evaluate the integral with the Gauss-Legendre method
    //  Note that we initialize the sum
    double int_gauss = 0;
    //  six-double loops
    for (int i=0;i<N;i++){
	      for (int j = 0;j<N;j++){
	          for (int k = 0;k<N;k++){
    	          for (int l = 0;l<N;l++){
        	          for (int m = 0;m<N;m++){
        	              for (int s = 0;s<N;s++){
                            int_gauss += w[i]*w[j]*w[k]*w[l]*w[m]*w[s]
                            *func_cart(x[i],x[j],x[k],x[l],x[m],x[s]);
                            //printf("%f\n", int_gauss);
                     		}
                    }
                }
            }
        }    
	  }
    delete [] x;
    delete [] w;
    return int_gauss;
}


int main(int argc, char *argv[]) {
    double a = atoi(argv[1]), b = atoi(argv[2]);
    int N = atoi(argv[3]);
    printf("%0.10f \n", gauss_legendre(a, b, N));
    return 0;
}