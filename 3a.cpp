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

//The integrand function in Cartesian Coordinates
double func_cart(double x1, double y1, double z1, double x2, double y2, double z2){
  if  ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) != 0)
    return  exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))) / 
            sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  else 
    return 0;
}

// Main function that find the apporximation to the integral using the Gauss-Legendre method
double gauss_legendre(double a, double b, int N) {
    double *x = new double [N];
    double *w = new double [N];

    //  Set up the mesh points and weights
    gauleg(a,b,x,w, N);

    //  Finds an apporximation to the integral with the Gauss-Legendre method
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
    // Takes in arguments for the command line
    double a = -atoi(argv[1]), b = atoi(argv[1]);
    int N = atoi(argv[2]);
    double gau;
    printf("Exact=\t%.8f\t\n", 5*M_PI*M_PI/256);
    clock_t start1, start2, finish1, finish2;   //  declare start and final time for each exponent to test the time of the algorithm
    start1 = clock();
    gau = gauss_legendre(a, b, N);
    finish1 = clock();
    printf("Time for Gauss-Legendre N=%d: %f     Result=%0.10f\n", N, ((finish1 - start1)/(double) CLOCKS_PER_SEC ), gau);  // Prints results
    return 0;
}