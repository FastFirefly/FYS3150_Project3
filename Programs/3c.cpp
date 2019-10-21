#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
#include "lib.h"
#include <omp.h>
#include <random>
using namespace std;

// The integrand function in Cartesian Coordinates
double func_cart(double x1, double y1, double z1, double x2, double y2, double z2){
	if  ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) != 0)
		return exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))) 
		          / sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	else 
		return 0;
}

// Find a random number between 0 and 1, uniform distribution
inline double ran(){
	return ((double) rand()) / RAND_MAX;
}

// Main function that find the apporximation to the integral using the Brute force Monte Carlo method
void Brute_MonteCarlo(int n, double a, double b, double  &integral, double  &var, double  &std) { 	
	double * x = new double [n];
	double x1, x2, y1, y2, z1, z2, f;
	
	double mc = 0.0;
	double sig = 0.0;
	double ba = (b-a);

	double jacob = pow(ba,6);		// Define the jacobi determinant 

	// Begin finding an approximation to the integral
	for (int i = 0; i < n; i++){
		x1=ran()*ba + a;
		x2=ran()*ba + a;
		y1=ran()*ba + a;
		y2=ran()*ba + a;
		z1=ran()*ba + a;
		z2=ran()*ba + a;
		f=func_cart(x1, x2, y1, y2, z1, z2);
		mc += f;
		x[i] = f;
	}
	mc = mc/((double) n );
	// Begin finding sigma
	for (int i = 0; i < n; i++){
		sig += (x[i] - mc)*(x[i] - mc); 
	}
	var = sig*jacob/((double) n );
	std = sqrt(var)/sqrt(n);		// Calculate the Standard Deviation
	integral = mc*jacob;			// Multiply with Jacobi determinant
	delete [] x;
}

int main(int argc, char *argv[]) {
    double a = -atoi(argv[1]), b = atoi(argv[1]);
    int n_mc = atoi(argv[2]);
    
    // Unit test for func_cart 
    if (func_cart(2, 3, 1, 2, 3, 1) != 0) {
        printf("Unit test failed, func_cart did not return a 0 \n");
        printf("Exiting program...\n");
        exit(1);
    }

    srand(time(NULL));		// Set seed for random number generation

    printf("Exact =\t%.8f\t\n", 5*M_PI*M_PI/256);

    double brute_mc, brute_var, brute_std;
	clock_t start1, start2, finish1, finish2;  //  declare start and final time for each exponent to test the time of the algorithm
	start1 = clock();
	Brute_MonteCarlo(n_mc, a, b, brute_mc, brute_var, brute_std);
	finish1 = clock();
	printf("Brute Monte Carlo, Time: %f    Results: %.16f    Standard Deviation: %.16f\n", ((finish1 - start1)/(double) CLOCKS_PER_SEC ), brute_mc, brute_std);
	//clock_t start1, start2, finish1, finish2;  //  declare start and final time for each exponent to test the time of the algorithm
	/*for (int i=0; i<1000; i++) {
		//start1 = clock();
		Brute_MonteCarlo(n_mc, a, b, brute_mc, brute_var, brute_std);
		//finish1 = clock();
		//printf("%.16f		%.16f		%.16f\n", brute_polar_mc, brute_polar_var, brute_polar_std);
		printf("%.16f\n", brute_mc);
	}*/
}