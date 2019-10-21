#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string>   
#include <sstream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <random>


#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10


using namespace std;

minstd_rand0 generator;

// The Integrand function in Polar Coordinates specifically for the Importance sampling Monte Carlo method
double func_polar_mc(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
}

// Find a random number between 0 and 1, uses generator to work with parallelization
inline double ran(){
	return ((double) generator())/2147483647;
}

// Main function that find the apporximation to the integral using Importance sampling in the Monte Carlo method
void Importance_MonteCarlo_Polar(int n, double &integral, double &var, double &std){
	double *x = new double [n];
	double r1, r2, rr1,rr2;
	double t1, t2;
	double p1, p2, f;

	double mc = 0.0;
	double sig = 0.0;

	double jacob = 4*pow(M_PI,4)/16;	// Define the jacobi determinant in polar coordinates
	
	// Begin parallelization for finding an approximation to the integral
	#pragma omp parallel 
	{
		#pragma omp for reduction(+:mc) private (r1, r2, t1, t2, p1, p2, rr1, rr2, f)
		for (int i = 0; i < n; i++){
			rr1=ran();
			rr2=ran();
			r1=-0.25*log(1-rr1);
			r2=-0.25*log(1-rr2);
			t1=ran()*M_PI;
			t2=ran()*M_PI;
			p1=ran()*2*M_PI;
			p2=ran()*2*M_PI;
			f=func_polar_mc(r1, t1, p1, r2, t2, p2);
			mc += f;
			x[i] = f;
		}
	}
	mc = mc/((double) n);
	// Begin parallelization for finding sigma
	#pragma omp parallel 
	{
		#pragma omp for reduction(+:sig)
		for (int i = 0; i < n; i++){
			sig += (x[i] - mc)*(x[i] - mc); 
		}
	}
	var = sig*jacob/((double) n );
	std = sqrt(var)/sqrt(n);		// Calculate the Standard Deviation
	integral = mc*jacob;			// Multiply with Jacobi determinant
	delete [] x;
}

int main(int argc, char *argv[]) {
    int n_mc = atoi(argv[1]);

	generator.seed(time(NULL));		// Set seed for random number generation

	printf("EXACT RESULT:\t%.16f\t\n", 5*M_PI*M_PI/256);

	double polar_imp_mc, polar_imp_var, polar_imp_std;
	// Time measurements for OpenMP version
	double start_time = omp_get_wtime();
	Importance_MonteCarlo_Polar(n_mc, polar_imp_mc, polar_imp_var, polar_imp_std);
	double time = omp_get_wtime() - start_time;
	printf("Time: %f 		 %.16f		%.16f\n", time, polar_imp_mc, polar_imp_std);

	// Time measurements for serial version
	/*clock_t start1, start2, finish1, finish2;  //  declare start and final time for each exponent to test the time of the algorithms
	start1 = clock();
	Importance_MonteCarlo_Polar(n_mc, polar_imp_mc, polar_imp_var, polar_imp_std);
	finish1 = clock();
	printf("Time: %f 		 %.16f		%.16f\n", ((finish1 - start1)/(double) CLOCKS_PER_SEC ), polar_imp_mc, polar_imp_std);*/
}