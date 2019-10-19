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
#define   ZERO       1.0E-8
#define NUM_THREADS 8

using namespace std;

minstd_rand0 generator;

//THE INTEGRAND FUNCTION IN POLAR COORDINATES
double func_polar(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = exp(-4*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
}

inline double ran(){
	return ((double) rand()) / RAND_MAX;
}

//Monte Carlo with finite domain
void Brute_MonteCarlo_Polar(int n, double a, double  &integral, double  &var, double  &std){
	double * x = new double [n];
	double r1, r2, t1, t2, p1, p2, f;
	double mc = 0.0;
	double sig = 0.0;
	double jacob = a*a*4*pow(M_PI,4);
	int i;

	#pragma omp parallel for reduction(+:mc)  private (i, r1, r2, t1, t2, p1, p2, f)
	for (i = 0; i < n; i++){
		r1=ran()*a;
		r2=ran()*a;
		t1=ran()*M_PI;
		t2=ran()*M_PI;
		p1=ran()*2*M_PI;
		p2=ran()*2*M_PI;
		f=func_polar(r1, t1, p1, r2, t2, p2);
		mc += f;
		x[i] = f;
	}
	mc = mc/((double) n );
	#pragma omp parallel for reduction(+:sig)  private (i)
	for (i = 0; i < n; i++){
		sig += (x[i] - mc)*(x[i] - mc); 
	}
	var = sig*jacob/((double) n );
	std = sqrt(sig)/sqrt(n);
	integral = mc*jacob;
	delete [] x;
}

int main(int argc, char *argv[]) {
    printf("EXACT RESULT:\t%.8f\t\n", 5*M_PI*M_PI/256);

    int n_mc = atoi(argv[1]);
    srand(time(NULL));
	generator.seed(time(NULL));

	double brute_polar_mc, brute_polar_var, brute_polar_std;
	Brute_MonteCarlo_Polar(n_mc, 4, brute_polar_mc, brute_polar_var, brute_polar_std);
	printf("%.16f		%.16f		%.16f\n", brute_polar_mc, brute_polar_var, brute_polar_std);	
	//clock_t start1, start2, finish1, finish2;  //  declare start and final time for each exponent to test the time of the algorithm
	/*for (int i=0; i<10; i++) {
		//start1 = clock();
		Brute_MonteCarlo_Polar(n_mc, 4, brute_polar_mc, brute_polar_var, brute_polar_std);
		//finish1 = clock();
		printf("%.16f		%.16f		%.16f\n", brute_polar_mc, brute_polar_var, brute_polar_std);	
	}*/
}