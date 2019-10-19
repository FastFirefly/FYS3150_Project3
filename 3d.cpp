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

//THE INTEGRAND FUNCTION IN POLAR COORDINATES FOR THE IMPORTANCE SAMPLING MONTE CARLO
double func_polar_mc(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
}

inline double ran(){
	//return ((double) generator())/2147483647;
	return ((double) rand()) / RAND_MAX;
}

//Monte Carlo with polar coordinates and importance sampling
void Importance_MonteCarlo_Polar(int n, double &integral, double &var, double &std){
	double * x = new double [n];
	double r1, r2, t1, t2, p1, p2, f,rr1,rr2;
	double mc = 0.0;
	double sig = 0.0;
	double jacob = 4*pow(M_PI,4)/16;
	int i;

	#pragma omp parallel for reduction(+:mc)  private (i, r1, r2, t1, t2, p1, p2, rr1, rr2, f)
	for (i = 0; i < n; i++){
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
	mc = mc/((double) n);
	#pragma omp parallel for reduction(+:sig)  private (i)
	for (i = 0; i < n; i++){
		sig += (x[i] - mc)*(x[i] - mc); 
	}
	var = sig*jacob/((double) n );
	std = sqrt(var)/sqrt(n);
	integral = mc*jacob;
	delete [] x;
}

int main(int argc, char *argv[]) {
    printf("EXACT RESULT:\t%.16f\t\n", 5*M_PI*M_PI/256);

    int n_mc = atoi(argv[1]);
    srand(time(NULL));
	generator.seed(time(NULL));

	double polar_imp_mc, polar_imp_var, polar_imp_std;

	// double start_time = omp_get_wtime();
	//printf("Polar MC IS: \t%.8f\t%.8f\n", polar_imp_mc, polar_imp_var);
	// double time = omp_get_wtime() - start_time;

	Importance_MonteCarlo_Polar(n_mc, polar_imp_mc, polar_imp_var, polar_imp_std);
		//finish1 = clock();
	printf("%.16f		%.16f		%.16f\n", polar_imp_mc, polar_imp_var, polar_imp_std);	
/*	for (int i=0; i<10; i++) {
		//start1 = clock();
		Importance_MonteCarlo_Polar(n_mc, polar_imp_mc, polar_imp_var, polar_imp_std);
		//finish1 = clock();
		printf("%.16f		%.16f		%.16f\n", polar_imp_mc, polar_imp_var, polar_imp_std);	
	}*/
}