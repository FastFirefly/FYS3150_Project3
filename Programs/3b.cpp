#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
#include "lib.h"
#include <random>

#define EPS 3.0e-14
#define MAXIT 10

using namespace std;

double gammln(double);

// The integrand function in Polar Coordinates specifically for the Laguerre
double func_polar_lag(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = exp(-3*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
}

// Computes the weight for the Laguerre polynomial
double gammln( double xx) {
	double x,y,tmp,ser;
	
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// Calculates weight and abcissas for Gauss-Laguerre
void gaulag(double *x, double *w, int n) {
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	double alf = 0;
	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}

// Main function that find the apporximation to the integral using the Gauss-Laguerre method
double Gauss_Laguerre(int n_lag, int n_leg) {
	double *xlag = new double [n_lag + 1];
	double *xt = new double [n_leg];
	double *xp = new double [n_leg];
	
	double *wlag = new double [n_lag + 1];
	double *wt = new double [n_leg];
	double *wp = new double [n_leg];

	//  Set up the mesh points and weights
	gaulag(xlag,wlag,n_lag);
	gauleg(0,M_PI,xt,wt,n_leg);
	gauleg(0,2*M_PI,xp,wp,n_leg);
	double int_gauss = 0.0;

	//  Finds an apporximation to the integral with the Gauss-Laguerre method
	for (int r1 = 1;  r1 <= n_lag;  r1++){    //r1
		for (int t1 = 0;  t1 <  n_leg;  t1++){    //t1
			for (int p1 = 0;  p1 <  n_leg;  p1++){    //p1
				for (int r2 = 1;  r2 <= n_lag;  r2++){    //r2
					for (int t2 = 0;  t2 <  n_leg;  t2++){    //t2
						for (int p2 = 0;  p2 <  n_leg;  p2++){    //p2
							int_gauss += wlag[r1]*wlag[r2]*wt[t1]*wp[p1]*wt[t2]*wp[p2]
							*func_polar_lag(xlag[r1],xt[t1],xp[p1],xlag[r2],xt[t2],xp[p2]);
						}
					}
				}
			}
		}
	}
	delete [] xt;
	delete [] wt;
	delete [] xp;
	delete [] wp;
	delete [] xlag;
	delete [] wlag;
	return int_gauss;
}

int main(int argc, char *argv[]) {
    int N = atoi(argv[1]);

    // Unit test for func_polar_lag
    if (func_polar_lag(2, 3, 1, 2, 3, 1) != 0) {
        printf("Unit test failed, func_polar_lag did not return a 0 \n");
        printf("Exiting program...\n");
        exit(1);
    }

	clock_t start1, start2, finish1, finish2;  //  declare start and final time for each exponent to test the time of the algorithm

	printf("Exact =\t%.8f\t\n", 5*M_PI*M_PI/256);
	double int_leg;
	int n_lag;
	int n_leg;
	double gau;
	n_lag = N;
    n_leg = N;
    start1 = clock();
    gau = Gauss_Laguerre(n_lag, n_leg);
    finish1 = clock();
    printf("Time for Gauss-Laguerre N=%d: %f     Result=%0.10f\n", n_lag, ((finish1 - start1)/(double) CLOCKS_PER_SEC ), gau);

	
/*    for (int i=5; i<=40; i+=2) {
        n_lag = i;
        n_leg = i;
        start1 = clock();
        gau = Gauss_Laguerre(n_lag, n_leg);
        finish1 = clock();
        // printf("Time for Gauss-Laguerre N=%d: %f     Result=%0.10f\n", i, ((finish1 - start1)/(double) CLOCKS_PER_SEC ), gau);
        printf("%d      %.10f\n", i, gau);

    }*/
}

#undef EPS
#undef MAXIT