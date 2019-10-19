#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include "time.h"


#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-8

using namespace std;

double gammln(double);

//THE INTEGRAND FUNCTION IN POLAR COORDINATES REDUCED FOR GAUSSIAN LAGUERRE
double func_polar_lag(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = exp(-3*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
}

void gauleg(double x1, double x2, double x[], double w[], int n) {
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()


//FUNCTIONS FOR COMPUTING THE LAGUERRE POLYNOMIALS WEIGHTS
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

void gaulag(double *x, double *w, int n, double alf) {
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

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

double Gauss_Laguerre(int n_lag, int n_leg) {
	double *xlag = new double [n_lag + 1];
	double *wlag = new double [n_lag + 1];
	double *xt = new double [n_leg];
	double *wt = new double [n_leg];
	double *xp = new double [n_leg];
	double *wp = new double [n_leg];
	gaulag(xlag,wlag,n_lag,0.0);
	gauleg(0,M_PI,xt,wt,n_leg);
	gauleg(0,2*M_PI,xp,wp,n_leg);
	double int_gauss = 0.0;

	int i,j,k,l,f,t;
	// #pragma omp parallel for reduction(+:int_gauss)  private (i,j,k,l,f,t)
	for (i = 1;  i <= n_lag;  i++){    //r1
		for (j = 0;  j <  n_leg;  j++){    //t1
			for (k = 0;  k <  n_leg;  k++){    //p1
				for (l = 1;  l <= n_lag;  l++){    //r2
					for (f = 0;  f <  n_leg;  f++){    //t2
						for (t = 0;  t <  n_leg;  t++){    //p2
							int_gauss += wlag[i]*wlag[l]*wt[j]*wp[k]*wt[f]*wp[t]
							*func_polar_lag(xlag[i],xt[j],xp[k],xlag[l],xt[f],xp[t]);
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
    double a = -atoi(argv[1]), b = atoi(argv[1]);
    int N = atoi(argv[2]);

	clock_t start1, start2, finish1, finish2;  //  declare start and final time for each exponent to test the time of the algorithm

	printf("EXACT RESULT:\t%.8f\t\n", 5*M_PI*M_PI/256);
	double int_leg;
	int n_lag;
	int n_leg;
	double gau;
	n_lag = 30;
    n_leg = 30;
    start1 = clock();
    gau = Gauss_Laguerre(n_lag, n_leg);
    finish1 = clock();
    printf("Time for Gauss-Laguerre N=%d: %f     Result=%0.10f\n", n_lag, ((finish1 - start1)/(double) CLOCKS_PER_SEC ), gau);

	
    /*for (int i=5; i<=20; i+=5) {
        n_lag = i;
        n_leg = i;
        start1 = clock();
        gau = Gauss_Laguerre(n_lag, n_leg);
        finish1 = clock();
        printf("Time for Gauss-Laguerre N=%d: %f     Result=%0.10f\n", i, ((finish1 - start1)/(double) CLOCKS_PER_SEC ), gau);

    }*/
}

#undef EPS
#undef MAXIT