#include<stdio.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

static int ncalls;

void Rosenbrock_grad(gsl_vector* r, gsl_vector* Rr){
	ncalls++;
	double x = gsl_vector_get(r,0);
	double y = gsl_vector_get(r,1);
	double Rx = 2*(x-1)+400*(x*x-y)*x;
	double Ry = 200*(y-x*x);
	gsl_vector_set(Rr, 0, Rx);
	gsl_vector_set(Rr,1,Ry);
}

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);
int OdeDriver(void f(double,gsl_vector*,gsl_vector*),double a,double b,gsl_vector* ya,gsl_vector* yb,double h,double acc,double eps,char *program);
double e; // the energy

void schroedinger(double x, gsl_vector* y, gsl_vector* dydx){
gsl_vector_set(dydx,0,gsl_vector_get(y,1));
gsl_vector_set(dydx,1,-2*e*gsl_vector_get(y,0)-2.0/x*gsl_vector_get(y,0));
}



void schroedingerSolved(double x, gsl_vector* y, gsl_vector* dydx){

}


void solver(gsl_vector* x, gsl_vector* M){
	ncalls++;
	int n = 2;
	e = gsl_vector_get(x,0);
	gsl_vector* ya = gsl_vector_alloc(n);
	gsl_vector* yb = gsl_vector_alloc(n);
	double a = 1e-7;
	double b = 8;
	double h = 0.001;
	double acc = 1e-4;
	double eps = 1e-4;
	gsl_vector_set(ya,0,a-a*a);
	gsl_vector_set(ya,1,1-2*a);
char* path = "Hydrogen.txt";
	OdeDriver(&schroedinger, a, b, ya, yb, h, acc, eps, path);

	gsl_vector_set(M,0,gsl_vector_get(yb,0));

}







int main()
{
printf("- - - - Part A - - - - \n\n");
double eps = 0.01;

int d = 2;
gsl_vector* r = gsl_vector_alloc(d);
double x0 = -2;
double y0 = 8;
gsl_vector_set(r,0, x0);
   gsl_vector_set(r, 1, y0);

     ncalls = 0;
    newton(Rosenbrock_grad, r, eps);
    printf("We search at (x,y) = (%g,%g)\n", x0, y0);
    printf("The function is called %i times\n", ncalls);
    printf("The found extremum of Rosenbrock's valley function is at:\n");
    for(int i=0; i<d; i++){
	  printf("%g\n", gsl_vector_get(r,i));
       }
       gsl_vector_free(r);

printf("- - - - Part B - - - - \n\n");

//path = "Hydrogen.txt";
	gsl_vector* x = gsl_vector_alloc(1);
	double xx0 = -2;
	gsl_vector_set(x,0,xx0);
	ncalls =0;
	newton(solver, x, eps);

	printf("Lowest root of M(e)=0:\n");
	printf("Chosen rmax = %g\n", 8.0);
	printf("We start searching at %g\n", xx0);
	printf("Lowest found root = %g\n", gsl_vector_get(x,0));
	printf("Exact result : -1/2\n");
	printf("The result is plotted on the figure Hydrogen.png");
}
