#include<stdio.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

void printdata(double x,gsl_vector* y,FILE *path){
int n=y->size;
fprintf(path,"%6g ",x);
for(int i=0;i<n;i++)fprintf(path,"%6g ",gsl_vector_get(y,i));
fprintf(path,"\n");
}


void rkstep12(
	void f(double t,gsl_vector* y,gsl_vector* dydt),
	double t,              					// the current value of the variable
	gsl_vector* yt,   				        // the current value y(t) of the sought function
	double h,   					        // the step to be taken
	gsl_vector* yh,					        // output: y(t+h)
	gsl_vector* dy)
	{             				// output: error estimate

	int i, n=yt->size;
	gsl_vector* k0			= gsl_vector_alloc(n);
	gsl_vector* k12			= gsl_vector_alloc(n);
	gsl_vector* y_half	= gsl_vector_alloc(n);

	f(t,yt,k0);
	for(i=0;i<n;i++)
			{
					double	y_half_i=gsl_vector_get(yt,i)+gsl_vector_get(k0,i)*h/2;
					gsl_vector_set(y_half,i,y_half_i);
			}

	f(t+h/2,y_half,k12);
	for(i=0;i<n;i++)
			{
					double yh_i=gsl_vector_get(yt,i)+gsl_vector_get(k12,i)*h;
					double dy_i=(gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*h/2;
					gsl_vector_set(yh,i,yh_i);
					gsl_vector_set(dy,i,dy_i);
			}

gsl_vector_free(k0);
gsl_vector_free(k12);
gsl_vector_free(y_half);
}


int OdeDriver(
	void f(double,gsl_vector*,gsl_vector*), 		// right-hand-side of dy/dt=f(t,y)
	double a,			                     	// the start-point a
	double b,                     				// the end-point of the integration
	gsl_vector* ya,                    				// y(a)
	gsl_vector* yb,                     			// y(b) to be calculated
	double h,                     				// initial step-size
	double acc,                   				// absolute accuracy goal
	double eps,
	char *program){                    			// relative accuracy goal

	int k=0, n=ya->size;
	double err, norm_y, tol, xi=a;
	gsl_vector* dy=gsl_vector_alloc(n);
	gsl_vector_memcpy(yb,ya);
	FILE *path;
	path = fopen(program,"w");
	printdata(xi,yb,path);
	while(xi<b)
	{
		if(xi+h>b) h=b-xi;	//Tjekker om det oplyste første step er for stort
		rkstep12(f,xi,ya,h,yb,dy);
		err=gsl_blas_dnrm2(dy);		//Summer alle dy'er op og tage absolut værdi for at få et udtryk for error
		norm_y=gsl_blas_dnrm2(yb);	//Summer alle y værdier op da de er relevante for relative uncertainty(bruges i tolerance)
		tol=(norm_y*eps+acc)*sqrt(h/(b-a));
		if(err<tol){
			xi+=h;
			gsl_vector_memcpy(ya,yb);
			k++;
			printdata(xi,yb,path);}
		if(err>0) h*=pow(tol/err,0.25)*0.95 ; // Vi skal lige passe på ikke at dele med 0
		else h*=2;
	}


	fclose(path);
	gsl_vector_free(dy);
return k;
}
