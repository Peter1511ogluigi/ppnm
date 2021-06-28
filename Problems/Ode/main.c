#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <complex.h>
#include <gsl/gsl_cblas.h>

#define FMT "%9.3f %9.3f "
#define TIME 6
#define M1 1.0
#define M2 1.0
#define M3 1.0
#define GRAVITY 1


void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

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

void odecos(double t,gsl_vector* y, gsl_vector* dydt){
	double y1=gsl_vector_get(y,1);
	double y2=-gsl_vector_get(y,0);
	gsl_vector_set(dydt,0,y1);
	gsl_vector_set(dydt,1,y2);
}

void epidemic2(double t,gsl_vector* y, gsl_vector* dydt){
	double N=6000000, Tr=10, Tc=0.5;
	double y0=-(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc);
	double y1=(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc)-gsl_vector_get(y,1)/Tr;
	double y2=gsl_vector_get(y,1)/Tr;
	gsl_vector_set(dydt,0,y0);				//dS/dt=-(I*S)/(N*T_c)
	gsl_vector_set(dydt,1,y1);		//dI/dt=(I*S)/(N*T_c)-I/T_r
	gsl_vector_set(dydt,2,y2);								//dR/dt=I/T_r
}

void epidemic5(double t,gsl_vector* y, gsl_vector* dydt){
	double N=6000000, Tr=10, Tc=1.5;
	double y0=-(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc);
	double y1=(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc)-gsl_vector_get(y,1)/Tr;
	double y2=gsl_vector_get(y,1)/Tr;
	gsl_vector_set(dydt,0,y0);				//dS/dt=-(I*S)/(N*T_c)
	gsl_vector_set(dydt,1,y1);		//dI/dt=(I*S)/(N*T_c)-I/T_r
	gsl_vector_set(dydt,2,y2);								//dR/dt=I/T_r
}

void epidemic10(double t,gsl_vector* y, gsl_vector* dydt){
	double N=6000000, Tr=10, Tc=3;
	double y0=-(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc);
	double y1=(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc)-gsl_vector_get(y,1)/Tr;
	double y2=gsl_vector_get(y,1)/Tr;
	gsl_vector_set(dydt,0,y0);				//dS/dt=-(I*S)/(N*T_c)
	gsl_vector_set(dydt,1,y1);		//dI/dt=(I*S)/(N*T_c)-I/T_r
	gsl_vector_set(dydt,2,y2);								//dR/dt=I/T_r
}

void f(double t,gsl_vector* u,gsl_vector*dudt);


int main(){
	printf("- - - - Exercise A - - - - \n\n");
	double a=0, b=M_PI*2, h=0.0001, acc=0.001, eps=0.001;
		gsl_vector* ya=gsl_vector_alloc(2);
		gsl_vector* yb=gsl_vector_alloc(2);
		gsl_vector* exact=gsl_vector_alloc(2);
		gsl_vector_set(ya,0,1); gsl_vector_set(ya,1,0);
		gsl_vector_set(exact,0,0); gsl_vector_set(exact,1,-1);

		int k=OdeDriver(odecos,a,b,ya,yb,h,acc,eps,"cos.txt");
printf("We solve for u''=-u, with solution cos(x).\n Therefore we have y(0)=1.\n");
vector_print("For this, the start vector is ",ya);
vector_print("After 2 pi we get the following vector",yb);
vector_print("The analytical answer is",exact);
printf("This took %i steps\n\n", k);
printf("the plot can be seen on cos.png.");

double t0=0, tdays=75, he=0.1, acce=0.0000001, epse=0.0000001;
double N=6000000;
double II=3200;
double R=226000;
double S=N-R-II;
gsl_vector* yae=gsl_vector_alloc(3);
gsl_vector* ybe=gsl_vector_alloc(3);
gsl_vector_set(yae,0,S);
gsl_vector_set(yae,1,II);
gsl_vector_set(yae,2,R);

int ke2=OdeDriver(epidemic2,t0,tdays,yae,ybe,he,acce,epse,"Epidemic2.txt");
gsl_vector_set(yae,0,S);
gsl_vector_set(yae,1,II);
gsl_vector_set(yae,2,R);
OdeDriver(epidemic5,t0,tdays,yae,ybe,he,acce,epse,"Epidemic5.txt");
gsl_vector_set(yae,0,S);
gsl_vector_set(yae,1,II);
gsl_vector_set(yae,2,R);

OdeDriver(epidemic10,t0,tdays,yae,ybe,he,acce,epse,"Epidemic10.txt");

printf("Now er try to use the SIR model to plot the evolution of corona in Denmark.\n");
printf("We choose the following parameters N=%g I=%g R=%g and investigate the development over %g days\n",N,II,R,tdays);
printf("The plot is shown on the figure Epidemic.png, which resembles the plot from the SIR wiki page.");
printf("It took %i steps to attain the plot with Tc=1/2\n",ke2);
printf("when comparing different Tc's we see that the larger Tc the more the solution is\n");
printf("shifted to the right and the number of infected decrease.\n");
printf("this makes sense since given the definition of Tr/Tc");
printf("\n\n\n");
printf("- - - - Exercise B - - - - \n\n");
printf("{ti,y(ti)} is stoed in the given path in the driver\n\n");
printf("- - - - Exercise C - - - - \n\n");
printf("The stable three body orbit can be seen solved on threeBody.png");

gsl_vector* u = gsl_vector_alloc(12);
gsl_vector* out = gsl_vector_alloc(12);

gsl_vector_set(u,0,-0.97000436);
gsl_vector_set(u,1,0.24308753);
gsl_vector_set(u,2,0);
gsl_vector_set(u,3,0);
gsl_vector_set(u,4,0.97000436);
gsl_vector_set(u,5,-0.24308753);
gsl_vector_set(u,6,0.4662036850);
gsl_vector_set(u,7,0.4323657300);
gsl_vector_set(u,8,-0.93240737);
gsl_vector_set(u,9,-0.86473146);
gsl_vector_set(u,10,0.4662036850);
gsl_vector_set(u,11,0.4323657300);


double aa=0,bb=TIME,hh=0.1,aacc=1e-3,eeps=1e-3;
OdeDriver(f,aa,bb,u,out, hh,aacc,eeps,"threeBodyProblem.txt");

gsl_vector_free(u);
gsl_vector_free(out);

	return 0;
}
