#include<math.h>
#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include<assert.h>


double F(double* p){

  #define R 0.9
	double x=p[0],y=p[1];
	if(x*x+y*y<R*R)return 1;
	else return 0;
	}

double f(double* x)
{return 1/(M_PI*M_PI*M_PI)*1/(1-cos(x[0])*cos(x[1])*cos(x[2]));}

double g(double* x)
{return 1/(exp(1-x[0]));}





complex plainmc(int dim,double f(double* x),double* a,double* b,int N);
complex quasimc(int d, double f(double*), double a[], double b[], int N);


int main()
{




double a[] = {0,0,0};
double b[] = {M_PI,M_PI,M_PI};
int N = 1e7;
complex value = plainmc(3,f,a,b,N);
complex value2 = plainmc(1,g,a,b,N);

printf("- - - - Exercise A - - - - \n\n");
printf("we check integral of 1/exp(1-x) from 0 to pi\n");
printf("numerical result %g +- i%g\n",creal(value2),cimag(value2));
printf("This is the same as found on wolfram alpha.\n\n");

printf("We now check the sought after integral from the exercise and get\n");
printf("numerical result %g +- i%g\n",creal(value),cimag(value));
printf("corresponding to the result given in the exercise.");

//Part B
printf("\n\n");
printf("- - - - Exercise B - - - -\n\n");
printf("In this exercise we need to 'Implement a multidimensional Monte-Carlo integrator\n");
printf(" that uses low-discrepancy sequences.'\n");
printf("and\n");
printf("Compare the scaling of the error with pseudo-random Monte-Carlo integrator.\n");
printf("The comparison can be seen on plot.png where it is clear that the error\n");
printf("becomes smaller with greater number of points.");


  FILE* error = fopen("data.txt","w");
  double bb[] = {M_PI,M_PI};
for(int i =1; i<100;i++){
int NN=1000*i;
complex result_p=plainmc(2,F,a,bb,NN);
complex result_q=quasimc(2,F,a,bb,NN);
double error_p = cimag(result_p);
double error_q = cimag(result_q);
fprintf(error,"%d\t%f\t%f\n",NN,error_p,error_q);
}
fclose(error);
  return 0;
}
