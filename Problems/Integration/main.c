
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>

//Defining functions from ai.c
double recursive_integrate(double f(double), double a, double b, double acc, double eps, double f2, double f3, int limit);
double integrate(double f(double), double a, double b, double acc, double eps);
double clenshaw_curtis(double f(double), double a, double b, double acc, double eps);
double CC24(double f(double), double a, double b, double acc, double eps, double f2, double f3, int nrec);
double integ_a_inf(double f(double),double a,double acc,double eps );
int calls;
double f(double x){ calls++; return exp(-x); };
double g(double x){ calls++; return 1/x/x; };




double g1(double x){
	calls++;
	return sqrt(x);
};
double g2(double x){
	calls++;
	return 4*sqrt(1-x*x);
}

double g2_GSL(double x, void* params)
{
	calls++;
	double alpha= *(double*) params;
	return alpha*4*sqrt(1-x*x);
}
int main(){
	//Interval from [0,1] with acc = 0.001 and eps = 0.001.
	double a = 0;
	double b = 1;
	double acc = 0.0001;
	double eps = 0.0001;
	calls = 0;
	double Q1 = integrate(g1,a,b,acc,eps);
	double exact = 2./3;
printf("- - - - Exercise A - - - -\n\n");

	printf("For f = sqrt(x) in interval %g to %g we get\n",a,b);
	printf("Q = %g\n",Q1);
	printf("with exact value being %g\n",exact);
	printf("We have used %d calls\n",calls);
	printf("with an estimated error = %g\n",acc+fabs(Q1)*eps);
	printf("compared to an actual error = %g\n\n\n",fabs(Q1-exact));


	//Printing results with f = 4*(1-x^2).
	calls = 0;
	exact = M_PI;
	double Q2 = integrate(g2,a,b,acc,eps);

	printf("For f = 4(1-x^2) in interval %g to %g we get \n",a,b);
	printf("Q = %g\n",Q2);
	printf("with exact value being %g\n",exact);
	printf("using %d calls \n",calls);
	printf("wih an estimated error = %g\n",acc + fabs(Q2)*eps);
	printf("compared to an actual error = %g\n",fabs(Q2-exact));



	//---B---
	calls = 0;

	double Q3 = clenshaw_curtis(g2,a,b,acc,eps);
	exact = 2;
	printf("-----------------------------\n\n");
	printf("- - - - Exercise B - - - - \n\n");
	printf("We implement CLENSHAW_CURTIS and calculate the integral of f = 4(1-x^2) in same interval\n");
	printf("We get\n");
	printf("Q = %g\n",Q3);
	printf("which we note is the correct result with greater accuracy\n");
	printf("this uses %d calls\n",calls);
	printf("which is more than the previous integrator\n");
	printf("We also get\n");
	printf("estimated error = %g\n",acc+fabs(Q3)*eps);
	printf("actual error = %g\n",fabs(Q3-exact));
	printf("So in short we have greater accuracy, a few more steps with a larger actual error\n\n");
	//compare gsl
double alpha = 1.0;
int limit = 10000;
double result;
double error;
gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	calls=0;
		gsl_function G;
		G.function = &g2_GSL;
		G.params = &alpha;
		gsl_integration_qags(&G, a, b, acc, eps, limit, w, &result, &error);
		printf("GSL value = %0.25f\nGSL number of calls = %d\n", result, calls);

printf("\n\n");
printf("- - - - Exercise C - - - -\n\n");

printf("we test some infinite integrals and see whether they converge:\n");
			double aa,QQ,eexact,aacc=0.0001,eeps=0.0001;

			aa=0; calls=0;
			QQ=integ_a_inf(f,aa,aacc,eeps);
			eexact=exp(-aa);
			printf("∫ exp(-x) from %g to INFINITY \n",aa);
			printf("              Q = %g\n",QQ);
			printf("          exact = %g\n",eexact);
			printf("          calls = %d\n",calls);
			printf("estimated error = %g\n",aacc+fabs(QQ)*eeps);
			printf("   actual error = %g\n",fabs(QQ-eexact));

			aa=0.5; calls=0;
			QQ=integ_a_inf(g,aa,aacc,eeps);
			exact=2;
			printf("\n∫ 1/x^2 from %g to INFINITY \n",aa);
			printf("              Q = %g\n",QQ);
			printf("          exact = %d\n",2);
			printf("          calls = %d\n",calls);
			printf("estimated error = %g\n",aacc+fabs(QQ)*eeps);
			printf("   actual error = %g\n",fabs(QQ-eexact));





	return 0 ;
}
