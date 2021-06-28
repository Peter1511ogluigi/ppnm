#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include<math.h>
#include<float.h>


int main(){
printf("- - - Exercise 1.1 - - -\n");
int i=1;				//Calculate INT_MAX from while loop
while(i+1>i){i++;}
printf("my max int from 	while loop 	is 		%i\n",i);

int n=0;				// Calculate INT_MAX from for loop
for(n=0;n<n+1;i++){n=n+1;
}
printf("my max int from  	for loop 	is 		%i\n",n);

do{printf("my max int from  	do while loop   is 		%i\n",n);}
while(n+1>n);

int m = INT_MAX;
printf("INT_MAX is given as 					%i\n",m);

printf("------------\n");
printf("- - - Exercise 1.2 - - -\n");
int ii=0;				//Calculate INT_MAX from while loop
while(-ii-1<ii){ii++;}
printf("my min int from 	while loop 	is 		%i\n",ii);


int nn=0;				// Calculate INT_MAX from for loop
for(nn=0;nn>nn-1;nn--){nn=nn-1;
}
printf("my min int from  	for loop 	is 		%i\n",nn);

do{printf("my min int from  	do while loop   is 		%i\n",nn);}
while(nn-1<nn);


int mm = INT_MIN;
printf("INT_MAX is given as 					%i\n",mm);


printf("------------\n");












printf("- - - Exercise 1.3 - - -\n");

float x_float=1;
	while(1+x_float!=1){x_float/=2;}x_float*=2;
	printf("float x for while loop is %f\n",x_float);
	double x_double=1;
	while(1+x_double!=1){x_double/=2;} x_double*=2;
	printf("double x for while loop is %g\n",x_double);
	long double x_longdouble=1;
	while(1+x_longdouble!=1){x_longdouble/=2;}x_longdouble*=2;
	printf("long double x for while loop is %Lg\n",x_longdouble);

	//for loop
	float e_float; double e_double; long double e_longdouble;
	for(e_float=1;1+e_float!=1;e_float/=2){} e_float*=2;
	for(e_double=1;1+e_double!=1;e_double/=2){} e_double*=2;
	for(e_longdouble=1;1+e_longdouble!=1;e_longdouble/=2){ }e_longdouble*=2;

	printf("float e for for loop is %f\n",e_float);
	printf("double e for for loop is %g\n",e_double);
	printf("long double e for for loop is %Lg\n",e_longdouble);
printf("\n\n");
	printf("float check %f\n",FLT_EPSILON);
	printf("double check %g\n",DBL_EPSILON);
	printf("long double check %Lg\n",LDBL_EPSILON);

	//while do
	float  y_f=1;double y_d=1; long double y_ld=1;
	do {y_f/=2;} while(1+y_f!=1); y_f*=2;
	do {y_d/=2;} while(1+y_d!=1); y_d*=2;
	do {y_ld/=2;} while(1+y_ld!=1); y_ld*=2;
	printf("float y for do while is %f\n",y_f);
	printf("double y for do while is %g\n",y_d);
	printf("long double y for do while is %Lg\n",y_ld);


/*



printf("- - - Exercise 2.1 - - - \n");
int max=INT_MAX/2;
float sum=0;

for(i=1;i<=max;i++){sum +=  1/(float) i;}
printf("The 'up' sum of floats gives f*%f\n",sum);


for(i=max;(i=0);i++){sum +=  1/(float) i;}
printf("The 'down' sum of floats gives f*%f\n",sum);
printf("These two sums should give the same since they are the same sum.\n The sum should diverge as max-> inf since this would the become the harmonic series.\n");
printf("- - - Exercise 3 - - - \n");

double a=1; double b=2; double epsilon = 0.1; double tau=1.2;
double equal(a,b,tau,epsilon);

*/
	return 0;


}
