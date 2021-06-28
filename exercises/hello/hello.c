#include "stdio.h" // for html display only, you should use angular brackets
#include "math.h"  // for html display only, you should use angular brackets
#include <complex.h>




int main(void){
	double r=gamma(5);
	printf("we calculate, gamma(5)=%g\n\n",r);
	double q=j1(0.5);
        printf("we calculate, bessel(0.5)=%g\n\n",q);

 double pi = 4 * atan(1.0);
 double complex sqrtm2 = csqrt(-2);
	 printf("we calculate sqrt(-2)=%f + %f * i\n\n", creal(sqrtm2), cimag(sqrtm2));

           double complex z = cexp(I * pi);
           printf("we calculate exp(i*pi)=%f + %f * i\n\n", creal(z), cimag(z));
					 double complex zz = cexp(I);
           printf("we calculate exp(i*pi)=%f + %f * i\n\n", creal(zz), cimag(zz));

					 	double complex ie = cpow(I,exp(1));
					 	                   printf("we calculate i^e=%f + %f * i\n\n", creal(ie), cimag(ie));

	double complex ipi = cpow(I,I);
	                   printf("we calculate i^i=%f + %f * i\n\n", creal(ipi), cimag(ipi));
   	   return 0;


}
