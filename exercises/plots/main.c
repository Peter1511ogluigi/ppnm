
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double mygam(double x){
	if(x<0) return M_PI/sin(M_PI*x)/mygam(1-x);
	if(x<9) return mygam(x+1)/x;

	double lnGamma= x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
	return exp(lnGamma);
}


double myerf(double x){
	if(x<0) return -myerf(-x);
	double a[]={0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429};
	double t=1/(1+0.3275911*x);
	double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));
	return 1-sum*exp(-x*x);
}

int main() {
	double xmin=0.1, xmax=3;
	for(double x=xmin; x<=xmax; x+=1.0/8){
		printf("%10g %10g %10g %10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),myerf(x),tgamma(x),gsl_sf_gamma(x),mygam(x));
	}
return 0;
}
