#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main(int argc, char* arg[]){
int i=0;
	while(i < argc){i++;
double x = atof(arg[i]);
printf("x= %g,   sin(x)=%g,   cos(x)=%g\n",x,sin(x),cos(x));
}
return 0;
}

