#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include<float.h>

#include "ClasseAusiliaria.h"
int  main(){
int equal(double a, double b, double tau, double epsilon){
if  (abs(a-b) < tau){ printf("1\n");}
else
if (abs(a-b)/(abs(a)+abs(b)) < epsilon/2){ printf("0\n");}
}
}
