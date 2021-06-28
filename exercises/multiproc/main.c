#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
//#define RND (double)rand()/RAND_MAX

int Number_of_points = 100000000;
unsigned int seed;

void* point_generator(void* arg){
double* Points_in_circle = (double*) arg;
for(int i = 1;i  < Number_of_points; i++){
        double x = rand_r(&seed);
        x/= RAND_MAX;
        double y = rand_r(&seed);
        y/= RAND_MAX;
        double z = x*x+y*y;
        if(z*z<=1)
        *Points_in_circle = *Points_in_circle + 1;
}
return NULL;
}
int main(){
	double c1,c2 = 0;
pthread_t thread1,thread2;
pthread_create(&thread1,NULL,point_generator,(void*)&c1);
pthread_create(&thread2,NULL,point_generator,(void*)&c2);
void* returnvalue = NULL;

pthread_join(thread1,returnvalue);
pthread_join(thread2,returnvalue);

double pi = 4*(c1+c2)/(2*Number_of_points);
printf("- - - - Exercise A - - - - \n\n");
printf("pi is estimated after 1e8 iterations to be %g\n",pi);

return 0;
}
