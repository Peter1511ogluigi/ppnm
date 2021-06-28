#include"stdio.h"
//#include"gsl_matrix.h"
#include<gsl/gsl_matrix.h>


int print_half_00(gsl_matrix* m)
{
	double half = 0.5; //changed from 1/2
 int status = printf( "half m_{00} = %lg\n", gsl_matrix_get(m,0,0)*half );
	//gsl_matrix_free(m);
	return 0;
}

int main(void)
{
	gsl_matrix* m = gsl_matrix_alloc(1,1); //dimensions must be strictly greater than zero
	gsl_matrix_set(m,0,0,66);
	printf("half m_{00} (should be 33):\n");
	int status = print_half_00(m);
	if(status>0){
		printf("status=%i : SOMETHING WENT TERRIBLY WRONG (status>0)\n",status);
	}
	else
{
		printf("status=%i : everything went just fine (status=0)\n",status);
	}
	gsl_matrix_free(m);
return 0;
}
