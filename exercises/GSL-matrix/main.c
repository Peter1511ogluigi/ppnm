#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g ",gsl_vector_get(v,i));
	printf("\n");
}

int main(){
	int n=3;
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	gsl_vector* b=gsl_vector_alloc(n);
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);
printf("1");
double a00=6.13, a01=-2.90, a02=5.86, a10=8.08, a11=-6.31, a12=-3.89, a20=-4.36, a21=1.00, a22=0.19;
	double b0=6.23, b1=5.37, b2=2.29;

	gsl_matrix_set(A,0,0,a00);
	gsl_matrix_set(A,0,1,a01);
	gsl_matrix_set(A,0,2,a02);
	gsl_matrix_set(A,1,0,a10);
	gsl_matrix_set(A,1,1,a11);
	gsl_matrix_set(A,1,2,a12);
	gsl_matrix_set(A,2,0,a20);
	gsl_matrix_set(A,2,1,a21);
	gsl_matrix_set(A,2,2,a22);

	 gsl_matrix_memcpy(Acopy,A);

	gsl_vector_set(b,0,b0);
	gsl_vector_set(b,1,b1);
	gsl_vector_set(b,2,b2);

	gsl_linalg_HH_solve(Acopy,b,x);

	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);
printf("- - - - Exercise A - - - -\n\n");
printf("I)   When solving the system for x using GSL we get\n");
	vector_print("    ",x);
printf("II)   We check if it is correct by multiplying the given matrix A with our solution x ");
	vector_print("This gives that b is ",y);
printf("which is the correct result");

gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(y);

return 0;
}
