#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <complex.h>
#include <gsl/gsl_cblas.h>
#include "diffClock.h"

void printMatrix(gsl_matrix *A);

void SVD(gsl_matrix* A, gsl_matrix* U, gsl_matrix* D, gsl_matrix* V);

int main()
{
int n=5;

gsl_matrix* A = gsl_matrix_alloc(n,n);
gsl_matrix* U = gsl_matrix_alloc(n,n);
gsl_matrix* D = gsl_matrix_alloc(n,n);
gsl_matrix* V = gsl_matrix_alloc(n,n);
gsl_matrix* UD = gsl_matrix_alloc(n,n);
gsl_matrix* UUT = gsl_matrix_alloc(n,n);
gsl_matrix* VVT = gsl_matrix_alloc(n,n);



gsl_matrix* hopefully_A = gsl_matrix_alloc(n,n);




for(int i=0; i<n; i++){
  	for(int j=0; j<n;j++){
      double number = ((double) rand()/RAND_MAX);
  		gsl_matrix_set(A,i,j,number);
  	}
  }

printf("- - - - Exam 2021 - - - - \n\n");
printf("In this exercise we are to demonstrate making the Singular Value decomposition\n");
printf("using a real square matrix A.\n\n\n");
printf("We start  by defining the matrix A as:\n");
  printMatrix(A);
printf("We now apply the SVD algorithm and attain the following matrices\n\n");
gsl_matrix_set_identity(U);
gsl_matrix_set_identity(V);
gsl_matrix_set_identity(D);
SVD(A,U,D,V);

printf("The Matrix A'->AJ is given as:\n\n");
printMatrix(A);
printf("\n\n");
printf("The SVD is given as A=UDV^T, so we have the following matrices.\n\n");
printf("The matrix U is given as:\n\n");
printMatrix(U);
printf("Which we note is an orthogonal matrix as demanded, since UU^T is:\n\n");
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,U,U,0.0,UUT);
printMatrix(UUT);

printf("The matrix D is given as\n\n");
printMatrix(D);
printf("which we note is diagonal with non-negative elements as demanded\n\n");
printf("matrix V is given as\n\n");

printf("\n\n");
printMatrix(V);
printf("\n\n");
printf("which we note is orthogonal as demanded since VV^T is:\n\n");
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,V,V,0.0,VVT);
printMatrix(VVT);
printf("\n\n");
printf("At last we can check that UDV^T is given as:\n\n");

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,U,D,0.0,UD);
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,UD,V,0.0,hopefully_A);
printMatrix(hopefully_A);
printf("which is clearly the matrix A as was given in the beginning of the exercise.\n");
printf("Therefore we conclude that the SVD algorithm works.\n");
printf("On the plot timePlot.png one can see a comparison with our SVD compared to\n");
printf("the GSL version. As can be clearly seen, our version is far slower.");

gsl_matrix_free(VVT);
gsl_matrix_free(UUT);
gsl_matrix_free(UD);
gsl_matrix_free(A);
gsl_matrix_free(U);
gsl_matrix_free(D);
gsl_matrix_free(V);
gsl_matrix_free(hopefully_A);
//comparison
FILE* my_out_stream = fopen("timeData.txt","w");



for(int n = 1; n<15;n++){
gsl_matrix* A_time     = gsl_matrix_alloc(n,n);
gsl_matrix* U_time     = gsl_matrix_alloc(n,n);
gsl_matrix* V_time     = gsl_matrix_alloc(n,n);
gsl_matrix* D_time     = gsl_matrix_alloc(n,n);
gsl_vector* S          = gsl_vector_alloc(n);


clock_t start, end;
double cpu_time_used;
double time;

for(int i=0; i<n; i++){
	for(int j=0; j<n;j++){
		gsl_matrix_set(A_time,i,j,((double) rand()/RAND_MAX));
	}
}

start = clock();
SVD(A_time,U_time,V_time,D_time);
end = clock();

cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

start = clock();
gsl_linalg_SV_decomp(A_time, V_time, S, S);
end = clock();
    time = ((double) (end - start)) / CLOCKS_PER_SEC;

		 fprintf(my_out_stream,"%d\t %g\t %g\n",n,cpu_time_used, time);


gsl_matrix_free(A_time);
gsl_matrix_free(U_time);
gsl_matrix_free(D_time);
gsl_matrix_free(V_time);

}

fclose(my_out_stream);


return 0;
}
