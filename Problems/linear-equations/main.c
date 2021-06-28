#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <complex.h>
#include <gsl/gsl_cblas.h>
#include "diffClock.h"
  #include <time.h>
#define FMT "%7.3f"
#define RND (double)rand()/RAND_MAX


//-------- Print Matrix function -----------------//

void printMatrix(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0; c<A->size2;c++){
			printf(FMT,gsl_matrix_get(A,r,c));
		}
		printf("\n");
	}
}
//---------Print Vector function-----------------------//
void vector_print(char s[],gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++){
		printf("%10g ",gsl_vector_get(v,i));
	}
	printf("\n");
}

//----------------QR-decomposition----------------//
void GramSchmidtDecomposition(gsl_matrix* A, gsl_matrix* R){ //This function replaces A with Q and determine R and put into the argument R
	int m = A->size2; //size2=columns --> m=#columns in A
	for(int i=0;i<m;i++){
		gsl_vector_view e = gsl_matrix_column(A,i); //Seperates the i'th column in A
		double r = gsl_blas_dnrm2(&e.vector); //Norm of the i'th column in A
		gsl_matrix_set(R,i,i,r);
		gsl_vector_scale(&e.vector,1/r); //Normalize e_i
		for(int j=i+1;j<m;j++){
			gsl_vector_view q =gsl_matrix_column(A,j); //Seperates columns in A that is to the right of the i'th column
			double s=0;
			gsl_blas_ddot(&e.vector,&q.vector,&s); //Dotproduct of e^T and q, returning it to s
			gsl_blas_daxpy(-s,&e.vector,&q.vector); //-s*e+q-->q
			gsl_matrix_set(R,i,j,s);
			//gsl_matrix_set(R,j,i,0);
		}
	}
}
//--------------------------------------------------//



void QR_backsub(gsl_matrix *R, gsl_vector *x){ //leaves the result in x
	int m= R->size1;
	for(int i=m-1;i>=0;i--){
		double s=0;
		for(int k=i+1;k<m;k++){
			s+=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
		}
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-s)/gsl_matrix_get(R,i,i));
	}
}

void QR_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);  //Calculates Q^T*b and leaves it in x (because Rx=Q^T*b)
	QR_backsub(R,x); 											// Turns x into the solution of Ax=b
}

//--------Inverse Matrix -------------//
void QR_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B){
	int n = R->size1;
	gsl_vector *b = gsl_vector_calloc(n);
	gsl_vector *x = gsl_vector_calloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set_basis(b,i);
		QR_solve(Q,R,b,x);
		gsl_matrix_set_col(B,i,x);
	}
	gsl_vector_free(b);
	gsl_vector_free(x);
}


int main(int argc, char* argv[] ){


	int m = 5;
	int n = 10;

	gsl_matrix* A     = gsl_matrix_alloc(n, m);
	gsl_matrix* AA		= gsl_matrix_alloc(m, m);
	gsl_matrix* B_inv	= gsl_matrix_alloc(m, m);
	gsl_matrix* BB_inv= gsl_matrix_alloc(m, m);
	gsl_matrix* B_invB= gsl_matrix_alloc(m, m);
	gsl_matrix* R     = gsl_matrix_alloc(m, m);
	gsl_matrix* RR    = gsl_matrix_alloc(m, m);
	gsl_matrix* Q     = gsl_matrix_alloc(n, m);
	gsl_matrix* QQ    = gsl_matrix_alloc(m, m);
	gsl_matrix* QTQ	  = gsl_matrix_alloc(m, m);
	gsl_vector* b			= gsl_vector_alloc(m);
	gsl_vector* x			= gsl_vector_alloc(m);

//------ creating random matrix A-----------//
	for(int i=0; i< A->size1; i++){
			for(int j=0; j<A->size2; j++)
			{
			double Aij=RND;
			gsl_matrix_set(A,i,j,Aij);
			}
	}
//	gsl_matrix_memcpy(AA,A);            //saving the memory of A in A_save

	GramSchmidtDecomposition(A, R);					//using GS on A to get Q,R


//---------------------------------------------------//
printf("------------ Part A (1) --------------\n\n");
printf("The random generated matrix A is given as\n\n");
printMatrix(A);
printf("\n\n");
printf("The matrix R is given as\n\n");
printMatrix(R);
printf("\n\n");
printf("The matrix Q is given as\n\n");
printMatrix(A);
printf("\n\n");
printf("The matrix product QTQ is given as\n\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,QTQ);
printMatrix(QTQ);
printf("\n\n");
printf("We check QR=A, since QR is:\n\n");
printMatrix(A);
printf("\n\n");

//gsl_matrix *AA = gsl_matrix_calloc(m,m);
//gsl_matrix *QQ = gsl_matrix_alloc(m,m);
for(int i=0; i<m; i++){
	for(int j=0; j<m;j++){
		gsl_matrix_set(AA,i,j,((double) rand()/RAND_MAX));
	}
}
printf("------------ Part A (2) --------------\n\n");


gsl_matrix_memcpy(QQ,AA);
printf("Our random-generated matrix, A, is:\n\n");
printMatrix(AA);
printf("\n\n");

for(int i=0; i<m; i++){
	gsl_vector_set(b,i,((double) rand()/RAND_MAX));
}
printf("Our random-generated vector, b, is:\n\n");
gsl_vector_fprintf(stdout,b,FMT);
printf("We now wanna solve Ax=b for x.\n\n");
GramSchmidtDecomposition(QQ,RR);
QR_solve(QQ,RR,b,x);
printf("Using our implementation we get x to be:\n\n");
gsl_vector_fprintf(stdout,x,FMT);
gsl_blas_dgemv(CblasNoTrans,1.0,AA,x,0.0,b);
printf("Using GSL's implementation we get that Ax is equal to:\n\n");
gsl_vector_fprintf(stdout,b,FMT);
printf("which equals b, so our implementation works.\n\n");


printf("\n\n");
printf("-------------- Part B --------------\n\n");
printf("the inverse of A(QR) is :\n\n");
QR_inverse(QQ,RR,B_inv);
printMatrix(B_inv);
printf("\n\n");
printf("We check that B^-1B becomes unity:\n\n");
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,AA,B_inv,0.0,BB_inv);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B_inv,AA,0.0,B_invB);
printMatrix(BB_inv);
printf("\n\n");
printf("We check that BB^_1 becomes unity");
printf("\n\n");
printMatrix(B_invB);


FILE* my_out_stream = fopen("timeData.txt","w");



for(int n = 1; n<30;n++){
gsl_matrix* A_time     = gsl_matrix_alloc(n*3,n*3);
gsl_matrix* R_time     = gsl_matrix_alloc(n*3,n*3);
gsl_vector* tau        = gsl_vector_alloc(n*3);

clock_t start, end;
double cpu_time_used;
double time;

for(int i=0; i<n*3; i++){
	for(int j=0; j<n*3;j++){
		gsl_matrix_set(A_time,i,j,((double) rand()/RAND_MAX));
	}
}

start = clock();

GramSchmidtDecomposition(A_time, R_time);

end = clock();


     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

start = clock();
gsl_linalg_QR_decomp(A_time, tau);
end = clock();
    time = ((double) (end - start)) / CLOCKS_PER_SEC;

		 fprintf(my_out_stream,"%d\t %g\t %g\n",n*3,cpu_time_used, time);


gsl_matrix_free(A_time);
gsl_matrix_free(R_time);
gsl_vector_free(tau);

}

fclose(my_out_stream);








gsl_matrix_free(A);
gsl_matrix_free(AA);
gsl_matrix_free(B_inv);
gsl_matrix_free(BB_inv);
gsl_matrix_free(B_invB);
gsl_matrix_free(R);
gsl_matrix_free(RR);
gsl_matrix_free(Q);
gsl_matrix_free(QQ);
gsl_matrix_free(QTQ);
gsl_vector_free(b);
gsl_vector_free(x);



//gsl_matrix_free(Acopy);
return 0;
}
