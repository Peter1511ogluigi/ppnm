#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <complex.h>
#include <gsl/gsl_cblas.h>
#define FMT "%7.3f"
#define RND (double)rand()/RAND_MAX

void printMatrix(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0; c<A->size2;c++){
			printf(FMT,gsl_matrix_get(A,r,c));

    }
		printf("\n\n");
	}
}

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
		}
}

void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
int changed;
do{
	changed=0;
  int m = A->size2;
	for(int p=0;p<m-1;p++)
	for(int q=p+1;q<m;q++){
		double apq=gsl_matrix_get(A,p,q);
		double app=gsl_matrix_get(A,p,p);
		double aqq=gsl_matrix_get(A,q,q);
		double theta=0.5*atan2(2*apq,aqq-app);
		double c=cos(theta),s=sin(theta);
		double new_app=c*c*app-2*s*c*apq+s*s*aqq;
		double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
		if(new_app!=app || new_aqq!=aqq) // do rotation
			{
			changed=1;
			timesJ(A,p,q, theta);
			Jtimes(A,p,q,-theta); // A←J^T*A*J
			timesJ(V,p,q, theta); // V←V*J
			}
	}
}while(changed!=0);
}




int main(){


int m = 5;
//int p = 1;
//int q = 3;
//double theta = 0;
  gsl_matrix* A     = gsl_matrix_alloc(m, m);
  gsl_matrix* V     = gsl_matrix_alloc(m, m);
  gsl_matrix* DVT     = gsl_matrix_alloc(m, m);
  gsl_matrix* VDVT     = gsl_matrix_alloc(m, m);

  gsl_matrix* hopefully_I = gsl_matrix_alloc(m, m);

//creating diagonal matrix
  for(int i=0; i<m; i++){
  	for(int j=0; j<m;j++){
      double number = ((double) rand()/RAND_MAX);
  		gsl_matrix_set(A,i,j,number);
      gsl_matrix_set(A,j,i,number);
  	}
  }


gsl_matrix_set_identity(V);
printf("Matrix A is given as:\n\n");
printMatrix(A);
printf("\n\n");
jacobi_diag(A,V);



printf("J^T*A*J is given as:\n\n");
printMatrix(A);
printf("\n\n");
printf("The matrix V is given as;\n");
printMatrix(V);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,V,0.0,hopefully_I);
printf("\n\n");
printf("The matrix product VTV gives\n");
printMatrix(hopefully_I);
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,A,V,0.0,DVT);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,DVT,0.0,VDVT);
printf("Hopefully A\n");
printMatrix(VDVT);

// part B

printf("\n\n");
printf("------------- part B -------------\n\n");
printf("The eigenvalues are given as:\n");

//build Hamiltonian matrix
int N=100;
double s=1.0/(N+1);
gsl_matrix* H = gsl_matrix_alloc(N,N);
for(int i=0;i<N-1;i++){
  gsl_matrix_set(H,i,i,-2);
  gsl_matrix_set(H,i,i+1,1);
  gsl_matrix_set(H,i+1,i,1);
  }
gsl_matrix_set(H,N-1,N-1,-2);
gsl_matrix_scale(H,-1/s/s);
//printMatrix(H);
gsl_matrix* VV = gsl_matrix_alloc(N,N);
gsl_matrix_set_identity(VV);
jacobi_diag(H,VV);
/*
printf("\n\n");
for (int k=0; k < n/3; k++){
    double exact = M_PI*M_PI*(k+1)*(k+1);
    double calculated = gsl_matrix_get(H,k,k);
    printf("%i %g %g\n",k,calculated,exact);
}
*/
char* InputFileName = "data.txt";
//printf("matrix V is given as;\n\n");
//printMatrix(VV);
FILE* my_out_stream = fopen(InputFileName, "w");
printf("\n\n");
for (int k=0; k < N/10; k++){
    double exact = M_PI*M_PI*(k+1)*(k+1);
    double calculated = gsl_matrix_get(H,k,k);
    printf("n=%i;  calculated = %g;  exact = %g;\n",k,calculated,exact);
}
for (int i = 0; i < N; ++i) {
		 gsl_vector_view view = gsl_matrix_column(VV, i);
		 gsl_vector* vi = &view.vector;
		 double first = gsl_vector_get(vi,0);
		 if (first < 0) {
				gsl_vector_scale(vi, -1);
		 }
	}

double x;
int nplots = 3;
   double scale = sqrt(2.0/(N + 1));

	fprintf(my_out_stream, "0 ");
	for (int k = 0; k < nplots; ++k) {
		 fprintf(my_out_stream, "0 0 ");
	}
	fprintf(my_out_stream, "\n");

	for(int i = 0; i < N; i++) {
		 x = (i + 1.0)/(N + 1);
		 fprintf(my_out_stream, "%g ", x);
		 for (int k = 0; k < nplots; ++k) {
				fprintf(my_out_stream, "%g %g ", sin((k+1)*M_PI*x), gsl_matrix_get(VV, i, k)/scale);
		 }
		 fprintf(my_out_stream, "\n");
	}

	fprintf(my_out_stream, "1 ");
	for (int k = 0; k < nplots; ++k) {
		 fprintf(my_out_stream, "0 0 ");
	}
	fprintf(my_out_stream, "\n");

  /*fprintf(my_out_stream,"0\t0\t0\t0\n");
  for(int i=0;i<n;i++)
	fprintf(my_out_stream,"%g\t%lf\t%lf\t%lf\n",(i+1.0)/(n+1), gsl_matrix_get(H,i,0),gsl_matrix_get(H,i,1),gsl_matrix_get(H,i,2));
  fprintf(my_out_stream,"1\t0\t0\t0\n");*/
fclose(my_out_stream);

gsl_matrix_free(A);
gsl_matrix_free(V);
gsl_matrix_free(VV);

gsl_matrix_free(H);
gsl_matrix_free(DVT);
gsl_matrix_free(VDVT);

gsl_matrix_free(hopefully_I);
  return 0;
}
