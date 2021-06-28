#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <complex.h>
#include <assert.h>
#include <gsl/gsl_cblas.h>
#define FMT "%7.3f"
#define RND (double)rand()/RAND_MAX

void matrix_column_view_vector (gsl_matrix* A, int j, gsl_vector* v)
{
  assert (v->size == A->size1);
  for (int i=0; i<A->size1; i++)
  {
    gsl_vector_set (v, i, gsl_matrix_get (A, i, j));
  }
}


double dot(gsl_vector* x, gsl_vector* y){
	double x_dot_y;
	gsl_blas_ddot(x, y, &x_dot_y); //GSL function that calculates dot product
	return x_dot_y;
}

double norm(gsl_vector* x){
	return sqrt(dot(x,x));
}

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


void SVD(gsl_matrix* A, gsl_matrix* U, gsl_matrix* D, gsl_matrix* V){

  //allocating for matrices and vectors
int n = A -> size1;
gsl_vector* ap            = gsl_vector_alloc(n);
gsl_vector* aq            = gsl_vector_alloc(n);
gsl_vector* a             = gsl_vector_alloc(n);


int changed;
do{
	changed=0;
  int m = A->size2;
	for(int p=0;p<m-1;p++){

	for(int q=p+1;q<m;q++){
  gsl_matrix_get_col(aq,A,q);
  gsl_matrix_get_col(ap,A,p);


    double apTq = dot(ap,aq);
    double aqTq = dot(aq,aq);
    double apTp = dot(ap,ap);
    double theta = 0.5*atan2(2*apTq,aqTq-apTp);
		double c=cos(theta),s=sin(theta);

    double new_app = c*apTp - s*apTq;
    double new_aqq = s*apTq + c*aqTq;

		if(new_app!=apTp || new_aqq!=aqTq) // do rotation
			{
			changed=1;
			timesJ(A,p,q, theta);
			timesJ(V,p,q, theta); 
			}
	}
}
}while(changed!=0);

for(int i = 0; i< n; i++){
gsl_matrix_get_col(a,A,i);
double norm_a = norm(a);
gsl_matrix_set(D,i,i,norm_a);

gsl_vector_scale(a,1/norm_a);
gsl_matrix_set_col(U,i,a);
}



gsl_vector_free(a);
gsl_vector_free(ap);
gsl_vector_free(aq);
}
