#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
/////////////////////////////////////////////////////////////
//--------include the QR decomposition and solver----------//
/////////////////////////////////////////////////////////////

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



void newton(
	void f(gsl_vector* x, gsl_vector* fx),
	gsl_vector* x,
	double eps){

int n = x -> size;
int LIMIT = 10000;
int steps = 0;

gsl_matrix* J = gsl_matrix_alloc(n,n);
gsl_matrix* R = gsl_matrix_alloc(n,n);
gsl_vector* Dx = gsl_vector_alloc(n);
gsl_vector* fx = gsl_vector_alloc(n);
gsl_vector* df = gsl_vector_alloc(n);
gsl_vector* y = gsl_vector_alloc(n);
gsl_vector* fy = gsl_vector_alloc(n);

while(steps<LIMIT){
f(x,fx); // For the given x, calculate f(x) and store in fx

// We create the Jacobian:.
for(int k=0; k<n; k++){
	double xk = gsl_vector_get(x,k);
	gsl_vector_set(x, k, xk+sqrt(DBL_EPSILON)); // Our finite difference is square root of machine epsilon
	f(x,df); // now df=f(x+sqrt())
        gsl_vector_sub(df,fx); // df=f(x+sqrt())-f(x), as it should
	// Now we have what we need to calculate J
	for(int i=0; i<n; i++){
		double Jik = (gsl_vector_get(df,i)/sqrt(DBL_EPSILON)); // Finite difference
		gsl_matrix_set(J,i,k,Jik);}
	gsl_vector_set(x,k,xk);
}
// Next we solve the equation J*Dx=-f(x) using the routines from the homework "Linear Equations"
GramSchmidtDecomposition(J,R);
QR_solve(J, R, fx, Dx); // solution of J*Dx=f(x) now stored in Dx. We want the solution of J*Dx=-f(x):
gsl_vector_scale(Dx,-1.0);
// Then, we introduce lambda and use the conservative step lambda*Dx
double lambda = 1.0; // initial value
while(lambda>1./64){
gsl_vector_memcpy(y,x); // copy of x
gsl_vector_add(y,Dx); // z=x+Dx
f(y, fy); // f(x+Dx)
if( gsl_blas_dnrm2(fy)<(1-lambda/2)*gsl_blas_dnrm2(fx)) break;
lambda*=0.5;
gsl_vector_scale(Dx,0.5);
}

// Put final result into fx
gsl_vector_memcpy(x,y);
gsl_vector_memcpy(fx,fy);
if(gsl_blas_dnrm2(Dx)<sqrt(DBL_EPSILON) || gsl_blas_dnrm2(fx)<eps ) break;
steps++;
}

gsl_matrix_free(J);
gsl_matrix_free(R);
gsl_vector_free(Dx);
gsl_vector_free(fx);
gsl_vector_free(df);
gsl_vector_free(y);
gsl_vector_free(fy);

}





////////////////////////////////////////////////////
//-------------implement newtons method-----------//
////////////////////////////////////////////////////
/*
void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps)
{
  int n = x -> size;
  int limits = 1000;
  int steps  = 0;

  gsl_matrix* J = gsl_matrix_alloc(n,n);
  gsl_matrix* R = gsl_matrix_alloc(n,n);
  gsl_vector* Dx = gsl_vector_alloc(n);
  gsl_vector* fx = gsl_vector_alloc(n);
  gsl_vector* df = gsl_vector_alloc(n);
  gsl_vector* y = gsl_vector_alloc(n);
  gsl_vector* fy = gsl_vector_alloc(n);

  while(steps < limits)
  {
      f(x,fx);
      for(int j=0; j<n; j++)
        {
          double xj = gsl_vector_get(x,j);
          gsl_vector_set(x,j,xj + sqrt(DBL_EPSILON));
          f(x,df);
          gsl_vector_sub(df,fx);
          for(int i=0; i< n; i++)
          {
            double Jij = (gsl_vector_get(df,i)/sqrt(DBL_EPSILON));
		        gsl_matrix_set(J,i,j,Jij);
          }
          gsl_vector_set(x,j,xj);
        }
  GramSchmidtDecomposition(J,R);
  QR_solve(J,R,fx,Dx);
  gsl_vector_scale(Dx,-1.0);



double s =  1;
 while(s>1.64){


gsl_vector_memcpy(y,x);
gsl_vector_add(y,Dx);
f(y,fy);
if(gsl_blas_dnrm2(fy)<(1-s*0.5)*gsl_blas_dnrm2(fx)) break;
s *= 0.5;
gsl_vector_scale(Dx,0.5);

}
gsl_vector_memcpy(x,y);
gsl_vector_memcpy(fx,fy);
if(gsl_blas_dnrm2(Dx)<sqrt(DBL_EPSILON) || gsl_blas_dnrm2(fx)<eps) break;
steps++;
}




gsl_matrix_free(J);
gsl_matrix_free(R);
gsl_vector_free(Dx);
gsl_vector_free(fx);
gsl_vector_free(df);
gsl_vector_free(y);
gsl_vector_free(fy);
}
*/
