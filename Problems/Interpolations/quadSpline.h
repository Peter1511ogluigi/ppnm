#ifndef HAVE_QUADSPLINE_H
#define HAVE_QUADSPLINE_H

#include <stdlib.h>
#include <assert.h>

typedef struct { int      numOfPts    ; /* n   */
                 double*  pts         ; /* x's */
                 double*  funcVals    ; /* y's */
                 double*  firstCoeff  ; /* b_i */
                 double*  secondCoeff ; /* c_i */ } quadSpline;

quadSpline* quadSpline_init( int numOfPts, double* pts, double* funcVals );
double quadSpline_eval( quadSpline* spline, double evalPt );
double quadSpline_eval_deriv( quadSpline* spline, double evalPt );
double quadSpline_eval_integ( quadSpline* spline, double evalPt );
void quadSpline_free( quadSpline* spline );

#endif
