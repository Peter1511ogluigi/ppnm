#ifndef HAVE_CUBICSPLINE_H
#define HAVE_CUBICSPLINE_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct { int      numOfPts    ; /* n   */
                 double*  pts         ; /* x's */
                 double*  funcVals    ; /* y's */
                 double*  firstCoeff  ; /* b_i */
                 double*  secondCoeff ; /* c_i */
                 double*  thirdCoeff  ; /* d_i */ } cubicSpline;

cubicSpline* cubicSpline_init( int numOfPts, double* pts, double* funcVals );
double cubicSpline_eval( cubicSpline* spline, double evalPt );
double cubicSpline_eval_integ( cubicSpline* spline, double evalPt );
double cubicSpline_eval_deriv( cubicSpline* spline, double evalPt );
void cubicSpline_free(cubicSpline* spline);

#endif
