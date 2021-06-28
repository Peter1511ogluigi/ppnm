
#ifndef HAVE_LINSPLINE_H
#define HAVE_LINSPLINE_H

double linSpline      ( int numOfPts, double* pts, double* funcVals, double evalPt );
double linSpline_integ( int numOfPts, double* pts, double* funcVals, double evalPt );

#endif
