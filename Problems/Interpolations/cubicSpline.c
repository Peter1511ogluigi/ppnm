
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "utilities.h"
#include "cubicSpline.h"


cubicSpline* cubicSpline_init( int numOfPts, double* pts, double* funcVals ){

  int numOfEqs           =  numOfPts - 1;
  cubicSpline* spline    =  (cubicSpline*)malloc( sizeof(cubicSpline) );
  spline -> pts          =  (double*)malloc( numOfPts*sizeof(double) );
  spline -> funcVals     =  (double*)malloc( numOfPts*sizeof(double) );
  spline -> firstCoeff   =  (double*)malloc( numOfPts*sizeof(double) );
  spline -> secondCoeff  =  (double*)malloc( numOfEqs*sizeof(double) );
  spline -> thirdCoeff   =  (double*)malloc( numOfEqs*sizeof(double) );
  spline -> numOfPts     =  numOfPts;

  // Fill field values from input data
  for( int it = 0; it < numOfPts; it++ ){
    spline -> pts[it]       =  pts[it];
    spline -> funcVals[it]  =  funcVals[it];
  }

  // Set up slope array (p_i and h_i) variables
  double ptsDiff[numOfEqs];
  double slope[numOfEqs];
  for(int it = 0; it < numOfEqs; it++){
    ptsDiff[it] = pts[it + 1] - pts[it];
    assert( ptsDiff[it] > 0 );

    slope[it] = (funcVals[it + 1] - funcVals[it]) / ptsDiff[it];
  }

  double LHS_matrixDiag[numOfPts];
  double LHS_matrixAboveDiag[numOfPts - 1];
  double RHS_vec[numOfPts];

  LHS_matrixDiag[0]       = 2             ;
  LHS_matrixAboveDiag[0]  = 1             ;
  RHS_vec[0]              = 3 * slope[0]  ;

  for( int it = 0; it < numOfEqs - 1; it++ ) {
    LHS_matrixDiag[it + 1]        =  2 * ptsDiff[it] / ptsDiff[it + 1] + 2;
    LHS_matrixAboveDiag[it + 1]   =      ptsDiff[it] / ptsDiff[it + 1];
    RHS_vec[it + 1]               =  3 * (slope[it] + slope[it + 1] * ptsDiff[it] / ptsDiff[it + 1]);
  }
  LHS_matrixDiag[numOfPts - 1]  =  2                      ;
  RHS_vec[numOfPts - 1]         =  3 * slope[numOfPts - 2];

  for(int it = 1; it < numOfPts; it++){
    LHS_matrixDiag[it]  -=  LHS_matrixAboveDiag[it - 1]  / LHS_matrixDiag[it - 1];
    RHS_vec[it]         -=  RHS_vec[it - 1]              / LHS_matrixDiag[it - 1];
  }

  spline -> firstCoeff[numOfEqs] = RHS_vec[numOfPts - 1] / LHS_matrixDiag[numOfPts - 1];
  for( int it = numOfEqs - 1; it >= 0; it--){
    spline -> firstCoeff[it] = (RHS_vec[it] - LHS_matrixAboveDiag[it]*(spline -> firstCoeff[it + 1])) / LHS_matrixDiag[it];
  }
                       
  for( int it = 0; it < numOfEqs; it++ ){
    spline -> secondCoeff[it]  =  (-2 * (spline -> firstCoeff[it]) - (spline -> firstCoeff[it + 1]) + 3*slope[it]) / ptsDiff[it];
    spline -> thirdCoeff[it]   =  (     (spline -> firstCoeff[it]) + (spline -> firstCoeff[it + 1]) - 2*slope[it]) / ptsDiff[it] /ptsDiff[it];
  }

  return spline;
}

double cubicSpline_eval( cubicSpline* spline, double evalPt ){

  assert( (evalPt >=  (spline -> pts[0])) && (evalPt <= (spline -> pts[spline -> numOfPts -1])) );

  //  Use a binary search to determine which subinterval evalPt is in
  int whichInterval  =  binarySearch(spline -> numOfPts, spline -> pts, evalPt);

  double ptsDiff     =  evalPt  - ( spline -> pts[whichInterval]                      );
  double thirdDiff   =  ptsDiff * ( spline -> thirdCoeff[whichInterval]               );
  double secondDiff  =  ptsDiff * ((spline -> secondCoeff[whichInterval]) + thirdDiff );
  double firstDiff   =  ptsDiff * ((spline -> firstCoeff[whichInterval])  + secondDiff);

  double interpVal = (spline -> funcVals[whichInterval]) + firstDiff ;
  return interpVal;
}

double cubicSpline_eval_deriv( cubicSpline* spline, double evalPt ){

  assert( (evalPt >=  (spline -> pts[0])) && (evalPt <= (spline -> pts[spline -> numOfPts -1])) );

  //  Use a binary search to determine which subinterval evalPt is in
  int whichInterval = binarySearch(spline -> numOfPts, spline -> pts, evalPt);
  double ptsDiff = evalPt - (spline -> pts[whichInterval]);

  double interpVal = (spline -> firstCoeff[whichInterval]) + 2*ptsDiff*(spline -> secondCoeff[whichInterval]) + 3*ptsDiff*ptsDiff*(spline -> thirdCoeff[whichInterval]) ;
  return interpVal;
}


double cubicSpline_eval_integ( cubicSpline* spline, double evalPt ){

  int whichInterval   =  binarySearch((spline -> numOfPts), (spline -> pts), evalPt);

  double integral = 0;
  double ptsDiff;
  for ( int Id = 0; Id <= whichInterval; Id++ ){

    if ( Id < whichInterval ){
      ptsDiff = ((spline -> pts[Id + 1]) - (spline -> pts[Id]));
    }
    else {
      ptsDiff = ( evalPt - (spline -> pts[Id]));
    }
    integral += (spline -> funcVals[Id]) * ptsDiff + (spline -> firstCoeff[Id]) * ptsDiff * ptsDiff / 2 + (spline -> secondCoeff[Id]) * ptsDiff * ptsDiff * ptsDiff / 3 + (spline -> thirdCoeff[Id]) * ptsDiff * ptsDiff * ptsDiff * ptsDiff/ 4;
  }

  return integral;
}


void cubicSpline_free(cubicSpline* spline){

  free(spline -> pts);
  free(spline -> funcVals);
  free(spline -> firstCoeff);
  free(spline -> secondCoeff);
  free(spline -> thirdCoeff);
  free(spline);
}
