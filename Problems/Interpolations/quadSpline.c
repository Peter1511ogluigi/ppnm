#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "utilities.h"
#include "quadSpline.h"

quadSpline* quadSpline_init( int numOfPts, double* pts, double* funcVals ){
  //  ------------------------------------------------------------------------------
  /*  quadSpline constructor to initiallize a quadratic spline from a quadSpline
      struct, by filling the various field values from respective function inputs.
      ¤ int       numOfPts  : The number of points to interpolate in.
      ¤ double*   pts       : A pointer to an array of doubles,
                              the known points {x_i}.
      ¤ double*   funcVals  : A pointer to an array of doubles,
                              the corresponding function values {f(x_i)}
      Returns: An inittialized quadSpline* struct                                    */
  //  ------------------------------------------------------------------------------
  quadSpline* spline = (quadSpline*)malloc( sizeof(quadSpline) );

  // NumOfEqs is the number of needed equations for the quadratic spline interpolation,
  // or the number of intervals to interpolate in (between a numOfPts amount of points)
  int numOfEqs  =  numOfPts - 1;

  // Fill out field values by first allocating memory to be initiallized below.
  spline -> firstCoeff   =  (double*)malloc( numOfEqs*sizeof( double ) );
  spline -> secondCoeff  =  (double*)malloc( numOfEqs*sizeof( double ) );
  spline -> pts          =  (double*)malloc( numOfPts*sizeof( double ) );
  spline -> funcVals     =  (double*)malloc( numOfPts*sizeof( double ) );
  spline -> numOfPts     =  numOfPts;

  // Fill out pts[] and funcVals[] arrays at fields
  for ( int it = 0; it < numOfPts; it++ ){
    spline -> pts[it]       =  pts[it];
    spline -> funcVals[it]  =  funcVals[it];
  }

  // Compute slopes at each intervals between individual pts[i], pts[i+1]
  double ptsDiff[numOfPts - 1];
  double slope[numOfPts - 1];
  for ( int it = 0; it < numOfPts - 1; it++ ){
    ptsDiff[it]  =  pts[it + 1] - pts[it];
    slope[it]    =  (funcVals[it + 1] - funcVals[it]) / ptsDiff[it];
  }

  // Forward recursion to compute the second quadratic spline coefficient
  spline -> secondCoeff[0] = 0;
  for ( int it = 0; it < numOfPts - 2; it++){
    spline -> secondCoeff[it + 1] = (slope[it + 1] - slope[it] - (spline->secondCoeff[it])*ptsDiff[it])/ptsDiff[it + 1];
  }

  // Backward recursion, starting from last value from forward recursion
  spline -> secondCoeff[numOfPts - 2] /=2;
  for ( int it = numOfPts - 3; it >= 0; it--){
    spline -> secondCoeff[it] = (slope[it + 1] - slope[it] - (spline->secondCoeff[it + 1]) * ptsDiff[it + 1]) / ptsDiff[it];
  }

  // Finally compute first coefficient
  for ( int it = 0; it < numOfPts - 1; it++ ){
    spline -> firstCoeff[it] = slope[it] - (spline->secondCoeff[it])*ptsDiff[it];
  }

  return spline;
}


double quadSpline_eval( quadSpline* spline, double evalPt ){
  //  ------------------------------------------------------------------------------
  /*  Do quadratic spline interpolation, using an already initiallized quadSpline*
      struct. The quadSpline* struct may be initiallized using quadSpline_init().
      The quadratic interpolant is computed at the evaluation pt. evalPt.
      ¤ quadSpline*       : A pointer to an initiallized quadSpline struct,
                            initiallized using quadSpline_init()
      ¤ double*   evalPt  : Point at which to evaluate interpolant at
      Returns: The function value of the interpolant at evalPt.                      */
  //  ------------------------------------------------------------------------------
  assert( (evalPt >=  (spline -> pts[0])) && (evalPt <= (spline -> pts[spline -> numOfPts -1])) );

  //  Use a binary search to determine which subinterval evalPt is in
  int whichInterval = binarySearch(spline -> numOfPts, spline -> pts, evalPt);
  double ptsDiff = evalPt - (spline -> pts[whichInterval]);

  double interpVal = (spline -> funcVals[whichInterval]) + ptsDiff*((spline -> firstCoeff[whichInterval]) + ptsDiff*(spline -> secondCoeff[whichInterval]));
  return interpVal;
}


double quadSpline_eval_deriv( quadSpline* spline, double evalPt ){
  //  ------------------------------------------------------------------------------
  /*  Do quadratic spline interpolation og function derivative, using an already
      initiallized quadSpline* struct. The quadSpline* struct may be initiallized
      using quadSpline_init(). The quadratic interpolant is computed at the evaluation
      pointt evalPt.
      ¤ quadSpline*       : A pointer to an initiallized quadSpline struct,
                            initiallized using quadSpline_init()
      ¤ double*   evalPt  : Point at which to evaluate interpolant at
      Returns: The function value of the interpolant derivative at evalPt.                      */
  //  ------------------------------------------------------------------------------
  assert( (evalPt >=  (spline -> pts[0])) && (evalPt <= (spline -> pts[spline -> numOfPts -1])) );

  //  Use a binary search to determine which subinterval evalPt is in
  int whichInterval = binarySearch(spline -> numOfPts, spline -> pts, evalPt);
  double ptsDiff = evalPt - (spline -> pts[whichInterval]);

  double interpVal = (spline -> firstCoeff[whichInterval]) + 2*ptsDiff*(spline -> secondCoeff[whichInterval]);
  return interpVal;
}


double quadSpline_eval_integ( quadSpline* spline, double evalPt ){
  //  ------------------------------------------------------------------------------
  /*  Integrate quadratic spline interpolation from a set of known points (x_i, f(x_i))
      using a binary search method. This asumes the dataset is an ordered array.
      Otherwise use sorting before passing to the function.
      ¤ int       numOfPts  : The number of points to interpolate in.
      ¤ double*   pts       : A pointer to an array of doubles,
                              the known points {x_i}.
      ¤ double*   funcVals  : A pointer to an array of doubles,
                              the corresponding function values {f(x_i)}
      ¤ double    evalPt    : The query point, double values at which to evaluate
                              the interpolant.
      Returns: A double P(evalPt), where P() is the interpolant function integrated  */
  //  ------------------------------------------------------------------------------

  //  Use a binary search to determine which subinterval evalPt is in
  int whichInterval   =  binarySearch((spline -> numOfPts), (spline -> pts), evalPt);

  double integral = 0;
  double ptsDiff;
  for ( int Id = 0; Id <= whichInterval; Id++ ){

    // Note: pow(x, 2) is slower than x*x in principle, but they are definitively comparable below 1e6 elements
    if ( Id < whichInterval ){
      ptsDiff = ((spline -> pts[Id + 1]) - (spline -> pts[Id]));
    }
    else {
      ptsDiff = ( evalPt - (spline -> pts[Id]) );
    }
    integral += (spline -> funcVals[Id]) * ptsDiff + (spline -> firstCoeff[Id]) * ptsDiff * ptsDiff / 2 + (spline -> secondCoeff[Id]) * ptsDiff * ptsDiff * ptsDiff / 3;
  }

  return integral;
}


void quadSpline_free( quadSpline* spline ){
  //  ------------------------------------------------------------------------------
  /*  quadSpline destructor, to free allocated memory.
      ¤ quadSpline*       : A pointer to an initiallized quadSpline struct,
                            initiallized using quadSpline_init()                     */
  //  ------------------------------------------------------------------------------
  free(spline -> pts);
  free(spline -> funcVals);
  free(spline -> firstCoeff);
  free(spline -> secondCoeff);
  free(spline);
}
