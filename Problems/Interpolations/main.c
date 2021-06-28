
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>

#include "utilities.h"
#include "integrateFunction.h"
#include "linSpline.h"
#include "quadSpline.h"
#include "cubicSpline.h"

double gsl_cos (double x, void * params) {
  return cos(x);
}

void defIntegral( int numOfPts,
                  double* xData,
                  double* yData,
                  double lowerLimit,
                  double upperLimit,
                  double absError,
                  double relError,
                  size_t iterationLimit ){

    double integralVal = linSpline_integ( numOfPts, xData, yData, upperLimit );

    gsl_function gslFuncCos;
    gslFuncCos.function    =  &gsl_cos;
    gslFuncCos.params      =  NULL;
    double integralValComp = integrateFunction( lowerLimit, upperLimit, &gslFuncCos, absError, relError, iterationLimit );

    printf("The integral of the interpolant            from x[0] = %g to z = %g is         = %g\n",       lowerLimit, upperLimit, integralVal     );
    printf("The integral of cos(x) (from math.h)       from        %g to     %g using GSL  = %g\n",       lowerLimit, upperLimit, integralValComp );
    printf("The integral of cos(x) (from WolframAlpha) from        %g to     %g            = -0.99999\n", lowerLimit, upperLimit                  );
}

int main( int argc, char* argv[]) {
  if ( argc < 2){	// Check that we have passed any arguments
    fprintf(stderr, "Error, no arguments were passed.\n"); // Else print to stderr
    exit(-1);
  }

  // __ Hyperparameters ________________________________________________________
  int numOfPts      =  20;
  int numOfSamples  =  (int)1e3;

  char* inputFilename  =  argv[1];
  FILE* outFileStream_lin    =  fopen(argv[2], "w");
  FILE* outFileStream_quad   =  fopen(argv[3], "w");
  FILE* outFileStream_cubic  =  fopen(argv[4], "w");
  double* xData  =  malloc( numOfPts*sizeof(double) );
  double* yData  =  malloc( numOfPts*sizeof(double) );

	inputToArray( xData, yData, inputFilename );

  double lowerLimit      =   xData[0] ;
  double upperLimit      =   11       ;
  double absError        =   1e-6     ;
  double relError        =   1e-6     ;
  size_t iterationLimit  =   999      ;
  // ___________________________________________________________________________
printf("- - - - Exercise A - - - -\n\n");
printf("We note that our integrator works, since\n\n");
  // __ Definite integral ______________________________________________________
  // Compute the definite integral and compare with GSL routines to make sure data is sensical
  defIntegral( numOfPts, xData, yData, lowerLimit, upperLimit, absError, relError, iterationLimit );
  // ___________________________________________________________________________
printf("\n");
printf("this can also be seen from the plot given in lin_spline_plot.png.\n\n");
printf("- - - - Exercise B - - - - \n\n");
printf("the comparison of quadratic spline can be seen on figure quad_spline_plot.png\n\n");
printf("- - - - Exercise C - - - - \n\n");
printf("The cubic spline implementation can be seen on figure cubic_spline_plot.png");


  // __ Initiallize GSL interpolation __________________________________________
  gsl_interp* gslInterp_lin    =  gsl_interp_alloc(gsl_interp_linear,     numOfPts);
  gsl_interp* gslInterp_quad   =  gsl_interp_alloc(gsl_interp_polynomial, numOfPts);
  gsl_interp* gslInterp_cubic  =  gsl_interp_alloc(gsl_interp_cspline,    numOfPts);
  gsl_interp_init(gslInterp_lin,   xData, yData, numOfPts);
  gsl_interp_init(gslInterp_quad,  xData, yData, numOfPts);
  gsl_interp_init(gslInterp_cubic, xData, yData, numOfPts);
  // ___________________________________________________________________________

  // __ LINEAR SPLINE INTERPOLATION ____________________________________________
  double resolution = fabs(xData[numOfPts-1] - xData[0]) / numOfSamples;
  for ( double evalPt = xData[0]; evalPt < xData[numOfPts]; evalPt += resolution ){

    double interp_tmp            =  linSpline(              numOfPts,       xData, yData,           evalPt       );
    double interp_gsl_tmp        =  gsl_interp_eval(        gslInterp_lin,  xData, yData,           evalPt, NULL );
    double interp_integ_tmp      =  linSpline_integ(        numOfPts,       xData, yData,           evalPt       );
    double interp_integ_gsl_tmp  =  gsl_interp_eval_integ(  gslInterp_lin,  xData, yData, xData[0], evalPt, NULL );

    fprintf(outFileStream_lin, "%g\t%g\t%g\t%g\t%g\n", evalPt, interp_tmp, interp_gsl_tmp, interp_integ_tmp, interp_integ_gsl_tmp);
  }
  // ___________________________________________________________________________

  // __ QUADRATIC SPLINE INTERPOLATION _________________________________________
  quadSpline* spline_quad = quadSpline_init( numOfPts, xData, yData );

  for ( double evalPt = xData[0]; evalPt < xData[numOfPts]; evalPt += resolution ){

    double interp_tmp            =  quadSpline_eval(        spline_quad, evalPt );
    double interp_gsl_tmp        =  gsl_interp_eval(        gslInterp_quad,  xData, yData,           evalPt, NULL );
    double interp_integ_tmp      =  quadSpline_eval_integ(  spline_quad, evalPt       );
    double interp_integ_gsl_tmp  =  gsl_interp_eval_integ(  gslInterp_quad,   xData, yData, xData[0], evalPt, NULL );
    double interp_deriv_tmp      =  quadSpline_eval_deriv(  spline_quad, evalPt       );
    double interp_deriv_gsl_tmp  =  gsl_interp_eval_deriv(  gslInterp_quad,   xData, yData, evalPt, NULL );

    fprintf(outFileStream_quad, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", evalPt, interp_tmp, interp_gsl_tmp, interp_integ_tmp, interp_integ_gsl_tmp, interp_deriv_tmp, interp_deriv_gsl_tmp);
  }
  // ___________________________________________________________________________

  // __ CUBIC SPLINE INTERPOLATION _____________________________________________
  cubicSpline* spline_cubic = cubicSpline_init( numOfPts, xData, yData );

  for ( double evalPt = xData[0]; evalPt < xData[numOfPts]; evalPt += resolution ){

    double interp_tmp             =  cubicSpline_eval( spline_cubic, evalPt );
    double interp_gsl_tmp         =  gsl_interp_eval(gslInterp_cubic,  xData, yData,           evalPt, NULL );
    double interp_integ_tmp       =  cubicSpline_eval_integ(  spline_cubic, evalPt       );
    double interp_integ_gsl_tmp   =  gsl_interp_eval_integ(  gslInterp_cubic,   xData, yData, xData[0], evalPt, NULL );
    double interp_deriv_tmp       =  cubicSpline_eval_deriv(  spline_cubic, evalPt       );
    double interp_deriv_gsl_tmp   =  gsl_interp_eval_deriv(  gslInterp_cubic,   xData, yData, evalPt, NULL );

    fprintf(outFileStream_cubic, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", evalPt, interp_tmp, interp_gsl_tmp, interp_integ_tmp, interp_integ_gsl_tmp, interp_deriv_tmp, interp_deriv_gsl_tmp);
  }
  // ___________________________________________________________________________


  // __ Close files and free dynamic memory ____________________________________
  fclose(outFileStream_lin);
  fclose(outFileStream_quad);

  quadSpline_free(spline_quad);
  cubicSpline_free(spline_cubic);
  gsl_interp_free(gslInterp_lin);
  gsl_interp_free(gslInterp_quad);
  // ___________________________________________________________________________

  return 0;
}
