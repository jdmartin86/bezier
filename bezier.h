/**
 * bezier.h
 *
 * interface file for the bezier library
 *
 * author: john martin jr. 
 * e-mail: jdmartin86@gmail.com
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// degree of bezier spline segment 
#define DEGREE 4

/**
 * curve_t
 * 
 * input curve structure to be interpolated
 */
typedef struct curve_
{
  double* x;
  unsigned int len;
} curve_t;

/**
 * spline_t
 *
 * output structure for all the spline coefficients
 */
typedef struct spline_
{
  curve_t* curve;
  double** coeff;
  unsigned int len;
  unsigned int num_segs;
} spline_t;

/**
 * create_curve
 *
 * create curve structure
 */
curve_t* create_curve( unsigned int len );

/**
 * free_curve
 *
 * releases memory of curve structure
 */
void free_curve( curve_t* curve );

/**
 * create_spline
 *
 * create spline structure
 */
spline_t* create_spline( curve_t* curve );

/**
 * free_spline
 *
 * releases memory of spline structure
 */
void free_spline( spline_t* spline );

/**
 * curve_to_file
 *
 * write the curve to the file "curve.txt"
 */
void curve_to_file( curve_t* curve );

/**
 * spline_to_file
 *
 * write the spline to the file "spline.txt". each row represents the
 * spline coefficients.
 */
void spline_to_file( spline_t* spline );

/**
 * gmatrix
 *
 * computes the geometry matrix for every point in a curve
 */
void gmatrix( spline_t* G , curve_t* curve );

/**
 * segment
 *
 * compute coefficients for one segment
 */
void segment( double S[DEGREE+1] , double G[DEGREE+1] );

/**
 * coefficients
 *
 * computes the bezier spline coefficients for an input geometry
 * matrix G. output is written to spline. 
 */
void coefficients( spline_t* spline , spline_t* G );

/**
 * polynomial
 *
 * evaluate a quartic polynomial with specified coefficients for the 
 * parameter value t
 */
double polynomial( double coeff[DEGREE+1] , double t );

/**
 * spline
 * 
 * computes the bezier spline coefficients for every point on a 
 * curve
 */
void spline( curve_t* curve , spline_t* spline );

/**
 * evaluate
 *
 * given spline coefficients, and a properly-sized solution curve,
 * this function will evaluate n = solution size / number of segments
 * equally-spaced points along each segment. evaluation starts at the
 * parameter value zero. 
 */
void evaluate( spline_t* spline , curve_t* solution );
