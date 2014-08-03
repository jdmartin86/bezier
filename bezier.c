/**
 * bezier.c
 *
 * This library computes the Bezier spline coefficients for an input 
 * curve and evaluates the spline at a given resolution.
 *
 * author: john martin jr.
 * e-mail: jdmartin86@gmail.com
 */
#include "bezier.h"

static const double PI = 3.1415926535897932384626;

// quartic Bezier basis matrix
static const double Mcr[DEGREE+1][DEGREE+1] = 
  { {  1.0 ,  -4.0 ,   6.0 , 4.0 , 1.0 },
    { -4.0 ,  12.0 , -12.0 , 4.0 , 0.0 },
    {  6.0 , -12.0 ,   6.0 , 0.0 , 0.0 },
    { -4.0 ,   4.0 ,   0.0 , 0.0 , 0.0 },
    {  1.0 ,   0.0 ,   0.0 , 0.0 , 0.0 } };

/**
 * create_curve
 *
 * create curve structure
 */
curve_t* create_curve( unsigned int len )
{
  curve_t* curve = (curve_t*) malloc( sizeof(curve_t) );
  
  // allocate and initialize curve
  curve->x = (double*) calloc( len , sizeof(double) );
  curve->len = len;
  return curve;
}

/**
 * free_curve
 *
 * releases memory of curve structure
 */
void free_curve( curve_t* curve )
{
  free( curve->x );
  curve->x = NULL;
  curve->len = 0;
}

/**
 * create_spline
 *
 * create spline structure
 */
spline_t* create_spline( curve_t* curve )
{
  int num_cpts = DEGREE + 1, 
      num_segs = curve->len - DEGREE,
      num_coeff = num_cpts*num_segs;

  spline_t* spline = (spline_t*) malloc( sizeof(spline_t) );
  
  // allocate and initialize spline
  spline->curve = curve;
  
  //
  // one big block for coefficients
  //
  // allocate a pointer to num_coeff element array for each segment
  spline->coeff = (double**) malloc( num_segs*sizeof(double*) ); 
  
  // allocate (DEGREE+1)*curve->len coefficients in total
  double* c     = (double*)  calloc( num_coeff , sizeof(double) );
  
  // arrange the coefficients contiguously w.r.t access
  for( int i = 0 ; i < num_segs ; i++ )
    spline->coeff[i] = c + i*num_cpts;
  
  spline->len = num_coeff;
  spline->num_segs = num_segs;
  
  return spline;
}

/**
 * free_spline
 *
 * releases memory of spline structure
 */
void free_spline( spline_t* spline )
{  
  spline->curve = NULL;

  free( *(spline->coeff) );
  free( spline->coeff );
  spline->coeff = NULL;

  spline->len = 0;
  spline->num_segs = 0;
}

/**
 * curve_to_file
 *
 * write the curve to the file "curve.txt"
 */
void curve_to_file( curve_t* curve )
{
  FILE* out;
  out = fopen( "curve.txt" , "w" );
  
  if ( out != NULL )
  {
    for( int i = 0 ; i < curve->len ; i++ )
      fprintf( out , "%g\n" , curve->x[i] );
  }
  else
    printf( "bezier: failed to write curve to file\n" );
}

/**
 * spline_to_file
 *
 * write the spline to the file "spline.txt". each row represents the
 * spline coefficients; it's only four wide now.
 */
void spline_to_file( spline_t* spline )
{
  FILE* out;
  out = fopen( "spline.txt" , "w" );
  
  if ( out != NULL )
  {
    for( int i = 0 ; i < spline->num_segs ; i++ )
      fprintf( out ,
	       "%g,%g,%g,%g\n",
	       spline->coeff[i][0],
	       spline->coeff[i][1],
	       spline->coeff[i][2],
	       spline->coeff[i][3] );
  }
  else
    printf( "bezier: failed to write spline to file\n" );
}

/**
 * gmatrix
 *
 * computes the geometry matrix for every point on a curve
 */
void gmatrix( spline_t* G , curve_t* curve )
{
  assert( G     != NULL );
  assert( curve != NULL );

  int num_segs = G->curve->len - DEGREE;

  for( int i = 0 , k = 0 ; i < num_segs ; i++ , k++ )
    for ( int j = 0 ; j < DEGREE + 1 ; j++ )
      G->coeff[i][j] = curve->x[j+k];
}

/**
 * segment
 *
 * compute coefficients for one segment
 */
void segment( double S[DEGREE+1] , double G[DEGREE+1] )
{
  assert( S != NULL );
  assert( G != NULL );

  memset( S , 0.0 , (DEGREE+1)*sizeof(double) );
  
  for( int i = 0 ; i < DEGREE + 1 ; i++ )
    for( int j = 0 ; j < DEGREE + 1 ; j++ ) 
      S[i] += Mcr[i][j]*G[j];     
}

/**
 * coefficients
 *
 * computes the bezier spline coefficients for an input geometry
 * matrix G. output is written to spline. 
 */
void coefficients( spline_t* spline , spline_t* G )
{
  assert( spline != NULL );
  assert( G      != NULL );

  for( int i = 0 ; i < spline->num_segs ; i++ )
    segment( spline->coeff[i] , G->coeff[i] );
}

/**
 * polynomial
 *
 * evaluate a quartic polynomial with specified coefficients for the 
 * parameter value t
 */
double polynomial( double coeff[DEGREE+1] , double t )
{
  assert( coeff != NULL );
  assert( t >= 0.0 || t <= 1.0 );

  return coeff[0]*t*t*t + coeff[1]*t*t + coeff[2]*t + coeff[3];
}


/**
 * spline
 * 
 * computes the bezier spline coefficients for every point on a 
 * curve. the first and last point will not be connected to the spline.
 */
void spline( curve_t* curve , spline_t* spline )
{
  assert( curve  != NULL );
  assert( spline != NULL );

  // allocate the geometry matrices
  spline_t* G = create_spline( curve );

  // compute the geometry matrix for every point along the curve
  gmatrix( G , curve );

  // compute polynomial coefficients
  coefficients( spline , G );
  
  // discard temporary geometry matrices
  free_spline( G );
}

/**
 * evaluate
 *
 * given spline coefficients, and a properly-sized solution curve,
 * this function will evaluate n = solution size / number of segments
 * equally-spaced points along each segment. evaluation starts at the
 * parameter value zero. 
 */
void evaluate( spline_t* spline , curve_t* solution )
{
  assert( spline   != NULL );
  assert( solution != NULL );

  // interpolated trajectory must be an integer multiple larger than 
  // the number of segments
  assert( spline->num_segs <= solution->len );
  
  int rem = (int) fmod( solution->len , spline->num_segs );
 
  assert( rem == 0 );
  // compute the interpolation resolution
  int res = solution->len/spline->num_segs;
  
  // evaluate the polynomial
  for( int i = 0 , k = 0 ; i < spline->num_segs ; i++ , k+=res )
    for( int j = 0 ; j < res ; j++ )
      solution->x[j+k] = polynomial( spline->coeff[i] , (double) j/res );
}

int main( void )
{
  // TODO: add test code
  return 1;
}
