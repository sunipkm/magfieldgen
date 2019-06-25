#ifndef magfieldgen_h
#define magfieldgen_h

/* Headers */
#include <stdio.h> //standard library
#include <math.h> //math library
#include <stdlib.h> //allocations
#include <fftw3.h> //fastest fourier transform in the west
#include <time.h> //time
#include <omp.h> //multithreading
#include <dc.h> //multithreaded mersenne twister random number generator
#include <limits.h> //limits of C variables
#include <dcmt_interface.h> //wrapper around dc mersenne twister

/* End Headers */

/*  Global Names  */
double * lxaxis , * lyaxis , * lzaxis ;
double * xaxis , * yaxis , * zaxis ;
/* End Global Names */

typedef struct solparamlist
{
	unsigned int nx , ny , nz ; double lx , ly , lz , solenoidrad ;
} solparam ;

typedef struct randparamlist
{
	unsigned int nx , ny , nz ; double lx , ly , lz , ran , var ;
} randparam ;

#define rnd(id) random_mt(id)

int magfieldgen_init ( int nthreads , int prime_id , randparam param ) ; //Inputs: Num_Threads, ID of Mersenne Prime , parameters for axes
/*
p is the exponent of the period. The period
     should be 2^p-1, but p must be an Mersenne exponent, i.e.,
     2^p-1 should be a prime. The list of usable p are as follows.
       521   607  1279  2203
       2281  3217  4253  4423
       9689  9941 11213 19937
       21701 23209 44497

       Any one of these can be selected using 1D array index in C style
*/
void linspace ( double * out , unsigned int n , double low , double high ) ; //Inputs: output array pointer, Nx , x_min , x_max
//Outputs the axis points, including the endpoint. Takes care of low > high, low == high and n == 0.
void fftfreq ( double * out , unsigned int n , double l ) ; //Inputs: output array pointer, Nx , X_max - X_min
//Outputs the Fourier sampling frequencies (Non-negative)

/* Distribution Function for the Random Field */
extern inline double dist ( double l , double ran )
{
	//currently Rayleigh Distribution
	double p = l * l ;
	return p * exp ( - p * ran ) ;
}

extern inline double powspec ( double l , double var , unsigned int npix , double psum , double ran ) ; //returns power spectrum

double psumcalc ( randparam param , int nthreads ) ;

int randomfield ( randparam param , int num_thread , double * b1 , double * b2 , double * b3 ) ;
// returns random field in contiguous arrays, with maximum magnitude normalized to 1
int magfieldgen ( int nx , int ny , int nz , double lx , double ly , double lz , double ran , double var , int nthread , double * bx , double * by , double * bz ) ;
#endif //magfieldgen_h