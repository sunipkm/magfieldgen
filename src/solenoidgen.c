#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h> //fastest fourier transform in the west
#include <time.h> //time
#include <omp.h> //multithreading
//#include <dc.h> //multithreaded mersenne twister random number generator
#include <limits.h> //limits of C variables

#define NTHREADS 2

double * lxaxis = NULL , * lyaxis = NULL , * lzaxis = NULL ; //global names, beware
double * xaxis = NULL , * yaxis = NULL , * zaxis = NULL ;

void fftfreq ( double * in , unsigned int n , double d ) //generates the fft sampling frequencies. refer to numpy fft fftfreq (python)
{
	//printf ( "FFTFREQ called, with n = %d and d = %lf.\n" , n , d ) ;
	unsigned int i ;
	d *= n ;
	for ( i = 0 ; i < n ; i++ )
		in [ i ] = i / d ; 

	return ;
}

void linspace ( double * in , double low , double high , unsigned int n ) //ref numpy linspace, this is default mode endpoint=True
{
	unsigned int i = n-1 ;
	if ( low == high )
	{
		printf ( "Linspace: Low == High, returning.\n" ) ;
		return ;
	}
	if ( n == 0 )
	{
		printf ( "Linspace: n = 0, returning.\n" ) ;
		return ;
	}
	if ( high < low )
	{
		double temp = low ;
		low = high ;
		high = temp ;
	}
	double d = ( high - low ) / ( n - 1 ) ;
	in [ i ] = high ;
	/*while ( i )
		in [ i ] = in [ i-- ] - d ;*/
	while ( i-- )
	{
		in [ i ] = in [ i + 1 ] - d ;
	}
	return ;
}

int solenoidfield ( unsigned int nx , unsigned int ny , unsigned int nz , double lx , double ly , double lz , double soleniodrad , double * b1 , double * b2 , double * b3 , int num_threads )
{
	omp_set_num_threads ( num_threads ) ;

	unsigned int npix = nx * ny * nz ;

	xaxis = ( double * ) malloc ( nx * sizeof ( double ) ) ;
	yaxis = ( double * ) malloc ( ny * sizeof ( double ) ) ;
	zaxis = ( double * ) malloc ( nz * sizeof ( double ) ) ;

	linspace ( xaxis , - lx , lx , nx ) ;
	linspace ( yaxis , - ly , ly , ny ) ;
	linspace ( zaxis , - lz , lz , nz ) ;

	fftw_complex * ax = NULL , * ay = NULL , * az = NULL ;
	ax = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	ay = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	az = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;

	unsigned int offsetx = nz * ny , offsety = nz ;
	//#pragma omp parallel for //generate AX for solenoidal field
	for ( unsigned int x = 0 ; x < nx ; x++ )
	{
		for ( unsigned int y = 0 ; y < ny ; y++ )
		{
			double rho1 = sqrt ( xaxis [ x ] * xaxis [ x ] + yaxis [ y ] * yaxis [ y ] ) ;
			double rho = rho1 > solenoidrad ? 0 : rho1 ; 
			double theta = atan2 ( yaxis [ y ] , xaxis [ x ] ) ;
			unsigned int coord = offsetx * x + nz * y ;
			#pragma omp parallel for
			for ( unsigned int z = 0 ; z < nz ; z++ )
			{
				unsigned int coord2 = coord + z ;
				ax [ coord2 ] [ 0 ] = - rho * sin ( theta ) ;
				ay [ coord2 ] [ 0 ] = rho * cos ( theta ) ;
				az [ coord2 ] [ 0 ] = 0 ;
				ax [ coord2 ] [ 1 ] = 0;
				ay [ coord2 ] [ 1 ] = 0 ;
				az [ coord2 ] [ 1 ] = 0 ;
			}
		}
	}

	free ( xaxis ) ; free ( yaxis ) ; free ( zaxis ) ;

	lxaxis = ( double * ) malloc ( nx * sizeof ( double ) ) ;
	lyaxis = ( double * ) malloc ( ny * sizeof ( double ) ) ;
	lzaxis = ( double * ) malloc ( nz * sizeof ( double ) ) ;

	fftfreq ( lxaxis , nx , 2 * lx / nx ) ;
	fftfreq ( lyaxis , ny , 2 * ly / ny ) ;
	fftfreq ( lzaxis , nz , 2 * lz / nz ) ;

	fftw_complex * fax , * fay , * faz ;
	fax = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	fay = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	faz = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;

	
	fftw_plan plan ;
	plan = fftw_plan_dft_3d ( nx , ny , nz , ax , fax , FFTW_FORWARD , FFTW_ESTIMATE ) ;
	fftw_execute ( plan ) ;
	fftw_free ( ax ) ;
	fftw_destroy_plan ( plan ) ;

	plan = fftw_plan_dft_3d ( nx , ny , nz , ay , fay , FFTW_FORWARD , FFTW_ESTIMATE ) ;
	fftw_execute ( plan ) ;
	fftw_free ( ay ) ;
	fftw_destroy_plan ( plan ) ;

	plan = fftw_plan_dft_3d ( nx , ny , nz , az , faz , FFTW_FORWARD , FFTW_ESTIMATE ) ;
	fftw_execute ( plan ) ;
	fftw_free ( az ) ;
	fftw_destroy_plan ( plan ) ;


	//cross product
	fftw_complex * axx , * ayy , * azz ;
	axx = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	ayy = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	azz = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	#pragma omp parallel for
	for ( unsigned int coord = 0 ; coord < npix ; coord++ )
	//for ( int x = 0 ; x < nx ; x++ )
	{
		//for ( int y = 0 ; y < ny ; y++ )
		//{
			//for ( int z = 0 ; z < nz ; z++ )
			//{
				unsigned int z = coord % offsety ;
				unsigned int x = coord / offsetx ;
				unsigned int y = ( coord - offsetx * x ) / offsety ;
				//int coord = nz * ny * x + nz * y + z ;
				axx [ coord ] [ 0 ] = lzaxis [ z ] * fay [ coord ] [ 1 ] - lyaxis [ y ] * faz [ coord ] [ 1 ] ;
				axx [ coord ] [ 1 ] = - lzaxis [ z ] * fay [ coord ] [ 0 ] + lyaxis [ y ] * faz [ coord ] [ 0 ] ;

				ayy [ coord ] [ 0 ] = lxaxis [ x ] * faz [ coord ] [ 1 ] - lzaxis [ z ] * fax [ coord ] [ 1 ] ;
				ayy [ coord ] [ 1 ] = - lxaxis [ x ] * faz [ coord ] [ 0 ] + lzaxis [ z ] * fax [ coord ] [ 0 ];

				azz [ coord ] [ 0 ] = lyaxis [ y ] * fax [ coord ] [ 1 ] - lxaxis [ x ] * fay [ coord ] [ 1 ] ;
				azz [ coord ] [ 1 ] = - lyaxis [ y ] * fax [ coord ] [ 0 ] + lxaxis [ x ] * fay [ coord ] [ 0 ] ;
			//}
		//}
	}

	free ( lxaxis ) ; free ( lyaxis ) ; free ( lzaxis ) ;
	fftw_free ( fax ) ; fftw_free ( fay ) ; fftw_free ( faz ) ; 
 
	//fftw_plan plan ;
	fftw_complex * bx ;
	bx = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	plan = fftw_plan_dft_3d ( nx , ny , nz , axx , bx , FFTW_BACKWARD, FFTW_ESTIMATE ) ; //printf ( "planning done\n") ;
	fftw_execute ( plan ) ; //printf("fft\n") ;
	fftw_free ( axx ) ; fftw_destroy_plan ( plan ) ;
	printf ( "Done FFT of X component, freed the potential and destroyed plan..\n" ) ;

	fftw_complex * by ;
	by = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	plan = fftw_plan_dft_3d ( nx , ny , nz , ayy , by , FFTW_BACKWARD, FFTW_ESTIMATE ) ; //printf ( "planning done\n") ;
	fftw_execute ( plan ) ; //printf("fft\n") ;
	fftw_free ( ayy ) ; fftw_destroy_plan ( plan ) ;
	printf ( "Done FFT of Y component, freed the potential and destroyed plan..\n" ) ;
	fftw_complex * bz ;
	bz = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	plan = fftw_plan_dft_3d ( nx , ny , nz , azz , bz , FFTW_BACKWARD, FFTW_ESTIMATE ) ; //printf ( "planning done\n") ;
	fftw_execute ( plan ) ; //printf("fft\n") ;
	fftw_free ( azz ) ; fftw_destroy_plan ( plan ) ;

	#pragma omp parallel for
	for ( unsigned int i = 0 ; i < npix ; i++ )
	{
		b1 [ i ] = bx [ i ] [ 0 ] / npix ; //real part, normalized
		b2 [ i ] = by [ i ] [ 0 ] / npix ;
		b3 [ i ] = bz [ i ] [ 0 ] / npix ;
	}
	fftw_free ( bx ) ; fftw_free ( by ) ; fftw_free ( bz ) ;

	return 0 ;
}