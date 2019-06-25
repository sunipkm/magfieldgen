#include <stdio.h> //standard library
#include <math.h> //math library
#include <stdlib.h> //allocations
#include <fftw3.h> //fastest fourier transform in the west
#include <time.h> //time
#include <omp.h> //multithreading
#include <dc.h> //multithreaded mersenne twister random number generator
#include <limits.h> //limits of C variables

typedef struct paramlist
{
	int nx , ny , nz ; double lx , ly , lz , ran , var ;
} plist ;

#include <magfieldgen.h> //inline header



double * lxaxis = NULL , * lyaxis = NULL , * lzaxis = NULL ; //global names, beware
double * xaxis = NULL , * yaxis = NULL , * zaxis = NULL ;

/*random number generator setting up*/
mt_struct ** mts = NULL ; //global mts struct

int generator ( int num_thread ) //thread safe random states with periodicity 2^521
{
	//printf ( "Dynamic creator called, creating %d states...\n" , num_thread ) ;
	mts = ( mt_struct ** ) malloc ( num_thread * sizeof ( mt_struct * ) ) ;
	for ( int i = 0 ; i < num_thread ; i++ )
	{
		mts [ i ] = get_mt_parameter_id_st ( 32 , 521 , rand ( ) % 65536 , 4172 ) ;

		if ( mts [ i ] == NULL )
			return ( -100 ) ;
		else
			sgenrand_mt ( 3241 + 40 * i , mts [ i ] ) ;
	}
	return 0 ;
}

double rnd ( int id ) //random number between 0 and 1
{
	return ( ( double ) genrand_mt ( mts [ id ] ) / 4294967295.0 ) ;
}

/*void fftfreq ( double * in , int n , double d ) //generates the fft sampling frequencies. refer to numpy fft fftfreq (python)
{
	printf ( "FFTFREQ called, with n = %d and d = %lf.\n" , n , d ) ;
	int i ;
	d *= n ;
	if ( n % 2 == 1 )
	{
		for ( i = 0 ; i <= n / 2 ; i++ )
		{
			in [ i ] = i / d ;
		}
		int j = 1 ;
		for ( ; i < n ; i++ )
		{
			in [ i ] = - in [ i - j ] ;
			j += 2 ;
		}
	}
	else
	{
		for ( i = 0 ; i < n / 2 ; i++ )
			in [ i ] = i / d ;
		int j = 0 ;
		for ( ; i < n ; i++ )
		{
			in [ i ] = ( j - i ) / d ;
			j += 2 ;
		}
	}
	return ;
}*/

void fftfreq ( double * in , int n , double d ) //generates the fft sampling frequencies. refer to numpy fft fftfreq (python)
{
	//printf ( "FFTFREQ called, with n = %d and d = %lf.\n" , n , d ) ;
	int i ;
	d *= n ;
	// if ( n % 2 == 0 )
	// {
	// 	for ( i = 0 ; i < n / 2 ; i++ )
	// 		in [ i ] = i / d ;
	// 	int j = 0 ;
	// 	for ( ; i < n ; i++ )
	// 	{
	// 		in [ i ] = ( j - i ) / d ;
	// 		j += 2 ;
	// 	}
	// }
	// else
	// {
	// 	for ( i = 0 ; i <= n / 2 ; i++ )
	// 	{
	// 		in [ i ] = i / d ;
	// 	}
	// 	int j = 1 ;
	// 	for ( ; i < n ; i++ )
	// 	{
	// 		in [ i ] = - in [ i - j ] ;
	// 		j += 2 ;
	// 	}
	// }
	for ( i = 0 ; i < n ; i++ )
		in [ i ] = i / d ; 

	return ;
}

double rayleigh ( double l , double ran ) //rayleigh distribution sigma
{
	double p = l * l ;
	return p * exp ( - p * ran ) ;
}

double powspec ( double l , double var , int npix , double psum , double ran )
{
	if ( l <= 0.0 )
		return 0.0 ;
	else
		return var * npix * 0.50 * rayleigh ( l , ran ) / psum ;
}

double psumcalc ( plist param , int num_thread )
{
	int nx = param . nx , ny = param . ny , nz = param . nz ;
	//double lx = param . lx , ly = param . ly , lz = param . lz ;
	double ran = param . ran ;

	//printf ( "PSUMCALC called, parameters: nx = %d , ny = %d , nz = %d , lx = %lf , ly = %lf , lz = %lf , ran = %lf, threads = %d.\n" , nx , ny , nz , lx , ly , lz , ran , num_thread ) ;

	/* WE ARE HOPING THAT THE FUNCTION THAT CALLS PSUMCALC WILL CREATE THE K AXES
		BEFORE PSUMCALC IS CALLED. THAT SAVES COMPUTATION TIME. THE THREE POINTERS
		ARE THEREFORE MADE GLOBAL AND ALLOCATED IN THE MAIN CALLING FUNCTION MAGFIELD
		BEWARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


	//lx /= nx ; ly /= ny ; lz /= nz ;
	//double * lxaxis = ( double * ) malloc ( nx * sizeof ( double ) ) ,
	//* lyaxis = ( double * ) malloc ( ny * sizeof ( double ) ) ,
	//* lzaxis = ( double * ) malloc ( nz * sizeof ( double ) ) ;

	//fftfreq ( lxaxis , nx , resx ) ;
	//fftfreq ( lyaxis , ny , resy ) ;
	//fftfreq ( lzaxis , nz , resz ) ;

	omp_set_num_threads ( num_thread ) ;

	double sum = 0 ;
	#pragma omp parallel for reduction (+:sum)
	for ( int i = 0 ; i < nx ; i++ )
	{
		for ( int j = 0 ; j < ny ; j++ )
		{
			for ( int k = 0 ; k < nz ; k++ )
			{
				sum += rayleigh ( sqrt ( lxaxis [ i ] * lxaxis [ i ] + lyaxis [ j ] * lyaxis [ j ] + lzaxis [ k ] * lzaxis [ k ] ) , ran ) ;
			}
		}
	}
	//free ( lxaxis ) ; free ( lyaxis ) ; free ( lzaxis ) ;
	return sum ;
}

int randomfield ( plist param , int num_thread , double * b1 , double * b2 , double * b3 )
{
	//retrieving variables from parameter
	int status = 0 ;
	printf ( "Randomfield called in C!\n" ) ;
	printf ( "Supplied parameters: " ) ;
	int nx = param . nx , ny = param . ny , nz = param . nz ;
	double lx = param . lx , ly = param . ly , lz = param . lz ;
	double ran = param . ran , var = param . var ;
	printf ( "nx = %d , ny = %d , nz = %d , lx = %lf , ly = %lf , lz = %lf , ran = %lf , var = %lf , threads = %d\n" , nx , ny , nz , lx , ly , lz , ran , var , num_thread ) ;

	//allocating the global pointers and creating the sampling frequency array
	lx /= nx ; ly /= ny ; lz /= nz ;
	lxaxis = ( double * ) malloc ( nx * sizeof ( double ) ) ;
	lyaxis = ( double * ) malloc ( ny * sizeof ( double ) ) ;
	lzaxis = ( double * ) malloc ( nz * sizeof ( double ) ) ;
	if ( lxaxis == NULL )
		return ( -55 ) ;
	/*else
		printf ( "Malloc for axes failed!\n" ) ;*/

	fftfreq ( lxaxis , nx , lx ) ;
	fftfreq ( lyaxis , ny , ly ) ;
	fftfreq ( lzaxis , nz , lz ) ;

	double psum = psumcalc ( param , num_thread ) ;
	int npix = nx * ny * nz ;

	//setting up parallel stuff
	status = generator ( num_thread ) ;
	if ( status != 0 )
		return status ;
	omp_set_num_threads ( num_thread ) ;

	//setting up the ax , ay , az as row major 3D arrays of complex data
	fftw_complex * ax = NULL , * ay = NULL , * az = NULL ;
	ax = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	ay = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	az = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;

	if ( ax == NULL )
		return ( -45 ) ;
	else
		printf ( "Magnetic potential x allocation successful!\n" ) ;
	if ( ay == NULL )
		return ( -45 ) ;
	else
		printf ( "Magnetic potential y allocation successful!\n" ) ;
	if ( az == NULL )
		return ( -45 ) ;
	else
		printf ( "Magnetic potential z allocation successful!\n" ) ;


	//double u , v ;
	printf ( "Generating AX...\n" ) ;
	int offsetx = nz * ny , offsety = nz ;
	#pragma omp parallel for// private(u,v)
	for ( int coord = 0 ; coord < npix ; coord++ )
	//for ( int x = 0 ; x < nx ; x++ )
	{
		//for ( int y = 0 ; y < ny ; y++ )
		//{
			//for ( int z = 0 ; z < nz ; z++ )
			//{
				int z = coord % offsety ;
				int x = coord / offsetx ;
				int y = ( coord - offsetx * x ) / offsety ;
				int id = omp_get_thread_num ( ) ;
				double l = sqrt ( lxaxis [ x ] * lxaxis [ x ] + lyaxis [ y ] * lyaxis [ y ] + lzaxis [ z ] * lzaxis [ z ] ) ;
				double sigma = sqrt ( powspec ( l , var , npix , psum , ran ) ) ;
				double s ;
				double u , v , z1 , z2 , fac ;
				//int coord = nz * ny * x + nz * y + z ;


				do
				{
					u = - 1 + 2 * rnd ( id ) ;
					v = - 1 + 2 * rnd ( id ) ;
					s = u * u + v * v ;
				} while ( s > 1.0 ) ;

				fac = sqrt ( - 2.0 * log ( s ) / s ) ;

				z1 = u * fac * sigma ;
				z2 = v * fac * sigma ;

				ax [ coord ] [ 0 ] = z1 ;
				ax [ coord ] [ 1 ] = z2 ;

				do
				{
					u = - 1 + 2 * rnd ( id ) ;
					v = - 1 + 2 * rnd ( id ) ;
					s = u * u + v * v ;
				} while ( s > 1.0 ) ;

				fac = sqrt ( - 2.0 * log ( s ) / s ) ;

				z1 = u * fac * sigma ;
				z2 = v * fac * sigma ;

				ay [ coord ] [ 0 ] = z1 ;
				ay [ coord ] [ 1 ] = z2 ;

				do
				{
					u = - 1 + 2 * rnd ( id ) ;
					v = - 1 + 2 * rnd ( id ) ;
					s = u * u + v * v ;
				} while ( s > 1.0 ) ;

				fac = sqrt ( - 2.0 * log ( s ) / s ) ;

				z1 = u * fac * sigma ;
				z2 = v * fac * sigma ;

				az [ coord ] [ 0 ] = z1 ;
				az [ coord ] [ 1 ] = z2 ;
			//}
		//}
	}
	printf ( "AX, AY, AZ generated.\n" ) ;
	/*#pragma omp parallel for
	for ( int x = 0 ; x < nx ; x++ )
	{
		for ( int y = 0 ; y < ny ; y++ )
		{
			for ( int z = 0 ; z < nz ; z++ )
			{
				int id = omp_get_thread_num ( ) ;
				double l = sqrt ( lxaxis [ x ] * lxaxis [ x ] + lyaxis [ y ] * lyaxis [ y ] + lzaxis [ z ] * lzaxis [ z ] ) ;
				double sigma = sqrt ( powspec ( l , var , npix , psum , ran ) ) ;
				double s , u , v ;
				do
				{
					u = - 1 + 2 * rnd ( id ) ;
					v = - 1 + 2 * rnd ( id ) ;
					s = u * u + v * v ;
				} while ( s > 1.0 ) ;

				double fac = sqrt ( - 2.0 * log ( s ) / s ) ;

				double z1 = u * fac * sigma ;
				double z2 = v * fac * sigma ;
				int coord = 2 * nz * ny * x + 2 * nz * y + 2 * z ;
				ay [ coord ] = z1 ;
				ay [ coord++ ] = z2 ;
			}
		}
	}

	printf ( "AY generated.\nGenerating AZ...\n" ) ;
	#pragma omp parallel for
	for ( int x = 0 ; x < nx ; x++ )
	{
		for ( int y = 0 ; y < ny ; y++ )
		{
			for ( int z = 0 ; z < nz ; z++ )
			{
				int id = omp_get_thread_num ( ) ;
				double l = sqrt ( lxaxis [ x ] * lxaxis [ x ] + lyaxis [ y ] * lyaxis [ y ] + lzaxis [ z ] * lzaxis [ z ] ) ;
				double sigma = sqrt ( powspec ( l , var , npix , psum , ran ) ) ;
				double s , u , v ;
				do
				{
					u = - 1 + 2 * rnd ( id ) ;
					v = - 1 + 2 * rnd ( id ) ;
					s = u * u + v * v ;
				} while ( s > 1.0 ) ;

				double fac = sqrt ( - 2.0 * log ( s ) / s ) ;

				double z1 = u * fac * sigma ;
				double z2 = v * fac * sigma ;
				int coord = 2 * nz * ny * x + 2 * nz * y + 2 * z ;
				az [ coord ] = z1 ;
				az [ coord++ ] = z2 ;
			}
		}
	}*/
	//printf ( "Trying to free the global axes..." ) ;
	//printf ( "Freed the global axes.\nPreparing to allocate memory to the FFTW arrays for cross product...\n" ) ;
	//cross product! >.<
	fftw_complex * axx , * ayy , * azz ;
	axx = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	ayy = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	azz = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	if ( axx == NULL )
		return ( -47 ) ;
	//printf ( "FFTW Malloc for the cross product arrays successful!\n" ) ;
	/*else
		printf ( "FFTW Malloc for the cross product arrays unsuccessful.\n" ) ;
	printf ( "Doing the cross product...\n" ) ;*/
	#pragma omp parallel for
	for ( int coord = 0 ; coord < npix ; coord++ )
	//for ( int x = 0 ; x < nx ; x++ )
	{
		//for ( int y = 0 ; y < ny ; y++ )
		//{
			//for ( int z = 0 ; z < nz ; z++ )
			//{
				int z = coord % offsety ;
				int x = coord / offsetx ;
				int y = ( coord - offsetx * x ) / offsety ;
				//int coord = nz * ny * x + nz * y + z ;
				axx [ coord ] [ 0 ] = lzaxis [ z ] * ay [ coord ] [ 1 ] - lyaxis [ y ] * az [ coord ] [ 1 ] ;
				axx [ coord ] [ 1 ] = - lzaxis [ z ] * ay [ coord ] [ 0 ] + lyaxis [ y ] * az [ coord ] [ 0 ] ;

				ayy [ coord ] [ 0 ] = lxaxis [ x ] * az [ coord ] [ 1 ] - lzaxis [ z ] * ax [ coord ] [ 1 ] ;
				ayy [ coord ] [ 1 ] = - lxaxis [ x ] * az [ coord ] [ 0 ] + lzaxis [ z ] * ax [ coord ] [ 0 ];

				azz [ coord ] [ 0 ] = lyaxis [ y ] * ax [ coord ] [ 1 ] - lxaxis [ x ] * ay [ coord ] [ 1 ] ;
				azz [ coord ] [ 1 ] = - lyaxis [ y ] * ax [ coord ] [ 0 ] + lxaxis [ x ] * ay [ coord ] [ 0 ] ;
			//}
		//}
	}
	printf ( "Cross product done.\n" ) ;
	free ( lxaxis ) ;
	printf("Lxaxis freed\n") ;
	free ( lyaxis ) ; printf("Lyaxis freed\n") ;
	free ( lzaxis ) ; printf("Lzaxis freed\n") ;
	fftw_free ( ax ) ; printf("AX freed\n") ;
	fftw_free ( ay ) ; printf("AY freed\n") ;
	fftw_free ( az ) ; printf("AZ freed\n") ;
	printf ( "Freed the magnetic potentials.\n" ) ;
	//fft
	fftw_plan plan ;
	fftw_complex * bx ;
	bx = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	plan = fftw_plan_dft_3d ( nx , ny , nz , axx , bx , FFTW_BACKWARD, FFTW_ESTIMATE ) ; //printf ( "planning done\n") ;
	fftw_execute ( plan ) ; //printf("fft\n") ;
	fftw_free ( axx ) ; fftw_destroy_plan ( plan ) ;
	//printf ( "Done FFT of X component, freed the potential and destroyed plan..\n" ) ;

	fftw_complex * by ;
	by = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	plan = fftw_plan_dft_3d ( nx , ny , nz , ayy , by , FFTW_BACKWARD, FFTW_ESTIMATE ) ; //printf ( "planning done\n") ;
	fftw_execute ( plan ) ; //printf("fft\n") ;
	fftw_free ( ayy ) ; fftw_destroy_plan ( plan ) ;
	//printf ( "Done FFT of Y component, freed the potential and destroyed plan..\n" ) ;
	fftw_complex * bz ;
	bz = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	plan = fftw_plan_dft_3d ( nx , ny , nz , azz , bz , FFTW_BACKWARD, FFTW_ESTIMATE ) ; //printf ( "planning done\n") ;
	fftw_execute ( plan ) ; //printf("fft\n") ;
	fftw_free ( azz ) ; fftw_destroy_plan ( plan ) ;
	//printf ( "Done FFT of Z component, freed the potential and destroyed plan..\n" ) ;
	//int totaldim = nx * ny * nz ;
	//printf ( "Saving output...\n" ) ;
	double bmax = - INFINITY ;
	for ( int i = 0 ; i < npix ; i++ )
	{
		//b1 [ i ] = bx [ i ] [ 0 ] / npix ; //real part
		//b2 [ i ] = by [ i ] [ 0 ] / npix ;
		//b3 [ i ] = bz [ i ] [ 0 ] / npix ;
		double try = pow ( bx [ i ] [ 0 ] , 2 ) + pow ( by [ i ] [ 0 ] , 2 ) + pow ( bz [ i ] [ 0 ] , 2 ) ;
		if ( try > bmax )
			bmax = try ;
	}
	bmax = sqrt ( bmax ) ;
	#pragma omp parallel for
	for ( int i = 0 ; i < npix ; i++ )
	{
		b1 [ i ] = bx [ i ] [ 0 ] / bmax ; //real part
		b2 [ i ] = by [ i ] [ 0 ] / bmax ;
		b3 [ i ] = bz [ i ] [ 0 ] / bmax ;
	}
	fftw_free ( bx ) ; fftw_free ( by ) ; fftw_free ( bz ) ;
	//printf ( "All done! Exiting!\n" ) ;
	return status ;
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

int randomfieldgen ( int nx , int ny , int nz , double lx , double ly , double lz , double ran , double var , int nthread , double * bx , double * by , double * bz )
{
	printf ( "In C: Magfieldgen\n" ) ;
	//time_t start = time ( NULL ) ;
	printf ( "Supplied pointers: %p, %p, %p\n" , bx , by , bz ) ;
	plist params = { nx , ny , nz , lx , ly , lz , ran , var } ;
	int status = randomfield ( params , nthread , bx , by , bz ) ;
	//printf ( "%lf seconds to generate the normalized magnetic field.\n" , ( double ) ( time ( NULL ) - start ) ) ;
	return status ;
}

int mixfieldgen ( unsigned int nx , unsigned int ny , unsigned int nz , double lx , double ly , double lz , double ran , double var , double solenrad , double fracsol , double fracran , int nthread , double * bx , double * by , double * bz )
{
	unsigned int npix = nx * ny * nz ;
	omp_set_num_threads ( nthread ) ;
	double * solbx = ( double * ) malloc ( npix * sizeof ( double ) ) ;
	double * solby = ( double * ) malloc ( npix * sizeof ( double ) ) ;
	double * solbz = ( double * ) malloc ( npix * sizeof ( double ) ) ;

	double * ranbx = ( double * ) malloc ( npix * sizeof ( double ) ) ;
	double * ranby = ( double * ) malloc ( npix * sizeof ( double ) ) ;
	double * ranbz = ( double * ) malloc ( npix * sizeof ( double ) ) ;

	int stat = solenoidfield ( nx , ny , nz , lx , ly , lz , solenrad , solbx , solby , solbz , nthread ) ;
	stat = randomfieldgen ( nx , ny , nz , lx , ly , lz , ran , var , nthread , nthread , ranbx , ranby , ranbz ) ;

	#pragma omp parallel for
	for ( unsigned int i = 0 ; i < npix ; i++ )
	{
		bx [ i ] = fracsol * solbx [ i ] + fracran * ranbx [ i ] ;
		by [ i ] = fracsol * solby [ i ] + fracran * ranby [ i ] ;
		bz [ i ] = fracsol * solbz [ i ] + fracran * ranbz [ i ] ;
	}
	free ( solbx ) ; free ( solby ) ; free ( solbz ) ;
	free ( ranbx ) ; free ( ranby ) ; free ( ranbz ) ;  
	return 0 ;
}
int main ( void )
{
	double * bx , * by , * bz ;
	int nx = 100 , ny = 100 , nz = 100 ;
	double lx = 300 , ly = 300 , lz = 300 ;
	double ran = 100 , var = 1 ;
	int nthread = 2 ;
	bx = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) ) ;
	by = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) ) ;
	bz = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) ) ;
	plist params = { nx , ny , nz , lx , ly , lz , ran , var } ;
	int status = randomfield ( params , nthread , bx , by , bz ) ;
	/*for ( int i = 0 ; i < nx * ny * nz ; i++ )
	{
		printf ( "%lf\t%lf\t%lf\n" , bx [ i ] , by [ i ] , bz [ i ] ) ;
	}*/
	printf ( "Status : %s\n" , ( status == 0 ) ? "Success!" : "Failure!" ) ;
	free ( bx ) ; free ( by ) ; free ( bz ) ;
	return 0 ;
}
