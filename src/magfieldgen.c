//magfieldgen.c

#include <magfieldgen.h>

int magfieldgen_init ( int nthreads , int prime_id , randparam param ) //initialize DCMT Functions
{
	int status = init_dcmt ( nthreads ) ;
	status = status | generator_dcmt ( prime_id ) ;

	unsigned int nx = param . nx , ny = param . ny , nz = param . nz ;
	double lx = param . lx , ly = param . ly , lz = param . lz ;

	lxaxis = ( double * ) malloc ( nx * sizeof ( double ) ) ;
	lyaxis = ( double * ) malloc ( ny * sizeof ( double ) ) ;
	lzaxis = ( double * ) malloc ( nz * sizeof ( double ) ) ;

	xaxis = ( double * ) malloc ( nx * sizeof ( double ) ) ;
	yaxis = ( double * ) malloc ( ny * sizeof ( double ) ) ;
	zaxis = ( double * ) malloc ( nz * sizeof ( double ) ) ;

	fftfreq (lxaxis,nx,2*lx) ;
	fftfreq (lyaxis,ny,2*ly) ;
	fftfreq (lzaxis,nz,2*lz) ;

	linspace (xaxis,nx,-lx,lx) ;
	linspace (yaxis,ny,-ly,ly) ;
	linspace (zaxis,nz,-lz,lz) ;

	return status ;
}

void fftfreq ( double * out , unsigned int n , double l ) 
{
	while ( n-- )
		out [ n ] = n / l ;
	return ;
}

void linspace ( double * out , unsigned int n , double low , double high )
{
	n-- ;
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
	double d = ( high - low ) / n ;
	out [ n ] = high ;
	while ( n-- )
		out [ n ] = out [ n + 1 ] - d ;
	return ;
}

extern inline double powspec ( double l , double var , unsigned int npix , double psum , double ran )
{
	if ( l <= 0.0 )
		return 0.0 ;
	else
		return var * npix * 0.50 * dist ( l , ran ) / psum ;
}

double psumcalc ( randparam param , int nthreads )
{
	/* WE ARE HOPING THAT THE FUNCTION THAT CALLS PSUMCALC WILL CREATE THE K AXES
		BEFORE PSUMCALC IS CALLED. THAT SAVES COMPUTATION TIME. THE THREE POINTERS
		ARE THEREFORE MADE GLOBAL AND ALLOCATED IN THE MAIN CALLING FUNCTION MAGFIELD
		BEWARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	//unpacking
	unsigned int nx = param . nx , ny = param . ny , nz = param . nz ;
	//double lx = param . lx , ly = param . ly , lz = param . lz ;
	double ran = param . ran ;

	omp_set_num_threads ( nthreads ) ;

	double sum = 0 ;
	#pragma omp parallel for reduction (+:sum)
	for ( unsigned int i = 0 ; i < nx ; i++ )
	{
		for ( unsigned int j = 0 ; j < ny ; j++ )
		{
			for ( unsigned int k = 0 ; k < nz ; k++ )
			{
				sum += dist ( sqrt ( xaxis [ i ] * xaxis [ i ] + yaxis [ j ] * yaxis [ j ] + zaxis [ k ] * zaxis [ k ] ) , ran ) ;
			}
		}
	}
	//free ( lxaxis ) ; free ( lyaxis ) ; free ( lzaxis ) ;
	return sum ;
}

int randomfield ( randparam param , int num_thread , double * b1 , double * b2 , double * b3 )
{
	//retrieving variables from parameter
	int status = 0 ;
	printf ( "Randomfield called in C!\n" ) ;
	printf ( "Supplied parameters: " ) ;
	unsigned int nx = param . nx , ny = param . ny , nz = param . nz ;
	double lx = param . lx , ly = param . ly , lz = param . lz ;
	double ran = param . ran , var = param . var ;
	printf ( "nx = %u , ny = %u , nz = %u , lx = %lf , ly = %lf , lz = %lf , ran = %lf , var = %lf , threads = %d\n" , nx , ny , nz , lx , ly , lz , ran , var , num_thread ) ;

	double psum = psumcalc ( param , num_thread ) ;
	int npix = nx * ny * nz ;

	//setting up parallel stuff
	omp_set_num_threads ( num_thread ) ;

	fftw_plan plan ;
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
	//printf ( "Trying to free the global axes..." ) ;
	//printf ( "Freed the global axes.\nPreparing to allocate memory to the FFTW arrays for cross product...\n" ) ;
	//cross product! >.<

	/*fftw_complex * fax , * fay , * faz ;

	fax = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	fay = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;
	faz = ( fftw_complex * ) fftw_malloc ( npix * sizeof ( fftw_complex ) ) ;

	plan ;
	printf ("Rand: FFT Set 1, A(r) -> A(k)\n") ;

	plan = fftw_plan_dft_3d ( nx , ny , nz , ax , fax , FFTW_FORWARD , FFTW_ESTIMATE ) ;
	fftw_execute ( plan ) ;
	fftw_destroy_plan ( plan ) ;
	fftw_free ( ax ) ;
	ax = fax ;

	plan = fftw_plan_dft_3d ( nx , ny , nz , ay , fay , FFTW_FORWARD , FFTW_ESTIMATE ) ;
	fftw_execute ( plan ) ;
	fftw_destroy_plan ( plan ) ;
	fftw_free ( ay ) ;
	ay = fay ;

	plan = fftw_plan_dft_3d ( nx , ny , nz , az , faz , FFTW_FORWARD , FFTW_ESTIMATE ) ;
	fftw_execute ( plan ) ;
	fftw_destroy_plan ( plan ) ;
	fftw_free ( az ) ;
	az = faz ;
	printf("FFT done, freed A(r).\n");*/

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
	fftw_free ( ax ) ; printf("AX freed\n") ;
	fftw_free ( ay ) ; printf("AY freed\n") ;
	fftw_free ( az ) ; printf("AZ freed\n") ;
	printf ( "Freed A(k).\n" ) ;
	//fft
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

/*
int solenoidfield ( solparam param , double * b1 , double * b2 , double * b3 , int num_threads )
{
	omp_set_num_threads ( num_threads ) ;

	unsigned int nx = param . nx , ny = param . ny , nz = param . nz ;
	double lx = param . lx , ly = param . ly , lz = param . lz ,
	solenoidrad = param . solenoidrad ;

	unsigned int npix = nx * ny * nz ;

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
}*/

void destroy_magfield ( void )
{
	free ( lxaxis ) ; free ( lyaxis ) ; free ( lzaxis ) ;
	free ( xaxis ) ; free ( yaxis ) ; free ( zaxis ) ;
	destroy_dcmt ( ) ; 
}

int magfieldgen ( int nx , int ny , int nz , double lx , double ly , double lz , double ran , double var , int nthread , double * bx , double * by , double * bz )
{
	printf ( "In C: Magfieldgen\n" ) ;
	//time_t start = time ( NULL ) ;
	printf ( "Supplied pointers: %p, %p, %p\n" , bx , by , bz ) ;
	randparam params = { nx , ny , nz , lx , ly , lz , ran , var } ;
	magfieldgen_init ( nthread , 3 , params ) ;
	int status = randomfield ( params , nthread , bx , by , bz ) ;
	//printf ( "%lf seconds to generate the normalized magnetic field.\n" , ( double ) ( time ( NULL ) - start ) ) ;
	destroy_magfield ( ) ;
	return status ;
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
	randparam params = { nx , ny , nz , lx , ly , lz , ran , var } ;
	int status = randomfield ( params , nthread , bx , by , bz ) ;
	/*for ( int i = 0 ; i < nx * ny * nz ; i++ )
	{
		printf ( "%lf\t%lf\t%lf\n" , bx [ i ] , by [ i ] , bz [ i ] ) ;
	}*/
	printf ( "Status : %s\n" , ( status == 0 ) ? "Success!" : "Failure!" ) ;
	free ( bx ) ; free ( by ) ; free ( bz ) ;
	return 0 ;
}