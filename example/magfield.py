"""
sys and numpy are required packages. Absence of these packages must not be circumvented.
Use of sys can be circumvented if 'main' program is not used.

Code for simulation of a random magnetic field, and some analyses of it.

Written by: Sunip K. Mukherjee, under guidance of Prof. Ritaban Chatterjee
			Presidency University, Kolkata
			
			Feb 03 2016

Bug report: sunipkmukherjee@gmail.com
"""
#-------------- COLOR CHARACTERS ---------------------------------------------------------
try :
	import os
	os_success = True
except ImportError :
	os_success = False
	
if os_success :
	if os . name == 'posix' :
		os_win = False
	else :
		os_win = True
else :
	os_win = True
if not os_win :
	HEAD = '\033[95m'#magenta
	OKB = '\033[94m'#blue
	OKG = '\033[92m'#green
	WARN = '\033[93m'#yellow ochre
	FAIL = '\033[91m'#red
	ENDC = '\033[0m'#ends an effect. So, use like HEAD+ULINE+"Hello"+ENDC+ENDC
	BOLD = '\033[1m'
	ULINE = '\033[4m'
else :
	HEAD = '' #blanks
	OKB = ''
	OKG = ''
	WARN = ''
	FAIL = ''
	ENDC = ''
	BOLD = ''
	ULINE = ''
#-----------------------------------------------------------------------------------------

#-------------- IMPORTING NUMPY FOR THE FUNCTIONS ----------------------------------------
try :
	import numpy as np
except ImportError :
	print FAIL+BOLD+"The most important package 'numpy' is missing. You can not proceed. Install numpy through 'pip' etc and try again."+ENDC+ENDC
	raise ImportError ( "No module named numpy" )
#-----------------------------------------------------------------------------------------

#-------------- IMPORTING SCIPY RELATED PARTS --------------------------------------------
scipypresence = True
try :
	import scipy
except ImportError :
	print WARN+BOLD+"Scipy is missing. You can not plot the wake or the light curve."+\
																				ENDC+ENDC
	scipypresence = False
if scipypresence :
	#from scipy . special import kv #Bessel function of the second kind
	from scipy . integrate import quad #Integrator
#-----------------------------------------------------------------------------------------

#-------------- IMPORTING INTEGRATION TABLE FOR KV ---------------------------------------
tab_int = np . load ( 'table_int.npy' )
#-----------------------------------------------------------------------------------------

#-------------- NUMERICAL PART OF CODE ---------------------------------------------------

#-------------- CRITICAL MAGFIELD PART OF CODE -------------------------------------------
try :
	from magfieldgen import *
except ImportError :
	print WARN+"The C backend library for computing the random magnetic field is either absent or not properly configured. The default Python implementation will be used instead, which is almost 100x slower. You may consider reconfiguration before executing this script. However, the functionalities provided remain the same."+ENDC
#-------------- RAYLEIGH SIGMA -----------------------------------------------------------
	def rayleigh ( l , ran ) : #rayleigh distribution
		return l * l * np . exp ( - l * l * ran )
#-----------------------------------------------------------------------------------------

#-------------- POWER SPECTRUM (NORMALIZED) ----------------------------------------------
	def powspec ( L , var , npix , psum , ran ) : #returning a sigma to use
		if L <= 0 :
			return 0.
		else :
			return var * npix * npix * 0.5 * rayleigh ( L , ran ) / psum
#-----------------------------------------------------------------------------------------

#-------------- CALCULATE TOTAL SUM FOR NORMALIZATION ------------------------------------
	def psumcalc ( nx , ny , nz , lx , ly , lz , pow ) : #calculates the normalization, pow is equivalent to range
		
		resx = float ( lx ) / float ( nx )
		resy = float ( ly ) / float ( ny )
		resz = float ( lz  ) / float ( nz )
		
		lxaxis = np . fft . fftfreq ( nx , resx )
		lyaxis = np . fft . fftfreq ( ny , resy )
		lzaxis = np . fft . fftfreq ( nz , resz )
		
		summ = 0
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				for k in xrange ( nz ) :
					summ += rayleigh ( ( np . sqrt ( lxaxis [ i ] * lxaxis [ i ] + lyaxis [ j ] * lyaxis [ j ] + lzaxis [ k ] * lzaxis 	[ k ] ) ) , pow )
	
		return summ
#-----------------------------------------------------------------------------------------

#-------------- MAGNETIC FIELD GENERATOR FUNCTION ----------------------------------------	
	def magfield ( nx , ny , nz , lx , ly , lz , pow , var ) : #magnetic field generator
		#resolutions
		resx = float ( lx ) / float ( nx )
		resy = float ( ly ) / float ( ny )
		resz = float ( lz  ) / float ( nz )
		#k space axes
		lxaxis = np . fft . fftfreq ( nx , resx )
		lyaxis = np . fft . fftfreq ( ny , resy )
		lzaxis = np . fft . fftfreq ( nz , resz )
		#normalization
		psum = psumcalc ( nx , ny , nz , lx , ly , lz , pow )
		#imaginary unit
		jj = 0. + 1j
		#magnetic potential components
		ax = np . zeros ( ( nx , ny , nz ) , dtype = complex )
		ay = np . zeros ( ( nx , ny , nz ) , dtype = complex )
		az = np . zeros ( ( nx , ny , nz ) , dtype = complex )
		#generating ax
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				for k in xrange ( nz ) :
					i1 = i - nx / 2 #defining coordinates
					j1 = j - ny / 2
					k1 = k - nz / 2
					
					#i2 = - ( i1 + nx / 2 ) #conjugates
					#j2 = - ( j1 + ny / 2 )
					#k2 = - ( k1 + nz / 2 )
					
					l = ( np . sqrt ( lxaxis [ i ] * lxaxis [ i ] + lyaxis [ j ] * lyaxis [ j ] + lzaxis [ k ] * lzaxis [ k ] ) ) #k = sqrt (kx**2 + ky**2 + kz**2 )
					
					sigma = np . sqrt ( powspec ( l , var , nx * ny * nz , psum , pow ) ) #width 
					
					s = 1.1 #just getting in the while loop
				
					while s > 1 : #box muller
						u = np . random . uniform ( - 1. , 1. ) 
						v = np . random . uniform ( - 1. , 1. )
					
						s = u * u + v * v
					
					fac = np . sqrt ( - 2. * np . log ( s ) / s )
				
					z1 = u * fac * sigma #real
					z2 = v * fac * sigma #imaginary
				
					ax [ i1 ] [ j1 ] [ k1 ] = z1 + jj * z2
					
			
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				for k in xrange ( nz ) :
					i1 = i - nx / 2
					j1 = j - ny / 2
					k1 = k - nz / 2

					l = ( np . sqrt ( lxaxis [ i ] * lxaxis [ i ] + lyaxis [ j ] * lyaxis [ j ] + lzaxis [ k ] * lzaxis [ k ] ) )
				
					sigma = np . sqrt ( powspec ( l , var , nx * ny * nz , psum , pow ) )
				
					s = 1.1
				
					while s > 1 :
						u = np . random . uniform ( - 1. , 1. )
						v = np . random . uniform ( - 1. , 1. )
						
						s = u * u + v * v
					
					fac = np . sqrt ( - 2. * np . log ( s ) / s )
				
					z1 = u * fac * sigma
					z2 = v * fac * sigma
				
					ay [ i1 ] [ j1 ] [ k1 ] = z1 + jj * z2
					
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				for k in xrange ( nz ) :
					i1 = i - nx / 2
					j1 = j - ny / 2
					k1 = k - nz / 2
				
					l = ( np . sqrt ( lxaxis [ i ] * lxaxis [ i ] + lyaxis [ j ] * lyaxis [ j ] + lzaxis [ k ] * lzaxis [ k ] ) )
					
					sigma = np . sqrt ( powspec ( l , var , nx * ny * nz , psum , pow ) )
				
					s = 1.1
					
					while s > 1 :
						u = np . random . uniform ( - 1. , 1. )
						v = np . random . uniform ( - 1. , 1. )
						
						s = u * u + v * v
					
					fac = np . sqrt ( - 2. * np . log ( s ) / s )
					
					z1 = u * fac * sigma
					z2 = v * fac * sigma

					az [ i1 ] [ j1 ] [ k1 ] = z1 + jj * z2
			
		#initializing the ik x A components
		ax1 = np . zeros ( ( nx , ny , nz ) , dtype = complex )
		ay1 = np . zeros ( ( nx , ny , nz ) , dtype = complex )
		az1 = np . zeros ( ( nx , ny , nz ) , dtype = complex )
		#cross product
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				for k in xrange ( nz ) :
					i1 = i# - nx / 2 
					j1 = j# - ny / 2
					k1 = k# - nz / 2
				
					x = lxaxis [ i ]
					y = lyaxis [ j ]
					z = lzaxis [ k ]
				
					ax1 [ i1 ] [ j1 ] [ k1 ] = jj * ( y * az [ i1 ] [ j1 ] [ k1 ] - z * ay [ i1 ] [ j1 ] [ k1 ] )
					ay1 [ i1 ] [ j1 ] [ k1 ] = jj * ( - x * az [ i1 ] [ j1 ] [ k1 ] + z * ax [ i1 ] [ j1 ] [ k1 ] )
					az1 [ i1 ] [ j1 ] [ k1 ] = jj * ( x * ay [ i1 ] [ j1 ] [ k1 ] - y * ax [ i1 ] [ j1 ] [ k1 ] )
		#doing the inverse fourier transform			
		bx1 = np . fft . ifftn ( ax1 )
		by1 = np . fft . ifftn ( ay1 )
		bz1 = np . fft . ifftn ( az1 )
		#shifting to get 0,0 in proper place and taking the real part only
		bx = np . fft . ifftshift ( bx1 ) . real # np . abs ( np . fft . ifftshift ( bx ) ) #. real
		by = np . fft . ifftshift ( by1 ) . real # np . abs ( np . fft . ifftshift ( by ) ) #. real
		bz = np . fft . ifftshift ( bz1 ) . real # np . abs ( np . fft . ifftshift ( bz ) ) #. real
		#getting modulus
		b = np . sqrt ( bx * bx + by * by + bz * bz )
		#normalizing the maximum field strength to 1
		bmax = np . max ( b )
		bx /= bmax
		by /= bmax
		bz /= bmax
		#removing any mean
		#bx -= np . mean ( bx )
		#by -= np . mean ( by )
		#bz -= np . mean ( bz )	
		#return the field components
		return ( bx , by , bz )
#-----------------------------------------------------------------------------------------

#-------------- END CRITICAL MAGFIELD PART OF CODE ---------------------------------------

#-------------- OTHER NUMERICAL FUNCTIONS ------------------------------------------------

#-------------- KV INTEGRATOR ------------------------------------------------------------
"""
	This function integrates Bessel function of the second kind of order 5/3 
	from some lower limit to infinity using analytic forms and integral tables.
"""
def kvint ( lowlim ) :
	global tab_int
	
	if lowlim <= 0.001 : #analytic
		return 2.145 * ( ( lowlim ** ( - 0.66667 ) ) - ( 0.001 ** ( - 0.66667 ) ) ) + \
				213.125823
	elif lowlim < 40 : #table interpolation
		return np . interp ( lowlim , tab_int [ 0 ] , tab_int [ 1 ] )
	else : #vanishing
		return 0.	
#-----------------------------------------------------------------------------------------

#-------------- ENERGY FUNCTION ----------------------------------------------------------
def efunc ( x , f , t , s , B ) : #gamma, freq, time, power of power law energy dist ,
									# magnitude of magfield
	k1 = 4.2 * ( 10 ** 6 ) * B
	k2 = 1.3 * ( 10 ** - 9 ) * B * B
	lowlim = f / ( k1 * x * x )
	
	return ( ( x ** ( - s - 2 ) ) * ( ( 1 - x * k2 * t ) ** ( s - 2  ) ) * f *
																kvint ( lowlim ) / k1 )
#-----------------------------------------------------------------------------------------

#-------------- INTEGRATOR WRAPPER -------------------------------------------------------
def efuncint ( freq , n0 , lowlim , uplim , t , B , s = 2.5 ) :
	if t >= 0 :
		if lowlim + 0.000001 < uplim :
			result , _ = quad ( efunc , lowlim , uplim , args = ( freq , t , s , B , ) )
		else :
			result = 0.
	else :
		result = 0.
	return n0 * result
#-----------------------------------------------------------------------------------------

#-------------- END NUMERICAL PART OF CODE -----------------------------------------------

#-------------- POST PROCESSING FUNCTIONS ------------------------------------------------

#-------------- DIVERGENCE CALCULATOR ----------------------------------------------------
def divcalc ( bx , by , bz , nx , ny , nz , lx , ly , lz ) :
	"""
	if np . shape ( bz ) != np . shape ( by ) :
		raise TypeError
	resx = float ( lx ) / float ( nx )
	resy = float ( ly ) / float ( ny )
	resz = float ( lz ) / float ( nz )
	print OKG+"Calculating divergence..."+ENDC
	b1 = np . zeros ( ( nx , ny , nz ) , dtype = float )
	b2 = np . zeros ( ( nx , ny , nz ) , dtype = float )
	b3 = np . zeros ( ( nx , ny , nz ) , dtype = float )
	# solenoidal check
	for i in xrange ( 1 , nx - 1 ) :
		for j in xrange ( 1 , ny - 1 ) :
			for k in xrange ( 1 , nz - 1 ) :
				b1 [ i ] [ j ] [ k ] = bx [ i + 1 ] [ j ] [ k ] - bx [ i - 1 ] [ j ] [ k ]
				b1 [ i ] [ j ] [ k ] /= 2 * resx
		
	for i in xrange ( 1 , nx - 1 ) :
		for j in xrange ( 1 , ny - 1 ) :
			for k in xrange ( 1 , nz - 1 ) :
				b2 [ i ] [ j ] [ k ] = bx [ i ] [ j + 1 ] [ k ] - bx [ i ] [ j - 1 ] [ k ]
				b2 [ i ] [ j ] [ k ] /= 2 * resy
		
	for i in xrange ( 1 , nx - 1 ) :
		for j in xrange ( 1 , ny - 1 ) :
			for k in xrange ( 1 , nz - 1 ) :
				b3 [ i ] [ j ] [ k ] = bx [ i ] [ j ] [ k + 1 ] - bx [ i ] [ j ] [ k - 1 ]
				b3 [ i ] [ j ] [ k ] /= 2 * resz
					
	b = b1 + b2 + b3
	print np . max ( np . abs ( b ) )
	plt . figure ( )
	print np . sum ( b )
	print ULINE+"Divergence calculation is done."+ENDC
	#for i in xrange ( nx ) :
	#	plt . imshow ( b [ i ] )
	#	plt . colorbar ( ) 
	#	plt . show ( )
	plt . imshow ( b [ nx / 5 ] )
	plt . colorbar ( )
	plt . show ( )
	plt . figure ( )
	plt . imshow ( b [ nx / 2 ] )
	plt . colorbar ( )
	plt . show ( )
	"""
	lxaxis = np.array([ i / ( 2 * lx ) for i in xrange(nx)])
	lyaxis = np.array([ i / ( 2 * ly ) for i in xrange(ny)])
	lzaxis = np.array([ i / ( 2 * lz ) for i in xrange(nz)])

	kx , ky , kz = np.meshgrid(lxaxis,lyaxis,lzaxis,indexing='ij')

	divergence = (np.fft.ifftn((np.fft.fftn(bx)*kx*1j))+np.fft.ifftn((np.fft.fftn(by)*ky*1j))+np.fft.ifftn((np.fft.fftn(bz)*kz*1j))).real
	p3.contour3d(divergence,contours=10,transparent=True)
	p3.show()
	print np.max(divergence),np.min(divergence)
	return None
#-----------------------------------------------------------------------------------------

#-------------- REDUCTION CALCULATOR -----------------------------------------------------
def reduction_ ( bx , by , bz , reduction ) :
	print OKG+"Reducing..."+ENDC
	if np . shape ( bz ) != np . shape ( by ) :
		raise TypeError
	nx , ny , nz = np . shape ( bx )
	bxx = np . zeros ( ( nx / reduction , ny / reduction , nz / reduction ) )
	byy = np . zeros ( ( nx / reduction , ny / reduction , nz / reduction ) )
	bzz = np . zeros ( ( nx / reduction , ny / reduction , nz / reduction ) )
		
	for i in xrange ( nx / reduction ) :
		for j in xrange ( ny / reduction ) :
			for k in xrange ( nz / reduction ) :
				for ii in xrange ( i , i + reduction ) :
					for jj in xrange ( j , j + reduction ) :
						for kk in xrange ( k , k + reduction ) :
							bxx [ i ] [ j ] [ k ] += bx [ ii ] [ jj ] [ kk ]
							byy [ i ] [ j ] [ k ] += by [ ii ] [ jj ] [ kk ]
							bzz [ i ] [ j ] [ k ] += bz [ ii ] [ jj ] [ kk ]
			#bxx /= reduction ** 3
			#byy /= reduction ** 3
			#bzz /= reduction ** 3
	print ULINE+"Reduction complete."+ENDC
	return ( bxx , byy , bzz )	
#-----------------------------------------------------------------------------------------

#-------------- IMAGE CREATOR (FACE ONLY, SLICES WILL BE IMPLEMENTED) --------------------	
def image_calc ( bx , by , bz , name ) :
	print name
	nx , ny , nz = np . shape ( bx )
	check = raw_input ( "Do you want to plot the images? ( Y / N ): " )
	while ( check != 'n' and check != 'y' and check != 'N' and check != 'Y' ):
		check = raw_input ( "Enter again: " )
	if check == 'n' or check == 'N':
		return None
	elif check == 'y' or check == 'Y':
		n = input ( "Enter spectral index = " )
		#Z axis:
		img = np . zeros ( ( nx , ny ) )
		#front Z axis:
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				img [ i ] [ j ] = np . sqrt ( bx [ i ] [ j ] [ 0 ] * bx [ i ] [ j ] [ 0 ] + by [ i ] [ j ] [ 0 ] * by [ i ] [ j ] [ 0 ] ) ** ( n + 1 )
		plt . figure ( "Z axis back" )
		plt . imshow ( img )
		plt . colorbar ( )
		plt . show ( )
		#back Z axis:
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				img [ i ] [ j ] = np . sqrt ( ( bx [ i ] [ j ] [ nz - 1 ] * bx [ i ] [ j ] [ nz - 1 ] + by [ i ] [ j ] [ nz - 1 ] * by [ i ] [ j ] [ nz - 1 ] ) ) ** ( n + 1 )
		plt . figure ( "Z axis front" )
		plt . imshow ( img )
		plt . colorbar ( )
		plt . show ( )
		#Y axis:
		img = np . zeros ( ( nx , nz ) )
		#front :
		for i in xrange ( nx ) :
			for j in xrange ( nz ) :
				img [ i ] [ j ] = np . sqrt ( bx [ i ] [ 0 ] [ j ] * bx [ i ] [ 0 ] [ j ] + bz [ i ] [ 0 ] [ j ] * bz [ i ] [ 0 ] [ j ] ) ** ( n + 1 )
		plt . figure ( "Y axis back" )
		plt . imshow ( img )
		plt . colorbar ( )
		plt . show ( )
		#back :
		for i in xrange ( nx ) :
			for j in xrange ( nz ) :
				img [ i ] [ j ] = np . sqrt ( bx [ i ] [ ny - 1 ] [ j ] * bx [ i ] [ ny - 1 ] [ j ] + bz [ i ] [ ny - 1 ] [ j ] * bz [ i ] [ ny - 1 ] [ j ] ) ** ( n + 1 )
		plt . figure ( "Y axis front" )
		plt . imshow ( img )
		plt . colorbar ( )
		plt . show ( )
		#X axis:
		img = np . zeros ( ( ny , nz ) )
		#front :
		for i in xrange ( ny ) :
			for j in xrange ( nz ) :
				img [ i ] [ j ] = np . sqrt ( bx [ 0 ] [ i ] [ j ] * bx [ 0 ] [ i ] [ j ] + bz [ 0 ] [ i ] [ j ] * bz [ 0 ] [ i ] [ j ] ) ** ( n + 1 )
		plt . figure ( "X axis back" )
		plt . imshow ( img )
		plt . colorbar ( )
		plt . show ( )
		#back :
		for i in xrange ( ny ) :
			for j in xrange ( nz ) :
				img [ i ] [ j ] = np . sqrt ( bx [ nx - 1 ] [ i ] [ j ] * bx [ nx - 1 ] [ i ] [ j ] + bz [ nx - 1 ] [ i ] [ j ] * bz [ nx - 1 ] [ i ] [ j ] ) ** ( n + 1 )
		plt . figure ( "X axis front" )
		plt . imshow ( img )
		plt . colorbar ( )
		plt . show ( )
		return None
#-----------------------------------------------------------------------------------------		

#-------------- LINEAR MAGNITUDE CALCULATOR ----------------------------------------------
def lincalc ( bx , by , bz , nx , ny , nz , lx , ly , lz ) :
	print BOLD+"Axis-wise variation of the magnetic field"+ENDC
	#lxaxis = np . linspace ( - lx / 2 , lx / 2 , nx )
	#lyaxis = np . linspace ( - ly / 2 , ly / 2 , ny )
	#lzaxis = np . linspace ( - lz / 2 , lz / 2 , nz )
	
	a = np . zeros ( nz )
	#bxx = np . zeros ( ( 1 , 1 , nz ) )
	#byy = np . zeros ( ( 1 , 1 , nz ) )
	#bzz = np . zeros ( ( 1 , 1 , nz ) )
	for i in xrange ( nx ) :
		for j in xrange ( ny ) :
			for k in xrange ( nz ) :
				#bxx [ 0 ] [ 0 ] [ k ] += bx [ i ] [ j ] [ k ]
				#byy [ 0 ] [ 0 ] [ k ] += by [ i ] [ j ] [ k ]
				#bzz [ 0 ] [ 0 ] [ k ] += bz [ i ] [ j ] [ k ]
				a [ k ] += np . sqrt ( bx [ i ] [ j ] [ k ] * bx [ i ] [ j ] [ k ] + by [ i ] [ j ] [ k ] * by [ i ] [ j ] [ k ] )
	#for i in xrange ( nz ) :
		#a [ i ] = ( np . sqrt ( bxx * bxx + byy * byy ) ) [ 0 ] [ 0 ] [ i ]
	plt . figure ( "Z axis" )
	plt . plot ( lzaxis , a )
	plt . show ( )
	
	a = np . zeros ( ny )
	#bxx = np . zeros ( ( 1 , ny , 1 ) )
	#byy = np . zeros ( ( 1 , ny , 1 ) )
	#bzz = np . zeros ( ( 1 , ny , 1 ) )
	for i in xrange ( nx ) :
		for j in xrange ( ny ) :
			for k in xrange ( nz ) :
				#bxx [ 0 ] [ j ] [ 0 ] += bx [ i ] [ j ] [ k ]
				#byy [ 0 ] [ j ] [ 0 ] += by [ i ] [ j ] [ k ]
				#bzz [ 0 ] [ j ] [ 0 ] += bz [ i ] [ j ] [ k ]
				a [ j ] += np . sqrt ( bx [ i ] [ j ] [ k ] * bx [ i ] [ j ] [ k ] + bz [ i ] [ j ] [ k ] * bz [ i ] [ j ] [ k ] )
	#for i in xrange ( ny ) :			
		#a [ i ] = ( np . sqrt ( bxx * bxx + bzz * bzz ) ) [ 0 ] [ i ] [ 0 ]
	plt . figure ( "Y axis" )
	plt . plot ( lyaxis , a )
	plt . show ( )
	
	a = np . zeros ( nx )
	#bxx = np . zeros ( ( nx , 1 , 1 ) )
	#byy = np . zeros ( ( nx , 1 , 1 ) )
	#bzz = np . zeros ( ( nx , 1 , 1 ) )
	for i in xrange ( nx ) :
		for j in xrange ( ny ) :
			for k in xrange ( nz ) :
				#bxx [ i ] [ 0 ] [ 0 ] += bx [ i ] [ j ] [ k ]
				#byy [ i ] [ 0 ] [ 0 ] += by [ i ] [ j ] [ k ]
				#bzz [ i ] [ 0 ] [ 0 ] += bz [ i ] [ j ] [ k ]
				a [ i ] += np . sqrt ( bz [ i ] [ j ] [ k ] * bz [ i ] [ j ] [ k ] + by [ i ] [ j ] [ k ] * by [ i ] [ j ] [ k ] )
	#a = np . sqrt ( bzz * bzz + byy * byy )
	#for i in xrange ( nx ) :
		#a [ i ] = ( np . sqrt ( bzz * bzz + byy * byy ) ) [ i ] [ 0 ] [ 0 ]
		
	plt . figure ( "X axis" )
	plt . plot ( lxaxis , a )
	plt . show ( )
#-----------------------------------------------------------------------------------------

#-------------- NEW WAKE -----------------------------------------------------------------
class FWake :
	def __init__ ( self , bx , by , bz , fieldpar , energypar , wwidth , wstrength = 1 , cprob = 1 , dt = 1 ) :
		self . b10 = bx
		self . b20 = by
		self . b30 = bz
		nx , ny , nz , lx , ly , lz , pow , var = fieldpar
		self . sizepar = ( nx , ny , nz )
		self . realpar = ( lx , ly , lz )
		self . fieldpar = ( nx , ny , nz , lx , ly , lz , pow , var )
		self . energypar = energypar
		#energypar = n0_min, n0_delta , n0_tau , gamma_min , gamma_max , gamma_cut , gamma_tau , freq
		self . wpos = 0
		self . wwidth = wwidth
		self . b11 = np . zeros ( self . sizepar )
		self . b21 = np . zeros ( self . sizepar )
		self . b31 = np . zeros ( self . sizepar )
		self . timestamp = np . zeros ( self . sizepar ) + - 1
		self . wstrength = wstrength
		self . cprob = cprob
		self . mastertime = 0
		self . dt = dt
		
	def randfieldupd ( self ) : #time evolution of the random magnetic field
		nx , ny , nz , lx , ly , lz , pow , var = self . fieldpar
		b11 , b21 , b31 = magfield ( nx , ny , nz , lx , ly , lz , pow , var )
		self . b10 += b11 * self . cprob
		self . b20 += b21 * self . cprob
		self . b30 += b31 * self . cprob
		b = np . sqrt ( self . b10 ** 2 + self . b20 ** 2 + self . b30 ** 2 )
		maxb = np . max ( b )
		self . b10 /= maxb
		self . b20 /= maxb
		self . b30 /= maxb
		return None
		
	def wakefieldinit ( self ) :
		self . wpos = 0
		nx , ny , nz = self . sizepar
		self . timestamp = np . zeros ( self . sizepar ) + - 1
		for k in xrange ( 1 ) :
			for i in xrange ( nx ) :
				for j in xrange ( ny ) :
					self . b11 [ i ] [ j ] [ k ] = self . wstrength * np . exp ( - k / self . wwidth ) + np . random . uniform ( - 0.2 , 0.2 )
					self . b21 [ i ] [ j ] [ k ] = self . wstrength * np . exp ( - k / self . wwidth ) + np . random . uniform ( - 0.2 , 0.2 )
					self . timestamp [ i ] [ j ] [ 0 ] = self . mastertime
		return None
	
	def wakefieldupd ( self ) :
		nx , ny , nz = self . sizepar
		for k in xrange ( self . wpos + 1 ) :
			for i in xrange ( nx ) :
				for j in xrange ( ny ) :
					if k < nz :
						self . b11 [ i ] [ j ] [ k ] = self . wstrength * np . exp ( ( k - self . wpos ) / self . wwidth ) + np . random . uniform ( - 0.2 , 0.2 )
						self . b21 [ i ] [ j ] [ k ] = self . wstrength * np . exp ( ( k - self . wpos ) / self . wwidth ) + np . random . uniform ( - 0.2 , 0.2 )
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				if self . wpos < nz :
					self . timestamp [ i ] [ j ] [ self . wpos ] = self . mastertime			
		return None
		
	def fluxcalc ( self ) :
		nx , ny , nz = self . sizepar
		n0_min , n0_del , n0_tau , g_min , g_max , g_cut , g_tau , freq = self . energypar
		bx = self . b10 + self . b11
		by = self . b20 + self . b21
		bz = self . b30 + self . b31
		
		b = np . sqrt ( bx * bx + by * by + bz * bz )
		
		site_flux = 0
		for k in xrange ( nz ) :
			if self . timestamp [ 0 ] [ 0 ] [ k ] >= 0 :
				t = self . mastertime - self . timestamp [ 0 ] [ 0 ] [ k ]
				n0 = n0_min + n0_del * np . exp ( - t / n0_tau )
				g_up = g_cut + ( g_max - g_cut ) * np . exp ( -t / g_tau )
				for i in xrange ( nx ) :
					for j in xrange ( ny ) :
						site_flux += efuncint ( freq , n0 , g_min , g_up , t , b [ i ] [ j ] [ k ] , s = 2.5 ) 
		
		return site_flux
		
	def wake_time ( self ) :
		return self . mastertime
	
	def wake_pos ( self ) :
		return self . wpos
		
	def wsys_update ( self ) :
		self . mastertime += self . dt
		self . wpos += 1
		self . wakefieldupd ( )
		#self . randfieldupd ( )
		return None
						
	
		
#-----------------------------------------------------------------------------------------

#-------------- WAKE ANIMATOR ------------------------------------------------------------
class WakeField :
	"""
	Field with wake class, initial state is the unperturbed magnetic field.
	"""
	def __init__ ( self , bx , by , bz , nx , ny , nz , lx , ly , lz , wwidth , wstrength ) :
		self . b10 = bx
		self . b20 = by
		self . b30 = bz
		self . sizepar = ( nx , ny , nz )
		self . realpar = ( lx , ly , lz )
		self . wpos = 0
		self . wwidth = wwidth
		self . b11 = bx
		self . b21 = by
		self . b31 = bz
		self . wstrength = wstrength
		
	def wake_ ( self ) :
		nx , ny , nz = self . sizepar
		temp = np . zeros ( self . sizepar )
		tot = float ( nx * ny )
		for x in xrange ( nx ) :
			for y in xrange ( ny ) :
				for z in xrange ( nz ) :
					if z > self . wpos :
						temp [ x ] [ y ] [ z ] += 0
					else :
						temp [ x ] [ y ] [ z ] += self . wstrength * np . exp ( ( - self . wpos + z ) / self . wwidth ) / tot

		return temp
		
	def lincalc_ ( self ) :
		_ , _ , nz = self . sizepar
		_ , _ , lz = self . realpar
		a = np . zeros ( nz )
		x = np . linspace ( - lz / 2 , lz / 2 , nz )
		bxx = np . zeros ( ( 1 , 1 , nz ) )
		byy = np . zeros ( ( 1 , 1 , nz ) )
		bzz = np . zeros ( ( 1 , 1 , nz ) )
		bx = self . b11
		by = self . b21
		bz = self . b31
		tot = float ( nx * ny )
		for i in xrange ( nx ) :
			for j in xrange ( ny ) :
				for k in xrange ( nz ) :
					bxx [ 0 ] [ 0 ] [ k ] += bx [ i ] [ j ] [ k ] / tot
					byy [ 0 ] [ 0 ] [ k ] += by [ i ] [ j ] [ k ] / tot
					bzz [ 0 ] [ 0 ] [ k ] += bz [ i ] [ j ] [ k ] / tot
					#a [ k ] += np . sqrt ( self . b11 [ i ] [ j ] [ k ] * self . b11 [ i ] [ j ] [ k ] + self . b21 [ i ] [ j ] [ k ] * self . b21 [ i ] [ j ] [ k ] )
		for i in xrange ( nz ) :
			a [ i ] = np . sqrt ( bxx * bxx + byy * byy ) [ 0 ] [ 0 ] [ i ]
		return ( x , a )
		
	def update ( self ) :
		self . wpos += 1
		self . b11 = self . b10 + self . wake_ ( )
		self . b21 = self . b20 + self . wake_ ( ) 
		self . b31 = self . b30
		return ( self . lincalc_ ( ) )
		
def init_anim ( ) :
	global Wake
	x , y = Wake . update ( )
	line . set_data ( x , y )
	line . axes . set_xlim ( np . min ( x ) , np . max ( x ) )
	line . axes . set_ylim ( np . min ( y ) , np . max ( y ) )
	return line ,
	
def animate ( i ) :
	global Wake
	x , y = Wake . update ( )
	line . set_data ( x , y )
	line . axes . set_xlim ( np . min ( x ) , np . max ( x ) )
	line . axes . set_ylim ( np . min ( y ) , np . max ( y ) )
	return line ,
#-----------------------------------------------------------------------------------------

#-------------- YES/NO STATEMENT VALIDATOR -----------------------------------------------
def yesno ( statement ) :
	check = raw_input ( statement )
	while ( check != 'n' and check != 'y' and check != 'N' and check != 'Y' ):
		check = raw_input ( FAIL+"Enter again: "+ENDC )
	if check == 'n' or check == 'N':
		return False
	elif check == 'y' or check == 'Y':
		return True
#-----------------------------------------------------------------------------------------

#-------------- END POST PROCESSING FUNCTIONS --------------------------------------------

#-------------- MAIN PROGRAM -------------------------------------------------------------			
if __name__ == "__main__" : #check if the code is invoked from terminal and then execute
	#ascerting Mayavi and Matplotlib imports
	mayavipresence = True
	matplotlibpresence = True
	try : 
		import sys
	except ImportError :
		print FAIL+BOLD+"Basic package 'sys' is missing. Terminating."+ENDC+ENDC
		raise ImportError ( "No module named sys" )
	try :
		import mayavi . mlab as p3
	except ImportError :
		print WARN+"Important 3D plotting system"+ENDC,
		print OKG+BOLD+"Mayavi"+ENDC+ENDC,
		print WARN+"is missing. You may continue executing the program, but you will not be able to"+ENDC,
		print FAIL+"see"+ENDC,
		print WARN+"the plots of the magnetic fields. You may only calculate the magnetic field, save it, reduce it, calculate the divergence and see the synchrotron radiation pattern."+ENDC
		mayavipresence = False
	try :
		import matplotlib . pyplot as plt
		import matplotlib . animation as matanim
	except ImportError :
		if mayavipresence == False :
			print OKG+BOLD+"Mayavi"+ENDC+ENDC,
			print FAIL+"and"+ENDC,
			print OKG+BOLD+"matplotlib"+ENDC+ENDC,
			print FAIL+"both are missing. You shall not be able to generate"+ENDC,
			print OKB+BOLD+ULINE+"ANY plots."+ENDC+ENDC+ENDC,
			print FAIL+"You can only create the magnetic field data and save it. Reduction will also be disabled as it is of no use without these packages."+ENDC
			print
			print OKG+"If you are using Anaconda, use the following command to install Mayavi:"+ENDC
			print BOLD+ULINE+"conda install -c https://conda.anaconda.org/anaconda mayavi\n"+ENDC+ENDC
			print OKB+"If you are using Enthought Canopy, just search for Mayavi in the Canopy Package Manager."+ENDC
			print WARN+"Warning: "+ENDC+"Mayavi for Anaconda installation seems to break matplotlib, as matplotlib is usually installed with Numpy 1.10, whereas installation of Mayavi requires a downgrade to 1.09. Proceed with care."
		else :
			print OKG+BOLD+"matplotlib"+ENDC+ENDC,
			print WARN+"is missing."+ENDC,
			print "You shall not obtain anything but the magnetic field plots. But as you may access the reduced magnetic field and plot it, reduction will be enabled."
		matplotlibpresence = False
		
	if mayavipresence == False and matplotlibpresence == False :
		reductionpos = False
	else :
		reductionpos = True
	#various usage scenarios
	# 1. Usage / Help message.
	if ( len ( sys . argv ) != 11 ) and ( len ( sys . argv ) != 9 ) : #usage message
		#Help statement
		print "\n"
		print OKG+BOLD+"Mayavi"+ENDC+ENDC,
		print "is required to plot the magnetic fields.\n"
		print OKG+"To generate new magnetic fields:"+ENDC
		print OKB+"Usage:\n<script> <pixel x> <pixel y> <pixel z> <length x> <length y> <length z> <reduction> <site width> <variance> <output file name>"+ENDC
		print "(No need to specify file extension, or the parameters. They will automatically be appended to the file name. The extension will be .dat as is the practice.)\n"
		print BOLD+"OR:"+ENDC
		print OKB+"To read from old magnetic field data:"+ENDC
		print OKG+"Usage:\n<script> <pixel x> <pixel y> <pixel z> <length x> <length y> <length z> <reduction> <input file name>\n\n"+ENDC
		sys . exit ( 0 )
	# 2. Production mode.
	elif len ( sys . argv ) == 11 : #if passed
		print BOLD+"Entering production mode."+ENDC
		nx = int ( sys . argv [ 1 ] )
		ny = int ( sys . argv [ 2 ] )
		nz = int ( sys . argv [ 3 ] )
		
		lx = float ( sys.argv [ 4 ] )
		ly = float ( sys.argv [ 5 ] )
		lz = float ( sys.argv [ 6 ] )
		
		reduction = int ( sys . argv [ 7 ] )
		
		if ( nx % reduction != 0 or ny % reduction != 0 or nz % reduction != 0 ) :
			print BOLD+FAIL+"Fatal: Use xpixel, ypixel, zpixel as multiple of reduction."+ENDC+ENDC
			sys . exit ( -1 )
		
		pown = float ( sys.argv [ 8 ] )
		var = float ( sys.argv [ 9 ] )
		
		fname = sys . argv [ 10 ] + "_" + str ( nx ) + "_" + str ( ny ) + "_" + str ( nz ) + "_" + str ( lx ) + "_" + str ( ly ) + "_" + str ( lz ) + "_" + str ( pown ) + "_" + str ( var ) + ".dat"
		
		resx = float ( lx ) / float ( nx - 1 )
		resy = float ( ly ) / float ( ny - 1 )
		resz = float ( lz  ) / float ( nz - 1 )
	
		lxaxis = np . linspace ( - lx / 2 , lx / 2 , nx )
		lyaxis = np . linspace ( - ly / 2 , ly / 2 , ny )
		lzaxis = np . linspace ( - lz / 2 , lz / 2 , nz )
		print OKB+"Generating the field..."+ENDC
		( bx , by , bz ) = magfield ( nx , ny , nz , lx , ly , lz , pown , var )
		print BOLD+"Field generated. "+ENDC,
		if reductionpos == True :
			bxx , byy , bzz = reduction_ ( bx , by , bz , reduction )
		
		#output to file
		if yesno ( BOLD+"Do you want to save the generated magnetic field? ( Y / N ): "+ENDC ) :
			print "Writing to file..."
			ofile = open ( fname , 'w' )
			for i in xrange ( nx ) :
				for j in xrange ( ny ) :
					for k in xrange ( nz ) :
						str1 = str ( i ) + " " + str ( j ) + " " + str ( k ) + " " + str ( bx [ i ] [ j ] [ k ] ) + " " + str ( by [ i ] [ j ] [ k ] ) + " " + str ( bz [ i ] [ j ] [ k ] ) + str ( lxaxis [ i ] ) + " " + str ( lyaxis [ j ] ) + " " + str ( lzaxis [ k ] ) + "\n"
						ofile . write ( str1 )
			ofile . close ( )
			print "Done."
			
		#plotting magnetic field
		if mayavipresence == True :
			if yesno ( OKB+"Do you want to plot the generated magnetic field? ( Y / N ): "+ENDC ) :
				p3 . quiver3d ( bx , by , bz )
				p3 . show ( )
		
		if matplotlibpresence == True :
			if yesno ( WARN+"Do you want calculate the divergence? ( Y / N ): "+ENDC ) :
				divcalc ( bx , by , bz , nx , ny , nz , lx , ly , lz )
		
		if mayavipresence == True :
			if yesno ( OKB+"Do you want to plot the reduced magnetic field? ( Y / N ): "+ENDC ) :
				p3 . quiver3d ( bxx , byy , bzz )
				p3 . show ( )
		if matplotlibpresence == True :	
			image_calc ( bx , by , bz , BOLD+"Total Field"+ENDC )
		
			image_calc ( bxx , byy , bzz , BOLD+"Reduced Field"+ENDC )
			
			if yesno ( BOLD + OKB + "Do you want to plot the linear axis-wise variation? ( Y / N ): " + ENDC + ENDC ) :
				lincalc ( bx , by , bz , nx , ny , nz , lx , ly , lz )
		"""# The famous wake
		fig = plt . figure ( )
		ax = fig . add_subplot ( 111 , autoscale_on = False )
		ax . grid ( )
		line , = ax . plot ( [ ] , [ ] )
		Wake = WakeField ( bx , by , bz , nx , ny , nz , lx , ly , lz , 4 , 0.01 )
		ani = matanim . FuncAnimation ( fig , animate , frames = nz - 2 , interval = 50 , blit = False , init_func = init_anim )
		plt . show ( )
		
		Wake = WakeField ( bx , by , bz , nx , ny , nz , lx , ly , lz , 4 , 0.5 )
		plt . figure ( )
		for i in xrange ( nz - 1 ) :
			x , y = Wake . update ( )
			plt . plot ( x , y )
			plt . show ( )
		"""
		fieldpar = ( nx , ny , nz , lx , ly , lz , pown , var )
		energypar = ( (1.05**(-33)) , 6. ** ( -33 ) , 0.2 , 50 , 3 * 10 ** 5 , 70 , 0.1 , 5 * 10 ** 14 )
		Wake = FWake ( bx , by , bz , fieldpar , energypar , 3 )
		Wake . wakefieldinit ( )
		time = [ ]
		flux = [ ]
		time . append ( Wake . mastertime )
		flux . append ( Wake . fluxcalc ( ) )
		while Wake . wpos < 2 * nz - 1 :
			print Wake . wpos
			Wake . wsys_update ( )
			time . append ( Wake . mastertime )
			flux . append ( Wake . fluxcalc ( ) )
		#Wake . wakefieldinit ( )
		#while Wake . wpos < 2 * nz - 1 :
		#	Wake . wsys_update ( )
		#	time . append ( Wake . mastertime )
		#	flux . append ( Wake . fluxcalc ( ) )
		plt . figure ( )
		plt . plot ( np . asarray ( time ) , ( np . asarray ( flux ) ) )
		plt . show ( )
		
		sys . exit ( 0 )
	# 3. Recovery mode.	
	elif len ( sys . argv ) == 9 :
		print BOLD+"Entering recovery mode."+ENDC
		
		nx = int ( sys . argv [ 1 ] )
		ny = int ( sys . argv [ 2 ] )
		nz = int ( sys . argv [ 3 ] )
		
		lx = float ( sys.argv [ 4 ] )
		ly = float ( sys.argv [ 5 ] )
		lz = float ( sys.argv [ 6 ] ) 
		
		resx = float ( lx ) / float ( nx - 1 )
		resy = float ( ly ) / float ( ny - 1 )
		resz = float ( lz  ) / float ( nz - 1 )
		
		reduction = int ( sys . argv [ 7 ] )
		
		if ( nx % reduction != 0 or ny % reduction != 0 or nz % reduction != 0 ) :
			print BOLD+FAIL+"Fatal: Use xpixel, ypixel, zpixel as multiple of reduction."+ENDC+ENDC
			sys . exit ( -1 )
					
		fname = sys . argv [ 8 ]
		
		ifile = open ( fname , 'r' )
		
		bx = np . zeros ( ( nx , ny , nz ) )
		by = np . zeros ( ( nx , ny , nz ) )
		bz = np . zeros ( ( nx , ny , nz ) )
		print OKB+"Reading from file..."+ENDC
		for line in ifile:
			word = line . rstrip ( '\n' )
			temp = word . split ( )
			i = int ( temp [ 0 ] )
			j = int ( temp [ 1 ] )
			k = int ( temp [ 2 ] )
			bx [ i ] [ j ] [ k ] = float ( temp [ 3 ] )
			by [ i ] [ j ] [ k ] = float ( temp [ 4 ] )
			bz [ i ] [ j ] [ k ] = float ( temp [ 5 ] )
		ifile . close ( )
		print BOLD+"File read."+ENDC
		#reduction
		if reductionpos == True :
			bxx , byy , bzz = reduction_ ( bx , by , bz , reduction )
			
		#plotting magnetic field
		if mayavipresence == True :
			if yesno ( OKB+"Do you want to plot the retrieved magnetic field? ( Y / N ): "+ENDC ) :
				p3 . quiver3d ( bx , by , bz )
				p3 . show ( )
		
		if matplotlibpresence == True :
			if yesno ( WARN+"Do you want calculate the divergence? ( Y / N ): "+ENDC ) :
				divcalc ( bx , by , bz , nx , ny , nz , lx , ly , lz )
		
		if mayavipresence == True :
			if yesno ( OKB+"Do you want to plot the reduced magnetic field? ( Y / N ): "+ENDC ) :
				p3 . quiver3d ( bxx , byy , bzz )
				p3 . show ( )
		if matplotlibpresence == True :	
			image_calc ( bx , by , bz , BOLD+"Total Field"+ENDC )
		
			image_calc ( bxx , byy , bzz , BOLD+"Reduced Field"+ENDC )
			
			if yesno ( BOLD + OKB + "Do you want to plot the linear axis-wise variation? ( Y / N ): " + ENDC + ENDC ) :
				lincalc ( bx , by , bz , nx , ny , nz , lx , ly , lz )
		
		sys . exit ( 0 )
#-----------------------------------------------------------------------------------------

#-------------- CHANGELOG ----------------------------------------------------------------
"""
Errors:
1. File reading has some problems. There is a grid-distortion of the magnetic field. (Fixed)
Changelog:
1. Changed the breaks, so that all of the field is generated. Previously only two diagonal
octants were being generated.
2. Fixed the cross product coordinate addressing.
3. Fixed the cross product itself.

Change 1 eliminated the regularity of the averaging on one axis, and changed it into an
irregular thing.
Change 2+3 did nothing more than reducing the irregular change of Y axis by 10^-6

4. Using abs instead of real/imaginary seems to take out the issues, but all components
are now strictly positive. Not allowed.
5. Removing ANY breaks and conjugate symmetry. Solved all the issues. Code working
perfectly.

5a. Further modularization of the code.

5b. Error 1 fixed.

5c. Color incompatibility in Windows is circumvented, colors in Windows removed.

6. In the linearity module, instead of acting on a net magnetic field, B * sin ( theta ) 
at all points are now being calculated and added.

7. Wake is realistically implemented, with flux information being obtained using
numerical model.

From __future__ import :
1. Plot disturbance caused by some wake. (With partial support, animation implementation
in progress)

2. Incorporate the new numerical wake system, and also a radiation from all other parts of
the box, incl. light travel time if possible. (Multi-state n all, got to be complicated)
DONE.

minor 3. Implement independent integration algorithm to circumvent the use of scipy quad.
"""
#-----------------------------------------------------------------------------------------
