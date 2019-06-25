import ctypes
from ctypes import *
import ctypes.util
from numpy.ctypeslib import ndpointer
import numpy as np
#Find and load the library

_lib = ctypes.CDLL('libsolenoid.dylib')
c_magfieldgen = _lib.magfield
"""
class param ( ctypes.Structure ) :
	_fields_ = [ ( "nx" , ctypes.c_int ) ,
					( "ny" , ctypes.c_int ) ,
					( "nz" , ctypes.c_int ) ,
					( "lx" , ctypes.c_double ) ,
					( "ly" , ctypes.c_double ) ,
					( "lz" , ctypes.c_double ) ,
					( "ran" , ctypes.c_double ) ,
					( "var" , ctypes.c_double ) ]
"""
def magfield ( nx , ny , nz , lx , ly , lz ):
	bx = np.zeros ( nx * ny * nz , dtype = np.double )
	by = np.zeros ( nx * ny * nz , dtype = np.double )
	bz = np.zeros ( nx * ny * nz , dtype = np.double )
	ax = np.zeros ( nx * ny * nz , dtype = np.double )
	ay = np.zeros ( nx * ny * nz , dtype = np.double )
	az = np.zeros ( nx * ny * nz , dtype = np.double )
	#params = param ( nx , ny , nz , lx , ly , lz , ran , var ) print params
	#c_magfieldgen ( params , ctypes.c_int ( nthreads ) ,
	#bx.ctypes.data_as(ctypes.POINTER(c_double)) ,
	#by.ctypes.data_as(ctypes.POINTER(c_double)) ,
	#bz.ctypes.data_as(ctypes.POINTER(c_double)) ) c_magfieldgen ( params , ctypes.c_int
	#( nthreads ) , c_void_p(bx.ctypes.data) , c_void_p(by.ctypes.data) ,
	#c_void_p(bz.ctypes.data) )
	status = c_magfieldgen ( c_int ( nx ) , c_int ( ny ) , c_int ( nz ) , c_double ( lx ) ,
c_double ( ly ) , c_double ( lz ) , c_void_p ( bx . ctypes . data ) , c_void_p ( by . ctypes . data ) , c_void_p ( bz . ctypes . data ) , c_void_p ( ax . ctypes . data ) , c_void_p ( ay . ctypes . data ) , c_void_p ( az . ctypes . data ) )
	if ( status != 0 ) :
		print "Failed to generate the magnetic field."
	return np . reshape ( bx , ( nx , ny , nz ) ) , np . reshape ( by , ( nx , ny , nz )) , np . reshape ( bz , ( nx , ny , nz ) ) , np . reshape ( ax , ( nx , ny , nz ) ) , np . reshape ( ay , ( nx , ny , nz )) , np . reshape ( az , ( nx , ny , nz ) )

if __name__ == '__main__' :
    mayavipresence = True
    try :
        import mayavi . mlab as p3
    except ImportError :
    	print "No Mayavi"
        mayavipresence = False
    bx , by , bz , ax , ay , az = magfield ( 2**8 , 2**8 , 2**3 , 32 , 32 , 1 )
    print np.max(ax),np.max(ay),np.max(az)
    print np.max(bx),np.max(by),np.max(bz)
    if mayavipresence :
    	p3.quiver3d(ax,ay,az)
    	p3.show()
        p3 . quiver3d ( bx , by , bz )
        p3 . show ( )
