import os, os.path, platform
import ctypes
from ctypes import *
import ctypes.util
from numpy.ctypeslib import ndpointer
import numpy as np
#Find and load the library
_lib = None
if platform.system()=='Darwin':
    _libraryname= 'libmagfieldgen.dylib'
else:
    _libraryname= 'libmagfieldgen.so'
_libname = ctypes.util.find_library(_libraryname)
#if _libname:
#    _lib = ctypes.cdll.LoadLibrary(_libname)
if _libname:
    _lib = ctypes.CDLL(_libname)
if _lib is None: #Hack
    p = os.path.join(TEMPLATE_LIBRARY_PATH,_libraryname)
    if os.path.exists(p):
        _lib = ctypes.CDLL(p)
if _lib is None:
    raise IOError(_libraryname+' library not found')


c_magfieldgen = _lib.magfieldgen
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
def magfield ( nx , ny , nz , lx , ly , lz , ran , var , nthreads = 2 ):
	bx = np.zeros ( nx * ny * nz , dtype = np.double )
	by = np.zeros ( nx * ny * nz , dtype = np.double )
	bz = np.zeros ( nx * ny * nz , dtype = np.double )
	#params = param ( nx , ny , nz , lx , ly , lz , ran , var ) print params
	#c_magfieldgen ( params , ctypes.c_int ( nthreads ) ,
	#bx.ctypes.data_as(ctypes.POINTER(c_double)) ,
	#by.ctypes.data_as(ctypes.POINTER(c_double)) ,
	#bz.ctypes.data_as(ctypes.POINTER(c_double)) ) c_magfieldgen ( params , ctypes.c_int
	#( nthreads ) , c_void_p(bx.ctypes.data) , c_void_p(by.ctypes.data) ,
	#c_void_p(bz.ctypes.data) )
	status = c_magfieldgen ( c_int ( nx ) , c_int ( ny ) , c_int ( nz ) , c_double ( lx ) ,
c_double ( ly ) , c_double ( lz ) , c_double ( ran ) , c_double ( var ) , c_int ( nthreads )
, c_void_p ( bx . ctypes . data ) , c_void_p ( by . ctypes . data ) , c_void_p ( bz . ctypes
. data ) )
	if ( status != 0 ) :
		print "Failed to generate the magnetic field."
	return np . reshape ( bx , ( nx , ny , nz ) ) , np . reshape ( by , ( nx , ny , nz )) , np . reshape ( bz , ( nx , ny , nz ) )

if __name__ == '__main__' :
    mayavipresence = True
    try :
        import mayavi . mlab as p3
    except ImportError :
        mayavipresence = False
    bx , by , bz = magfield ( 100 , 100 , 100 ,10 , 10 , 10 , 100 , 1 , 2 )
    print np.max(bx),np.max(by),np.max(bz)
    if mayavipresence :
        p3 . quiver3d ( bx , by , bz )
        p3 . show ( )
