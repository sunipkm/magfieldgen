#ifndef magfieldgen_h
#define magfieldgen_h
int generator ( int num_thread ) ;
double rnd ( int id ) ;
void fftfreq ( double * in , int n , double d ) ;
double rayleigh ( double l , double ran ) ;
double powspec ( double l , double var , int npix , double psum , double ran ) ;
double psumcalc ( plist param , int num_thread ) ;
int magfield ( plist , int num_thread , double * b1 , double * b2 , double * b3 ) ;
int magfieldgen ( int nx , int ny , int nz , double lx , double ly , double lz , double ran , double var , int nthread , double * bx , double * by , double * bz );
#endif