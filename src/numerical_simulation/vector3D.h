#ifndef VECTOR3D_H
#define VECTOR3D_H
/* This header file contains the definitions of 2D and 3D vectors and of the classes particle_3D and particle_2D, containing all
     the information about single particles and managing all the operations performed over them */

using namespace std;
#include <cmath>
#include <sstream>

template < typename T >
struct posizione_2D ;

template < typename T >
struct posizione_3D {
  T x, y, z ;

  friend ostream & operator<< ( ostream &out, const struct posizione_3D &pos ) {
    out << pos.x << " " << pos.y << " " << pos.z ;
    return out ;
  };
  inline void operator=( const posizione_3D<T> &other ) ;
  inline bool operator==( const posizione_3D<T> &other ) ;
  inline bool operator!=( const posizione_3D<T> &other ) ;
  inline double operator*( const posizione_3D<T> &other ) ;
  inline posizione_3D<T> operator*( double factor ) ;
  inline void operator*=( double factor ) ;
  inline posizione_3D<T> operator/( double factor ) ;
  inline posizione_3D<T> operator/( posizione_3D<int> int_vec ) ;
  inline posizione_3D<T> operator/( posizione_3D<double> vec ) ;
  inline void operator/=( double factor ) ;
  inline posizione_3D<T> operator-( const posizione_3D<T> &other ) ;
  inline void operator-=( const posizione_3D<T> &other ) ;
  inline void operator-=( double disp ) ;
  inline posizione_3D<T> operator+( const posizione_3D<T> &other ) ;
  inline void operator+=( const posizione_3D<T> &other ) ;
  inline void operator+=( double delta ) ;
  inline posizione_3D<T> vectorial_product( const posizione_3D<T> &other ) ;
  inline posizione_3D<T> vectorial_product( const posizione_2D<T> &other ) ;
  inline double square_norm( void ) ;
  inline double norm( void ) ;
  inline double distance( const posizione_3D<T> *other ) ;
  inline double distance( const posizione_3D<T> &other ) ;
  inline double zenit( void ) ;
  inline double azimut( void ) ;
  inline posizione_3D<int> convert2int_floor( void ) ;
  inline posizione_3D<int> convert2int_ceil( void ) ;
  inline void rndPos_generate( double max ) ;
  inline void rndPos_generate( double x_max , double y_max , double z_max ) ;
  inline void rndPos_generate( const struct posizione_3D<double> &max ) ;
  inline void clear( void ) ;
};


#include "vector3D.tpp"

#endif
