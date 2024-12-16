/* This header file contains the definitions of 3D vectors */
#include "vector2D.h"

/********************************************************************************/

template < typename T >
inline void posizione_3D<T>::operator=( const posizione_3D<T> &other ) {
  x = other.x ;
  y = other.y ;
  z = other.z ;
};
/********************************************************************************/

template < typename T >
inline bool posizione_3D<T>::operator==( const posizione_3D<T> &other ) {
  return x == other.x && y == other.y && z == other.z ;
};
/********************************************************************************/

template < typename T >
inline bool posizione_3D<T>::operator!=( const posizione_3D<T> &other ) {
  return x != other.x || y != other.y || z != other.z ;
};
/********************************************************************************/

template < typename T >
inline double posizione_3D<T>::operator*( const posizione_3D<T> &other ) {
  return x * other.x + y * other.y + z * other.z ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::operator*( double factor ) {
  posizione_3D<T> product ;
  product.x = x * factor ;
  product.y = y * factor ;
  product.z = z * factor ;
  return product ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::operator/( double factor ) {
  posizione_3D<T> product ;
  product.x = x / factor ;
  product.y = y / factor ;
  product.z = z / factor ;
  return product ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::operator/( posizione_3D<int> int_vec ) {
  posizione_3D<T> division ;
  division.x = ( x / (double)int_vec.x ) ;
  division.y = ( y / (double)int_vec.y ) ;
  division.z = ( z / (double)int_vec.z ) ;
  return division ;
}
/********************************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::operator/( posizione_3D<double> vec ) {
  posizione_3D<T> division ;
  division.x = ( x / vec.x ) ;
  division.y = ( y / vec.y ) ;
  division.z = ( z / vec.z ) ;
  return division ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_3D<T>::operator/=(double factor) {
  x /= factor ;
  y /= factor ;
  z /= factor ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_3D<T>::operator*=( double factor ) {
  x = x * factor ;
  y = y * factor ;
  z = z * factor ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::operator-( const posizione_3D<T> &other ) {
  posizione_3D<T> difference ;
  difference.x = x - other.x ;
  difference.y = y - other.y ;
  difference.z = z - other.z ;
  return difference ;
};
/********************************************************************************/

template < typename T >
inline void posizione_3D<T>::operator-=( const posizione_3D<T> &other ) {
  x = x - other.x ;
  y = y - other.y ;
  z = z - other.z ;
};
/********************************************************************************/

template < typename T >
inline void posizione_3D<T>::operator-=( double disp ) {
  x = x - disp ;
  y = y - disp ;
  z = z - disp ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::operator+( const posizione_3D<T> &other ) {
  posizione_3D<T> sum ;
  sum.x = x + other.x ;
  sum.y = y + other.y ;
  sum.z = z + other.z ;
  return sum ;
};
/********************************************************************************/

template < typename T >
inline void posizione_3D<T>::operator+=( const posizione_3D<T> &other ) {
  x = x + other.x ;
  y = y + other.y ;
  z = z + other.z ;
};
/********************************************************************************/

template < typename T >
inline void posizione_3D<T>::operator+=( double delta ) {
  x = x + delta ;
  y = y + delta ;
  z = z + delta ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::vectorial_product(const posizione_3D<T> &other) {
  posizione_3D<T> vector_product ;
  vector_product.x = y * other.z - z * other.y ;
  vector_product.y = z * other.x - x * other.z ;
  vector_product.z = x * other.y - y * other.x ;

  return vector_product;
}
/********************************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_3D<T>::vectorial_product(const posizione_2D<T> &other) {
  posizione_3D<T> vector_product ;
  vector_product.x = - z * other.y ;
  vector_product.y = z * other.x ;
  vector_product.z = x * other.y - y * other.x ;

  return vector_product;
}
/********************************************************************************************/

template < typename T >
inline double posizione_3D<T>::square_norm( void ) {
  return x*x + y*y + z*z ;
};
/********************************************************************************/

template < typename T >
inline double posizione_3D<T>::norm( void ) {
  return sqrt( x*x + y*y + z*z ) ;
};
/********************************************************************************/

template < typename T >
inline double posizione_3D<T>::distance( const posizione_3D<T> *other ) {
  return sqrt( (x - other->x)*(x - other->x) + (y - other->y)*(y - other->y) + (z - other->z)*(z - other->z) );
}
/********************************************************************************************/

template < typename T >
inline double posizione_3D<T>::distance( const posizione_3D<T> &other ) {
  return sqrt( (x - other.x)*(x - other.x) + (y - other.y)*(y - other.y) + (z - other.z)*(z - other.z) );
}
/********************************************************************************************/

template < typename T >
inline double posizione_3D<T>::zenit( void ) {
  return acos( z / norm() ) ;
};
/********************************************************************************/

template < typename T >
inline double posizione_3D<T>::azimut( void ) {
  //  -0.5*pi < theta <= 1.5*pi
  if ( x > 0 ) return atan( y/x ) ;
  else if ( x < 0 ) return atan( y/x ) + M_PI ;
  else {
    if ( y > 0 ) return 0.5 * M_PI ;
    else if ( y < 0 ) return 1.5 * M_PI ;
    else return nan("") ;
  }
};
/********************************************************************************/

template < typename T >
inline struct posizione_3D<int> posizione_3D<T>::convert2int_floor( void ) {
  struct posizione_3D<int> ret_val ;
  ret_val.x = (int)floor(x) ;
  ret_val.y = (int)floor(y) ;
  ret_val.z = (int)floor(z) ;

  return ret_val ;
}
/********************************************************************************************/

template < typename T >
inline struct posizione_3D<int> posizione_3D<T>::convert2int_ceil( void ) {
  struct posizione_3D<int> ret_val ;
  ret_val.x = (int)ceil(x) ;
  ret_val.y = (int)ceil(y) ;
  ret_val.z = (int)ceil(z) ;

  return ret_val ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_3D<T>::rndPos_generate( double max ) {  // generate a random position for the particle in the specified box
  x = drand48() * max ;
  y = drand48() * max ;
  z = drand48() * max ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_3D<T>::rndPos_generate( double x_max, double y_max, double z_max ) {  // generate a random position for the particle in the specified box
  x = drand48()*x_max;
  y = drand48()*y_max;
  z = drand48()*z_max;
}
/********************************************************************************************/

template < typename T >
inline void posizione_3D<T>::rndPos_generate( const struct posizione_3D<double> &max ) {  // generate a random position for the particle in the specified box
  x = drand48() * max.x ;
  y = drand48() * max.y ;
  z = drand48() * max.z ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_3D<T>::clear( void ) {  x = 0 , y = 0 , z = 0 ;  };
/********************************************************************************/
