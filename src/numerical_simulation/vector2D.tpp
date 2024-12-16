/* This header file contains the definitions of 2D vectors */
#include "vector3D.h"

/********************************************************************************/

template < typename T >
inline void posizione_2D<T>::operator=( const posizione_2D<T> &other ) {
  x = other.x ;
  y = other.y ;
};
/********************************************************************************/

template < typename T >
inline bool posizione_2D<T>::operator==( const posizione_2D<T> &other ) {
  return x == other.x && y == other.y ;
};
/********************************************************************************/

template < typename T >
inline bool posizione_2D<T>::operator!=( const posizione_2D<T> &other ) {
  return x != other.x || y != other.y ;
};
/********************************************************************************/

template < typename T >
inline double posizione_2D<T>::operator*( const posizione_2D<T> &other ) {
  return x * other.x + y * other.y ;
};
/********************************************************************************/

template < typename T >
inline posizione_2D<T> posizione_2D<T>::operator*( double factor ) {
  posizione_2D<T> product ;
  product.x = x * factor ;
  product.y = y * factor ;
  return product ;
};
/********************************************************************************/

template < typename T >
inline posizione_2D<T> posizione_2D<T>::operator/( double factor ) {
  posizione_2D<T> product ;
  product.x = x / factor ;
  product.y = y / factor ;
  return product ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_2D<T>::operator/( posizione_3D<int> int_vec ) {
  posizione_3D<T> division ;
  division.x = ( x / (double)int_vec.x ) ;
  division.y = ( y / (double)int_vec.y ) ;
  division.z = 0.0 ;
  return division ;
}
/********************************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_2D<T>::operator/( posizione_3D<double> vec ) {
  posizione_3D<T> division ;
  division.x = ( x / vec.x ) ;
  division.y = ( y / vec.y ) ;
  division.z = 0.0 ;
  return division ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_2D<T>::operator/=( double factor ) {
  x /= factor ;
  y /= factor ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_2D<T>::operator*=( double factor ) {
  x = x * factor ;
  y = y * factor ;
};
/********************************************************************************/

template < typename T >
inline posizione_2D<T> posizione_2D<T>::operator-( const posizione_2D<T> &other ) {
  posizione_2D<T> difference ;
  difference.x = x - other.x ;
  difference.y = y - other.y ;
  return difference ;
};
/********************************************************************************/

template < typename T >
inline void posizione_2D<T>::operator-=( const posizione_2D<T> &other ) {
  x = x - other.x ;
  y = y - other.y ;
};
/********************************************************************************/

template < typename T >
inline void posizione_2D<T>::operator-=( double disp ) {
  x = x - disp ;
  y = y - disp ;
};
/********************************************************************************/

template < typename T >
inline posizione_2D<T> posizione_2D<T>::operator+( const posizione_2D<T> &other ) {
  posizione_2D<T> sum ;
  sum.x = x + other.x ;
  sum.y = y + other.y ;
  return sum ;
};
/********************************************************************************/

template < typename T >
inline void posizione_2D<T>::operator+=( const posizione_2D<T> &other ) {
  x = x + other.x ;
  y = y + other.y ;
};
/********************************************************************************/

template < typename T >
inline void posizione_2D<T>::operator+=( double delta ) {
  x = x + delta ;
  y = y + delta ;
};
/********************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_2D<T>::vectorial_product(const posizione_2D<T> &other) {
  posizione_3D<T> vector_product;
  vector_product.x = 0 ;
  vector_product.y = 0 ;
  vector_product.z = x * other.y - y * other.x ;

  return vector_product;
}
/********************************************************************************************/

template < typename T >
inline posizione_3D<T> posizione_2D<T>::vectorial_product(const posizione_3D<T> &other) {
  posizione_3D<T> vector_product;
  vector_product.x = y * other.z ;
  vector_product.y = - x * other.z ;
  vector_product.z = x * other.y - y * other.x ;

  return vector_product;
}
/********************************************************************************************/

template < typename T >
inline double posizione_2D<T>::norm( void ) {
  return sqrt( x*x + y*y ) ;
};
/********************************************************************************/

template < typename T >
inline double posizione_2D<T>::square_norm( void ) {
  return x*x + y*y ;
};
/********************************************************************************/

template < typename T >
inline double posizione_2D<T>::distance(const posizione_2D<T> *other) {
  return sqrt( (x - other->x)*(x - other->x) + (y - other->y)*(y - other->y) ) ;
}
/********************************************************************************************/

template < typename T >
inline double posizione_2D<T>::distance(const posizione_2D<T> &other) {
  return sqrt( (x - other.x)*(x - other.x) + (y - other.y)*(y - other.y) ) ;
}
/********************************************************************************************/

template < typename T >
inline double posizione_2D<T>::theta( void ) {
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
inline struct posizione_3D<int> posizione_2D<T>::convert2int_floor( void ) {
  struct posizione_3D<int> ret_val ;
  ret_val.x = (int)floor(x) ;
  ret_val.y = (int)floor(y) ;
  ret_val.z = 0 ;

  return ret_val ;
}
/********************************************************************************************/

template < typename T >
inline struct posizione_3D<int> posizione_2D<T>::convert2int_ceil( void ) {
  struct posizione_3D<int> ret_val ;
  ret_val.x = (int)ceil(x) ;
  ret_val.y = (int)ceil(y) ;
  ret_val.z = 0 ;

  return ret_val ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_2D<T>::rndPos_generate( double max ) {  // generate a random position for the particle in the specified box
  x = drand48() * max ;
  y = drand48() * max ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_2D<T>::rndPos_generate( double x_max, double y_max ) {  // generate a random position for the particle in the specified box
  x = drand48()*x_max;
  y = drand48()*y_max;
}
/********************************************************************************************/

template < typename T >
inline void posizione_2D<T>::rndPos_generate( const struct posizione_2D<double> &max ) {  // generate a random position for the particle in the specified box
  x = drand48() * max.x ;
  y = drand48() * max.y ;
}
/********************************************************************************************/

template < typename T >
inline void posizione_2D<T>::clear( void ) {  x = 0 , y = 0 ;  };
/********************************************************************************/
