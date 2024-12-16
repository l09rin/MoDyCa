#ifndef MATRIX_H
#define MATRIX_H
/* This library contains the class matrix, that allows to easily do matrix calculations */
#include <cstddef>
#include <iostream>

using namespace std;

class col_vector {
 public:
  int D = 0 ;
  double *element = NULL ;

  col_vector( void ) {
    element = NULL ;
    D = 0 ;
  };
  col_vector( int N ) {
    element = new double[N] ;
    D = N ;
  };
  col_vector( const col_vector &other ) {
    D = other.D ;
    element = new double[D] ;
    for( int i=0 ; i<D ; i++ ) element[i] = other.element[i] ;
  };
  ~col_vector() {
    delete [] element ;
    element = NULL ;
    D = 0 ;
  };

  void cartesian_base( int e ) ;
  void operator=( const col_vector &other ) ;
  void clear( void ) {    for( int i=0 ; i<D ; i++ ) element[i] = 0 ;  };
  friend ostream & operator<< ( ostream &out, const col_vector &vec ) ;
};



class matrix {
 public:
  int rows , columns ;
  double **element = NULL ;

  matrix( void ) ;
  matrix( int rowsN, int columnsN ) ;
  matrix( const matrix &other_mat ) ;
  ~matrix() ;

  void rebuild( int rowsN, int columnsN ) ;
  double determinant( void ) ;
  double minor_determinant( int i , int j ) ;
  matrix *construct_minor( int i , int j ) ;
  matrix compute_inverse( void ) ;
  col_vector solve_cramer( const col_vector constant_vec ) ;
  friend ostream & operator<< ( ostream &out, const col_vector &vec ) ;
  matrix dot_prod( const matrix &B ) ;
};

#endif
