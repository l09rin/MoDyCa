#include "matrix.h"

/*============================================**
                 CLASS COL_VECTOR
**============================================*/

/**************************************
     CARTESIAN_BASE
  set the vector as unitary in the direction e
**************************************/
void col_vector::cartesian_base( int e ) {
  for( int i=0 ; i<D ; i++ ) {
    if( i == e ) element[i] = 1 ;
    else element[i] = 0 ;
  }
};


/**************************************
     OPERATOR=
  assignment as a copy
**************************************/
void col_vector::operator=( const col_vector &other ) {
  if( D != other.D ) {
    delete [] element ;
    D = other.D ;
    element = new double[D] ;
  } else {
    for( int i=0 ; i<D ; i++ ) element[i] = other.element[i] ;
  }
};


/**************************************
     OPERATOR<<
  print on output stream
**************************************/
ostream & operator<< ( ostream &out, const col_vector &vec ) {
  out << "(\t" ;
  for( int i=0 ; i<vec.D ; i++ ) out << vec.element[i] << "\t" ;
  out << ")" ;
  return out ;
};



/*============================================**
                 CLASS MATRIX
**============================================*/

/**************************************
     CONSTRUCTOR : DEFAULT
  unallocated matrix
**************************************/
matrix::matrix( void ) {
  rows = 0 ;
  columns = 0 ;
  element = NULL ;
}

/**************************************
     CONSTRUCTOR 2
  allocates matrix with rowsN rows and columnsN columns
**************************************/
matrix::matrix( int rowsN , int columnsN ) {
  rows = rowsN ;
  columns = columnsN ;
  if( rows != 0 ) {
    element = new double* [rows] ;
    for( int i=0 ; i<rows ; i++ ) element[i] = new double [columns] ;
  }else{
    element = NULL ;
  }
}

/**************************************
     CONSTRUCTOR of COPY
**************************************/
matrix::matrix( const matrix &other_mat ) {
  rows = other_mat.rows ;
  columns = other_mat.columns ;
  if( rows != 0 ) {
    element = new double* [rows] ;
    for( int i=0; i<rows; i++ ) {
      element[i] = new double [columns] ;
      for( int j=0; j<columns; j++ ) element[i][j] = other_mat.element[i][j] ;
    }
  }else{
    element = NULL ;
  }
}

/**************************************
     DESTRUCTOR : DEFAULT
**************************************/
matrix::~matrix() {
  if( rows != 0 ) {
    for( int i=0; i<rows; i++ ) delete[] element[i] ;
  }
  delete[] element ;
}


/**************************************
     REBUILD
  clear all the memory and reallocate a new matrix
**************************************/
void matrix::rebuild( int rowsN, int columnsN ) {
  if( rows != 0 ) {
    for( int i=0; i<rows; i++ ) delete[] element[i] ;
  }
  delete[] element ;
  rows = rowsN ;
  columns = columnsN ;
  if( rows != 0 ) {
    element = new double* [rows] ;
    for( int i=0 ; i<rows ; i++ ) element[i] = new double [columns] ;
  }else{
    element = NULL ;
  }
}


/**************************************
     DETERMINANT
  calculation of the determinant with the
  usual recursive method of minors
**************************************/
double matrix::determinant( void ) {
  double determinante = 0 ;
  if( rows != columns ) return 0 ;
  if( rows == 1 ) {
    return element[0][0] ;
  }else{
    for( int i=0; i<rows; i++ ) {
      matrix *minor = NULL ;
      minor = construct_minor( i , 0 ) ;
      determinante += (1-2*(i%2)) * element[i][0] * minor->determinant() ;
      delete minor ;
    }
    return determinante ;
  }
}


/**************************************
     CONSTRUCT_MINOR
  returns a given minor of the matrix
  (eliminates the row i and column j)
**************************************/
matrix *matrix::construct_minor( int i, int j ) {
  matrix *minor = NULL ;
  minor = new matrix( rows-1 , columns-1 ) ;
  for( int l=0; l<i; l++ ) {
    for( int k=0; k<j; k++ ) {
      minor->element[l][k] = element[l][k] ;
    }
    for( int k=j+1; k<columns; k++ ) {
      minor->element[l][k-1] = element[l][k] ;
    }
  }
  for( int l=i+1; l<rows; l++ ) {
    for( int k=0; k<j; k++ ) {
      minor->element[l-1][k] = element[l][k] ;
    }
    for( int k=j+1; k<columns; k++ ) {
      minor->element[l-1][k-1] = element[l][k] ;
    }
  }

  return minor ;
}


/**************************************
     MINOR_DETERMINANT
  wraps the construct_minor function and returns
  the determinant of a given minor
**************************************/
double matrix::minor_determinant( int i, int j ) {
  matrix *minor = NULL ;
  minor = construct_minor( i, j ) ;
  double det = minor->determinant() ;
  delete minor ;
  return det ;
}


/**************************************
     SOLVE_CRAMER
  solve the linear system M * x = constant_vec
  with the cramer method
**************************************/
col_vector matrix::solve_cramer( const col_vector constant_vec ) {
  double det = determinant() ;
  col_vector solution( constant_vec.D ) ;
  if( columns != rows || rows != constant_vec.D || det == 0 ) {
    cout << "The system is not solvable !" ;
    exit( EXIT_FAILURE ) ;
  }
  matrix appoggio( *this ) ;
  for( int j=0 ; j<constant_vec.D ; j++ ) {
    for( int i=0 ; i<constant_vec.D ; i++ ) {
      appoggio.element[i][j] = constant_vec.element[i] ;
    }
    solution.element[j] = appoggio.determinant() / det ;
    for( int i=0 ; i<constant_vec.D ; i++ ) {
      appoggio.element[i][j] = element[i][j] ;
    }
  }

  return solution ;
}


/**************************************
     COMPUTE_INVERSE
  returns the inverse matrix of the matrix
**************************************/
matrix matrix::compute_inverse( void ) {
  double det = determinant() ;
  if( columns != rows || det == 0 ) {
    cout << "The inverse matrix cannot be computed !" ;
    exit( EXIT_FAILURE ) ;
  }
  matrix inverse( rows , columns ) ;
  col_vector solution( rows ) , cartesian_vec( rows ) ;

  for( int j=0 ; j<columns ; j++ ) {
    cartesian_vec.cartesian_base( j ) ;
    solution = solve_cramer( cartesian_vec ) ;
    for( int i=0 ; i<rows ; i++ ) inverse.element[i][j] = solution.element[i] ;
  }

  return inverse ;
}


/**************************************
     OPERATOR<<
  print the matrix on the out stream
**************************************/
ostream & operator<< ( ostream &out, const matrix &mat ) {
  for( int i=0 ; i<mat.rows ; i++ ) {
    out << "|\t" ;
    for( int j=0 ; j<mat.columns ; j++ ) out << mat.element[i][j] << "\t" ;
    out << "|" << endl ;
  }
  return out ;
};


/**************************************
     DOT_PROD
  returns the matrix resulting by the
  rows by columns product with matrix B
**************************************/
matrix matrix::dot_prod( const matrix &B ) {
  if( columns != B.rows ) {
    cout << "The outer product cannot be computed, the number of rows and columns of the two matrices does not match !" << endl ;
    exit( EXIT_FAILURE ) ;
  } else {
    matrix result( rows , B.columns ) ;
    for( int i=0 ; i<result.rows ; i++ ) {
      for( int j=0 ; j<result.columns ; j++ ) {
	result.element[i][j] = 0 ;
	for( int k=0 ; k<columns ; k++ ) result.element[i][j] += ( element[i][k] * B.element[k][j] ) ;
      }
    }

    return result ;
  }
}
