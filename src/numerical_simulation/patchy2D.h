#ifndef PATCHY2D_H
#define PATCHY2D_H

class particle_2D ;

class patchy_2D : public particle_2D {
public:
  double rotation ;
  int span_round_angles = 0 ;
  int Npatches = 0 ;
  struct posizione_2D<double> *patches = NULL ;
  patchy_2D **bonds = NULL ;
  vlist_element<patchy_2D> **list_ptr ;

  void set_patches( int N_patches ) ;
  inline void operator=( const patchy_2D &part ) ;
  inline void operator=( const particle_2D &part ) ;
  inline void rndPos_generate( double max ) ;
  inline void rndPos_generate( const struct posizione_2D<double> &max ) ;
  inline void rand_rotate( double maxangle ) ;
  inline struct posizione_2D<double> patch( int i ) ; // this gives the actual orientation of patch i
  inline string dump_position_variables( int Ndigits = 12 , const char *format = NULL , double *parameters = NULL ) ;

  patchy_2D( double x=0 , double y=0 , int val=0 , int id_number=0 , int N_patches = 0 ) ;
  patchy_2D( int val , int N_patches = 0 ) ;
  patchy_2D( const patchy_2D &copy ) ;
  ~patchy_2D() ;
};


#include "patchy2D.ipp"

#endif
