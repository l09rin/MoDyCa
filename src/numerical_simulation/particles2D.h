#ifndef PARTICLES2D_H
#define PARTICLES2D_H

using namespace std;
#include <cmath>
#include <sstream>
#include "vector2D.h"

class particle_3D ;
template <typename particle> struct vlist_element ;

class particle_2D {

public:
  struct posizione_2D<double> position ;
  struct posizione_2D<int> periodic_box ;
  int id , type , mol ;
  double charge ;
  double mass ;
  double radius ;
  int valence ;
  int *bonded_monomer ;
  vlist_element<particle_2D> **list_ptr = NULL ;

  inline void set_position( const posizione_2D<double> &vec ) ;
  inline void operator=( const particle_2D &part ) ;
  inline void operator=( const particle_3D &part ) ;
  //  inline void operator%( double max_size ) ;
  //  inline void operator%( const particle_2D &max_size ) ;
  inline void rndPos_generate( double max ) ;
  inline void rndPos_generate( const struct posizione_2D<double> &max ) ;
  inline void Gaussian_particle( void ) ;
  inline void random_unit_vector( void ) ;
  void clear( void ) ;
  inline string dump_position_variables( int Ndigits = 12 , const char *format = NULL , double *parameters = NULL ) ;
  friend ostream & operator<< ( ostream &out, const particle_2D &part ) ;
  inline void read_by_string( char *line ) ;
  inline void read_by_string( char **words ) ;
  inline void read_by_string_pbc( char *line ) ;
  inline void read_by_string_pbc( char **words ) ;
  inline double max_comp( void ) ;
  inline double min_comp( void ) ;
  inline double avg_comp( void ) ;
  inline void set_equal_comp( double coord ) ;
  inline double volume( void ) ;
  inline particle_2D PBC_NearNeigh2( const particle_2D &part0 , const particle_2D &box_side ) ;
  inline posizione_2D<double> unwrapped( posizione_2D<double> box_sides ) ;
  void set_patches( int n ) ;

  particle_2D( double x=0 , double y=0 , int val=0 , int id_number=0 ) ;
  particle_2D( int val ) ;
  particle_2D( const particle_2D &copy ) ;
  ~particle_2D() ;
protected:

};

#include "particles2D.ipp"

#endif
