#ifndef PARTICLES3D_H
#define PARTICLES3D_H

using namespace std;
#include <cmath>
#include <sstream>
#include "vector3D.h"

class particle_2D ;
template <typename particle> struct vlist_element ;

class particle_3D {

public:
  struct posizione_3D<double> position ;
  struct posizione_3D<int> periodic_box ;
  int id , type , mol ;
  double charge ;
  double mass ;
  double radius ;
  int valence ;
  int *bonded_monomer ;
  vlist_element<particle_3D> **list_ptr = NULL ;

  inline void set_position( const posizione_3D<double> &vec ) ;
  inline void operator=( const particle_3D &part ) ;
  inline void operator=( const particle_2D &part ) ;
  //  inline void operator%( double max_size ) ;
  //  inline void operator%( const particle_3D &max_size ) ;
  inline particle_2D reduce_z( void ) ;
  inline void rndPos_generate( double max ) ;
  inline void rndPos_generate( const struct posizione_3D<double> &max ) ;
  inline void Gaussian_particle( void ) ;
  inline void random_unit_vector( void ) ;
  void clear( void ) ;
  inline string dump_position_variables( int Ndigits = 12 , const char *format = NULL , double *parameters = NULL ) ;
  friend ostream & operator<< ( ostream &out, const particle_3D &part ) ;
  //  friend ofstream & operator<< ( ofstream &fout, const particle_3D &part ) ;
  inline void read_by_string( char* line ) ;
  inline void read_by_string( char** words ) ;
  inline void read_by_string_pbc( char* line ) ;
  inline void read_by_string_pbc( char** words ) ;
  inline double max_comp( void ) ;
  inline double min_comp( void ) ;
  inline double avg_comp( void ) ;
  inline void set_equal_comp( double coord ) ;
  inline double volume( void ) ;
  inline particle_3D PBC_NearNeigh2( const particle_3D &part0 , const particle_3D &box_side ) ;
  inline posizione_3D<double> unwrapped( posizione_3D<double> box_sides ) ;
  void set_patches( int n ) ;

  particle_3D( double x=0 , double y=0 , double z=0 , int val=0 , int id_number=0 ) ;
  particle_3D( int val ) ;
  particle_3D( const particle_3D &copy ) ;
  ~particle_3D() ;
protected:

};

#include "particles3D.ipp"

#endif
