#include "particles3D.h"

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
particle_2D::particle_2D( double x , double y , int val , int id_number ) {
  position.x = x ;
  position.y = y ;
  periodic_box.x = 0 ;
  periodic_box.y = 0 ;

  valence = val ;
  if( valence == 0 ) {
    bonded_monomer = NULL ;
  }else{
    bonded_monomer = new int [valence] ;
    for( int i=0; i<valence; i++ ) bonded_monomer[i] = -1 ;
  }
  charge = 0 ;
  id = id_number ;
  type = 0 ;
  mol = 0 ;
  mass = 1.0 ;
  radius = 1.0 ;
  list_ptr = NULL ;
}
/********************************************************************************************/

particle_2D::particle_2D( int val ) {
  position.x = 0 ;
  position.y = 0 ;
  periodic_box.x = 0 ;
  periodic_box.y = 0 ;

  valence = val ;
  if( valence == 0 ) {
    bonded_monomer = NULL ;
  }else{
    bonded_monomer = new int [valence] ;
    for( int i=0; i<valence; i++ ) bonded_monomer[i] = -1 ;
  }
  charge = 0 ;
  id = 0 ;
  type = 0 ;
  mol = 0 ;
  mass = 1.0 ;
  radius = 1.0 ;
  list_ptr = NULL ;
}
/********************************************************************************************/

particle_2D::particle_2D( const particle_2D &copy ) {
  position.x = copy.position.x ;
  position.y = copy.position.y ;
  periodic_box.x = copy.periodic_box.x ;
  periodic_box.y = copy.periodic_box.y ;

  valence = copy.valence ;
  if( valence == 0 ) {
    bonded_monomer = NULL ;
  }else{
    bonded_monomer = new int [valence] ;
    for( int i=0; i<valence; i++ ) bonded_monomer[i] = copy.bonded_monomer[i] ;
  }
  charge = copy.charge ;
  id = copy.id ;
  type = copy.type ;
  mol = copy.mol ;
  mass = copy.mass ;
  radius = copy.radius ;
  list_ptr = NULL ;
}
/********************************************************************************************/

particle_2D::~particle_2D() {
  delete[] bonded_monomer ;
  if( list_ptr != NULL ) free(list_ptr) ;
  list_ptr = NULL ;
}
/********************************************************************************************/
/********************************************************************************************/

inline void particle_2D::set_position( const posizione_2D<double> &vec ) {
  position.x = vec.x ;
  position.y = vec.y ;
}
/********************************************************************************************/

void particle_2D::clear(void) {         // sets to 0 all coordinates
  position.x = 0 ;
  position.y = 0 ;
}
/********************************************************************************************/

inline void particle_2D::rndPos_generate( double max ) {  // generate a random position for the particle in the specified box
  position.x = drand48() * max ;
  position.y = drand48() * max ;
}
/********************************************************************************************/

inline void particle_2D::rndPos_generate( const struct posizione_2D<double> &max ) {  // generate a random position for the particle in the specified box
  position.x = drand48() * max.x ;
  position.y = drand48() * max.y ;
}
/********************************************************************************************/

inline void particle_2D::Gaussian_particle( void ) {                                                  // generate a position distributed according a 2D gaussian with 0 mean and unitary variance
  double aux1 = drand48() ;
  double aux2 = drand48() ;

  position.x = -sqrt(-2.*log(aux1))*cos(2.*M_PI*aux2) ;
  position.y = -sqrt(-2.*log(aux1))*sin(2.*M_PI*aux2) ;
}
/********************************************************************************************/

inline void particle_2D::random_unit_vector( void ) {                                                 // generate a position uniformly distributed on the surface of a sphere of unit radius
  double auxAngle = drand48() ;

  position.x = cos(2.*M_PI*auxAngle) ;
  position.y = sin(2.*M_PI*auxAngle) ;
}
/********************************************************************************************/

/* inline void particle_2D::operator%( double max_size ) { */
/*   while( position.x < 0 ) (position.x) += max_size ; */
/*   while( position.x > max_size ) (position.x) -= max_size ; */
/*   while( position.y < 0 ) (position.y) += max_size ; */
/*   while( position.y > max_size ) (position.y) -= max_size ; */
/* } */
/* /\********************************************************************************************\/ */

/* inline void particle_2D::operator%(const particle_2D &max_size) { */
/*   while( position.x < 0 ) (position.x) += max_size.position.x ; */
/*   while( position.x > max_size.position.x  ) (position.x) -= max_size.position.x ; */
/*   while( position.y < 0 ) (position.y) += max_size.position.y ; */
/*   while( position.y > max_size.position.y  ) (position.y) -= max_size.position.y ; */
/* } */
/* /\********************************************************************************************\/ */

inline void particle_2D::operator=(const particle_2D &part) {
  position.x = part.position.x ;
  position.y = part.position.y ;
}
/********************************************************************************************/

inline void particle_2D::operator=(const particle_3D &part) {
  position.x = part.position.x ;
  position.y = part.position.y ;
}
/********************************************************************************************/

inline string particle_2D::dump_position_variables( int Ndigits , const char *format , double *parameters ) {
  ostringstream linestream ;
  linestream << setprecision(Ndigits) << position.x << " " << position.y ;
  string line = linestream.str() ;

  return line ;
}
/********************************************************************************************/

ostream & operator<< ( ostream &out, const particle_2D &part ) {
  out << "(" << part.position.x << " , " << part.position.y << ")" << endl ;
  out << "id: " << part.id << " , charge: " << part.charge << " , mass: " << part.mass << " , valence: " << part.valence << endl ;
  for( int i=0 ; i<part.valence ; i++ ) out << part.bonded_monomer[i] << " " ;
  out << endl ;

  return out ;
}
/********************************************************************************************/

/* ofstream & operator<< ( ofstream &fout, const particle_2D &part ) { */
/*   fout << part.position.x << " " << part.position.y ; */

/*   return fout ; */
/* } */
/********************************************************************************************/

inline void particle_2D::read_by_string( char *line ) {
  sscanf( line , "%lf %lf" , &(position.x) , &(position.y) ) ;
}
/********************************************************************************************/

inline void particle_2D::read_by_string( char **words ) {
  sscanf( words[0] , "%lf" , &(position.x) ) ;
  sscanf( words[1] , "%lf" , &(position.y) ) ;
}
/********************************************************************************************/

inline void particle_2D::read_by_string_pbc( char *line ) {
  sscanf( line , "%lf %lf %d %d" , &(position.x) , &(position.y) , &(periodic_box.x) , &(periodic_box.y) ) ;
}
/********************************************************************************************/

inline void particle_2D::read_by_string_pbc( char **words ) {
  sscanf( words[0] , "%lf" , &(position.x) ) ;
  sscanf( words[1] , "%lf" , &(position.y) ) ;
  sscanf( words[2] , "%d" , &(periodic_box.x) ) ;
  sscanf( words[3] , "%d" , &(periodic_box.y) ) ;
}
/********************************************************************************************/

inline double particle_2D::max_comp( void ) {
  double max = (position.x>position.y) ? position.x : position.y ;

  return max ;
}
/********************************************************************************************/

inline double particle_2D::min_comp( void ) {
  double min = (position.x<position.y) ? position.x : position.y ;

  return min ;
}
/********************************************************************************************/

inline double particle_2D::avg_comp( void ) {
  double avg = ( position.x + position.y ) / 2.0 ;

  return avg ;
}
/********************************************************************************************/

inline void particle_2D::set_equal_comp( double coord ) {
  position.x = coord ;
  position.y = coord ;
}
/********************************************************************************************/

inline double particle_2D::volume( void ) {
  double vol = position.x * position.y ;

  return vol ;
}
/********************************************************************************************/

inline particle_2D particle_2D::PBC_NearNeigh2( const particle_2D &part0 , const particle_2D &box_side ) {
  particle_2D nearest_neighbour , midside = box_side ;
  midside.position *= 0.5 ;

  if( fabs( position.x - part0.position.x ) < midside.position.x ) {   // Do i-particle interact with j-particle rather than one of its copies?  x axis...
    nearest_neighbour.position.x = position.x ;
  }else if( ( position.x - part0.position.x ) > 0 ) {
    nearest_neighbour.position.x = position.x - box_side.position.x * (int)( ( (position.x - part0.position.x) / midside.position.x + 1 ) / 2 ) ;
  }else{
    nearest_neighbour.position.x = position.x - box_side.position.x * (int)( ( (position.x - part0.position.x) / midside.position.x - 1 ) / 2 ) ;
  }
  if( fabs( position.y - part0.position.y ) < midside.position.y ) {   // y axis...
    nearest_neighbour.position.y = position.y ;
  }else if( ( position.y - part0.position.y ) > 0 ) {
    nearest_neighbour.position.y = position.y - box_side.position.y * (int)( ( (position.y - part0.position.y) / midside.position.y + 1 ) / 2 ) ;
  }else{
    nearest_neighbour.position.y = position.y - box_side.position.y * (int)( ( (position.y - part0.position.y) / midside.position.y - 1 ) / 2 ) ;
  }

  return nearest_neighbour ;
}
/********************************************************************************************/

inline posizione_2D<double> particle_2D::unwrapped( posizione_2D<double> box_sides ) {
  posizione_2D<double> upos ;
  upos.x = position.x + box_sides.x * periodic_box.x ;
  upos.y = position.y + box_sides.y * periodic_box.y ;

  return upos ;
}
/********************************************************************************************/

void particle_2D::set_patches( int n ) {
  cout << "*** ERROR : This particle class is not compatible with patchy particles !!" << endl ;
  exit(EXIT_FAILURE) ;
}
/********************************************************************************************/
