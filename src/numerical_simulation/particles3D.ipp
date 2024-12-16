#include "particles2D.h"

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
particle_3D::particle_3D( double x , double y , double z , int val , int id_number ) {
  position.x = x ;
  position.y = y ;
  position.z = z ;
  periodic_box.x = 0 ;
  periodic_box.y = 0 ;
  periodic_box.z = 0 ;

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

particle_3D::particle_3D( int val ) {
  position.x = 0 ;
  position.y = 0 ;
  position.z = 0 ;
  periodic_box.x = 0 ;
  periodic_box.y = 0 ;
  periodic_box.z = 0 ;

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

particle_3D::particle_3D( const particle_3D &copy ) {
  position.x = copy.position.x;
  position.y = copy.position.y;
  position.z = copy.position.z;
  periodic_box.x = copy.periodic_box.x;
  periodic_box.y = copy.periodic_box.y;
  periodic_box.z = copy.periodic_box.z;

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

particle_3D::~particle_3D() {
  delete [] bonded_monomer ;
  if( list_ptr != NULL ) free(list_ptr) ;
  list_ptr = NULL ;
}
/********************************************************************************************/
/********************************************************************************************/

inline void particle_3D::set_position( const posizione_3D<double> &vec ) {
  position.x = vec.x ;
  position.y = vec.y ;
  position.z = vec.z ;
}
/********************************************************************************************/

void particle_3D::clear( void ) {         // sets to 0 all coordinates
  position.x = 0;
  position.y = 0;
  position.z = 0;
}
/********************************************************************************************/

inline void particle_3D::rndPos_generate( double max ) {  // generate a random position for the particle in the specified box
  position.x = drand48() * max ;
  position.y = drand48() * max ;
  position.z = drand48() * max ;
}
/********************************************************************************************/

inline void particle_3D::rndPos_generate( const struct posizione_3D<double> &max ) {  // generate a random position for the particle in the specified box
  position.x = drand48() * max.x ;
  position.y = drand48() * max.y ;
  position.z = drand48() * max.z ;
}
/********************************************************************************************/

inline void particle_3D::Gaussian_particle( void ) {                                                  // generate a position distributed according a 3D gaussian with 0 mean and unitary variance
  double aux1 = drand48();
  double aux2 = drand48();
  double aux3 = drand48();
  double aux4 = drand48();

  position.x = -sqrt(-2.*log(aux1))*cos(2.*M_PI*aux2);
  position.y = -sqrt(-2.*log(aux1))*sin(2.*M_PI*aux2);
  position.z = -sqrt(-2.*log(aux3))*cos(2.*M_PI*aux4);
}
/********************************************************************************************/

inline void particle_3D::random_unit_vector( void ) {                                                 // generate a position uniformly distributed on the surface of a sphere of unit radius
  double aux1, aux2, R = 2 ;

  while( R > 1. ) {
    aux1 = 1.0 - 2.0 * drand48() ;
    aux2 = 1.0 - 2.0 * drand48() ;
    R = aux1 * aux1 + aux2 * aux2 ;
  }

  double aux3 = 2. * sqrt( 1. - R ) ;
  position.x = aux1 * aux3 ;
  position.y = aux2 * aux3 ;
  position.z = 1. - 2. * R ;
}
/********************************************************************************************/

/* inline void particle_3D::operator%(double max_size) { */
/*   while( position.x < 0 ) (position.x) += max_size; */
/*   while( position.x > max_size ) (position.x) -= max_size; */
/*   while( position.y < 0 ) (position.y) += max_size; */
/*   while( position.y > max_size ) (position.y) -= max_size; */
/*   while( position.z < 0 ) (position.z) += max_size; */
/*   while( position.z > max_size ) (position.z) -= max_size; */
/* } */
/* /\********************************************************************************************\/ */

/* inline void particle_3D::operator%(const particle_3D &max_size) { */
/*   while( position.x < 0 ) (position.x) += max_size.position.x ; */
/*   while( position.x > max_size.position.x  ) (position.x) -= max_size.position.x ; */
/*   while( position.y < 0 ) (position.y) += max_size.position.y ; */
/*   while( position.y > max_size.position.y  ) (position.y) -= max_size.position.y ; */
/*   while( position.z < 0 ) (position.z) += max_size.position.z ; */
/*   while( position.z > max_size.position.z  ) (position.z) -= max_size.position.z ; */
/* } */
/* /\********************************************************************************************\/ */

inline void particle_3D::operator=(const particle_3D &part) {
  position.x = part.position.x;
  position.y = part.position.y;
  position.z = part.position.z;
}
/********************************************************************************************/

inline void particle_3D::operator=(const particle_2D &part) {
  position.x = part.position.x;
  position.y = part.position.y;
  position.z = 0.0;
}
/********************************************************************************************/

inline string particle_3D::dump_position_variables( int Ndigits , const char *format , double *parameters ) {
  ostringstream linestream ;
  linestream << setprecision(Ndigits) << position.x << " " << position.y << " " << position.z ;
  string line = linestream.str() ;

  return line ;
}
/********************************************************************************************/

ostream & operator<< ( ostream &out, const particle_3D &part ) {
  out << "(" << part.position.x << " , " << part.position.y << " , " << part.position.z << ")" << endl ;
  out << "id: " << part.id << " , charge: " << part.charge << " , mass: " << part.mass << " , valence: " << part.valence << endl ;
  for( int i=0 ; i<part.valence ; i++ ) out << part.bonded_monomer[i] << " " ;
  out << endl ;

  return out ;
}
/********************************************************************************************/

/* ofstream & operator<< ( ofstream &fout, const particle_3D &part ) { */
/*   fout << part.position.x << " " << part.position.y << " " << part.position.z ; */

/*   return fout ; */
/* } */
/********************************************************************************************/

inline void particle_3D::read_by_string( char *line ) {
  sscanf( line , "%lf %lf %lf" , &(position.x) , &(position.y) , &(position.z) ) ;
}
/********************************************************************************************/

inline void particle_3D::read_by_string( char **words ) {
  sscanf( words[0] , "%lf" , &(position.x) ) ;
  sscanf( words[1] , "%lf" , &(position.y) ) ;
  sscanf( words[2] , "%lf" , &(position.z) ) ;
}
/********************************************************************************************/

inline void particle_3D::read_by_string_pbc( char *line ) {
  sscanf( line , "%lf %lf %lf %d %d %d" , &(position.x) , &(position.y) , &(position.z) , &(periodic_box.x) , &(periodic_box.y) , &(periodic_box.z) ) ;
}
/********************************************************************************************/

inline void particle_3D::read_by_string_pbc( char **words ) {
  sscanf( words[0] , "%lf" , &(position.x) ) ;
  sscanf( words[1] , "%lf" , &(position.y) ) ;
  sscanf( words[2] , "%lf" , &(position.z) ) ;
  sscanf( words[3] , "%d" , &(periodic_box.x) ) ;
  sscanf( words[4] , "%d" , &(periodic_box.y) ) ;
  sscanf( words[5] , "%d" , &(periodic_box.z) ) ;
}
/********************************************************************************************/

inline particle_2D particle_3D::reduce_z( void ) {
  particle_2D reduced_part = particle_2D( position.x, position.y, valence ) ;
  if( valence > 0 ) {
    reduced_part.bonded_monomer = new int [valence] ;
    for( int i=0; i<valence; i++ ) reduced_part.bonded_monomer[i] = bonded_monomer[i] ;
  }
  reduced_part.charge = charge ;
  reduced_part.mass = mass ;

  return reduced_part ;
}
/********************************************************************************************/

inline double particle_3D::max_comp( void ) {
  double max = (position.x>position.y) ? ( (position.x>position.z) ? position.x : position.z ) : ( (position.y>position.z) ? position.y : position.z ) ;

  return max ;
}
/********************************************************************************************/

inline double particle_3D::min_comp( void ) {
  double min = (position.x<position.y) ? ( (position.x<position.z) ? position.x : position.z ) : ( (position.y<position.z) ? position.y : position.z ) ;

  return min ;
}
/********************************************************************************************/

inline double particle_3D::avg_comp( void ) {
  double avg = ( position.x + position.y + position.z ) / 3.0 ;

  return avg ;
}
/********************************************************************************************/

inline void particle_3D::set_equal_comp( double coord ) {
  position.x = coord ;
  position.y = coord ;
  position.z = coord ;
}
/********************************************************************************************/

inline double particle_3D::volume( void ) {
  double vol = position.x * position.y * position.z ;

  return vol ;
}
/********************************************************************************************/

inline particle_3D particle_3D::PBC_NearNeigh2( const particle_3D &part0 , const particle_3D &box_side ) {
  particle_3D nearest_neighbour , midside = box_side ;
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
  if( fabs( position.z - part0.position.z ) < midside.position.z ) {    // z axis...
    nearest_neighbour.position.z = position.z ;
  }else if( ( position.z - part0.position.z ) > 0 ) {
    nearest_neighbour.position.z = position.z - box_side.position.z * (int)( ( (position.z - part0.position.z) / midside.position.z + 1 ) / 2 ) ;
  }else{
    nearest_neighbour.position.z = position.z - box_side.position.z * (int)( ( (position.z - part0.position.z) / midside.position.z - 1 ) / 2 ) ;
  }

  return nearest_neighbour ;
}
/********************************************************************************************/

inline posizione_3D<double> particle_3D::unwrapped( posizione_3D<double> box_sides ) {
  posizione_3D<double> upos ;
  upos.x = position.x + box_sides.x * periodic_box.x ;
  upos.y = position.y + box_sides.y * periodic_box.y ;
  upos.z = position.z + box_sides.z * periodic_box.z ;

  return upos ;
}
/********************************************************************************************/

void particle_3D::set_patches( int n ) {
  cout << "*** ERROR : This particle class is not compatible with patchy particles !!" << endl ;
  exit(EXIT_FAILURE) ;
}
/********************************************************************************************/
