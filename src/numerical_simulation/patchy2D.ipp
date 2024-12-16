#include "particles2D.h"

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
patchy_2D::patchy_2D( double x , double y , int val , int id_number , int N_patches ) :
  particle_2D::particle_2D( x , y , val , id_number ) {
  rotation = 0.0 ;
  span_round_angles = 0 ;
  Npatches = 0 ;
  patches = NULL ;
  bonds = NULL ;
  set_patches( N_patches ) ;
  list_ptr = NULL ;
}
/********************************************************************************************/

patchy_2D::patchy_2D( int val , int N_patches ) :
  particle_2D::particle_2D( val ) {
  rotation = 0.0 ;
  span_round_angles = 0 ;
  Npatches = 0 ;
  patches = NULL ;
  bonds = NULL ;
  set_patches( N_patches) ;
  list_ptr = NULL ;
}
/********************************************************************************************/

patchy_2D::patchy_2D( const patchy_2D &copy ) :
  particle_2D::particle_2D( copy ) {
  rotation = copy.rotation ;
  span_round_angles = copy.span_round_angles ;
  Npatches = 0 ;
  patches = NULL ;
  bonds = NULL ;
  set_patches( copy.Npatches) ;
  for( int i=0; i<Npatches; i++ ) patches[i] = copy.patches[i] ;
  list_ptr = NULL ;
}
/********************************************************************************************/

patchy_2D::~patchy_2D() {
  delete[] patches ;
  patches = NULL ;
  delete[] bonds ;
  bonds = NULL ;
  Npatches = 0 ;
  if( list_ptr != NULL ) free(list_ptr) ;
  list_ptr = NULL ;
}
/********************************************************************************************/
/********************************************************************************************/

void patchy_2D::set_patches( int N_patches ) {
  if ( Npatches != 0 ) delete [] patches ;
  Npatches = N_patches ;
  if ( Npatches == 0 ) patches = NULL ;
  else {
    double angle = 2 * M_PI / Npatches ;
    patches = new struct posizione_2D<double> [Npatches] ;
    for( int i=0; i<Npatches; i++ ) {
      patches[i].x = cos(i*angle) ;
      patches[i].y = sin(i*angle) ;
    }
  }
}
/********************************************************************************************/

inline void patchy_2D::operator=(const patchy_2D &part) {
  position.x = part.position.x ;
  position.y = part.position.y ;

  rotation = part.rotation ;
  span_round_angles = part.span_round_angles ;
  for( int i=0; i<Npatches; i++ ) {
    patches[i].x = part.patches[i].x ;
    patches[i].y = part.patches[i].y ;
  }
}
/********************************************************************************************/

inline void patchy_2D::operator=(const particle_2D &part) {
  position.x = part.position.x ;
  position.y = part.position.y ;
}
/********************************************************************************************/

inline void patchy_2D::rndPos_generate( double max ) {  // generate a random position for the particle in the specified box
  particle_2D::rndPos_generate( max ) ;
  rotation = ( M_PI*2 ) * (drand48() - 0.5) ;
}
/********************************************************************************************/

inline void patchy_2D::rndPos_generate( const struct posizione_2D<double> &max ) {  // generate a random position for the particle in the specified box
  particle_2D::rndPos_generate( max ) ;
  rotation = ( M_PI*2 ) * (drand48() - 0.5) ;
}
/********************************************************************************************/

inline void patchy_2D::rand_rotate( double maxangle ) {
  rotation += maxangle * (drand48() - 0.5) ;
}
/********************************************************************************************/

inline struct posizione_2D<double> patchy_2D::patch( int i ) { // this gives the actual orientation of patch i
  struct posizione_2D<double> up2date_patch ;
  // forse conviene lavorare solo con gli angoli? con pbc? faster?
  up2date_patch.x = patches[i].x * cos(rotation) - patches[i].y * sin(rotation) ;
  up2date_patch.y = patches[i].x * sin(rotation) + patches[i].y * cos(rotation) ;
  return up2date_patch ;
}
/********************************************************************************************/

inline string patchy_2D::dump_position_variables( int Ndigits , const char *format , double *parameters ) {
  ostringstream linestream ;
  // to avoid loss of precision on the angle in very long simulations:
  if( rotation > M_PI ) {
    span_round_angles += floor( ( rotation / M_PI + 1 ) / 2 ) ;
    rotation -= floor( ( rotation / M_PI + 1 ) / 2 ) * M_PI * 2 ;
  } else if( rotation < -M_PI ) {
    span_round_angles += ceil( ( rotation / M_PI - 1 ) / 2 ) ;
    rotation -= ceil( ( rotation / M_PI - 1 ) / 2 ) * M_PI * 2 ;
  }
  if( format == NULL )  linestream << setprecision(Ndigits) << position.x << " " << position.y ;
  else {
    if( strcmp( format , "all" ) == 0 )  linestream << setprecision(Ndigits) << position.x << " " << position.y << " " << rotation + M_PI*2*span_round_angles ;
    else if( strcmp( format , "ptc" ) == 0 || strcmp( format , "patch" ) == 0 ) {
      // here parameters are the patch opening and the maximum interaction range of the patches
      linestream << setprecision(Ndigits) << (char)('a' + type-1) << " " << position.x << " " << position.y << " 0.0 " << radius << " " << cos(parameters[1]) << " " << 2.0*(parameters[0] + parameters[2]) ;
      // 3D rotation matrix expanded along the rows
      linestream << " " << cos(rotation) << " " << sin(rotation) << " 0.0" ;
      linestream << " " << -sin(rotation) << " " << cos(rotation) << " 0.0" ;
      linestream << " 0.0 0.0 1.0" ;
      // bonded neighbours ids :
      for( int ip=0; ip<Npatches; ip++ ) linestream << " " << ((bonds[ip] == NULL) ? -1 : bonds[ip]->id) ;
    } else  linestream << setprecision(Ndigits) << position.x << " " << position.y ;
  }
  string line = linestream.str() ;

  return line ;
}
/********************************************************************************************/
