/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
molecule<particle>::molecule( void ) {
  gyration = 0 ;
  gyration_square = 0 ;
  monomers_diameter = 0 ;
  Natoms = 0 ;
  atom = NULL ;
}
/********************************************************************************************/

template <typename particle>
molecule<particle>::molecule( const molecule &other ) {
  cm = other.cm ;
  vcm = other.vcm ;
  fcm = other.fcm ;
  omega = other.omega ;
  gyration = other.gyration ;
  gyration_square = other.gyration_square ;
  monomers_diameter = other.monomers_diameter ;
  Natoms = other.Natoms ;
  if( Natoms != 0 ) {
    atom = (particle **)calloc( Natoms , sizeof(particle *) ) ;
    if( atom == NULL ) {
      cout << "\n Allocation error in molecule constructor" ;
      exit(EXIT_FAILURE) ;
    }
    for( int i=0 ; i<Natoms ; i++ ) atom[i] = other.atom[i] ;
  } else {
    atom = NULL ;
  }
}
/********************************************************************************************/

template <typename particle>
molecule<particle>::~molecule() {
  for( int i=0 ; i<Natoms ; i++ ) delete atom[i] ;
  free( atom ) ;
  atom = NULL ;
  Natoms = 0 ;
}
/********************************************************************************************/
/********************************************************************************************/

template <typename particle>
inline void molecule<particle>::compute_CoM( particle box_sides ) {
  cm.position.clear() ;
  cm.periodic_box.clear() ;
  for( int i=0 ; i<Natoms ; i++ ) cm.position += atom[i]->unwrapped( box_sides.position ) ;
  cm.position /= (double)Natoms ;
}
/********************************************************************************************/

template <typename particle>
inline void molecule<particle>::compute_gyration_square( particle box_sides ) {
  gyration_square = 0 ;
  for( int i=0 ; i<Natoms ; i++ ) gyration_square += ( atom[i]->unwrapped( box_sides.position ) - cm.unwrapped( box_sides.position ) ).square_norm() ;
  gyration_square /= (double)Natoms ;
}
/********************************************************************************************/

template <typename particle>
inline void molecule<particle>::compute_gyration( void ) {
  gyration = sqrt( gyration_square ) ;
}
/********************************************************************************************/
