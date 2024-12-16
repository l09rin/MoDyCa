/*****************************************************************************************/

template <typename particle>
class hs_square_well : public interaction<particle> {
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::name ;
  using interaction<particle>::idx ;

public:
  double CORE_RADIUS ;   // hard core radius
  double SQUARE_WELL_RADIUS ;   // radial extension of the attractive/repulsive well
  double SQUARE_WELL_ENERGY ;   // square well potential energy
  hs_square_well( void ) ;
  hs_square_well( double core , double sqwell , double energy ) ;

private:
  double distance ;
  double CORE_DIAMETER = 0.0 , SQUARE_WELL_DIAMETER = 0.0;
  inline double pair_potential( particle *i_part , particle *j_part ) ;
  inline particle pair_force( particle *i_part , particle *j_part ) ;

 public:
  double return_pair_potential( double distance ) ;
  double return_pair_potential( particle *part1 , particle *part2 ) ;
  particle return_pair_force( particle *i_part , particle *j_part ) ;
  void generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) ;

  double compute_energy_contribution( void ) ;
  inline double compute_energy_contribution( particle *part ) ;
  double delta_energy_constrained_CoM( particle *disp ) ;
  double lambda_derivative_perpart( void ) ;
  void compute_forces_contribution( void ) ;
  double compute_virial_contribution( void ) ;
  void compute_per_part_energy_contribution( double *envec ) ;
  void compute_per_part_virial_contribution( double *virvec ) ;

  string ostr( void ) ;
  void build_bonds_network( double bonds_cutoff = 0.0 ) ;
  double *dump_parameters( void ) ;
  hs_square_well *copy( void ) ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
hs_square_well<particle>::hs_square_well( void ) {
  CORE_RADIUS = 0.5 ;
  CORE_DIAMETER = 1.0 ;
  SQUARE_WELL_RADIUS = 0.1 ;
  SQUARE_WELL_DIAMETER = 1.2 ;
  SQUARE_WELL_ENERGY = -1.0 ;
  cutoff = SQUARE_WELL_DIAMETER ;
  sprintf( name , "%s" , "hs_squarewell\0" ) ;
};
/********************************************************************************************/

template <typename particle>
hs_square_well<particle>::hs_square_well( double core , double sqwell , double energy ) {
  CORE_RADIUS = core ;
  CORE_DIAMETER = 2.0 * CORE_RADIUS ;
  SQUARE_WELL_RADIUS = sqwell ;
  SQUARE_WELL_DIAMETER = 2.0 * SQUARE_WELL_RADIUS ;
  SQUARE_WELL_ENERGY = energy ;
  cutoff = SQUARE_WELL_DIAMETER ;
  sprintf( name , "%s" , "hs_squarewell\0" ) ;
};
/********************************************************************************************/

template <typename particle>
void hs_square_well<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
  SHIFT = 0 ;
  double verlet_radius , verlet_skin ;
  sscanf( file_rec->word[1] , "%lf" , &CORE_RADIUS ) ;
  sscanf( file_rec->word[2] , "%lf" , &SQUARE_WELL_RADIUS ) ;
  sscanf( file_rec->word[3] , "%lf" , &SQUARE_WELL_ENERGY ) ;
  sscanf( file_rec->word[4] , "%lf" , &verlet_radius ) ;
  CORE_DIAMETER = 2.0 * CORE_RADIUS ;
  SQUARE_WELL_DIAMETER = 2.0 * SQUARE_WELL_RADIUS ;
  cutoff = SQUARE_WELL_DIAMETER ;
  verlet_skin = verlet_radius - cutoff ;
  if( verlet_skin <= 0 ) {
    cout << "*** ERROR : verlet cells must be bigger than the potential's cutoff!" << endl ;
    exit(EXIT_FAILURE) ;
  }
      // generation of the verlet list
  verlet_list = new verlet<particle>( conf_ptr , parts_number , verlet_radius , verlet_skin ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 1 ;
	              // initialization of the verlet list
  verlet_list->initialize_cellist() ;
  verlet_list->initialize_particles_list( "ORDERED" ) ;
  verlet_list->verlet_creation_bycell( simulation_kind ) ; // verlet list calculation
      // setting up the particles' radii
  for( int i=0 ; i<verlet_list->particles_in_list ; i++ ) verlet_list->element[i].ptr->radius = CORE_RADIUS ;
};
/********************************************************************************************/

template <typename particle>
inline double hs_square_well<particle>::pair_potential( particle *i_part , particle *j_part ) {
    // energy in units of EPSILON
    /***********   HARD SPHERE + SQUARE WELL POTENTIAL          **********/
  distance = ( j_part->position - i_part->position ).norm() ;
  if( distance < CORE_DIAMETER ) return INFINITY ;
  else if( distance < SQUARE_WELL_DIAMETER ) return SQUARE_WELL_ENERGY ;
  else return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
inline particle hs_square_well<particle>::pair_force( particle *i_part , particle *j_part ) {
    // returns the force j_part applies on i_part, in units of EPSILON
  cout << "This potential has an Hard Core functional shape, cannot be used for force computations !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
double hs_square_well<particle>::return_pair_potential( double distance ) {
  if( distance < CORE_DIAMETER ) return INFINITY ;
  else if( distance < SQUARE_WELL_DIAMETER ) return SQUARE_WELL_ENERGY ;
  else return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double hs_square_well<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  double energy ;

  verlet_list->god->right_copy( part2 , part1 , apparent_near ) ;
  energy = pair_potential( part1 , apparent_near ) ;

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
particle hs_square_well<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  cout << "*** ERROR : HARD SPHERE + SQUARE WELL potential can only be used with MonteCarlo integration !!" << endl ;
  particle force ;
  force.position = force.position * nanl("") ;
  return force ;
};
/********************************************************************************************/

template <typename particle>
double hs_square_well<particle>::compute_energy_contribution( void ) {
  double energy = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;
      energy += pair_potential( element[i].ptr , apparent_near ) ;
      walker = walker->next ;
    }
  }

  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) return 0.5 * energy ;
  else return energy ;
};
/********************************************************************************************/

template <typename particle>
inline double hs_square_well<particle>::compute_energy_contribution( particle *part ) {
  //  verlet_list->DOUBLE_PAIRS  HAS TO BE 1 !!
  // This function computes the energy of ONE particle only
  static double energy ;
  static neighbour_list<particle> *walker = NULL ;
  static particle apparent_near ;
  energy = 0 ;

  if( part->list_ptr[idx] ) {
    walker = part->list_ptr[idx]->neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , part , &apparent_near ) ;
      energy += pair_potential( part , &apparent_near ) ;
      walker = walker->next ;
    }
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
double hs_square_well<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double hs_square_well<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from ideal lattice sites with the external harmonic spring potential
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
void hs_square_well<particle>::compute_forces_contribution(void) {
  cout << "*** ERROR : I cannot properly compute forces for step-wise discontinuous potentials !!" << endl ;
  exit(EXIT_FAILURE) ;
};
/********************************************************************************************/

template <typename particle>
double hs_square_well<particle>::compute_virial_contribution( void ) {
  cout << "*** ERROR : I cannot properly compute the virial for step-wise discontinuous potentials !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
string hs_square_well<particle>::ostr( void ) {
  ostringstream st ;
  st <<  "HARD SPHERE + SQUARE WELL potential parameters: ( CORE_RADIUS=" << CORE_RADIUS << " , SQ-WELL_RADIUS=" << SQUARE_WELL_RADIUS << " , SQUARE_WELL_ENERGY=" << SQUARE_WELL_ENERGY << " )" << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
hs_square_well<particle> *hs_square_well<particle>::copy( void ) {
  hs_square_well *copy = new hs_square_well( CORE_RADIUS , SQUARE_WELL_RADIUS , SQUARE_WELL_ENERGY ) ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void hs_square_well<particle>::compute_per_part_energy_contribution( double *envec ) {
  double energy = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) envec[i] = 0 ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;
      energy = 0.5 * pair_potential( element[i].ptr , apparent_near ) ;
      envec[element[i].ptr->id] += energy ;
      envec[walker->ptr->id] += energy ;
      walker = walker->next ;
    }
  }

  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) envec[i] *= 0.5 ;
};
/********************************************************************************************/

template <typename particle>
void hs_square_well<particle>::compute_per_part_virial_contribution( double *virvec ) {
  cout << "*** ERROR : I cannot properly compute the virial for step-wise discontinuous potentials !!" << endl ;
};
/********************************************************************************************/

template <typename particle>
void hs_square_well<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * hs_square_well<particle>::dump_parameters( void ) {
  double *params = new double [3] ;
  params[0] = CORE_RADIUS ;
  params[1] = SQUARE_WELL_RADIUS ;
  params[2] = SQUARE_WELL_ENERGY ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
