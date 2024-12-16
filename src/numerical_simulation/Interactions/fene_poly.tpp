/*****************************************************************************************/

template <typename particle>
class fene_poly : public interaction<particle> {     // MD  // energy in units of EPSILON
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::idx ;

 public:
  double k_F , R_f , two_R_f , energy_scale , two_energy_scale ;   // spring constant and maximum extension factor for the polydisperse FENE_POLY potential (bonding interaction)
  fene_poly( double kf , double max_ext_factor ) ;
 private:
  particle force ;
  double multiplier ;  //auxiliary variable                       // max_extension = max_ext_factor * monomers_diameter
  inline double pair_potential( particle *i_part , particle *j_part , double distance , double max_extension ) ;
  inline particle pair_force( particle *i_part , particle *j_part , double distance , double max_extension ) ;
 public:
  double return_pair_potential( double distance ) ;
  double return_pair_potential( double distance , double max_extension ) ;
  double return_pair_potential( particle *part1 , particle *part2 ) ;
  particle return_pair_force( particle *i_part , particle *j_part ) ;
  particle return_pair_force( particle *i_part , particle *j_part , double max_extension ) ;

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
  fene_poly *copy( void ) ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
fene_poly<particle>::fene_poly( double kf , double max_ext_factor ) {
  k_F = kf ;
  R_f = max_ext_factor ;
  energy_scale = k_F * R_f * R_f ;
  two_R_f = 2.0 * R_f ;
  two_energy_scale = 2.0 * energy_scale ;
};
/********************************************************************************************/

template <typename particle>
void fene_poly<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
	              // CUTOFF
  cutoff = R_f ;
	              // SHIFT
  SHIFT = 0 ;
	              // generation of the verlet list
  verlet_list = new verlet<particle>( conf_ptr , parts_number, cutoff , 0 , "NO_DISPLACEMENT_VECTOR" ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 0 ; // Indeed it is the default value
	              // initialization of the verlet list
  verlet_list->initialize_particles_list( "BONDS" ) ;
  verlet_list->bonding_verlet_creation( simulation_kind ) ; // verlet list calculation
};
/********************************************************************************************/

template <typename particle>
inline double fene_poly<particle>::pair_potential( particle *i_part , particle *j_part , double distance , double max_extension ) {
    /***********   FENE_POLY POTENTIAL          **********/
  multiplier = distance/max_extension ;
  return - energy_scale * log( 1.0 - (multiplier*multiplier) ) ;
};
/********************************************************************************************/

template <typename particle>
inline particle fene_poly<particle>::pair_force( particle *i_part , particle *j_part , double distance , double max_extension ) {
    /***********   FENE_POLY FORCE    **********/                                           // forces are in units of EPSILON/SIGMA
  multiplier = two_energy_scale / ( max_extension*max_extension - distance*distance ) ;
  force.position = j_part->position - i_part->position ;
  force.position *= multiplier ;

  return force ;
};
/********************************************************************************************/

template <typename particle>
double fene_poly<particle>::return_pair_potential( double distance ) {
  cout << "To return the value of the polydisperse fene potential the dimer's maximum extension is needed !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
double fene_poly<particle>::return_pair_potential( double distance , double max_extension ) {
  return pair_potential( NULL , NULL , distance , max_extension ) + SHIFT ;
};
/********************************************************************************************/

template <typename particle>
double fene_poly<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  double energy = 0 ;
  double dist = 0 ;
  double actual_max_ext = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  verlet_list->god->right_copy( part2 , part1 , apparent_near ) ;
  actual_max_ext = two_R_f * part2->radius ;
  if( ( dist = part1->position.distance( apparent_near->position ) ) < actual_max_ext ) {
    energy = ( pair_potential( part1 , apparent_near , dist , actual_max_ext ) ) ;  // NO SHIFT HERE !
  }

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
particle fene_poly<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  cout << "To return the value of the polydisperse fene force the dimer's maximum extension is needed !!" << endl ;
  particle force ;
  force.position = force.position * nanl("") ;
  return force ;
};
/********************************************************************************************/

template <typename particle>
particle fene_poly<particle>::return_pair_force( particle *i_part , particle *j_part , double max_extension ) {
  double dist = i_part->position.distance( j_part->position ) ;
  return pair_force( i_part , j_part , dist , max_extension ) ;
};
/********************************************************************************************/

template <typename particle>
double fene_poly<particle>::compute_energy_contribution( void ) {
  double energy = 0 ;
  double dist = 0 ;
  double actual_max_ext = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
      actual_max_ext = two_R_f * walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < actual_max_ext ) {
	energy += ( pair_potential( element[i].ptr , apparent_near , dist , actual_max_ext ) ) ;  // NO SHIFT HERE !
      }
      walker = walker->next ;
    }
  }

  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) return 0.5 * energy ;
  else return energy ;
};
/********************************************************************************************/

template <typename particle>
inline double fene_poly<particle>::compute_energy_contribution( particle *part ) {
  // This function computes the energy of ONE particle only
  static double actual_max_ext ;
  static double energy ;
  static double dist ;
  static neighbour_list<particle> *walker = NULL ;
  static particle apparent_near ;
  energy = 0 ;
  dist = 0 ;
  actual_max_ext = 0 ;

  if( part->list_ptr[idx] ) {
    walker = part->list_ptr[idx]->neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , part , &apparent_near ) ;
      actual_max_ext = two_R_f * walker->ptr->radius ;
      if( ( dist = part->position.distance( apparent_near.position ) ) < actual_max_ext ) {
	energy += ( pair_potential( part , &apparent_near , dist , actual_max_ext ) ) ;  // NO SHIFT HERE !
      }
      walker = walker->next ;
    }
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
double fene_poly<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double fene_poly<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from ideal lattice sites with the external harmonic spring potential
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
void fene_poly<particle>::compute_forces_contribution(void) {
  double dist = 0 ;
  double actual_max_ext = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL , deltaF , **forces_t2 = verlet_list->god->forces_t2 ;
  apparent_near = new particle ;
  if( verlet_list->DOUBLE_PAIRS ) {
    cout << " *** ERROR : The forces computation with lists built in MC mode has not been implemented yet ! " << endl ;
    exit( EXIT_FAILURE ) ;
  }

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
      actual_max_ext = two_R_f * walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < actual_max_ext ) {
	deltaF = pair_force( element[i].ptr , apparent_near , dist , actual_max_ext ) ;
	forces_t2[ element[i].ptr->id ]->position += deltaF.position ;         // in verlet list i store i-j connection only for i or j
	forces_t2[ walker->ptr->id ]->position -= deltaF.position ;
      }
      walker = walker->next ;
    }
  }
  delete apparent_near ;
};
/********************************************************************************************/

template <typename particle>
double fene_poly<particle>::compute_virial_contribution(void) {
  double dist = 0 ;
  double actual_max_ext = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL , deltaF ;
  apparent_near = new particle ;
  double virial = 0 ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
      actual_max_ext = two_R_f * walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < actual_max_ext ) {
	deltaF = pair_force( element[i].ptr , apparent_near , dist , actual_max_ext ) ;
	// in verlet list i store i-j connection only for i or j
	virial += ( deltaF.position * ( element[i].ptr->position - apparent_near->position ) ) ;
      }
      walker = walker->next ;
    }
  }
  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) return 0.5 * virial ;
  else return virial ;
};
/********************************************************************************************/

template <typename particle>
string fene_poly<particle>::ostr( void ) {
  ostringstream st ;
  st << "FENE_POLY POTENTIAL ( force constant , maximum extension factor ) : " << k_F << " , " << R_f << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
fene_poly<particle> *fene_poly<particle>::copy( void ) {
  fene_poly *copy = new fene_poly( k_F , R_f ) ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void fene_poly<particle>::compute_per_part_energy_contribution( double *envec ) {
  double energy = 0 ;
  double dist = 0 ;
  double actual_max_ext = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) envec[i] = 0 ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
      actual_max_ext = two_R_f * walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < actual_max_ext ) {
	energy = 0.5 * ( pair_potential( element[i].ptr , apparent_near , dist , actual_max_ext ) ) ;  // NO SHIFT HERE !
	envec[element[i].ptr->id] += energy ;
	envec[walker->ptr->id] += energy ;
      }
      walker = walker->next ;
    }
  }

  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) envec[i] *= 0.5 ;
};
/********************************************************************************************/

template <typename particle>
void fene_poly<particle>::compute_per_part_virial_contribution( double *virvec ) {
  double dist = 0 ;
  double actual_max_ext = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL , deltaF ;
  apparent_near = new particle ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) virvec[i] = 0 ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
      actual_max_ext = two_R_f * walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < actual_max_ext ) {
	deltaF = pair_force( element[i].ptr , apparent_near , dist , actual_max_ext ) ;
	// in verlet list i store i-j connection only for i or j
	virvec[ element[i].ptr->id ] += ( deltaF.position * element[i].ptr->position ) ;
	virvec[ walker->ptr->id ] -= ( deltaF.position * apparent_near->position ) ;
      }
      walker = walker->next ;
    }
  }
  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) virvec[i] *= 0.5 ;
};
/********************************************************************************************/

template <typename particle>
void fene_poly<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * fene_poly<particle>::dump_parameters( void ) {
  double *params = new double [3] ;
  params[0] = k_F ;
  params[1] = R_f ;
  params[2] = energy_scale ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
