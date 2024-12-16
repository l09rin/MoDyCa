/*****************************************************************************************/

template <typename particle>
class lennard_jones_poly : public interaction<particle> { // MD  // energy in units of EPSILON
  double cutoff_factor = 0 ;    // The actual cutoff is cutoff_factor * ( radius_atom[i] + radius_atom[j] )
  using interaction<particle>::verlet_list ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::idx ;

 private:
  particle force ;
  double multiplier = 0 , r6_inv = 0 , r2_inv = 0 ;   // auxiliary variables
  inline double pair_potential( particle *i_part , particle *j_part , double distance , double diameter ) ;
  inline particle pair_force( particle *i_part , particle *j_part , double distance , double diameter ) ;

 public:
  double return_pair_potential( double distance ) ;
  double return_pair_potential( double distance , double diameter ) ;
  double return_pair_potential( particle *part1 , particle *part2 ) ;
  particle return_pair_force( particle *i_part , particle *j_part ) ;
  particle return_pair_force( particle *i_part , particle *j_part , double diameter ) ;
  void generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) ;

 public:
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
  lennard_jones_poly *copy( void ) ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

/********************************************************************************************/

template <typename particle>
void lennard_jones_poly<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
	              // CUTOFF
  sscanf( file_rec->word[1] , "%lf" , &cutoff_factor ) ;
	              // SHIFT
  if( strcmp(file_rec->word[2], "auto") == 0 ) {
    SHIFT = 0 ;
    SHIFT -= pair_potential( NULL , NULL , cutoff_factor , 1.0 ) ;
  }else{
    sscanf( file_rec->word[2] , "%lf" , &SHIFT ) ;
  }
	              // generation of the verlet list
  double verlet_ray = 0 , max_radius = 0 ;
  sscanf( file_rec->word[3] , "%lf" , &verlet_ray ) ;   // computation of the verlet radius from the maximum diameter of atoms
  // This cycle has to be modified in case the lj does not act on all the particles in the system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111
  for( int i=0 ; i<conf_ptr->numberOfPart ; i++ ) max_radius = ( conf_ptr->particles[i]->radius > max_radius ) ? conf_ptr->particles[i]->radius : max_radius ;
  verlet_list = new verlet<particle>( conf_ptr , parts_number , 2.0 * max_radius + verlet_ray , verlet_ray ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 1 ; // Indeed it has been already done in the constr
	              // initialization of the verlet list
  verlet_list->initialize_cellist() ;
  verlet_list->initialize_particles_list( "ORDERED" ) ;
  verlet_list->verlet_creation_bycell( simulation_kind ) ; // verlet list calculation
};
/********************************************************************************************/

template <typename particle>
inline double lennard_jones_poly<particle>::pair_potential( particle *i_part , particle *j_part , double distance , double diameter ) {
    /***********   LENNARD-JONES POTENTIAL     offset energy Vlj(Rc=3)    **********/
  r2_inv = diameter * diameter / ( distance*distance ) ;
  r6_inv = r2_inv * r2_inv * r2_inv ;
  return 4.0 * ( r6_inv * r6_inv - r6_inv ) ;
};
/********************************************************************************************/

template <typename particle>
inline particle lennard_jones_poly<particle>::pair_force( particle *i_part , particle *j_part , double distance , double diameter ) {
    /***********   LENNARD-JONES POTENTIAL     **********/                                           // forces are in units of EPSILON/SIGMA
  force.position = i_part->position - j_part->position ;
  multiplier = diameter * diameter ;
  r2_inv = multiplier / ( distance*distance ) ;
  r6_inv = r2_inv * r2_inv * r2_inv ;
  multiplier = ( 48.0 * r6_inv - 24.0 ) * r6_inv * r2_inv / multiplier ;
  force.position *= multiplier ;

  return force ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_poly<particle>::return_pair_potential( double distance ) {
  cout << "To return the value of the polydisperse lennard jones potential the dimer's diameter is needed !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_poly<particle>::return_pair_potential( double distance , double diameter ) {
  return pair_potential( NULL , NULL , distance , diameter ) + SHIFT ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_poly<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  double energy = 0 ;
  double dist = 0 ;
  double actual_diameter = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  verlet_list->god->right_copy( part2 , part1 , apparent_near ) ;
  actual_diameter = part1->radius + part2->radius ;
  if( ( dist = part1->position.distance( apparent_near->position ) ) < cutoff_factor * actual_diameter ) {
    energy = ( pair_potential( part1 , apparent_near , dist , actual_diameter ) + SHIFT ) ;
  }

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
particle lennard_jones_poly<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  cout << "To return the value of the polydisperse lennard jones force the dimer's diameter is needed !!" << endl ;
  particle force ;
  force.position = force.position * nanl("") ;
  return force ;
};
/********************************************************************************************/

template <typename particle>
particle lennard_jones_poly<particle>::return_pair_force( particle *i_part , particle *j_part , double diameter ) {
  double dist = i_part->position.distance( j_part->position ) ;
  return pair_force( i_part , j_part , dist , diameter ) ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_poly<particle>::compute_energy_contribution( void ) {
  double energy = 0 ;
  double dist = 0 ;
  double actual_diameter = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
      actual_diameter = element[i].ptr->radius + walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff_factor * actual_diameter ) {
	energy += ( pair_potential( element[i].ptr , apparent_near , dist , actual_diameter ) + SHIFT ) ;
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
inline double lennard_jones_poly<particle>::compute_energy_contribution( particle *part ) {
  // This function computes the energy of ONE particle only
  static double dist ;
  static double energy ;
  static double actual_diameter ;
  static neighbour_list<particle> *walker = NULL ;
  static particle apparent_near ;
  energy = 0 ;

  if( part->list_ptr[idx] ) {
    walker = part->list_ptr[idx]->neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , part , &apparent_near ) ;
      actual_diameter = part->radius + walker->ptr->radius ;
      if( ( dist = part->position.distance( apparent_near.position ) ) < cutoff_factor * actual_diameter ) {
	energy += ( pair_potential( part , &apparent_near , dist , actual_diameter ) + SHIFT ) ;
      }
      walker = walker->next ;
    }
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_poly<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_poly<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from ideal lattice sites with the external harmonic spring potential
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
void lennard_jones_poly<particle>::compute_forces_contribution(void) {
  double dist = 0 ;
  double actual_diameter = 0 ;
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
      actual_diameter = element[i].ptr->radius + walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff_factor * actual_diameter ) {
	deltaF = pair_force( element[i].ptr , apparent_near , dist , actual_diameter ) ;
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
double lennard_jones_poly<particle>::compute_virial_contribution(void) {
  double dist = 0 ;
  double actual_diameter = 0 ;
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
      actual_diameter = element[i].ptr->radius + walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff_factor * actual_diameter ) {
	deltaF = pair_force( element[i].ptr , apparent_near , dist , actual_diameter ) ;
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
string lennard_jones_poly<particle>::ostr( void ) {
  ostringstream st ;
  st << "POLYDISPERSE LENNARD-JONES potential cutoff factor: " << cutoff_factor << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
lennard_jones_poly<particle> *lennard_jones_poly<particle>::copy( void ) {
  lennard_jones_poly *copy = new lennard_jones_poly() ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void lennard_jones_poly<particle>::compute_per_part_energy_contribution( double *envec ) {
  double energy = 0 ;
  double dist = 0 ;
  double actual_diameter = 0 ;
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
      actual_diameter = element[i].ptr->radius + walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff_factor * actual_diameter ) {
	energy = 0.5 * ( pair_potential( element[i].ptr , apparent_near , dist , actual_diameter ) + SHIFT ) ;
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
void lennard_jones_poly<particle>::compute_per_part_virial_contribution( double *virvec ) {
  double dist = 0 ;
  double actual_diameter = 0 ;
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
      actual_diameter = element[i].ptr->radius + walker->ptr->radius ;
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff_factor * actual_diameter ) {
	deltaF = pair_force( element[i].ptr , apparent_near , dist , actual_diameter ) ;
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
void lennard_jones_poly<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * lennard_jones_poly<particle>::dump_parameters( void ) {
  double *params = new double [1] ;
  params[0] = cutoff_factor ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
