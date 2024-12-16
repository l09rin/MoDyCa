/*****************************************************************************************
   This version of the potential considers the actual particle's charge in the computations .
   If you want a faster and simpler version, you must use the potential
       debye_huckel (it cosiders all the interacting particles as having charge = 1.0e-) .
*****************************************************************************************/

template <typename particle>
class debye_huckel_poly : public interaction<particle> {  // MD  // energy in units of EPSILON
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::idx ;

 public:
  double bjerrum_length , debye_length , amplitude ;
  debye_huckel_poly( double kT , double bl , double dl ) ;
 private:
  particle force ;
  double multiplier = 0 ;   // auxiliary variables
  inline double pair_potential( particle *i_part , particle *j_part , double distance ) ;
  inline particle pair_force( particle *i_part , particle *j_part , double distance ) ;
  inline double pair_potential( double i_charge , double j_charge , double distance ) ;
  inline particle pair_force( particle *i_part , particle *j_part , double i_charge , double j_charge , double distance ) ;
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
  debye_huckel_poly *copy( void ) ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
debye_huckel_poly<particle>::debye_huckel_poly( double kT , double bl , double dl ) {
  bjerrum_length = bl ;
  debye_length = dl ;
  amplitude = kT * bjerrum_length ;
};
/********************************************************************************************/

template <typename particle>
void debye_huckel_poly<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
	              // CUTOFF
  sscanf( file_rec->word[3] , "%lf" , &cutoff ) ;
  cutoff *= debye_length ;
	              // SHIFT
  if( strcmp(file_rec->word[4], "auto") == 0 ) {
    cout << " *** ERROR in interactions_array::debye_huckel_poly::generate_by_record()  : if the charge of particles is polydisperse You cannot shift the interaction potential to 0 for all pairs !" << endl ;
    exit( EXIT_FAILURE ) ;
  }else{
    sscanf( file_rec->word[4] , "%lf" , &SHIFT ) ;
  }
	              // generation of the verlet list
                      // I create the verlet list having to host all bare/screened coulombic interactions
  double verlet_ray ;
  sscanf( file_rec->word[5] , "%lf" , &verlet_ray ) ;
  // OLD, WRONG: verlet_list = new verlet<particle>( conf_ptr , parts_number , cutoff + verlet_ray , verlet_ray , "NO_DISPLACEMENT_VECTOR" ) ;
  verlet_list = new verlet<particle>( conf_ptr , parts_number , cutoff + verlet_ray , verlet_ray ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 1 ; // HERE it is NECESSARY, unless it set to default 0
	              // initialization of the verlet list
  verlet_list->initialize_cellist() ;
  verlet_list->initialize_particles_list( "CHARGES" ) ;    // This function attempts to create the cell list for charges and computes that verlet list for the first time
  verlet_list->verlet_creation_bycell( simulation_kind ) ; // verlet list calculation
};
/********************************************************************************************/

template <typename particle>
inline double debye_huckel_poly<particle>::pair_potential( particle *i_part , particle *j_part , double distance ) {
                  // energy in units of EPSILON
    /***********   DEBYE-HUCKEL POTENTIAL          **********/
  cout << " *** ERROR: the use of this function is not allowed ! Nearest neighbours pointers cannot access to charge particles information !" << endl ;
  exit( EXIT_FAILURE ) ;
};
/********************************************************************************************/

template <typename particle>
inline double debye_huckel_poly<particle>::pair_potential( double i_charge , double j_charge , double distance ) {
                  // energy in units of EPSILON
    /***********   DEBYE-HUCKEL POTENTIAL          **********/
  return amplitude * i_charge * j_charge / distance * exp( -distance / debye_length ) ;
};
/********************************************************************************************/

template <typename particle>
inline particle debye_huckel_poly<particle>::pair_force( particle *i_part , particle *j_part , double distance ) {
              // returns the force j_part applies on i_part, in units of EPSILON/SIGMA
    /***********   DEBYE-HUCKEL FORCE          **********/
  cout << " *** ERROR: the use of this function is not allowed ! Nearest neighbours pointers cannot access to charge particles information !" << endl ;
  exit( EXIT_FAILURE ) ;
};
/********************************************************************************************/

template <typename particle>
inline particle debye_huckel_poly<particle>::pair_force( particle *i_part , particle *j_part , double i_charge , double j_charge , double distance ) {
              // returns the force j_part applies on i_part, in units of EPSILON/SIGMA
    /***********   DEBYE-HUCKEL FORCE          **********/
  multiplier = ( 1.0/distance + 1.0/debye_length ) * amplitude * i_charge * j_charge / (distance*distance) * exp( -distance/debye_length ) ;
  force.position = i_part->position - j_part->position ;
  force.position *= multiplier ;

  return force ;
};
/********************************************************************************************/

template <typename particle>
double debye_huckel_poly<particle>::return_pair_potential( double distance ) {
  return pair_potential( 1.0 , 1.0 , distance ) + SHIFT ;
};
/********************************************************************************************/

template <typename particle>
double debye_huckel_poly<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  double energy = 0 ;
  double dist = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  if( ! verlet_list->DOUBLE_PAIRS ) {
    cout << " *** ERROR : The energy contribution of a single pair of particles cannot be computed if the verlet list is built in MD mode ! " << endl ;
    exit( EXIT_FAILURE ) ;
  }

  verlet_list->god->right_copy( part2 , part1 , apparent_near ) ;
  if( ( dist = part1->position.distance( apparent_near->position ) ) < cutoff ) {
    energy = ( pair_potential( part1->charge , part2->charge , dist ) + SHIFT ) ;
  }

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
particle debye_huckel_poly<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  double dist = i_part->position.distance( j_part->position ) ;
  return pair_force( i_part , j_part , 1.0 , 1.0 , dist ) ;
};
/********************************************************************************************/

template <typename particle>
double debye_huckel_poly<particle>::compute_energy_contribution( void ) {
  double energy = 0 ;
  double dist = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<particle> *element = verlet_list->element ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff ) {
	energy += ( pair_potential( element[i].ptr->charge , walker->ptr->charge , dist ) + SHIFT ) ;
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
inline double debye_huckel_poly<particle>::compute_energy_contribution( particle *part ) {
  // This function computes the energy of ONE particle only
  static double energy ;
  static double dist ;
  static neighbour_list<particle> *walker = NULL ;
  static particle apparent_near ;
  energy = 0 ;
  dist = 0 ;

  if( part->list_ptr[idx] ) {
    walker = part->list_ptr[idx]->neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , part , &apparent_near ) ;
      if( ( dist = part->position.distance( apparent_near.position ) ) < cutoff ) {
	energy += ( pair_potential( part->charge , walker->ptr->charge , dist ) + SHIFT ) ;
      }
      walker = walker->next ;
    }
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
double debye_huckel_poly<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double debye_huckel_poly<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from ideal lattice sites with the external harmonic spring potential
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
void debye_huckel_poly<particle>::compute_forces_contribution(void) {
  double dist = 0 ;
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
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff ) {
	deltaF = pair_force( element[i].ptr , apparent_near , element[i].ptr->charge , walker->ptr->charge , dist ) ;
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
double debye_huckel_poly<particle>::compute_virial_contribution(void) {
  double dist = 0 ;
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
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff ) {
	deltaF = pair_force( element[i].ptr , apparent_near , element[i].ptr->charge , walker->ptr->charge , dist ) ;
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
string debye_huckel_poly<particle>::ostr( void ) {
  ostringstream st ;
  st << "DEBYE HUCKEL CHARGE-POLYDISPERSE POTENTIAL (Bjerrum length , Amplitude , Debye length) : " << bjerrum_length << " , " << amplitude << " , " << debye_length << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
debye_huckel_poly<particle> *debye_huckel_poly<particle>::copy( void ) {
  debye_huckel_poly *copy = new debye_huckel_poly( amplitude/bjerrum_length , bjerrum_length , debye_length ) ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void debye_huckel_poly<particle>::compute_per_part_energy_contribution( double *envec ) {
  double energy = 0 ;
  double dist = 0 ;
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
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff ) {
	energy = 0.5 * ( pair_potential( element[i].ptr->charge , walker->ptr->charge , dist ) + SHIFT ) ;
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
void debye_huckel_poly<particle>::compute_per_part_virial_contribution( double *virvec ) {
  double dist = 0 ;
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
      if( ( dist = element[i].ptr->position.distance( apparent_near->position ) ) < cutoff ) {
	deltaF = pair_force( element[i].ptr , apparent_near , element[i].ptr->charge , walker->ptr->charge , dist ) ;
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
void debye_huckel_poly<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * debye_huckel_poly<particle>::dump_parameters( void ) {
  double *params = new double [3] ;
  params[0] = bjerrum_length ;
  params[1] = debye_length ;
  params[2] = amplitude ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
