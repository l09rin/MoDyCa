/*****************************************************************************************/
/* This is an external field for which each particle is bounded with an harmonic potential to its initial position */

template <typename particle>
class harmonic_spring_field_wignerseitz : public interaction<particle> { // energy in units of EPSILON
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::name ;
  using interaction<particle>::idx ;

 public:
  double lambda = 0 ; // spring constant
  harmonic_spring_field_wignerseitz( double Lambda = 0.0 ) ;
 private:
  particle *equilibrium_positions = NULL ;
  particle f ;
  double multiplier = 0 ;   // auxiliary variables
  int Nvac = 0 ;
  particle ***near_vacancies = NULL , *vacancy_sites = NULL ;
  inline double potential( particle *part , particle *part0 ) ;
  inline particle force( particle *part , particle *part0 ) ;

 public:
  double return_pair_potential( double distance ) ;
  double return_pair_potential( particle *part1 , particle *part2 ) ;
  particle return_pair_force( particle *part , particle *part0 ) ;
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
  harmonic_spring_field_wignerseitz *copy( void ) ;
  ~harmonic_spring_field_wignerseitz() ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
harmonic_spring_field_wignerseitz<particle>::harmonic_spring_field_wignerseitz( double Lambda ) {
  lambda = Lambda ;
  sprintf( name , "%s" , "sprWignerSeitz\0" ) ;
};
/********************************************************************************************/

template <typename particle>
harmonic_spring_field_wignerseitz<particle>::~harmonic_spring_field_wignerseitz( void ) {
  delete [] equilibrium_positions ;
  equilibrium_positions = NULL ;
  lambda = 0 , cutoff = 0 ;
  for( int i=0 ; i<verlet_list->particles_in_list ; i++ ) free( near_vacancies[i] ) ;
  delete [] near_vacancies ;
  near_vacancies = NULL ;
  if( Nvac > 0 ) {
    delete [] vacancy_sites ;
    vacancy_sites = NULL ;
  }
};
/********************************************************************************************/

template <typename particle>
void harmonic_spring_field_wignerseitz<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
  SHIFT = 0.0 ;    // SHIFT , not used here
  sscanf( file_rec->word[1] , "%lf" , &lambda ) ;
  sscanf( file_rec->word[2] , "%lf" , &cutoff ) ;   // CUTOFF to build neighbour list for WS check
	              // non - generation of the verlet list
  verlet_list = new verlet<particle>( conf_ptr , parts_number, cutoff , 0.0 , "NO_DISPLACEMENT_VECTOR" ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 0 ; // Indeed it has been already done in the constr
	              // initialization of the verlet list
  verlet_list->initialize_cellist() ;
  verlet_list->initialize_particles_list( "ORDERED" ) ;
  verlet_list->verlet_creation_bycell( simulation_kind ) ; // verlet list calculation
  // initialization of the array of equilibrium positions
  // here we consider only the case in which every particle in the system feels this potential
  if( parts_number != conf_ptr->numberOfPart ) {
    cout << "*** ERROR: if you want to apply the harmonic spring field potential you should change the code!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  equilibrium_positions = new particle [conf_ptr->numberOfPart] ;
  for( int i=0; i<conf_ptr->numberOfPart; i++ ) equilibrium_positions[i].position = conf_ptr->particles[i]->unwrapped( conf_ptr->box_sides.position ) ;

  // initialization of the array of vacancies positions, if any
  near_vacancies = new particle** [conf_ptr->numberOfPart] ;
  for( int i=0; i<conf_ptr->numberOfPart; i++ ) {
    near_vacancies[i] = (particle **)calloc( 1 , sizeof(particle *) ) ;
    near_vacancies[i][0] = NULL ;
  }
  if( file_rec->words_number > 3 ) {
    ifstream vacancies_file ;
    record *rec = new record ;
    int i = 0 ;
    if( strcmp( file_rec->word[3] , "vacancies") == 0 ) {
      vacancies_file.open( file_rec->word[4] ) ;
      if( vacancies_file.fail() ) {
	cout << "*** ERROR: Impossible to open the file " << file_rec->word[4] << " !" << endl << endl ;
	exit( EXIT_FAILURE ) ;
      } else {
	rec->getrecord( vacancies_file ) ;
	rec->split() ;
	if( rec->words_number > 0 ) {
	  if( rec->word[0][0] == '#' ) i++ ;
	  if( strcmp( rec->word[i] , "N" ) == 0 ) sscanf( rec->word[i+1] , "%d" , &Nvac ) ;
	}
	if( Nvac < 1 ) {
	  cout << "*** ERROR: Unrecognized format of file " << file_rec->word[4] << " !" << endl << endl ;
	  exit( EXIT_FAILURE ) ;
	} else {
	  vacancy_sites = new particle [Nvac] ;
	  int vac_Nneighatoms = 0 , *neigh_array = NULL , neigh_array_length = 0 , neigh_atom = 0 ;
	  int *neigh_Nvac = new int [conf_ptr->numberOfPart] ;
	  for( int i=0; i<Nvac; i++ ) {
	    rec->getrecord( vacancies_file ) ;
	    vacancy_sites[i].read_by_string( rec->line ) ;
	    vac_Nneighatoms = verlet_list->get_neighbour_indices( vacancy_sites+i , &neigh_array_length , &neigh_array ) ;
	    for( int j=0; j<vac_Nneighatoms; j++ ) {
	      neigh_atom = neigh_array[j] ;
	      near_vacancies[ neigh_atom ][ neigh_Nvac[neigh_atom] ] = vacancy_sites+i ;
	      neigh_Nvac[ neigh_atom ]++ ;
	      near_vacancies[ neigh_atom ] = (particle **)realloc( near_vacancies[ neigh_atom ] , (neigh_Nvac[neigh_atom]+1)*sizeof(particle *) ) ;
	      near_vacancies[ neigh_atom ][ neigh_Nvac[neigh_atom] ] = NULL ;
	    }
	  }
	  free( neigh_array ) ;
	  neigh_array = NULL ;
	  neigh_array_length = 0 ;
	  delete [] neigh_Nvac ;
	}
      }
    }
  }
};
/********************************************************************************************/

template <typename particle>
inline double harmonic_spring_field_wignerseitz<particle>::potential( particle *part , particle *part0 ) {
    /***********   HARMONIC EXTERNAL POTENTIAL     Vi = lambda * (ri-ri0)^2    **********/
  return lambda * ( part->position - part0->position ).square_norm() ;
};
/********************************************************************************************/

template <typename particle>
inline particle harmonic_spring_field_wignerseitz<particle>::force( particle *part , particle *part0 ) {
    /***********   HARMONIC POTENTIAL     **********/               // forces are in units of EPSILON/SIGMA
  f.position = part0->position - part->position ;
  f.position *= ( 2.0*lambda ) ;

  return f ;
};
/********************************************************************************************/

template <typename particle>
double harmonic_spring_field_wignerseitz<particle>::return_pair_potential( double distance ) {
  return lambda * distance * distance ;
};
/********************************************************************************************/

template <typename particle>
double harmonic_spring_field_wignerseitz<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  // This is an external field, it doesn't generate any pair contribution, hence gives 0
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
particle harmonic_spring_field_wignerseitz<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  // This is an external field, it doesn't generate any pair contribution, hence gives 0
  f.position *= 0.0 ;
  return f ;
};
/********************************************************************************************/

template <typename particle>
double harmonic_spring_field_wignerseitz<particle>::compute_energy_contribution( void ) {
  double energy = 0 ;
  static particle upart , apparent_near , **near_vacancy ;
  neighbour_list<particle> *walker = NULL ;
  vlist_element<particle> *element = verlet_list->element ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  particle *box_sides = &( verlet_list->god->box_sides ) ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    upart.position = particles[i]->unwrapped( box_sides->position ) ;
    energy += ( potential( &upart , &(equilibrium_positions[i]) ) + SHIFT ) ;
    // check if the particle is in its own Wigner-Seitz cell
    upart.position -= ( equilibrium_positions[i].position + verlet_list->god->CoM_displacement.position ) ;
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( &(equilibrium_positions[walker->ptr->id]) , &(equilibrium_positions[i]) , &apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour site
      apparent_near.position -= equilibrium_positions[i].position ;
      if( ( (upart.position * apparent_near.position) / apparent_near.position.square_norm() ) > 0.5 ) 	return INFINITY ;
      walker = walker->next ;
    }
    // check if the particle does not invade defects' Wigner-Seitz cells
    near_vacancy = near_vacancies[i] ;
    while( *near_vacancy ) {
      verlet_list->god->right_copy( *near_vacancy , &(equilibrium_positions[i]) , &apparent_near ) ;
      apparent_near.position -= equilibrium_positions[i].position ;
      if( ( (upart.position * apparent_near.position) / apparent_near.position.square_norm() ) > 0.5 ) 	return INFINITY ;
      near_vacancy++ ;
    }
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
inline double harmonic_spring_field_wignerseitz<particle>::compute_energy_contribution( particle *part ) {
  // This function computes the energy of ONE particle only
  static double energy ;
  static particle upart , apparent_near , **near_vacancy ;
  static neighbour_list<particle> *walker = NULL ;
  energy = 0 ;

  particle **particles = verlet_list->god->particles ;
  particle *box_sides = &( verlet_list->god->box_sides ) ;
  if( part->id == particles[part->id]->id ) {
    upart.position = particles[ part->id ]->unwrapped( box_sides->position ) ;
    energy = ( potential( &upart , &(equilibrium_positions[ part->id ]) ) + SHIFT ) ;
    // check if the particle is in its own Wigner-Seitz cell
    upart.position -= ( equilibrium_positions[ part->id ].position + verlet_list->god->CoM_displacement.position ) ;
    walker = verlet_list->element[ part->id ].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( &(equilibrium_positions[walker->ptr->id]) , &(equilibrium_positions[part->id]) , &apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour site
      apparent_near.position -= equilibrium_positions[part->id].position ;
      if( ( (upart.position * apparent_near.position) / apparent_near.position.square_norm() ) > 0.5 ) 	return INFINITY ;
      walker = walker->next ;
    }
    // check if the particle does not invade defects' Wigner-Seitz cells
    near_vacancy = near_vacancies[ part->id ] ;
    while( *near_vacancy ) {
      verlet_list->god->right_copy( *near_vacancy , &(equilibrium_positions[part->id]) , &apparent_near ) ;
      apparent_near.position -= equilibrium_positions[part->id].position ;
      if( ( (upart.position * apparent_near.position) / apparent_near.position.square_norm() ) > 0.5 ) 	return INFINITY ;
      near_vacancy++ ;
    }
  } else {
    cout << "ERROR: FIX THIS PROBLEM WITH NEIGHBOUR LISTS, PLEASE..." << endl ;
    exit( EXIT_FAILURE ) ;
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
double harmonic_spring_field_wignerseitz<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  return - lambda * ( disp->position * ( (verlet_list->god->CoM_displacement.position * 2.0) - (disp->position / verlet_list->god->numberOfPart) ) ) ;
};
/********************************************************************************************/

template <typename particle>
double harmonic_spring_field_wignerseitz<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from the ideal lattice sites with the constrained centre of mass
  double msd = 0 ;
  particle upart ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  particle *box_sides = &( verlet_list->god->box_sides ) ;
  particle CoM_displacement = verlet_list->god->CoM_displacement ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    upart.position = particles[i]->unwrapped( box_sides->position ) ;
    msd += ( upart.position - CoM_displacement.position - equilibrium_positions[i].position ).square_norm() ;
  }

  return msd / particles_in_list ;
};
/********************************************************************************************/

template <typename particle>
void harmonic_spring_field_wignerseitz<particle>::compute_forces_contribution(void) {
  particle upart ;
  particle deltaF , **forces_t2 = verlet_list->god->forces_t2 ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  particle *box_sides = &( verlet_list->god->box_sides ) ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    upart.position = particles[i]->unwrapped( box_sides->position ) ;
    deltaF = force( &upart , &(equilibrium_positions[i]) ) ;
    forces_t2[i]->position += deltaF.position ;         // in verlet list i store i-j connection only for i or j
  }
};
/********************************************************************************************/

template <typename particle>
double harmonic_spring_field_wignerseitz<particle>::compute_virial_contribution(void) {
  particle upart ;
  particle deltaF ;
  double virial = 0 ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  particle *box_sides = &( verlet_list->god->box_sides ) ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    upart.position = particles[i]->unwrapped( box_sides->position ) ;
    deltaF = force( &upart , &(equilibrium_positions[i]) ) ;
    virial += ( deltaF.position * upart.position ) ;
  }
  return virial ;
};
/********************************************************************************************/

template <typename particle>
string harmonic_spring_field_wignerseitz<particle>::ostr( void ) {
  ostringstream st ;
  st << "HARMONIC SPRING FIELD with WIGNER-SEITZ cells constraint (spring constant) : " << lambda << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
harmonic_spring_field_wignerseitz<particle> *harmonic_spring_field_wignerseitz<particle>::copy( void ) {
  harmonic_spring_field_wignerseitz *copy = new harmonic_spring_field_wignerseitz( lambda ) ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void harmonic_spring_field_wignerseitz<particle>::compute_per_part_energy_contribution( double *envec ) {
  particle upart ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  particle *box_sides = &( verlet_list->god->box_sides ) ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    upart.position = particles[i]->unwrapped( box_sides->position ) ;
    envec[i] += ( potential( &upart , &(equilibrium_positions[i]) ) + SHIFT ) ;
  }
};
/********************************************************************************************/

template <typename particle>
void harmonic_spring_field_wignerseitz<particle>::compute_per_part_virial_contribution( double *virvec ) {
  particle upart ;
  particle deltaF ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  particle *box_sides = &( verlet_list->god->box_sides ) ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    upart.position = particles[i]->unwrapped( box_sides->position ) ;
    deltaF = force( &upart , &(equilibrium_positions[i]) ) ;
    virvec[i] += ( deltaF.position * upart.position ) ;
  }
};
/********************************************************************************************/

template <typename particle>
void harmonic_spring_field_wignerseitz<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * harmonic_spring_field_wignerseitz<particle>::dump_parameters( void ) {
  double *params = new double [1] ;
  params[0] = lambda ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
