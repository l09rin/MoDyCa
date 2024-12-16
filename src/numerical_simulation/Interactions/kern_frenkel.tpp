/*****************************************************************************************/

template <typename particle>
class kern_frenkel : public interaction<particle> {
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::name ;
  using interaction<particle>::idx ;

public:
  double CORE_RADIUS ;   // hard core radius
  double EPSILON ; // patch to patch bonding energy
  double LAMBDA ;  // LAMBDA is the radial extension of the square well patch
  double DELTA ; // delta is the angular opening of patches: solid angle (3D) or half angle (2D)
  kern_frenkel( void ) ;
  kern_frenkel( double core , double epsilon , double lambda , double delta ) ;

private:
  particle Rij , Rji ;
  int bond_flag = 0 ;
  double distance ;
  double COSDELTA = 0 ;
  double CORE_DIAMETER = 0.0 ;
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
  kern_frenkel *copy( void ) ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
kern_frenkel<particle>::kern_frenkel( void ) {
  CORE_RADIUS = 1.0 ;
  CORE_DIAMETER = 2.0 * CORE_RADIUS ;
  EPSILON = 1.0 ;
  LAMBDA = 0.2 ;
  DELTA = 0.1 ;
  COSDELTA = cos(DELTA) ;
  cutoff = CORE_DIAMETER + 2.0 * LAMBDA ;
  sprintf( name , "%s" , "kernfrenkel\0" ) ;
};
/********************************************************************************************/

template <typename particle>
kern_frenkel<particle>::kern_frenkel( double core , double epsilon , double lambda , double delta ) {
  CORE_RADIUS = core ;
  CORE_DIAMETER = 2.0 * CORE_RADIUS ;
  EPSILON = epsilon ;
  LAMBDA = lambda ;
  DELTA = delta ;
  COSDELTA = cos(DELTA) ;
  cutoff = CORE_DIAMETER + 2.0 * LAMBDA ;
  sprintf( name , "%s" , "kernfrenkel\0" ) ;
};
/********************************************************************************************/

template <typename particle>
void kern_frenkel<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
  SHIFT = 0 ;
  double verlet_radius , verlet_skin ;
  sscanf( file_rec->word[1] , "%lf" , &CORE_RADIUS ) ;
  sscanf( file_rec->word[2] , "%lf" , &LAMBDA ) ;
  sscanf( file_rec->word[3] , "%lf" , &DELTA ) ;
  sscanf( file_rec->word[4] , "%lf" , &EPSILON ) ;
  sscanf( file_rec->word[5] , "%lf" , &verlet_radius ) ;
  CORE_DIAMETER = 2.0 * CORE_RADIUS ;
  cutoff = CORE_DIAMETER + 2.0 * LAMBDA ;
  verlet_skin = verlet_radius - cutoff ;
  COSDELTA = cos(DELTA) ;
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
inline double kern_frenkel<particle>::pair_potential( particle *i_part , particle *j_part ) {
  cout << "*** ERROR : KERN-FRENKEL potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
};
/********************************************************************************************/

template <>
inline double kern_frenkel<patchy_2D>::pair_potential( patchy_2D *i_part , patchy_2D *j_part ) {
    // energy in units of EPSILON
  //   ONE BOND PER PATCH IMPOSED !!
    /***********   KERN-FRENKEL PATCHY POTENTIAL          **********/
  Rij.position = j_part->position - i_part->position ;
  distance = Rij.position.norm() ;
  if( distance < CORE_DIAMETER ) return INFINITY ;
  else if( distance < cutoff ) {
    Rij.position /= distance ;
    Rji.position = Rij.position * (-1.0) ;
    bond_flag = 1 ;
    for( int i=0 ; i<i_part->Npatches && bond_flag ; i++ ) {
      if( i_part->patch(i) * Rij.position > COSDELTA ) {
	for( int j=0 ; j<j_part->Npatches && bond_flag ; j++ ) {
	  if( j_part->patch(j) * Rji.position > COSDELTA ) bond_flag = 0 ;
	}
      }
    }
    if( bond_flag == 1 ) return 0.0 ;
    else return -EPSILON ;
  } else return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
inline particle kern_frenkel<particle>::pair_force( particle *i_part , particle *j_part ) {
    // returns the force j_part applies on i_part, in units of EPSILON
  cout << "This potential has an Hard Core functional shape, cannot be used for force computations !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
double kern_frenkel<particle>::return_pair_potential( double distance ) {
  cout << "To return the value of the patchy potential the distance is not enough !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
double kern_frenkel<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  cout << "*** ERROR : KERN-FRENKEL potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
};
/********************************************************************************************/

template <>
double kern_frenkel<patchy_2D>::return_pair_potential( patchy_2D *part1 , patchy_2D *part2 ) {
  patchy_2D *apparent_near = NULL ;
  apparent_near = new patchy_2D ;
  double energy ;

  verlet_list->god->right_copy( part2 , part1 , apparent_near ) ;
  apparent_near->Npatches = part2->Npatches ;
  apparent_near->patches = part2->patches ;
  apparent_near->rotation = part2->rotation ;
  energy = pair_potential( part1 , apparent_near ) ;
  apparent_near->Npatches = 0 ;
  apparent_near->patches = NULL ;

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
particle kern_frenkel<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  cout << "*** ERROR : KERN-FRENKEL potential can only be used with MonteCarlo integration !!" << endl ;
  particle force ;
  force.position = force.position * nanl("") ;
  return force ;
};
/********************************************************************************************/

template <typename particle>
double kern_frenkel<particle>::compute_energy_contribution( void ) {
  cout << "*** ERROR : KERN-FRENKEL potential can only be used with the particle class patchy_2D !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <>
double kern_frenkel<patchy_2D>::compute_energy_contribution( void ) {
  double energy = 0 ;
  neighbour_list<patchy_2D> *walker = NULL ;
  patchy_2D *apparent_near = NULL ;
  apparent_near = new patchy_2D ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<patchy_2D> *element = verlet_list->element ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;
      apparent_near->Npatches = walker->ptr->Npatches ;
      apparent_near->patches = walker->ptr->patches ;
      apparent_near->rotation = walker->ptr->rotation ;
      energy += pair_potential( element[i].ptr , apparent_near ) ;
      walker = walker->next ;
    }
  }

  apparent_near->Npatches = 0 ;
  apparent_near->patches = NULL ;
  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) return 0.5 * energy ;
  else return energy ;
};
/********************************************************************************************/

template <typename particle>
inline double kern_frenkel<particle>::compute_energy_contribution( particle *part ) {
  cout << "*** ERROR : KERN-FRENKEL potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
};
/********************************************************************************************/

template <>
inline double kern_frenkel<patchy_2D>::compute_energy_contribution( patchy_2D *part ) {
  //  verlet_list->DOUBLE_PAIRS  HAS TO BE 1 !!
  // This function computes the energy of ONE particle only
  static double energy ;
  static neighbour_list<patchy_2D> *walker = NULL ;
  static patchy_2D apparent_near ;
  energy = 0 ;

  if( part->list_ptr[idx] ) {
    walker = part->list_ptr[idx]->neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , part , &apparent_near ) ;
      apparent_near.Npatches = walker->ptr->Npatches ;
      apparent_near.patches = walker->ptr->patches ;
      apparent_near.rotation = walker->ptr->rotation ;
      energy += pair_potential( part , &apparent_near ) ;
      walker = walker->next ;
    }
  }

  apparent_near.Npatches = 0 ;
  apparent_near.patches = NULL ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
double kern_frenkel<particle>::delta_energy_constrained_CoM( particle *disp ) {
  cout << "*** ERROR : KERN-FRENKEL potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
};
/********************************************************************************************/

template <>
double kern_frenkel<patchy_2D>::delta_energy_constrained_CoM( patchy_2D *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double kern_frenkel<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from ideal lattice sites with the external harmonic spring potential
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
void kern_frenkel<particle>::compute_forces_contribution(void) {
  cout << "*** ERROR : I cannot properly compute forces for step-wise discontinuous potentials !!" << endl ;
  exit(EXIT_FAILURE) ;
};
/********************************************************************************************/

template <typename particle>
double kern_frenkel<particle>::compute_virial_contribution( void ) {
  cout << "*** ERROR : I cannot properly compute the virial for step-wise discontinuous potentials !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
string kern_frenkel<particle>::ostr( void ) {
  ostringstream st ;
  st <<  "KERN-FRENKEL potential parameters: ( CORE_RADIUS=" << CORE_RADIUS << " ; EPSILON=" << EPSILON << " ; LAMBDA=" << LAMBDA << " ; DELTA=" << DELTA << " )" << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
kern_frenkel<particle> *kern_frenkel<particle>::copy( void ) {
  kern_frenkel *copy = new kern_frenkel( CORE_RADIUS , EPSILON , LAMBDA , DELTA ) ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void kern_frenkel<particle>::compute_per_part_energy_contribution( double *envec ) {
  cout << "*** ERROR : KERN-FRENKEL potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
};
/********************************************************************************************/

template <>
void kern_frenkel<patchy_2D>::compute_per_part_energy_contribution( double *envec ) {
  double energy = 0 ;
  neighbour_list<patchy_2D> *walker = NULL ;
  patchy_2D *apparent_near = NULL ;
  apparent_near = new patchy_2D ;

  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<patchy_2D> *element = verlet_list->element ;
  for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) envec[i] = 0 ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;
      apparent_near->Npatches = walker->ptr->Npatches ;
      apparent_near->patches = walker->ptr->patches ;
      apparent_near->rotation = walker->ptr->rotation ;
      energy = 0.5 * pair_potential( element[i].ptr , apparent_near ) ;
      envec[element[i].ptr->id] += energy ;
      envec[walker->ptr->id] += energy ;
      walker = walker->next ;
    }
  }

  apparent_near->Npatches = 0 ;
  apparent_near->patches = NULL ;
  delete apparent_near ;
  if( verlet_list->DOUBLE_PAIRS ) for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) envec[i] *= 0.5 ;
};
/********************************************************************************************/

template <typename particle>
void kern_frenkel<particle>::compute_per_part_virial_contribution( double *virvec ) {
  cout << "*** ERROR : I cannot properly compute the virial for step-wise discontinuous potentials !!" << endl ;
};
/********************************************************************************************/

template <typename particle>
void kern_frenkel<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <>
void kern_frenkel<patchy_2D>::build_bonds_network( double bonds_cutoff ) {
  neighbour_list<patchy_2D> *walker = NULL ;
  patchy_2D *apparent_near = NULL ;
  apparent_near = new patchy_2D ;
  int ip , jp ;
  bool bond_flag = 1 ;
  int particles_in_list = verlet_list->particles_in_list ;
  vlist_element<patchy_2D> *element = verlet_list->element ;

  for( int i=0 ; i<particles_in_list ; i++ ) {
    if( element[i].ptr->bonds == NULL )  element[i].ptr->bonds = new patchy_2D* [ element[i].ptr->Npatches ] ;
    for( int ip=0 ; ip<element[i].ptr->Npatches ; ip++ ) element[i].ptr->bonds[ip] = NULL ;
  }
  for( int i=0 ; i<particles_in_list ; i++ ) {
    walker = element[i].neighbours ;
    while( walker != NULL ) {
      verlet_list->god->right_copy( walker->ptr , element[i].ptr , apparent_near ) ;
      apparent_near->Npatches = walker->ptr->Npatches ;
      apparent_near->patches = walker->ptr->patches ;
      apparent_near->rotation = walker->ptr->rotation ;

      Rij.position = apparent_near->position - element[i].ptr->position ;
      distance = Rij.position.norm() ;
      if( distance < cutoff ) {
	Rij.position /= distance ;
	Rji.position = Rij.position * (-1.0) ;
	bond_flag = 1 ;
	for( ip=0 ; ip<element[i].ptr->Npatches && bond_flag ; ip++ ) {
	  if( element[i].ptr->patch(ip) * Rij.position > COSDELTA ) {
	    for( jp=0 ; jp<apparent_near->Npatches && bond_flag ; jp++ ) {
	      if( apparent_near->patch(jp) * Rji.position > COSDELTA ) {
		bond_flag = 0 ;
		element[i].ptr->bonds[ip] = walker->ptr ;
		walker->ptr->bonds[jp] = element[i].ptr ;
	      }
	    }
	  }
	}
      }

      walker = walker->next ;
    }
  }

  apparent_near->Npatches = 0 ;
  apparent_near->patches = NULL ;
  delete apparent_near ;
};
/********************************************************************************************/

template <typename particle>
double * kern_frenkel<particle>::dump_parameters( void ) {
  double *params = new double [4] ;
  params[0] = CORE_RADIUS ;
  params[1] = DELTA ;
  params[2] = LAMBDA ;
  params[3] = EPSILON ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
