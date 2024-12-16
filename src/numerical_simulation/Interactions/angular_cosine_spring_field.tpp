/*****************************************************************************************/
/* This is an external field for which each particle's orientation is (periodically) bonded with a sine square spring to its initial value */

template <typename particle>
class angular_c_spring_field : public interaction<particle> { // energy in units of EPSILON
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::name ;
  using interaction<particle>::idx ;

 public:
  double lambda = 0 , alpha_theta = 0 ; // spring constants, lambda is the thermodynamic integration variable, lambda*alpha_theta is the prefactor of the potential
  int n_minima = 1 ;                    // number of potential minima in the interval [-pi,pi] (periodicity: 2*pi/n_minima)
  angular_c_spring_field( double Lambda = 0.0 , double Alpha_theta = 1.0 , int N_minima = 1 ) ;
 private:
  double PERIOD , HALF_PERIOD ;   // angular periodicity of the potential
  double *equilibrium_positions = NULL ;
  particle f ;
  double multiplier = 0 ;   // auxiliary variables
  inline double potential( particle *part , double theta0 ) ;
  inline particle force( particle *part , double theta0 ) ;

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
  angular_c_spring_field *copy( void ) ;
  ~angular_c_spring_field() ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
angular_c_spring_field<particle>::angular_c_spring_field( double Lambda , double Alpha_theta , int N_minima ) {
  lambda = Lambda ;
  alpha_theta = Alpha_theta ;
  n_minima = N_minima ;
  HALF_PERIOD = M_PI / n_minima ;
  PERIOD = 2 * HALF_PERIOD ;
  sprintf( name , "%s" , "ang_c_spring\0" ) ;
};
/********************************************************************************************/

template <typename particle>
angular_c_spring_field<particle>::~angular_c_spring_field( void ) {
  delete [] equilibrium_positions ;
  equilibrium_positions = NULL ;
  lambda = 0 ;
  alpha_theta = 1 ;
  n_minima = 1 ;
  HALF_PERIOD = M_PI / n_minima ;
  PERIOD = 2 * HALF_PERIOD ;
  sprintf( name , "%s" , "ang_c_spring\0" ) ;
};
/********************************************************************************************/

template <typename particle>
void angular_c_spring_field<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
  cout << "*** ERROR : ANGULAR-COSINE-SPRING potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
}

template <>
void angular_c_spring_field<patchy_2D>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<patchy_2D> *conf_ptr , int parts_number ) {
  cutoff = 0.0 ;   // CUTOFF, not used here
  SHIFT = 0.0 ;    // SHIFT , not used here
  sscanf( file_rec->word[1] , "%lf" , &lambda ) ;
  sscanf( file_rec->word[2] , "%lf" , &alpha_theta ) ;
  sscanf( file_rec->word[3] , "%d" , &n_minima ) ;
  HALF_PERIOD = M_PI / n_minima ;
  PERIOD = 2 * HALF_PERIOD ;
	              // non - generation of the verlet list
  verlet_list = new verlet<patchy_2D>( conf_ptr , 0 , 0.0 , 0.0 ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 0 ; // Indeed it has been already done in the constr
	              // initialization of the verlet list
  // initialization of the array of equilibrium positions
  // here we consider only the case in which every particle in the system feels this potential
  if( parts_number != conf_ptr->numberOfPart ) {
    cout << "*** ERROR: if you want to apply the cosine spring field potential you should change the code!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  equilibrium_positions = new double [conf_ptr->numberOfPart] ;
  for( int i=0; i<conf_ptr->numberOfPart; i++ ) equilibrium_positions[i] = conf_ptr->particles[i]->rotation ;
};
/********************************************************************************************/

template <typename particle>
inline double angular_c_spring_field<particle>::potential( particle *part , double theta0 ) {
  cout << "*** ERROR : ANGULAR-COSINE-SPRING potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
}

template <>
inline double angular_c_spring_field<patchy_2D>::potential( patchy_2D *part , double theta0 ) {
  /***********   ANGULAR COSINE EXTERNAL POTENTIAL     Vi = lambda * alpha_theta * sin( 0.5*n_minima*(rotation-theta0) )^2 ;    **********/
  static double sine ;
  sine = sin( 0.5 * n_minima * ( part->rotation - theta0 ) ) ;
  return lambda * alpha_theta * sine * sine  ;
};
/********************************************************************************************/

template <typename particle>
inline particle angular_c_spring_field<particle>::force( particle *part , double theta0 ) {
    /***********   ANGULAR COSINE POTENTIAL     **********/               // forces are in units of EPSILON/SIGMA
  f.position *= 0 ;
  cout << "*** ERROR: ANGULAR-COSINE-SPRING potential generates no force, but a torque!" << endl;
  exit( EXIT_FAILURE ) ;

  return f ;
};
/********************************************************************************************/

template <typename particle>
double angular_c_spring_field<particle>::return_pair_potential( double distance ) {
  static double sine ;
  sine = sin( 0.5 * n_minima * distance ) ;
  return lambda * alpha_theta * sine * sine  ;
};
/********************************************************************************************/

template <typename particle>
double angular_c_spring_field<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  // This is an external field, it doesn't generate any pair contribution, hence gives 0
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
particle angular_c_spring_field<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  // This is an external field, it doesn't generate any pair contribution, hence gives 0
  f.position *= 0.0 ;
  return f ;
};
/********************************************************************************************/

template <typename particle>
double angular_c_spring_field<particle>::compute_energy_contribution( void ) {
  double energy = 0 ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    energy += ( potential( particles[i] , equilibrium_positions[i] ) + SHIFT ) ;
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
inline double angular_c_spring_field<particle>::compute_energy_contribution( particle *part ) {
  // This function computes the energy of ONE particle only
  double energy = 0 ;

  particle **particles = verlet_list->god->particles ;
  if( part->id == particles[part->id]->id ) {
    // can release this check...
    energy = ( potential( particles[part->id] , equilibrium_positions[part->id] ) + SHIFT ) ;
  } else {
    int particles_in_list = verlet_list->god->numberOfPart ;
    bool keep_running = 1 ;
    for( int i=0 ; i<particles_in_list && keep_running ; i++ ) {
      if( particles[i] == part ) {
	keep_running = 0 ;
	energy = ( potential( particles[i] , equilibrium_positions[i] ) + SHIFT ) ;
      }
    }
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
double angular_c_spring_field<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For an angular potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double angular_c_spring_field<particle>::lambda_derivative_perpart( void ) {
  cout << "*** ERROR : ANGULAR-COSINE-SPRING potential can only be used with the particle class patchy_2D !!" << endl ;
  exit(EXIT_FAILURE) ;
}

template <>
double angular_c_spring_field<patchy_2D>::lambda_derivative_perpart( void ) {
  // This function computes the PER PARTICLE derivative of the external field potential as a function of lambda,
  double der = 0 , sine ;

  int particles_in_list = verlet_list->god->numberOfPart ;
  patchy_2D **particles = verlet_list->god->particles ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    sine = sin( 0.5 * n_minima * ( particles[i]->rotation - equilibrium_positions[i] ) ) ;
    der += ( sine * sine ) ;
  }

  return alpha_theta * der / particles_in_list ;
};
/********************************************************************************************/

template <typename particle>
void angular_c_spring_field<particle>::compute_forces_contribution(void) {
  /******* ANGULAR-COSINE-SPRING potential generates no force, but a torque  ****/
  particle upart ;
};
/********************************************************************************************/

template <typename particle>
double angular_c_spring_field<particle>::compute_virial_contribution(void) {
  /******* ANGULAR-COSINE-SPRING potential generates no force, but a torque  ****/
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
string angular_c_spring_field<particle>::ostr( void ) {
  ostringstream st ;
  st << "ANGULAR COSINE SPRING FIELD (spring constant, prefactor, number of minima) : ( LAMBDA=" << lambda << " ; ALPHA_THETA=" << alpha_theta << " ; N_MINIMA=" << n_minima << " )" << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
angular_c_spring_field<particle> *angular_c_spring_field<particle>::copy( void ) {
  angular_c_spring_field *copy = new angular_c_spring_field( lambda , alpha_theta , n_minima ) ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void angular_c_spring_field<particle>::compute_per_part_energy_contribution( double *envec ) {
  int particles_in_list = verlet_list->god->numberOfPart ;
  particle **particles = verlet_list->god->particles ;
  for( int i=0 ; i<particles_in_list ; i++ )  envec[i] += ( potential( particles[i] , equilibrium_positions[i] ) + SHIFT ) ;
};
/********************************************************************************************/

template <typename particle>
void angular_c_spring_field<particle>::compute_per_part_virial_contribution( double *virvec ) {
  /******* ANGULAR-COSINE-SPRING potential generates no force, but a torque  ****/
  particle upart ;
};
/********************************************************************************************/

template <typename particle>
void angular_c_spring_field<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * angular_c_spring_field<particle>::dump_parameters( void ) {
  double *params = new double [3] ;
  params[0] = lambda ;
  params[1] = alpha_theta ;
  params[2] = n_minima ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
