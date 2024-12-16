/* This library contains the definition of the Nosè Hoover thermostat classes, a simple simplectic deterministic thermostat */

class NoseHooverChain {
 public:
  double *chain_coordinate ;      // 0<j<=M
  double *chain_momenta ;
  double *chain_forces ;          // [1] = kin_energy_particles - 3NkT , [j>1] = (p_eta_(j-1))^2/Q_(j-1) - kT
  double *chain_friction ;
  double *scale_factor_by_friction ;         // friction[ 0<j<=M ] = p_eta_j/Q_j   ==>   scale_factor_by_friction[j] = exp(-delta_alpha*friction[j]/8)     (Tuckerman's notation)

  NoseHooverChain( int elements_per_chain = 0 ) ;
  ~NoseHooverChain() ;

  void single_NHchain_copy( const NoseHooverChain &other_chain , int elements_per_chain ) ;
};

class NoseHooverThermostat {
 public:
  int thermostats_per_chain ;               // M
  double *Q ;                      // thermostat "masses"
  int RESPA_NHC_integration_steps ;
  int Yoshida_parameters_number ;
  double *Yoshida_parameters ;
  int number_of_chains ;
  NoseHooverChain **chain ;

  NoseHooverThermostat() ;
  NoseHooverThermostat( const NoseHooverThermostat &other_thermostat ) ;
  ~NoseHooverThermostat() ;
};


/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
NoseHooverChain::NoseHooverChain( int elements_per_chain ) {
  if( elements_per_chain == 0 ) {
    chain_coordinate = NULL ;
    chain_momenta = NULL ;
    chain_forces = NULL ;
    chain_friction = NULL ;
    scale_factor_by_friction = NULL ;
  }else{
    chain_coordinate = (double *)calloc( elements_per_chain , sizeof(double) ) ;
    chain_momenta = (double *)calloc( elements_per_chain , sizeof(double) ) ;
    chain_forces = (double *)calloc( elements_per_chain , sizeof(double) ) ;
    chain_friction = (double *)calloc( elements_per_chain , sizeof(double) ) ;
    scale_factor_by_friction = (double *)calloc( elements_per_chain , sizeof(double) ) ;
    if( chain_coordinate == NULL  ||  chain_momenta == NULL  ||  chain_forces == NULL  ||  scale_factor_by_friction == NULL  ||  chain_friction == NULL  ) {
      cout << "\n Allocation error in NoseHooverChain constructor" ;
      exit( EXIT_FAILURE ) ;
    }
  }
}
/********************************************************************************************/

void NoseHooverChain::single_NHchain_copy( const NoseHooverChain &other_chain , int elements_per_chain ) {
  if( elements_per_chain > 0 ) {
    for(int i=0; i<elements_per_chain; i++) {
      chain_coordinate[i] = other_chain.chain_coordinate[i] ;
      chain_momenta[i] = other_chain.chain_momenta[i] ;
      chain_forces[i] = other_chain.chain_forces[i] ;
      chain_friction[i] = other_chain.chain_friction[i] ;
      scale_factor_by_friction[i] = other_chain.scale_factor_by_friction[i] ;
    }
  }else{
    chain_coordinate = NULL ;
    chain_momenta = NULL ;
    chain_forces = NULL ;
    chain_friction = NULL ;
    scale_factor_by_friction = NULL ;
  }
}
/********************************************************************************************/

NoseHooverChain::~NoseHooverChain() {
  if( chain_coordinate != NULL ) free( chain_coordinate ) ;
  chain_coordinate = NULL ;
  if( chain_momenta != NULL ) free( chain_momenta ) ;
  chain_momenta = NULL ;
  if( chain_forces != NULL ) free( chain_forces ) ;
  chain_forces = NULL ;
  if( chain_friction != NULL ) free( chain_friction ) ;
  chain_friction = NULL ;
  if( scale_factor_by_friction != NULL ) free( scale_factor_by_friction ) ;
  scale_factor_by_friction = NULL ;
}

/********************************************************************************************/

NoseHooverThermostat::NoseHooverThermostat() {
  thermostats_per_chain = 0 ;
  Yoshida_parameters_number = 0 ;
  RESPA_NHC_integration_steps = 0 ;
  number_of_chains = 0 ;
  Q = NULL ;
  Yoshida_parameters = NULL ;
  chain = NULL ;
}
/********************************************************************************************/

NoseHooverThermostat::NoseHooverThermostat( const NoseHooverThermostat &other_thermostat ) {
  thermostats_per_chain = other_thermostat.thermostats_per_chain ;
  Yoshida_parameters_number = other_thermostat.Yoshida_parameters_number ;
  RESPA_NHC_integration_steps = other_thermostat.RESPA_NHC_integration_steps ;
  number_of_chains = other_thermostat.number_of_chains ;
  if( thermostats_per_chain > 0 ) {
    Q = (double *)calloc( thermostats_per_chain , sizeof(double) ) ;
    if( Q == NULL ) {
      cout << "\n Allocation error in NoseHooverThermostat copy constructor" ;
      exit( EXIT_FAILURE ) ;
    }
    for( int i=0 ; i<thermostats_per_chain ; i++ ) {
      Q[i] = other_thermostat.Q[i] ;
    }
  }

  if( Yoshida_parameters_number > 0 ) {
    Yoshida_parameters = (double *)calloc( Yoshida_parameters_number , sizeof(double) ) ;
    if( Yoshida_parameters == NULL ) {
      cout << "\n Allocation error in NoseHooverThermostat copy constructor" ;
      exit( EXIT_FAILURE ) ;
    }
    for( int i=0 ; i<Yoshida_parameters_number ; i++ ) {
      Yoshida_parameters[i] = other_thermostat.Yoshida_parameters[i] ;
    }
  }

  if( number_of_chains > 0 ) {
    chain = (NoseHooverChain **)calloc( number_of_chains , sizeof(NoseHooverChain *) ) ;
    if( chain == NULL ) {
      cout << "\n Allocation error in NoseHooverThermostat copy constructor" ;
      exit( EXIT_FAILURE ) ;
    }
    for( int j=0 ; j<number_of_chains ; j++ ) {
      chain[j] = new NoseHooverChain( thermostats_per_chain ) ;
      chain[j]->single_NHchain_copy( *(other_thermostat.chain[j]) , thermostats_per_chain ) ;
    }
  }
}
/********************************************************************************************/

NoseHooverThermostat::~NoseHooverThermostat() {
  if( Yoshida_parameters_number > 0 ) free( Yoshida_parameters ) ;
  if( number_of_chains > 0 ) {
    for( int j=0 ; j<number_of_chains ; j++ ) {
      delete chain[j] ;
      chain[j] = NULL ;
    }
    free( chain ) ;
  }
  if( thermostats_per_chain > 0 ) free( Q ) ;
  Yoshida_parameters = NULL ;
  chain = NULL ;
  Q = NULL ;
}
