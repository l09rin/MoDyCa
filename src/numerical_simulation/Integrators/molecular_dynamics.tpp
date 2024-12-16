/* methods to integrate the system in time with molecular dynamics */
#include "../../system_comm.h"

/**************************   BLOCK 3 - MD integrators   ********************************/
/****************************************************************************************/

template <typename particle>
void configuration<particle>::MD( const char *questions ) {
  char control, command_output ;
  int current_step = 1 ;
  int starting_time = time(0) ;
  // Variables needed for NVT integrators
  double half_timestep = 0.5 * time_step , N3kT = 3.0 * ((double)numberOfPart) * k_T ;
  // Initialization of NVT NHequation variables
  NoseHoover_degreesOFfreedom = 3 * numberOfPart - 7 ;
  double gkT = ((double)NoseHoover_degreesOFfreedom) * k_T ;
  NoseHooverMass = gkT * 0.01 * 0.01 ;
  NoseHooverFriction = 0 ;

  if( strcmp( simulation_mode , "MD" ) != 0 ) {
    cout << "ERROR: The simulation is not in the Molecular Dynamics mode !" << endl ;
    fflush( 0 ) ;
    exit( EXIT_FAILURE ) ;
  }
  if( strcmp( simulation_ensemble , "NVT_NHE" ) == 0 ) {
    compute_kinetic_energy() ;
    calculate_molecules_CoM() ;
  }else if( strcmp( simulation_ensemble , "NVT_NHC" ) == 0 ) {
    compute_kinetic_energy() ;             // I need initial kinetic energy to start NHC integration
    initialize_NHC_thermostat() ;
  }else if( strcmp( simulation_ensemble , "BROWNIAN" ) == 0 ) {
    sqrt_2kT_over_xidt = sqrt( 2.0 * k_T / langevin_xi / time_step ) ;
    cout << "  Brownian integrator constants:" << endl ;
    cout << "  _ xi = " << langevin_xi << endl ;
    cout << "  _ sqrt(2*kT/(xi*dt)) = " << sqrt_2kT_over_xidt << endl << endl ;
  }else if( strcmp( simulation_ensemble , "LANGEVIN" ) == 0 ) {
    if( particle_mass == 0 ) {cout << "*** ERROR : In the LANGEVIN algorithm you need to specify a PARTICLE_MASS value!" << endl ; exit(EXIT_FAILURE) ; };
    langevin_C0 = exp( -langevin_xi * time_step / particle_mass ) ;
    langevin_C1 = ( 1. - langevin_C0 ) / ( langevin_xi * time_step / particle_mass ) ;
    langevin_C2 = ( 1. - langevin_C1 ) / ( langevin_xi * time_step / particle_mass ) ;
    sigma_dr = sqrt(  k_T * time_step / langevin_xi * ( 2.0 - ( 3.0 - 4.0 * langevin_C0 + langevin_C0 * langevin_C0 ) / ( langevin_xi * time_step / particle_mass ) )  ) ;
    sigma_dv = sqrt(  k_T / particle_mass * ( 1.0 - langevin_C0 * langevin_C0 )  ) ;
    correlation_dv_dr = k_T / langevin_xi * ( 1.0 - langevin_C0 ) * ( 1.0 - langevin_C0 ) / sigma_dr / sigma_dv ;
    uncorrelation_dv_dr = sqrt( 1.0 - correlation_dv_dr*correlation_dv_dr ) ;
    cout << "  Langevin integrator constants:" << endl ;
    cout << "  _ xi = " << langevin_xi << endl ;
    cout << "  _ C0 = " << langevin_C0 << endl ;
    cout << "  _ C1 = " << langevin_C1 << endl ;
    cout << "  _ C2 = " << langevin_C2 << endl ;
    cout << "  _ sqrt(<(dr)^2>) = " << sigma_dr << endl ;
    cout << "  _ sqrt(<(dv)^2>) = " << sigma_dv << endl ;
    cout << "  _ <dr*dv>/sqrt(<(dv)^2>)/sqrt(<(dr)^2>) = " << correlation_dv_dr << endl << endl ;
  }
  set_saving_options( "EQUILIBRATION" ) ;
  save_simulation_conditions( "MD_simulation_conditions.dat" ) ;

  /*************                  equilibration stage                  *************/
  if( EQUILIBRATION_STEPS != 0 ) {
    if( SAVE.backup_interval == 0 ) SAVE.backup_interval = EQUILIBRATION_STEPS + 10 ;
    clear_runtime_calculations() ;
                 // generation of the sub-directory containing the equilibration data
    char ROOT[100] ;
    if( generate_equilibration_subdirectory( ROOT ) == NULL ) cout << " Impossible to open equilibration sub-directory\n";
    initialize_storage_file( -1 , seed ) ;
    //    save_data(0) ;
    //    save_configuration(0) ;
    // molecules positions update
    if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)

    printf( "\n  EQUILIBRATION STAGE . . .\n" ) ; fflush(0) ;
    int equilibration_steps = EQUILIBRATION_STEPS ;
    int total_step_div_10 = floor( equilibration_steps/10. ) ;
    if( total_step_div_10 == 0 ) total_step_div_10 = 1 ;
    int input_int = 0 ;
    check_action( 0 , starting_time , total_step_div_10 ) ;
    do{
      lists_rebuilds = 0 ;  // number of events of recalculation of the verlet lists
      equilibration_steps += input_int ;              //  if system needs another equilibration run I can specify n

      if( strcmp( simulation_ensemble , "NVE" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  VelocityVerletIntegrationStep() ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}
      }else if( strcmp( simulation_ensemble , "NVT_NHE" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  NoseHooverEquation_direct_integration() ;
	  if( current_step % vcm_reset_interval == 0 ) NHE_cm_reset() ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}
      }else if( strcmp( simulation_ensemble , "NVT_NHC" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  NHC_integration_step( half_timestep ) ;
	  VelocityVerletIntegrationStep() ;
	  KINETIC_ENERGY.time = -1 ;  // This is needed in order to force the function compute_kinetic_energy() to compute the value, otherwise it will believe it is up to date
	  // NHC needs of K to compute first thermostat's coupling force. Global thermostat case
	  NHC_thermostat->chain[0]->chain_forces[0] = 2.0 * compute_kinetic_energy() - N3kT ;
	  NHC_integration_step( half_timestep ) ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}
      }else if( strcmp( simulation_ensemble , "BROWNIAN" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  BROWN_integration_step( time_step ) ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}
      }else if( strcmp( simulation_ensemble , "LANGEVIN" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  //	  print_all_interactions() ;
	  LANGEVIN_integration_step( time_step ) ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}
      }

      cout << "\tNumber of events of verlet lists recalculation : " << lists_rebuilds << endl ;
      input_int = 0 ;
      if( strcmp( questions, "NO_STOP" ) != 0 ) {
	control = request( "\n Do you want to continue the equilibration stage?" , "y/n" , '/' ) ; // if system needs more time to equilibrate
	if( control == 'y' ) {
	  printf( "  Insert the number of further equilibration steps: _" ) ;
	  command_output = scanf( "%d%*c" , &input_int ) ;
	  if (command_output == 0) cout << "*** Warning: input steps: " << input_int << endl ;
	}
      }

    } while( input_int != 0 ) ;
    if( CHECKS.bonds_lists_number > 0 ) cout << endl << "  Maximum bonding distance after equilibration:  " << compute_max_bond_distance() << endl ;
    close_storage_file() ;
    if( strcpy( _DIRECTORY_NAME , ROOT ) == NULL ) cout << " Impossible to open simulation directory" ;
  }

  /*************                  production stage                  *************/
  starting_time = time(0) ;
  set_saving_options( "PRODUCTION" ) ;
  if( PRODUCTION_STEPS != 0 ) {
    int production_steps = PRODUCTION_STEPS ;
    int total_step_div_10 = floor( production_steps/10. ) ;
    if( total_step_div_10 == 0 ) total_step_div_10 = 1 ;
    if( strcmp( questions, "NO_STOP" ) != 0 ) {
      control = request( "\n Do you want to change the number of production steps?" , "y/n" , '/' ) ;
      if( control == 'y' ) {
	printf( "  Insert the number of production steps: _" ) ;
	command_output = scanf( "%d%*c" , &production_steps ) ;
	if (command_output == 0) cout << "*** Warning: input steps: " << production_steps << endl ;
      }
    }

    printf( "\n  PRODUCTION STAGE . . .\n" ) ; fflush(0) ;
    lists_rebuilds = 0 ; // number of events of recalculation of the verlet lists
    initialize_storage_file( production_steps , seed ) ;     // open the file in which configurations are periodically saved
    //    save_data(0) ;
    if( SAVE.backup_interval == 0 ) SAVE.backup_interval = PRODUCTION_STEPS + 10 ;
    clear_runtime_calculations() ;
    // molecules positions update
    if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)

    check_action( 0 , starting_time , total_step_div_10 ) ;
    if( strcmp( simulation_ensemble , "NVE" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	VelocityVerletIntegrationStep() ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }else if( strcmp( simulation_ensemble , "NVT_NHE" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	NoseHooverEquation_direct_integration() ;
	if( current_step % vcm_reset_interval == 0 ) NHE_cm_reset() ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }else if( strcmp( simulation_ensemble , "NVT_NHC" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	NHC_integration_step( half_timestep ) ;
	VelocityVerletIntegrationStep() ;
	// NHC needs of K to compute first thermostat's coupling force. Global thermostat case
	KINETIC_ENERGY.time = -1 ;  // This is needed in order to force the function compute_kinetic_energy() to compute the value, otherwise it will believe it is up to date
	NHC_thermostat->chain[0]->chain_forces[0] = 2.0 * compute_kinetic_energy() - N3kT ;
	NHC_integration_step( half_timestep ) ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }else if( strcmp( simulation_ensemble , "BROWNIAN" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	BROWN_integration_step( time_step ) ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }else if( strcmp( simulation_ensemble , "LANGEVIN" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	LANGEVIN_integration_step( time_step ) ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }
  
    cout << "\tNumber of events of verlet lists recalculation : " << lists_rebuilds << endl ;
    save_runtime_calculations( current_step ) ;
    if( CHECKS.bonds_lists_number > 0 ) cout << endl << "  Maximum bonding distance in last configuration:  " << compute_max_bond_distance() << endl ;

    cout << endl ;
    fflush( NULL ) ;
    close_storage_file() ;
  }
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::VelocityVerletIntegrationStep(void) {                        // SI PUO' PARALLELIZZARE?
  // positions update
  for( int i=0 ; i<numberOfPart ; i++ ) {
    particles[i]->position += ( ( velocities[i]->position*time_step ) + ( forces_t1[i]->position*(0.5*time_step*time_step) ) ) ;
    // pbc( particles[i] ) ;                   // PBC // has been put in the function updating verlet lists
  }
  // molecules positions update
  if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)
  check_verlet_update() ;
  // velocities update
  calculate_forces_t2() ;
  for( int i=0 ; i<numberOfPart ; i++ ) {
    velocities[i]->position += ( ( forces_t2[i]->position + forces_t1[i]->position ) * (0.5*time_step) ) ;
  }
  support_pp2particle = forces_t1 ;
  forces_t1 = forces_t2 ;
  forces_t2 = support_pp2particle ;
  global_time += time_step ;
}
/****************************************************************************************/

/**************************   BLOCK 4 - MD-NVT   *************************************/
/****************************************************************************************/

template <typename particle>
void configuration<particle>::initialize_NHC_thermostat() {
  double tau = 20.0 ;
  NHC_thermostat = new NoseHooverThermostat() ;
  NHC_thermostat->thermostats_per_chain = 7 ;          // Instert a number of thermostats_per_chain > 1 in order to avoid errors in the integration algorithm!
  NHC_thermostat->RESPA_NHC_integration_steps = 10 ;
  ifstream config_file ;                               // the function tries to load parameters from an input file
  config_file.open( "NHC_input_file.dat" );
  if( config_file.fail() ) {
    cout << "  Default parameters will be used for the NHC thermostat\n" ;
  }else{
    record *file_rec = NULL ;
    file_rec = new record ;

    while( file_rec->getrecord( config_file ) ) {
      file_rec->split() ;
      if( file_rec->words_number > 0 ) {
	if( file_rec->word[0][0] != '#' ) {
	  if( strcmp( file_rec->word[0], "TAU" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &tau ) ;
	  if( strcmp( file_rec->word[0], "THERMOSTATS_PER_CHAIN" ) == 0 ) {
	    sscanf( file_rec->word[1] , "%d" , &(NHC_thermostat->thermostats_per_chain) ) ;
	    if( NHC_thermostat->thermostats_per_chain < 2 ) {
	      cout << "  The number of thermostats per chain has to be more then 1\n" ;
	      exit( EXIT_FAILURE ) ;
	    }
	  }
	  if( strcmp( file_rec->word[0], "NUMBER_RESPA_STEPS" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &(NHC_thermostat->RESPA_NHC_integration_steps) ) ;
	}
      }
      delete file_rec ;
      file_rec = new record ;
    }
    delete file_rec ;
    file_rec = NULL ;
    config_file.close() ;
  }

  NHC_thermostat->Q = (double *)calloc( NHC_thermostat->thermostats_per_chain, sizeof(double) ) ;
  NHC_thermostat->Q[0] = 3.0 * ((double)numberOfPart) * k_T * tau*tau * time_step*time_step ;     // For the global thermostat
  for( int i=1; i<NHC_thermostat->thermostats_per_chain; i++ ) NHC_thermostat->Q[i] = k_T * tau*tau * time_step*time_step ;

  // In this case I implement the factorization of the NHC propagator of the 4th order in the Yoshida scheme (Tuckerman, p. 194)
  NHC_thermostat->Yoshida_parameters_number = 3 ;
  NHC_thermostat->Yoshida_parameters = (double *)calloc( NHC_thermostat->Yoshida_parameters_number, sizeof(double) ) ;
  NHC_thermostat->Yoshida_parameters[0] = 1.0 / ( 2.0 - pow(2.0, 1.0/3.0) ) ;
  NHC_thermostat->Yoshida_parameters[1] = 1.0 - 2.0 * NHC_thermostat->Yoshida_parameters[0] ;
  NHC_thermostat->Yoshida_parameters[2] = NHC_thermostat->Yoshida_parameters[0] ;

  // In this case I use a global thermostat
  NHC_thermostat->number_of_chains = 1 ;
  NHC_thermostat->chain = (NoseHooverChain **)calloc( NHC_thermostat->number_of_chains, sizeof(NoseHooverChain *) ) ;
  for(int i=0; i<NHC_thermostat->number_of_chains; i++) NHC_thermostat->chain[i] = new NoseHooverChain( NHC_thermostat->thermostats_per_chain ) ;

  double aux1 ;
  double aux2 ;
  for( int i=0; i<NHC_thermostat->number_of_chains; i++ ) {
    for( int j=0; j<NHC_thermostat->thermostats_per_chain; j++ ) {
  // Initialization of thermostats' coordinates and momenta to random gaussian values
      aux1 = drand48() ;
      aux2 = drand48() ;
      NHC_thermostat->chain[i]->chain_coordinate[j] = - sqrt(-2.0*log(aux1)) * cos(2.0*M_PI*aux2) ;
      NHC_thermostat->chain[i]->chain_momenta[j] = - sqrt(-2.0*log(aux1)) * sin(2.0*M_PI*aux2) ;
      // calculation of thermostats' forces. NHC integration doesn't need of calculation of frictions to start
      if( j == 0 ) {
	NHC_thermostat->chain[i]->chain_forces[j] = 2.0 * compute_kinetic_energy() - 3.0 * ((double)numberOfPart) * k_T ;
      }else{
	NHC_thermostat->chain[i]->chain_forces[j] = NHC_thermostat->chain[i]->chain_momenta[j-1] * NHC_thermostat->chain[i]->chain_momenta[j-1] / NHC_thermostat->Q[j-1] - k_T ;
      }
    }
  }
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::NHC_integration_step(double time_delta) {
  double RESPA_Yoshida_time_step ;
  for( int i=0 ; i<NHC_thermostat->Yoshida_parameters_number ; i++ ) {
    RESPA_Yoshida_time_step = time_delta * NHC_thermostat->Yoshida_parameters[i] / NHC_thermostat->RESPA_NHC_integration_steps ;
    for( int j=0 ; j<NHC_thermostat->RESPA_NHC_integration_steps ; j++ ) {
      Factorized_NHC_evolutionOperator( RESPA_Yoshida_time_step ) ;
    }
  }
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::Factorized_NHC_evolutionOperator(double time_delta) {
  double half_time_delta = 0.5 * time_delta ;                  // delta_alpha/4 in the notation used by Tuckerman
  int thermostats_per_chain = NHC_thermostat->thermostats_per_chain ;
  double *Q = NHC_thermostat->Q ;

  double *chain_momenta = NHC_thermostat->chain[0]->chain_momenta ;
  double *chain_forces = NHC_thermostat->chain[0]->chain_forces ;
  double *chain_friction = NHC_thermostat->chain[0]->chain_friction ;
  double *scale_factor_by_friction = NHC_thermostat->chain[0]->scale_factor_by_friction ;

  // Global thermostat
  // thermostats' momenta update for the first time
  chain_momenta[ thermostats_per_chain-1 ] +=   half_time_delta * chain_forces[ thermostats_per_chain-1 ] ;               // p_eta_M increment
  chain_friction[ thermostats_per_chain-1 ] =   chain_momenta[ thermostats_per_chain-1 ] / Q[ thermostats_per_chain-1 ] ;     // friction_M
  scale_factor_by_friction[ thermostats_per_chain-1 ] =   exp( -half_time_delta * chain_friction[ thermostats_per_chain-1 ] / 2.0 ) ;   // scale_factor_M

  for( int i=thermostats_per_chain-1 ; i>1 ; i-- ) {                       // p_eta_j update ( j=i-1 in Tuckerman's notation ) , j>1
    chain_momenta[i-1] =    scale_factor_by_friction[i] * scale_factor_by_friction[i] * chain_momenta[i-1] +
      half_time_delta * scale_factor_by_friction[i] * chain_forces[i-1] ;
    chain_friction[i-1] =   chain_momenta[i-1] / Q[i-1] ;                  // friction_j update (i-1=j in Tuckerman's notation)
    scale_factor_by_friction[i-1] =   exp( -half_time_delta * chain_friction[i-1] / 2.0 ) ;            // scale factor j
  }
  // p_eta_1 update
  chain_momenta[0] =   scale_factor_by_friction[1] * scale_factor_by_friction[1] * chain_momenta[0] +
    half_time_delta * scale_factor_by_friction[1] * chain_forces[0] ;
  chain_friction[0] =   chain_momenta[0] / Q[0] ;                  // friction_j update (i-1=j in Tuckerman's notation)
  scale_factor_by_friction[0] =   exp( -time_delta * chain_friction[0] ) ;  // scale factor 1. Definition is different from other scale factors, this as to multiply particles' velocity
  
  // coordinates update
  for( int i=0 ; i<thermostats_per_chain ; i++ ) NHC_thermostat->chain[0]->chain_coordinate[i] -=   time_delta * chain_friction[i] ;

  // particles' momenta rescaling
  for( int i=0 ; i<numberOfPart ; i++ ) velocities[i]->position *=   scale_factor_by_friction[0] ;
  // thermostats' momenta update for the second time
  // first thermostat force update
  KINETIC_ENERGY.time = -1 ;  // This is needed in order to force the function compute_kinetic_energy() to compute the value, otherwise it will believe it is up to date
  chain_forces[0] =   2.0 * compute_kinetic_energy() - 3.0 * ((double)numberOfPart) * k_T ;         // p_eta_1 increment
  chain_momenta[0] =   scale_factor_by_friction[1] * scale_factor_by_friction[1] * chain_momenta[0] +
    half_time_delta * scale_factor_by_friction[1] * chain_forces[0] ;
  for( int i=1 ; i<thermostats_per_chain-1 ; i++ ) {
    chain_forces[i] =   chain_momenta[i-1] * chain_momenta[i-1] / Q[i-1] - k_T ;              // p_eta_j increment (j<M)
    chain_momenta[i] =   scale_factor_by_friction[i+1] * scale_factor_by_friction[i+1] * chain_momenta[i] +
      half_time_delta * scale_factor_by_friction[i+1] * chain_forces[i] ;
  }
  chain_forces[ thermostats_per_chain-1 ] =
    chain_momenta[ thermostats_per_chain-2 ] * chain_momenta[ thermostats_per_chain-2 ] / Q[ thermostats_per_chain-2 ] - k_T ;                      // p_eta_M increment
  chain_momenta[ thermostats_per_chain-1 ] +=   half_time_delta * chain_forces[ thermostats_per_chain-1 ] ;
  
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::NoseHooverEquation_direct_integration( void ) {
          // positions update
  for( int i=0 ; i<numberOfPart ; i++ ) {
    particles[i]->position += ( ( velocities[i]->position * time_step ) + ( ( forces_t1[i]->position - ( velocities[i]->position * NoseHooverFriction ) ) * (0.5*time_step*time_step) ) ) ;
    // pbc( particles[i] ) ;                   // PBC // now in function updating verlet lists
  }
  // molecules positions update
  if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)
  check_verlet_update() ;
	  // velocities update
  calculate_forces_t2() ;
  for(int i=0; i<numberOfPart; i++) {
    velocities[i]->position += ( ( forces_t1[i]->position - ( velocities[i]->position * NoseHooverFriction ) ) * (0.5*time_step) ) ;
  }
  KINETIC_ENERGY.time = -1 ;  // This is needed in order to force the function compute_kinetic_energy() to compute the value, otherwise it will believe it is up to date
  NoseHooverFriction += ( time_step * ( 2.0 * compute_kinetic_energy() - ((double)NoseHoover_degreesOFfreedom)*k_T ) / NoseHooverMass ) ;
  for(int i=0; i<numberOfPart; i++) {
    velocities[i]->position = ( ( velocities[i]->position + ( forces_t2[i]->position * (0.5*time_step) ) ) / ( 1.0 + 0.5*time_step*NoseHooverFriction ) ) ;
  }
  KINETIC_ENERGY.time = -1 ;  // This is needed in order to force the function compute_kinetic_energy() to compute the value, otherwise it will believe it is up to date
  support_pp2particle = forces_t1;
  forces_t1 = forces_t2;
  forces_t2 = support_pp2particle;
  global_time += time_step ;
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::NHE_cm_reset( void ) {
  CoM_velocity_reset() ;
}
/*****************************************************************************************/

/**************************   BLOCK 5 - BROWNIAN DYNAMICS   *****************************/
/************ Eqn : v(t) = (xi^-1)*F(x,t) + sqrt(2kT/m^2/xi)*gauss(t;0,3) ***************/
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::BROWN_integration_step( double time_delta ) {
  calculate_forces_t2() ;
  for( int i=0 ; i<numberOfPart ; i++ ) {
    // velocities update
    brown_dv_gauss.Gaussian_particle() ;
    velocities[i]->position = ( ( brown_dv_gauss.position * sqrt_2kT_over_xidt ) + ( forces_t2[i]->position / langevin_xi ) ) ;
    // positions update
    particles[i]->position += ( velocities[i]->position * time_step ) ;
    // pbc( particles[i] ) ;                   // PBC // Now in function updating verlet lists
  }
  // molecules positions update
  if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)
  check_verlet_update() ;
  global_time += time_step ;
}
/****************************************************************************************/

/**************************   BLOCK 6 - LANGEVIN DYNAMICS   *****************************/
/************ Eqn : m*a(t) = -m*xi*v(t) + F(x,t) + sqrt(2*m*kT*xi)*R(t) *******************/
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::LANGEVIN_integration_step( double time_delta ) {
  // positions update
  for( int i=0 ; i<numberOfPart ; i++ ) {
    dr_gauss[i]->Gaussian_particle() ;
    dv_gauss[i]->Gaussian_particle() ;
    particles[i]->position += ( ( velocities[i]->position * (time_step*langevin_C1) ) + ( forces_t1[i]->position * (langevin_C2*time_step*time_step) ) + ( dr_gauss[i]->position * sigma_dr ) ) ;
    // pbc( particles[i] ) ;                   // PBC // now in function updating verlet lists
  }
  // molecules positions update
  if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)
  calculate_forces_t2() ;
    // velocities update
  for( int i=0 ; i<numberOfPart ; i++ ) {
    velocities[i]->position = ( ( velocities[i]->position * langevin_C0 ) + ( forces_t1[i]->position * (time_step*(langevin_C1-langevin_C2)) ) + ( forces_t2[i]->position * (time_step*langevin_C2) ) + ( ((dr_gauss[i]->position*correlation_dv_dr)+(dv_gauss[i]->position*uncorrelation_dv_dr))*sigma_dv ) ) ;
  }
  check_verlet_update() ;
  support_pp2particle = forces_t1 ;
  forces_t1 = forces_t2 ;
  forces_t2 = support_pp2particle ;
  global_time += time_step ;
}
/*****************************************************************************************/
