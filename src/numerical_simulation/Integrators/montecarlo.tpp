/* methods to integrate the system through MonteCarlo algorithms */
#include "../../system_comm.h"

/**************************   BLOCK 3 - MonteCarlo integrators   ************************/
/****************************************************************************************/

template <typename particle>
void configuration<particle>::MonteCarlo( const char *questions ) {
  char control, command_output ;
  int current_step = 1 ;
  int starting_time = time(0) ;
  particle old_position( *(particles[0]) ) , displacement ;  // needed if using patchy particles, you cannot define as static in the integrator function before knowing Npatches
  double *old_coordinates = new double[ numberOfPart ] ;

  if( strcmp( simulation_mode , "MC" ) != 0 ) {
    cout << "ERROR: The simulation is not in the Monte Carlo mode !" << endl ;
    fflush( 0 ) ;
    exit( EXIT_FAILURE ) ;
  }

  set_saving_options( "EQUILIBRATION" ) ;
  save_simulation_conditions( "MC_simulation_conditions.dat" ) ;
  // Check for right preservation of verlet lists, in case of long translational moves
  for( int k=0 ; k<INTERACTIONS.number ; k++ )
    if( INTERACTIONS.type[k]->verlet_list->must_be_updated  &&
	( INTERACTIONS.type[k]->verlet_list->delta_verlet <= sqrt( MC_maxdisplacement * MC_maxdisplacement * 3 ) ||
	  INTERACTIONS.type[k]->verlet_list->delta_verlet <= sqrt( MC_PARAMS_MAX.MC_maxdisplacement * MC_PARAMS_MAX.MC_maxdisplacement * 3 ) ) ) {
      cout << "*** ERROR: the maximum translation distance in the MC algorithm (" << sqrt( MC_maxdisplacement * MC_maxdisplacement * 3 ) <<
	" , " << sqrt( MC_PARAMS_MAX.MC_maxdisplacement * MC_PARAMS_MAX.MC_maxdisplacement * 3 ) <<
	") could be too large compared with the verlet skin (" << INTERACTIONS.type[k]->verlet_list->delta_verlet <<
	") for potential " << k << endl ;
      cout << "\t This would lead to a bad calculation of the verlet lists!" << endl ;
      exit( EXIT_FAILURE ) ;
    }

  compute_potential_energy() ;
  // Finite potential energy check, useful for hard core potentials
  if ( POTENTIAL_ENERGY.val == INFINITY ) {
    cout << "*** ERROR: The potential energy of the starting configuration is not bounded!" << endl <<
      "\t If you are using hard core potentials check the packing fraction." << endl ;
    exit( EXIT_FAILURE ) ;
  }
  if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
    bool same_mass = 1 ;
    for( int i=1; i<numberOfPart; i++ ) same_mass &= ( particles[i]->mass == particles[0]->mass ) ;
    if( ! same_mass ) {
      cout << "*** ERROR : This version of MC NVT integration with the constrained centre of mass cannot be used with systems of particles with different mass !" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    compute_centre_of_mass() ;
    CoM_displacement.clear() ;
    deltaEnergyOfCoM = 0.0 ;
  }

  /*************                  equilibration stage                  *************/
  if( EQUILIBRATION_STEPS != 0 ) {
    if( SAVE.backup_interval == 0 ) SAVE.backup_interval = EQUILIBRATION_STEPS + 10 ;
    clear_runtime_calculations() ;
                 // generation of the sub-directory containing the equilibration data
    char ROOT[100] ;
    if( generate_equilibration_subdirectory( ROOT ) == NULL ) cout << " Impossible to open equilibration sub-directory\n";
    initialize_storage_file( -1 , seed ) ;
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
      MonteCarlo_stats.clear() ;
      equilibration_steps += input_int ;              //  if system needs another equilibration run I can specify n

      if( strcmp( simulation_ensemble , "MC_SWAP_CHARGE" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  MonteCarloSwapChargeStep() ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}

      } else if( strcmp( simulation_ensemble , "MC_NVT" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  MonteCarlo_NVT( old_position ) ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}

      } else if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  MonteCarlo_NVT_CMfixed( old_position ) ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}

      } else if( strcmp( simulation_ensemble , "MC_NPT" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  MonteCarlo_NPT( old_position , old_coordinates ) ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}
      }

      input_int = 0 ;
      if( strcmp( questions, "NO_STOP" ) != 0 ) {
	control = request( "\n Do you want to continue the equilibration stage?" , "y/n" , '/' ) ; // if system needs more time to equilibrate
	if( control == 'y' ) {
	  printf( "  Insert the number of further equilibration steps: _" ) ;
	  command_output = scanf( "%d%*c" , &input_int ) ;
	  if (command_output == 0) cout << "*** ERROR: input value " << input_int << endl ;
	}
      } else {
	double min_acceptance_pos = 0.25 , max_acceptance_pos = 0.35 ;
	double min_acceptance_rot = 0.25 , max_acceptance_rot = 0.35 ;
	double min_acceptance_vol = 0.15 , max_acceptance_vol = 0.25 ;
	if( MonteCarlo_stats.pos_attempted != 0 ) {
	  double acceptance_rate = (double)MonteCarlo_stats.pos_accepted/MonteCarlo_stats.pos_attempted ;
	  if( acceptance_rate > max_acceptance_pos  &&  MC_maxdisplacement < MC_PARAMS_MAX.MC_maxdisplacement ) {
	    MC_maxdisplacement *= 1.1 ;
	    if( MC_maxdisplacement > MC_PARAMS_MAX.MC_maxdisplacement ) MC_maxdisplacement = MC_PARAMS_MAX.MC_maxdisplacement ;
	    input_int = EQUILIBRATION_STEPS ;
	  } else if( acceptance_rate < min_acceptance_pos  &&  MC_maxdisplacement > MC_PARAMS_MIN.MC_maxdisplacement ) {
	    MC_maxdisplacement *= 0.9 ;
	    if( MC_maxdisplacement < MC_PARAMS_MIN.MC_maxdisplacement ) MC_maxdisplacement = MC_PARAMS_MIN.MC_maxdisplacement ;
	    input_int = EQUILIBRATION_STEPS ;
	  }
	  if( input_int > 0 ) printf( "\n\tAcceptance rate for coordinates changes = %lf%% [%lf:%lf]\n\t\t(max displacement : %lf [%lf:%lf])" ,
				      (double)MonteCarlo_stats.pos_accepted/MonteCarlo_stats.pos_attempted*100. ,
				      min_acceptance_pos*100. , max_acceptance_pos*100. ,
				      MC_maxdisplacement , MC_PARAMS_MIN.MC_maxdisplacement , MC_PARAMS_MAX.MC_maxdisplacement ) ;
	}
	if( MonteCarlo_stats.vol_attempted != 0 ) {
	  double acceptance_rate = (double)MonteCarlo_stats.vol_accepted/MonteCarlo_stats.vol_attempted ;
	  if( acceptance_rate > max_acceptance_vol  &&  MC_maxboxdisp < MC_PARAMS_MAX.MC_maxboxdisp ) {
	    MC_maxboxdisp *= 1.1 ;
	    if( MC_maxboxdisp > MC_PARAMS_MAX.MC_maxboxdisp ) MC_maxboxdisp = MC_PARAMS_MAX.MC_maxboxdisp ;
	    input_int = EQUILIBRATION_STEPS ;
	  } else if( acceptance_rate < min_acceptance_vol  &&  MC_maxboxdisp > MC_PARAMS_MIN.MC_maxboxdisp ) {
	    MC_maxboxdisp *= 0.9 ;
	    if( MC_maxboxdisp < MC_PARAMS_MIN.MC_maxboxdisp ) MC_maxboxdisp = MC_PARAMS_MIN.MC_maxboxdisp ;
	    input_int = EQUILIBRATION_STEPS ;
	  }
	  if( input_int > 0 ) printf( "\n\tAcceptance rate for box sides changes = %lf%% [%lf:%lf]\n\t\t(max side variation : %lf [%lf:%lf])" ,
				      (double)MonteCarlo_stats.vol_accepted/MonteCarlo_stats.vol_attempted*100. ,
				      min_acceptance_vol*100. , max_acceptance_vol*100. ,
				      MC_maxboxdisp , MC_PARAMS_MIN.MC_maxboxdisp , MC_PARAMS_MAX.MC_maxboxdisp ) ;
	}
	if( MonteCarlo_stats.rot_attempted != 0 ) {
	  double acceptance_rate = (double)MonteCarlo_stats.rot_accepted/MonteCarlo_stats.rot_attempted ;
	  if( acceptance_rate > max_acceptance_rot  &&  MC_maxrotation < MC_PARAMS_MAX.MC_maxrotation ) {
	    MC_maxrotation *= 1.1 ;
	    if( MC_maxrotation > MC_PARAMS_MAX.MC_maxrotation ) MC_maxrotation = MC_PARAMS_MAX.MC_maxrotation ;
	    input_int = EQUILIBRATION_STEPS ;
	  } else if( acceptance_rate < min_acceptance_rot  &&  MC_maxrotation > MC_PARAMS_MIN.MC_maxrotation ) {
	    MC_maxrotation *= 0.9 ;
	    if( MC_maxrotation < MC_PARAMS_MIN.MC_maxrotation ) MC_maxrotation = MC_PARAMS_MIN.MC_maxrotation ;
	    input_int = EQUILIBRATION_STEPS ;
	  }
	  if( input_int > 0 ) printf( "\n\tAcceptance rate for orientation changes = %lf%% [%lf:%lf]\n\t\t(max rotation : %lf [%lf:%lf])" ,
				      (double)MonteCarlo_stats.rot_accepted/MonteCarlo_stats.rot_attempted*100. ,
				      min_acceptance_rot*100. , max_acceptance_rot*100. ,
				      MC_maxrotation , MC_PARAMS_MIN.MC_maxrotation , MC_PARAMS_MAX.MC_maxrotation ) ;
	}
	cout << endl << "\tNumber of events of verlet lists recalculation : " << lists_rebuilds << endl ;
	cout << endl ;
      }

    } while( input_int != 0 ) ;
    cout << endl ;
    if( CHECKS.bonds_lists_number > 0 ) cout << endl << "  Maximum bonding distance in last configuration:  " << compute_max_bond_distance() << endl ;
    cout << "## Attempted changes:\tN\tpos\tvol\trot" << endl ;
    cout << "\t\t" << MonteCarlo_stats.NoP_attempted << "\t" << MonteCarlo_stats.pos_attempted << "\t" << MonteCarlo_stats.vol_attempted << "\t" << MonteCarlo_stats.rot_attempted << endl ;
    if( MonteCarlo_stats.NoP_attempted != 0 ) printf( "\tAcceptance rate for changes in the number of particles = %lf%%\n" ,
						     (double)MonteCarlo_stats.NoP_accepted/MonteCarlo_stats.NoP_attempted*100. ) ;
    if( MonteCarlo_stats.pos_attempted != 0 ) printf( "\tAcceptance rate for coordinates changes = %lf%%\n\t\t(max displacement : %lf)\n" ,
						     (double)MonteCarlo_stats.pos_accepted/MonteCarlo_stats.pos_attempted*100. , MC_maxdisplacement ) ;
    if( MonteCarlo_stats.vol_attempted != 0 ) printf( "\tAcceptance rate for volume changes = %lf%%\n\t\t(max box side change : %lf)\n" ,
						     (double)MonteCarlo_stats.vol_accepted/MonteCarlo_stats.vol_attempted*100. , MC_maxboxdisp ) ;
    if( MonteCarlo_stats.chr_swap_attempted != 0 ) printf( "\tAcceptance rate for charge swap moves = %lf%%\n" ,
							  (double)MonteCarlo_stats.chr_swap_accepted/MonteCarlo_stats.chr_swap_attempted*100. ) ;
    if( MonteCarlo_stats.rot_attempted != 0 ) printf( "\tAcceptance rate for reorientation moves = %lf%%\n\t\t(max rotation : %lf)\n" ,
						     (double)MonteCarlo_stats.rot_accepted/MonteCarlo_stats.rot_attempted*100. , MC_maxrotation ) ;
    cout << "\tNumber of events of verlet lists recalculation : " << lists_rebuilds << endl ;
    cout << endl ;
    close_storage_file() ;
    if( strcpy( _DIRECTORY_NAME , ROOT ) == NULL ) cout << " Impossible to open the simulation directory" ;
  }

  /*************                  production stage                  *************/
  starting_time = time(0) ;
  set_saving_options( "PRODUCTION" ) ;
  if( PRODUCTION_STEPS != 0 ) {
    int production_steps = PRODUCTION_STEPS ;
    int total_step_div_10 = floor( production_steps/10. ) ;
    if( total_step_div_10 == 0 ) total_step_div_10 = 1 ;
    MonteCarlo_stats.clear() ;
    if( strcmp( questions, "NO_STOP" ) != 0 ) {
      control = request( "\n Do you want to change the number of production steps?" , "y/n" , '/' ) ;
      if( control == 'y' ) {
	printf( "  Insert the number of production steps: _" ) ;
	command_output = scanf( "%d%*c" , &production_steps ) ;
	  if (command_output == 0) cout << "*** ERROR: input value " << production_steps << endl ;
      }
    }

    printf( "\n  PRODUCTION STAGE . . .\n" ) ; fflush(0) ;
    lists_rebuilds = 0 ;  // number of events of recalculation of the verlet lists
    initialize_storage_file( production_steps , seed ) ;     // open the file in which configurations are periodically saved
    if( SAVE.backup_interval == 0 ) SAVE.backup_interval = PRODUCTION_STEPS + 10 ;
    clear_runtime_calculations() ;
    // molecules positions update
    if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)

    check_action( 0 , starting_time , total_step_div_10 ) ;
    if( strcmp( simulation_ensemble , "MC_SWAP_CHARGE" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	MonteCarloSwapChargeStep() ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }

    } else if( strcmp( simulation_ensemble , "MC_NVT" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	MonteCarlo_NVT( old_position ) ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }

    } else if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	MonteCarlo_NVT_CMfixed( old_position ) ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }

    } else if( strcmp( simulation_ensemble , "MC_NPT" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	MonteCarlo_NPT( old_position , old_coordinates ) ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }
  
    save_runtime_calculations( current_step ) ;
    cout << endl ;
    if( CHECKS.bonds_lists_number > 0 ) cout << endl << "  Maximum bonding distance in last configuration:  " << compute_max_bond_distance() << endl ;
    cout << "## Attempted changes:\tN\tpos\tvol\trot" << endl ;
    cout << "\t\t" << MonteCarlo_stats.NoP_attempted << "\t" << MonteCarlo_stats.pos_attempted << "\t" << MonteCarlo_stats.vol_attempted << "\t" << MonteCarlo_stats.rot_attempted << endl ;
    if( MonteCarlo_stats.NoP_attempted != 0 ) printf( "\tAcceptance rate for changes in the number of particles = %lf%%\n" ,
						     (double)MonteCarlo_stats.NoP_accepted/MonteCarlo_stats.NoP_attempted*100. ) ;
    if( MonteCarlo_stats.pos_attempted != 0 ) printf( "\tAcceptance rate for coordinates changes = %lf%%\n\t\t(max displacement : %lf)\n" ,
						     (double)MonteCarlo_stats.pos_accepted/MonteCarlo_stats.pos_attempted*100. , MC_maxdisplacement ) ;
    if( MonteCarlo_stats.vol_attempted != 0 ) printf( "\tAcceptance rate for volume changes = %lf%%\n\t\t(max box side change : %lf)\n" ,
						     (double)MonteCarlo_stats.vol_accepted/MonteCarlo_stats.vol_attempted*100. , MC_maxboxdisp ) ;
    if( MonteCarlo_stats.chr_swap_attempted != 0 ) printf( "\tAcceptance rate for charge swap moves = %lf%%\n" ,
							  (double)MonteCarlo_stats.chr_swap_accepted/MonteCarlo_stats.chr_swap_attempted*100. ) ;
    if( MonteCarlo_stats.rot_attempted != 0 ) printf( "\tAcceptance rate for reorientation moves = %lf%%\n\t\t(max rotation : %lf)\n" ,
						     (double)MonteCarlo_stats.rot_accepted/MonteCarlo_stats.rot_attempted*100. , MC_maxrotation ) ;
    cout << "\tNumber of events of verlet lists recalculation : " << lists_rebuilds << endl ;
    cout << endl ;
    fflush( NULL ) ;
    close_storage_file() ;
  }

  delete [] old_coordinates ;
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::MonteCarloSwapChargeStep( void ) {
  int candidate1 , candidate2 ;
  double charge1 , charge2 , old_energy , new_energy , random ;
  for( int i=0 ; i<MCS_LENGTH ; i++ ) {

    do candidate1 = (int)( (double)numberOfPart * drand48() ) ; while ( particles[candidate1]->charge == 0 ) ;
    charge1 = particles[candidate1]->charge ;
    do candidate2 = (int)( (double)numberOfPart * drand48() ) ; while ( particles[candidate2]->charge == 0 ) ;
    charge2 = particles[candidate2]->charge ;
    // When the fraction of charged particles is too low/high you could do it in a more efficient way . . .
    //                                   cout << ++ERR_FLAG << " " << candidate1 << " " << charge1 << " " << candidate2 << " " << charge2 << " " << POTENTIAL_ENERGY << endl ; fflush(0) ;
   if( charge2 == charge1 ) {
     while( particles[candidate2]->charge == charge1 || particles[candidate2]->charge == 0 ) candidate2 = (int)( (double)numberOfPart * drand48() ) ;
     charge2 = particles[candidate2]->charge ;
   }
   cout << "Consider crosslinkers!" << endl ; exit(1) ;

    // charge swap
    MonteCarlo_stats.chr_swap_attempted ++ ;
    old_energy = compute_potential_energy( particles[candidate1] ) + compute_potential_energy( particles[candidate2] ) - compute_potential_energy( particles[candidate1] , particles[candidate2] ) ;
    particles[candidate1]->charge = charge2 ;
    particles[candidate2]->charge = charge1 ;
    new_energy = compute_potential_energy( particles[candidate1] ) + compute_potential_energy( particles[candidate2] ) - compute_potential_energy( particles[candidate1] , particles[candidate2] ) ;
    // if random [0;1) number is lower than acceptance probability, then accept the swap move
    random = drand48() ;
    if( random < exp( (old_energy-new_energy)/k_T ) ) {
      MonteCarlo_stats.chr_swap_accepted ++ ;
      POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
    }else{
      particles[candidate1]->charge = charge1 ;
      particles[candidate2]->charge = charge2 ;
    }
    //    cout << MonteCarlo_stats.chr_swap_accepted << " " << POTENTIAL_ENERGY << " " << old_energy << " " << new_energy << " " << exp( (old_energy-new_energy)/k_T ) << " " << random << endl ;
  }
  global_time += time_step ;
  POTENTIAL_ENERGY.time = global_time ;
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::MonteCarlo_NVT( particle old_position ) {
  static int candidate ;
  static double old_energy , new_energy , random , acceptance_probability = 0 , halfmaxdisp ;
  static particle displacement ;
  halfmaxdisp = 0.5*MC_maxdisplacement ;
  for( int i=0 ; i<MCS_LENGTH ; i++ ) {
    candidate = (int)( (double)numberOfPart * drand48()) ;
    old_energy = compute_potential_energy( particles[candidate] ) ;
    old_position.position = particles[candidate]->position ;
    old_position.periodic_box = particles[candidate]->periodic_box ;

    MonteCarlo_stats.pos_attempted ++ ;
    displacement.position.rndPos_generate( MC_maxdisplacement ) ;
    displacement.position -= halfmaxdisp ;
    particles[candidate]->position += displacement.position ;
    check_verlet_update( particles[candidate] ) ;            // control if I need to update verlet list
    new_energy = compute_potential_energy( particles[candidate] ) ;
    if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
    else  acceptance_probability = exp( (old_energy - new_energy) / k_T ) ;
    random = drand48() ;
    if( random < acceptance_probability ) {
      MonteCarlo_stats.pos_accepted ++ ;
      POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
    } else {
      // To correctly assess pbc if verlet list rebuild event occurred
      if( particles[candidate]->periodic_box != old_position.periodic_box ) particles[candidate]->position -= displacement.position ;
      else particles[candidate]->position = old_position.position ;
    }

  }
  global_time += time_step ;
  POTENTIAL_ENERGY.time = global_time ;
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::MonteCarlo_NVT_CMfixed( particle old_position ) {
  static int candidate ;
  static double old_energy , new_energy , CoM_energy_correction , random , acceptance_probability = 0 , halfmaxdisp ;
  static particle displacement , old_CoM_disp ;
  halfmaxdisp = 0.5*MC_maxdisplacement ;
  for( int i=0 ; i<MCS_LENGTH ; i++ ) {
    candidate = (int)( (double)numberOfPart * drand48()) ;
    old_energy = compute_potential_energy( particles[candidate] ) ;
    old_position.position = particles[candidate]->position ;
    old_position.periodic_box = particles[candidate]->periodic_box ;
    old_CoM_disp.position = CoM_displacement.position ;

    MonteCarlo_stats.pos_attempted ++ ;
    displacement.position.rndPos_generate( MC_maxdisplacement ) ;
    displacement.position -= halfmaxdisp ;
    particles[candidate]->position += displacement.position ;
    CoM_displacement.position += ( displacement.position / (double)numberOfPart ) ;
    check_verlet_update( particles[candidate] ) ;            // check if I need to update verlet list
    new_energy = compute_potential_energy( particles[candidate] ) ;
    CoM_energy_correction = delta_energy_constrained_CoM( &displacement ) ;
    if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
    else  acceptance_probability = exp( (old_energy - new_energy - CoM_energy_correction) / k_T ) ;
    random = drand48() ;
    if( random < acceptance_probability ) {
      MonteCarlo_stats.pos_accepted ++ ;
      POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
      deltaEnergyOfCoM += CoM_energy_correction ;
    } else {
      // To correctly assess pbc if verlet list rebuild event occurred
      if( particles[candidate]->periodic_box != old_position.periodic_box ) particles[candidate]->position -= displacement.position ;
      else particles[candidate]->position = old_position.position ;
      CoM_displacement.position = old_CoM_disp.position ;
    }

  }
  global_time += time_step ;
  POTENTIAL_ENERGY.time = global_time ;
}
/*****************************************************************************************/

template <>
inline void configuration<patchy_2D>::MonteCarlo_NVT( patchy_2D old_position ) {
  static int candidate ;
  static double old_energy , new_energy , random , acceptance_probability = 0 , halfmaxdisp ;
  halfmaxdisp = 0.5*MC_maxdisplacement ;  // needed in the equilibration stage, when MC_maxdisplacement changes...
  static double rotjump_thresh = 1.0 - (1.0-MC_translation_probability) / 10 , rotsign_thresh = 1.0 - (1.0-MC_translation_probability) / 20 ;
  static patchy_2D displacement ;

  for( int i=0 ; i<MCS_LENGTH ; i++ ) {
    candidate = (int)( (double)numberOfPart * drand48()) ;
    old_energy = compute_potential_energy( particles[candidate] ) ;

    random = drand48() ;
    if ( random < MC_translation_probability ) {
      MonteCarlo_stats.pos_attempted ++ ;
      old_position.position = particles[candidate]->position ;
      old_position.periodic_box = particles[candidate]->periodic_box ;

      displacement.position.rndPos_generate( MC_maxdisplacement ) ;
      displacement.position -= halfmaxdisp ;
      particles[candidate]->position += displacement.position ;
      check_verlet_update( particles[candidate] ) ;
      new_energy = compute_potential_energy( particles[candidate] ) ;
      if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
      else {
	acceptance_probability = exp( (old_energy - new_energy) / k_T ) ;
	random = drand48() ;
      }
      if( random < acceptance_probability ) {
	MonteCarlo_stats.pos_accepted ++ ;
	POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
      } else {
	// To correctly assess pbc if verlet list rebuild event occurred
	if( particles[candidate]->periodic_box != old_position.periodic_box ) particles[candidate]->position -= displacement.position ;
	else particles[candidate]->position = old_position.position ;
      }

    } else {
      MonteCarlo_stats.rot_attempted ++ ;
      old_position.rotation = particles[candidate]->rotation ;

      particles[candidate]->rand_rotate( MC_maxrotation ) ;
      if( random > rotjump_thresh ) particles[candidate]->rotation += MC_rotationjump * ( 2*(random > rotsign_thresh) - 1 ) ;
      new_energy = compute_potential_energy( particles[candidate] ) ;
      acceptance_probability = exp( (old_energy - new_energy) / k_T ) ;
      random = drand48() ;
      if( random < acceptance_probability ) {
	MonteCarlo_stats.rot_accepted ++ ;
	POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
      } else {
	// particles[candidate]->position = old_position.position ;
	particles[candidate]->rotation = old_position.rotation ;
      }
    }

  }
  global_time += time_step ;
  POTENTIAL_ENERGY.time = global_time ;
}
/*****************************************************************************************/

template <>
inline void configuration<patchy_2D>::MonteCarlo_NVT_CMfixed( patchy_2D old_position ) {
  static int candidate ;
  static double old_energy , new_energy , CoM_energy_correction , random , acceptance_probability = 0 , halfmaxdisp ;
  static double rotjump_thresh = 1.0 - (1.0-MC_translation_probability) / 10 , rotsign_thresh = 1.0 - (1.0-MC_translation_probability) / 20 ;
  static patchy_2D displacement , old_CoM_disp ;
  halfmaxdisp = 0.5*MC_maxdisplacement ;
  for( int i=0 ; i<MCS_LENGTH ; i++ ) {
    candidate = (int)( (double)numberOfPart * drand48()) ;
    old_energy = compute_potential_energy( particles[candidate] ) ;

    random = drand48() ;
    if ( random < MC_translation_probability ) {
      MonteCarlo_stats.pos_attempted ++ ;
      old_position.position = particles[candidate]->position ;
      old_position.periodic_box = particles[candidate]->periodic_box ;
      old_CoM_disp.position = CoM_displacement.position ;

      displacement.position.rndPos_generate( MC_maxdisplacement ) ;
      displacement.position -= halfmaxdisp ;
      particles[candidate]->position += displacement.position ;
      CoM_displacement.position += ( displacement.position / (double)numberOfPart ) ;
      check_verlet_update( particles[candidate] ) ;
      new_energy = compute_potential_energy( particles[candidate] ) ;
      CoM_energy_correction = delta_energy_constrained_CoM( &displacement ) ;
      if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
      else  acceptance_probability = exp( (old_energy - new_energy - CoM_energy_correction) / k_T ) ;
      random = drand48() ;
      if( random < acceptance_probability ) {
	MonteCarlo_stats.pos_accepted ++ ;
	POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
	deltaEnergyOfCoM += CoM_energy_correction ;
      } else {
	// To correctly assess pbc if verlet list rebuild event occurred
	if( particles[candidate]->periodic_box != old_position.periodic_box ) particles[candidate]->position -= displacement.position ;
	else particles[candidate]->position = old_position.position ;
      }

    } else {
      MonteCarlo_stats.rot_attempted ++ ;
      old_position.rotation = particles[candidate]->rotation ;

      particles[candidate]->rand_rotate( MC_maxrotation ) ;
      if( random > rotjump_thresh ) particles[candidate]->rotation += MC_rotationjump * ( 2*(random > rotsign_thresh) - 1 ) ;
      new_energy = compute_potential_energy( particles[candidate] ) ;
      acceptance_probability = exp( (old_energy - new_energy) / k_T ) ;
      random = drand48() ;
      if( random < acceptance_probability ) {
	MonteCarlo_stats.rot_accepted ++ ;
	POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
      } else {
	// particles[candidate]->position = old_position.position ;
	particles[candidate]->rotation = old_position.rotation ;
      }
    }

  }
  global_time += time_step ;
  POTENTIAL_ENERGY.time = global_time ;
}
/*****************************************************************************************/



template <typename particle>
inline void configuration<particle>::MonteCarlo_NPT( particle old_position , double *old_coordinates ) {
  // eliminare displacement dagli argomenti, può essere definito qui static, perchè a lui non servono le patches, come per old_position... in realtà se ne potrebbe fare a meno, se disentanglo completamente rotazione e traslazione . . .
  static int candidate ;
  static double old_energy , new_energy , random , acceptance_probability = 0 , halfmaxdisp ;
  static particle displacement ;
  halfmaxdisp = 0.5*MC_maxdisplacement ;
  for( int i=0 ; i<MCS_LENGTH ; i++ ) {
    candidate = (int)( (double)(numberOfPart+DIMENSION) * drand48()) ;
    if( candidate < numberOfPart ) {
      old_energy = compute_potential_energy( particles[candidate] ) ;
      old_position.position = particles[candidate]->position ;
      old_position.periodic_box = particles[candidate]->periodic_box ;

      MonteCarlo_stats.pos_attempted ++ ;
      displacement.position.rndPos_generate( MC_maxdisplacement ) ;
      displacement.position -= halfmaxdisp ;
      particles[candidate]->position += displacement.position ;
      check_verlet_update( particles[candidate] ) ;
      new_energy = compute_potential_energy( particles[candidate] ) ;
      if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
      else {
	acceptance_probability = exp( (old_energy - new_energy) / k_T ) ;
	random = drand48() ;
      }
      if( random < acceptance_probability ) {
	MonteCarlo_stats.pos_accepted ++ ;
	POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
      } else {
	// To correctly assess pbc if verlet list rebuild event occurred
	if( particles[candidate]->periodic_box != old_position.periodic_box ) particles[candidate]->position -= displacement.position ;
	else particles[candidate]->position = old_position.position ;
      }

    } else {
      // volume_move( candidate - numberOfPart , old_coordinates ) ;
      volume_move_isotropic() ;

    }
  }
  global_time += time_step ;
  POTENTIAL_ENERGY.time = global_time ;
}
/*****************************************************************************************/

template <>
inline void configuration<patchy_2D>::MonteCarlo_NPT( patchy_2D old_position , double *old_coordinates ) {
  static int candidate ;
  static double old_energy , new_energy , random , acceptance_probability = 0 , halfmaxdisp ;
  static double rotjump_thresh = 1.0 - (1.0-MC_translation_probability) / 10 , rotsign_thresh = 1.0 - (1.0-MC_translation_probability) / 20 ;
  static patchy_2D displacement ;
  halfmaxdisp = 0.5*MC_maxdisplacement ;
  for( int i=0 ; i<MCS_LENGTH ; i++ ) {
    candidate = (int)( (double)(numberOfPart+DIMENSION) * drand48()) ;
    if( candidate < numberOfPart ) {
      old_energy = compute_potential_energy( particles[candidate] ) ;

      random = drand48() ;
      if ( random < MC_translation_probability ) {
	MonteCarlo_stats.pos_attempted ++ ;
	old_position.position = particles[candidate]->position ;
	old_position.periodic_box = particles[candidate]->periodic_box ;

	displacement.position.rndPos_generate( MC_maxdisplacement ) ;
	displacement.position -= halfmaxdisp ;
	particles[candidate]->position += displacement.position ;
	check_verlet_update( particles[candidate] ) ;
	new_energy = compute_potential_energy( particles[candidate] ) ;
	if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
	else {
	  acceptance_probability = exp( (old_energy - new_energy) / k_T ) ;
	  random = drand48() ;
	}
	if( random < acceptance_probability ) {
	  MonteCarlo_stats.pos_accepted ++ ;
	  POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
	} else {
	// To correctly assess pbc if verlet list rebuild event occurred
	  if( particles[candidate]->periodic_box != old_position.periodic_box ) particles[candidate]->position -= displacement.position ;
	  else particles[candidate]->position = old_position.position ;
	}

      } else {
	MonteCarlo_stats.rot_attempted ++ ;
	old_position.rotation = particles[candidate]->rotation ;

	particles[candidate]->rand_rotate( MC_maxrotation ) ;
	if( random > rotjump_thresh ) particles[candidate]->rotation += MC_rotationjump * ( 2*(random > rotsign_thresh) - 1 ) ;
	new_energy = compute_potential_energy( particles[candidate] ) ;
	acceptance_probability = exp( (old_energy - new_energy) / k_T ) ;
	random = drand48() ;
	if( random < acceptance_probability ) {
	  MonteCarlo_stats.rot_accepted ++ ;
	  POTENTIAL_ENERGY.val += ( new_energy - old_energy ) ;
	} else {
	  particles[candidate]->rotation = old_position.rotation ;
	}
      }
    } else {
      // volume_move( candidate - numberOfPart , old_coordinates ) ;
      volume_move_isotropic() ;
    }
  }
  global_time += time_step ;
  POTENTIAL_ENERGY.time = global_time ;
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::volume_move_isotropic( void ) {
  static double old_energy , new_energy , side_ratio , old_volume , new_volume , random , acceptance_probability = 0 ;
  static particle old_side ;
  static particle *old_coordinates = new particle [numberOfPart] ;
  bool list_rebuild ;

  MonteCarlo_stats.vol_attempted ++ ;
  // volume change
  old_side.position = box_sides.position ;
  old_volume = box_sides.volume() ;
  old_energy = POTENTIAL_ENERGY.val ;
  box_sides.position += ( (drand48() - 0.5) * MC_maxboxdisp ) ;
  midside.position = (box_sides.position * 0.5) ;
  side_ratio = box_sides.position.x / old_side.position.x ;
  for( int i=0 ; i<numberOfPart ; i++ ) {
    old_coordinates[i].position = particles[i]->position ;
    particles[i]->position *= side_ratio ;
  }
  list_rebuild = check_verlet_update() ;
  new_volume = box_sides.volume() ;
  new_energy = compute_potential_energy() ;
  if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
  else {
    acceptance_probability = pow(new_volume/old_volume , (double)numberOfPart) * exp(( NPT_pressure * (old_volume-new_volume) + old_energy - new_energy ) / k_T) ;
    random = drand48() ;
  }
  if( random < acceptance_probability ) {
    MonteCarlo_stats.vol_accepted ++ ;
  } else {
    box_sides.position = old_side.position ;
    midside.position = (box_sides.position * 0.5) ;
    POTENTIAL_ENERGY.val = old_energy ;
    if( list_rebuild ) for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position /= side_ratio ;
    else for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position = old_coordinates[i].position ;
    list_rebuild = check_verlet_update() ;            // check if I need to update verlet list after moving back
  }
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::volume_move( int direction , double *old_coordinates ) {
  static double old_energy , new_energy , old_side , side_ratio , old_volume , new_volume , random , acceptance_probability = 0 ;
  bool list_rebuild ;

  MonteCarlo_stats.vol_attempted ++ ;
  if( direction == 0 ) {  // Lx change
    old_side = box_sides.position.x ;
    old_volume = box_sides.volume() ;
    old_energy = POTENTIAL_ENERGY.val ;
    box_sides.position.x += (drand48() - 0.5) * MC_maxboxdisp ;
    midside.position = (box_sides.position * 0.5) ;
    side_ratio = box_sides.position.x / old_side ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      old_coordinates[i] = particles[i]->position.x ;
      particles[i]->position.x *= side_ratio ;
    }
    list_rebuild = check_verlet_update() ;            // check if I need to update verlet list
    new_volume = box_sides.volume() ;
    new_energy = compute_potential_energy() ;
    if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
    else {
      acceptance_probability = pow(new_volume/old_volume , (double)numberOfPart) * exp(( NPT_pressure * (old_volume-new_volume) + old_energy - new_energy ) / k_T) ;
      random = drand48() ;
    }
    if( random < acceptance_probability ) {
      MonteCarlo_stats.vol_accepted ++ ;
    } else {
      box_sides.position.x = old_side ;
      midside.position = (box_sides.position * 0.5) ;
      POTENTIAL_ENERGY.val = old_energy ;
      if( list_rebuild ) for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.x /= side_ratio ;
      else for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.x = old_coordinates[i] ;
      list_rebuild = check_verlet_update() ;            // check if I need to update verlet list after moving back
    }

  } else {  // Ly change
    old_side = box_sides.position.y ;
    old_volume = box_sides.volume() ;
    old_energy = POTENTIAL_ENERGY.val ;
    box_sides.position.y += (drand48() - 0.5) * MC_maxboxdisp ;
    midside.position = (box_sides.position * 0.5) ;
    side_ratio = box_sides.position.y / old_side ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      old_coordinates[i] = particles[i]->position.y ;
      particles[i]->position.y *= side_ratio ;
    }
    list_rebuild = check_verlet_update() ;            // check if I need to update verlet list
    new_volume = box_sides.volume() ;
    new_energy = compute_potential_energy() ;
    if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
    else {
      acceptance_probability = pow(new_volume/old_volume , (double)numberOfPart) * exp(( NPT_pressure * (old_volume-new_volume) + old_energy - new_energy ) / k_T) ;
      random = drand48() ;
    }
    if( random < acceptance_probability ) {
      MonteCarlo_stats.vol_accepted ++ ;
    } else {
      box_sides.position.y = old_side ;
      midside.position = (box_sides.position * 0.5) ;
      POTENTIAL_ENERGY.val = old_energy ;
      if( list_rebuild ) for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.y /= side_ratio ;
      else for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.y = old_coordinates[i] ;
      list_rebuild = check_verlet_update() ;            // check if I need to update verlet list after moving back
    }
  }
}
/*****************************************************************************************/

template <>
inline void configuration<particle_3D>::volume_move( int direction , double *old_coordinates ) {
  static double old_energy , new_energy , old_side , side_ratio , old_volume , new_volume , random , acceptance_probability = 0 ;
  bool list_rebuild ;

  MonteCarlo_stats.vol_attempted ++ ;
  if( direction == 0 ) {  // Lx change
    old_side = box_sides.position.x ;
    old_volume = box_sides.volume() ;
    old_energy = POTENTIAL_ENERGY.val ;
    box_sides.position.x += (drand48() - 0.5) * MC_maxboxdisp ;
    midside.position = (box_sides.position * 0.5) ;
    side_ratio = box_sides.position.x / old_side ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      old_coordinates[i] = particles[i]->position.x ;
      particles[i]->position.x *= side_ratio ;
    }
    list_rebuild = check_verlet_update() ;            // check if I need to update verlet list
    new_volume = box_sides.volume() ;
    new_energy = compute_potential_energy() ;
    if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
    else {
      acceptance_probability = pow(new_volume/old_volume , (double)numberOfPart) * exp(( NPT_pressure * (old_volume-new_volume) + old_energy - new_energy ) / k_T) ;
      random = drand48() ;
    }
    if( random < acceptance_probability ) {
      MonteCarlo_stats.vol_accepted ++ ;
    } else {
      box_sides.position.x = old_side ;
      midside.position = (box_sides.position * 0.5) ;
      POTENTIAL_ENERGY.val = old_energy ;
      if( list_rebuild ) for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.x /= side_ratio ;
      else for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.x = old_coordinates[i] ;
      list_rebuild = check_verlet_update() ;            // check if I need to update verlet list after moving back
    }

  } else if( direction == 1 ) {  // Ly change
    old_side = box_sides.position.y ;
    old_volume = box_sides.volume() ;
    old_energy = POTENTIAL_ENERGY.val ;
    box_sides.position.y += (drand48() - 0.5) * MC_maxboxdisp ;
    midside.position = (box_sides.position * 0.5) ;
    side_ratio = box_sides.position.y / old_side ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      old_coordinates[i] = particles[i]->position.y ;
      particles[i]->position.y *= side_ratio ;
    }
    list_rebuild = check_verlet_update() ;            // check if I need to update verlet list
    new_volume = box_sides.volume() ;
    new_energy = compute_potential_energy() ;
    if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
    else {
      acceptance_probability = pow(new_volume/old_volume , (double)numberOfPart) * exp(( NPT_pressure * (old_volume-new_volume) + old_energy - new_energy ) / k_T) ;
      random = drand48() ;
    }
    if( random < acceptance_probability ) {
      MonteCarlo_stats.vol_accepted ++ ;
    } else {
      box_sides.position.y = old_side ;
      midside.position = (box_sides.position * 0.5) ;
      POTENTIAL_ENERGY.val = old_energy ;
      if( list_rebuild ) for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.y /= side_ratio ;
      else for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.y = old_coordinates[i] ;
      list_rebuild = check_verlet_update() ;            // check if I need to update verlet list after moving back
    }

  } else {  // Lz change
    old_side = box_sides.position.z ;
    old_volume = box_sides.volume() ;
    old_energy = POTENTIAL_ENERGY.val ;
    box_sides.position.z += (drand48() - 0.5) * MC_maxboxdisp ;
    midside.position = (box_sides.position * 0.5) ;
    side_ratio = box_sides.position.z / old_side ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      old_coordinates[i] = particles[i]->position.z ;
      particles[i]->position.z *= side_ratio ;
    }
    list_rebuild = check_verlet_update() ;            // check if I need to update verlet list
    new_volume = box_sides.volume() ;
    new_energy = compute_potential_energy() ;
    if ( new_energy == INFINITY ) acceptance_probability = 0.0 ;
    else {
      acceptance_probability = pow(new_volume/old_volume , (double)numberOfPart) * exp(( NPT_pressure * (old_volume-new_volume) + old_energy - new_energy ) / k_T) ;
      random = drand48() ;
    }
    if( random < acceptance_probability ) {
      MonteCarlo_stats.vol_accepted ++ ;
    } else {
      box_sides.position.z = old_side ;
      midside.position = (box_sides.position * 0.5) ;
      POTENTIAL_ENERGY.val = old_energy ;
      if( list_rebuild ) for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.z /= side_ratio ;
      else for( int i=0 ; i<numberOfPart ; i++ ) particles[i]->position.z = old_coordinates[i] ;
      list_rebuild = check_verlet_update() ;            // check if I need to update verlet list after moving back
    }
  }
}
/*****************************************************************************************/
