/* this file contains methods to initialize output files and to save thermodynamic and per particle data along the simulations */
// MD VERSION
#include <limits>
#include "../system_comm.h"

template <typename particle>
void configuration<particle>::plot_equilibration_data(void) {
  cout << " ERROR: The function plot_equilibration_data() has still to be implemented!! " << endl ;
  exit( EXIT_FAILURE ) ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::save_options::copy_options( const save_options &other ) {
  PARTICLES = other.PARTICLES ;
  PARTICLES_IMG = other.PARTICLES_IMG ;
  VELOCITIES = other.VELOCITIES ;
  FORCES = other.FORCES ;
  PER_PART_ENERGY = other.PER_PART_ENERGY ;
  PER_PART_VIRIAL = other.PER_PART_VIRIAL ;
  FORCES_SPLIT = other.FORCES_SPLIT ;
  MOL_IDX = other.MOL_IDX ;
  IDX = other.IDX ;
  TYPE = other.TYPE ;
  CHARGE = other.CHARGE ;
  MOLECULES = other.MOLECULES ;
  MOLECULES_IMG = other.MOLECULES_IMG ;
  MOL_VELOCITIES = other.MOL_VELOCITIES ;
  MOL_FORCES = other.MOL_FORCES ;
  MOL_FORCES_SPLIT = other.MOL_FORCES_SPLIT ;
  POTENTIAL_ENERGY = other.POTENTIAL_ENERGY ;
  PARTIAL_POTENTIAL_ENERGY = other.PARTIAL_POTENTIAL_ENERGY ;
  KINETIC_ENERGY = other.KINETIC_ENERGY ;
  SECONDARY_POTENTIAL = other.SECONDARY_POTENTIAL ;
  TEMPERATURE = other.TEMPERATURE ;
  VOLUME = other.VOLUME ;
  MSD = other.MSD ;
  GYR = other.GYR ;
  RUNTIME_CALCULATIONS = other.RUNTIME_CALCULATIONS ;
  VIRIAL = other.VIRIAL ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::save_options::clear( void ) {
  PARTICLES = 0 ;
  PARTICLES_IMG = 0 ;
  VELOCITIES = 0 ;
  FORCES = 0 ;
  PER_PART_ENERGY = 0 ;
  PER_PART_VIRIAL = 0 ;
  FORCES_SPLIT = 0 ;
  MOL_IDX = 0 ;
  IDX = 0 ;
  TYPE = 0 ;
  CHARGE = 0 ;
  MOLECULES = 0 ;
  MOLECULES_IMG = 0 ;
  MOL_VELOCITIES = 0 ;
  MOL_FORCES = 0 ;
  MOL_FORCES_SPLIT = 0 ;
  POTENTIAL_ENERGY = 0 ;
  PARTIAL_POTENTIAL_ENERGY = 0 ;
  KINETIC_ENERGY = 0 ;
  SECONDARY_POTENTIAL = 0 ;
  TEMPERATURE = 0 ;
  VOLUME = 0 ;
  MSD = 0 ;
  GYR = 0 ;
  RUNTIME_CALCULATIONS = 0 ;
  VIRIAL = 0 ;
}
/****************************************************************************************/

template <typename particle>
configuration<particle>::save_options::save_options( void ) {
  PARTICLES = 0 ;
  PARTICLES_IMG = 0 ;
  VELOCITIES = 0 ;
  FORCES = 0 ;
  PER_PART_ENERGY = 0 ;
  PER_PART_VIRIAL = 0 ;
  FORCES_SPLIT = 0 ;
  MOL_IDX = 0 ;
  IDX = 0 ;
  TYPE = 0 ;
  CHARGE = 0 ;
  MOLECULES = 0 ;
  MOLECULES_IMG = 0 ;
  MOL_VELOCITIES = 0 ;
  MOL_FORCES = 0 ;
  MOL_FORCES_SPLIT = 0 ;
  POTENTIAL_ENERGY = 0 ;
  PARTIAL_POTENTIAL_ENERGY = 0 ;
  KINETIC_ENERGY = 0 ;
  SECONDARY_POTENTIAL = 0 ;
  TEMPERATURE = 0 ;
  VOLUME = 0 ;
  MSD = 0 ;
  GYR = 0 ;
  RUNTIME_CALCULATIONS = 0 ;
  VIRIAL = 0 ;
}
/****************************************************************************************/

template <typename particle>
char *configuration<particle>::generate_new_directory_for_data(void) {
  int command_output;
  if( strcmp( _DIRECTORY_NAME, "" ) == 0 ) {
    time_t current;
    time(&current);
    struct tm *date_hour = localtime(&current);
    strftime( _DIRECTORY_NAME , 300 , "simulation_%I:%M_%d-%m-%Y" , date_hour );
  }
  char dir_creation[400];
  strcpy( dir_creation , "mkdir " );
  if( strcat( dir_creation , _DIRECTORY_NAME ) == NULL) cout << " Problems in generate_new_directory_for_data()\n";
  command_output = system( dir_creation );
  if (command_output < 0) cout << "*** ERROR: " << dir_creation << endl << "in method generate_new_directory_for_data()" << endl ;

  return _DIRECTORY_NAME;
}
/****************************************************************************************/

template <typename particle>
char *configuration<particle>::generate_equilibration_subdirectory(char *root) {
  int command_output;
  if( strcpy( root , _DIRECTORY_NAME ) == NULL) cout << " Problems in generate_equilibration_subdirectory\n";
  if( strcat( _DIRECTORY_NAME , "/equilibration_run" ) == NULL) cout << " Problems in generate_equilibration_subdirectory()\n";
  char dir_creation[300];
  strcpy( dir_creation , "mkdir " );
  if( strcat( dir_creation , _DIRECTORY_NAME ) == NULL) cout << " Problems in generate_equilibration_subdirectory()\n";
  command_output = system( dir_creation );
  if (command_output < 0) cout << "*** ERROR: " << dir_creation << endl << "in method generate_equilibration_subdirectory()" << endl ;

  return root;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle> :: initialize_storage_file ( int mcsteps , double seed ) {
  char file_path[300] ;
  // initialization of the supporting variables to compute the contribution of the different kind of interactions
  if( SAVE.MOL_FORCES == 1 && SAVE.MOL_FORCES_SPLIT == 1 ) {
    SAVE.force_contributions = new particle ** [ 2*INTERACTIONS.number ] ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      SAVE.force_contributions[k] = new particle * [ numberOfPart ] ;
      for( int i=0 ; i<numberOfPart ; i++ ) SAVE.force_contributions[k][i] = new particle ;
    }
    for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) {
      SAVE.force_contributions[k] = new particle * [ numberOfMol ] ;
      for( int i=0 ; i<numberOfMol ; i++ ) SAVE.force_contributions[k][i] = new particle ;
    }
  } else if( SAVE.FORCES == 1 && SAVE.FORCES_SPLIT == 1 ) {
    SAVE.force_contributions = new particle ** [ INTERACTIONS.number ] ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      SAVE.force_contributions[k] = new particle * [ numberOfPart ] ;
      for( int i=0 ; i<numberOfPart ; i++ ) SAVE.force_contributions[k][i] = new particle ;
    }
  } else SAVE.force_contributions = NULL ;
  if( SAVE.PER_PART_ENERGY == 1 ) {
    SAVE.per_part_energy = new double * [ INTERACTIONS.number+1 ] ;   // the last is the total
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) SAVE.per_part_energy[k] = new double [ numberOfPart ] ;
  } else SAVE.per_part_energy = NULL ;
  if( SAVE.PER_PART_VIRIAL == 1 ) {
    SAVE.per_part_virial = new double * [ INTERACTIONS.number+1 ] ;   // the last is the total
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) SAVE.per_part_virial[k] = new double [ numberOfPart ] ;
  } else SAVE.per_part_virial = NULL ;

  strcpy( file_path , _DIRECTORY_NAME ) ;
  if( strcat( file_path , "/interactions_data.dat" ) == NULL) cout << " Problems in initialize_storage_file\n" ;
  ofstream _int_data ;
  _int_data.open( file_path ) ;
  _int_data << "total_steps " << mcsteps << " drand48_seed " << seed << endl ;
  for( int i=0 ; i<INTERACTIONS.number ; i++ ) {
    _int_data << endl << "INTERACTION_POTENTIAL " << i << " : cut-off " << INTERACTIONS.type[i]->cutoff << " ; shift " << INTERACTIONS.type[i]->SHIFT << " ; verlet_radii " << INTERACTIONS.type[i]->verlet_list->r_verlet << " , " << INTERACTIONS.type[i]->verlet_list->delta_verlet << " ; list_updated " << INTERACTIONS.type[i]->verlet_list->must_be_updated << " ; list_displacements " << INTERACTIONS.type[i]->verlet_list->disp_on << " ; DOUBLE_PAIRS " << INTERACTIONS.type[i]->verlet_list->DOUBLE_PAIRS << " ." << endl << INTERACTIONS.type[i]->ostr() << endl ;
  }
  _int_data.close() ;

  if( SAVE.PARTICLES == 1 || SAVE.VELOCITIES == 1 || SAVE.FORCES == 1 || SAVE.PER_PART_ENERGY == 1 || SAVE.PER_PART_VIRIAL == 1 ) {
    if( strcmp( SAVE.configuration_format , "xyz" ) == 0 ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      if( strcat( file_path , "/particles.dat" ) == NULL) cout << " Problems in initialize_storage_file\n" ;
      _particles_file.open( file_path ) ;
      _particles_file << setprecision(12) ;
      _particles_file << "# Columns_meaning:" ;
      if( SAVE.PARTICLES == 1 ) {
	_particles_file << " X/SIGMA Y" ;
	if( DIMENSION == 3 ) _particles_file << " Z" ;
	if( ClassName( particles[0] ).compare(string( "patchy_2D" )) == 0 ) _particles_file << " angle" ;
	if( SAVE.PARTICLES_IMG == 1 ) {
	  _particles_file << " imgX/box_side imgY" ;
	  if( DIMENSION == 3 ) _particles_file << " imgZ" ;
	}
      }
      if( SAVE.VELOCITIES == 1 ) {
	_particles_file << " Vx/sqrt(EPSILON/MASS) Vy" ;
	if( DIMENSION == 3 ) _particles_file << " Vz" ;
      }
      if( SAVE.FORCES == 1 ) {
	if( SAVE.FORCES_SPLIT == 1 ) {
	  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
	    _particles_file << " Fx_" << k << "/(EPSILON/SIGMA) Fy_" << k ;
	    if( DIMENSION == 3 ) _particles_file << " Fz_" << k ;
	  }
	} else if( SAVE.FORCES_SPLIT == 0 ) {
	  _particles_file << " Fx/(EPSILON/SIGMA) Fy" ;
	  if( DIMENSION == 3 ) _particles_file << " Fz" ;
	}
      }
      if( SAVE.PER_PART_ENERGY == 1 ) {
	_particles_file << " energy/EPSILON" ;
      }
      if( SAVE.PER_PART_VIRIAL == 1 ) {
	_particles_file << " virial/EPSILON" ;
      }
      _particles_file << endl << endl ;

    } else if( strcmp( SAVE.configuration_format , "patch" ) == 0 ||
	       strcmp( SAVE.configuration_format , "ptc" ) == 0 ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      if( strcat( file_path , "/particles.patch" ) == NULL ) cout << " Problems in initialize_storage_file\n" ;
      _particles_file.open( file_path ) ;
      _particles_file << setprecision(12) ;

    } else if( strcmp( SAVE.configuration_format , "sph" ) == 0 ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      if( strcat( file_path , "/particles.sph" ) == NULL ) cout << " Problems in initialize_storage_file\n" ;
      _particles_file.open( file_path ) ;
      _particles_file << setprecision(12) ;
    }

    if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      if( strcat( file_path , "/unconstrained_CoM.dat" ) == NULL ) cout << " Problems in initialize_storage_file\n" ;
      ofstream _uCoM_file ;
      _uCoM_file.open( file_path ) ;
      _uCoM_file << "# step time CoM_position(SIGMA)" << endl ;
      _uCoM_file.close() ;
    }
  }

  if( SAVE.MOLECULES == 1 || SAVE.MOL_VELOCITIES == 1 || SAVE.MOL_FORCES == 1 ) {
    strcpy( file_path , _DIRECTORY_NAME ) ;
    if( strcat( file_path , "/molecules.dat" ) == NULL) cout << " Problems in initialize_storage_file\n" ;
    _molecules_file.open( file_path ) ;
    _molecules_file << setprecision(12) ;
    _molecules_file << "# Columns_meaning:" ;
    if( SAVE.MOLECULES == 1 ) {
      _molecules_file << " X/SIGMA Y" ;
      if( DIMENSION == 3 ) _molecules_file << " Z" ;
      if( SAVE.MOLECULES_IMG == 1 ) {
	_molecules_file << " imgX/box_side imgY" ;
	if( DIMENSION == 3 ) _molecules_file << " imgZ" ;
      }
    }
    if( SAVE.MOL_VELOCITIES == 1 ) {
      _molecules_file << " Vx/sqrt(EPSILON/MASS) Vy" ;
      if( DIMENSION == 3 ) _molecules_file << " Vz" ;
    }
    if( SAVE.MOL_FORCES == 1 ) {
      if( SAVE.MOL_FORCES_SPLIT == 1 ) {
	for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
	  _molecules_file << " Fx_" << k << "/(EPSILON/SIGMA) Fy_" << k ;
	  if( DIMENSION == 3 ) _molecules_file << " Fz_" << k ;
	}
      } else if( SAVE.MOL_FORCES_SPLIT == 0 ) {
	_molecules_file << " Fx/(EPSILON/SIGMA) Fy" ;
	if( DIMENSION == 3 ) _molecules_file << " Fz" ;
      }
    }
    _molecules_file << endl << endl ;
  }

  if( SAVE.POTENTIAL_ENERGY == 1 ||
      SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ||
      SAVE.KINETIC_ENERGY == 1 ||
      SAVE.SECONDARY_POTENTIAL == 1 ||
      SAVE.TEMPERATURE == 1 || SAVE.GYR == 1 || SAVE.VOLUME == 1 || SAVE.VIRIAL == 1 ) {

    strcpy( file_path , _DIRECTORY_NAME ) ;
    if( strcat( file_path , "/thermoinfo.dat" ) == NULL) cout << " Problems in initialize_storage_file\n" ;
    _thermoinfo_file.open( file_path ) ;
    _thermoinfo_file << "# step time" ;

    if( SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ) {
      for( int i=0 ; i<INTERACTIONS.number ; i++ ) _thermoinfo_file << " potential_" << i << "(EPSILON)" ;
    }
    if( SAVE.POTENTIAL_ENERGY == 1 ) {
      _thermoinfo_file << " total_potential_en(EPSILON)" ;
      if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
	_thermoinfo_file << " unconstrainedCM_delta_en(EPSILON) constrainedCM_total_en(EPSILON) lambda_derivative_perpart" ;
      }
    }
    if( SAVE.KINETIC_ENERGY == 1 ) {
      _thermoinfo_file << " kinetic_en(EPSILON)" ;
    }
    if( SAVE.SECONDARY_POTENTIAL == 1 ) {
      _thermoinfo_file << " secondary_potential_energy(EPSILON)" ;
    }
    if( SAVE.TEMPERATURE == 1 ) {
      _thermoinfo_file << " temperature(EPSILON)" ;
    }
    if( SAVE.GYR == 1 ) {
      _thermoinfo_file << " gyration_radius(SIGMA)" ;
    }
    if( SAVE.VOLUME == 1 ) {
      _thermoinfo_file << " boxside_x/(SIGMA) boxside_y" ;
      if( DIMENSION == 3 ) _thermoinfo_file << " boxside_z" ;
      _thermoinfo_file << " volume(SIGMA^" << DIMENSION << ")" ;
    }
    if( SAVE.VIRIAL == 1 ) {
      _thermoinfo_file << " virial_pressure(EPSILON/SIGMA^" << DIMENSION << ")" ;
    }
    _thermoinfo_file << endl ;
  }
  if( SAVE.MSD == 1 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) MSD[i]->initialize_MSDfiles( _DIRECTORY_NAME ) ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle> :: check_action ( int current_step , int starting_time , int screen_print_interval ) {
  if( current_step % screen_print_interval == 0 ) {
    int hour , min , sec , elapsed_time = time(0)-starting_time ;
    sec = elapsed_time % 60 ;
    min = ( (elapsed_time-sec) / 60 ) % 60 ;
    hour = (elapsed_time-sec-60*min) / 3600 ;
    printf( "\r   step %d ; elapsed time: %d:%d:%d", current_step , hour , min , sec ) ;
    fflush(0) ;
  }
  if( current_step % SAVE.global_interval == 0 ) save_data( current_step ) ;
  if( current_step % SAVE.configuration_interval == 0 ) save_configuration( current_step ) ;
  if( current_step % SAVE.backup_interval == 0 ) write_backup_data( current_step ) ;
  if( current_step % SAVE.runtime_interval == 0 ) save_runtime_calculations( current_step ) ;
  if( SAVE.MSD == 1 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) if( current_step == MSD[i]->saving_step ) MSD[i]->save_msd_conf( global_time , box_sides ) ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::save_configuration( int step ) {
  // Preliminary computations :
  if( ( SAVE.FORCES == 1 && SAVE.FORCES_SPLIT == 1 ) || ( SAVE.MOL_FORCES == 1 && SAVE.MOL_FORCES_SPLIT == 1 ) ) {
    particle **true_forces_t2 = forces_t2 ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      for( int i=0 ; i<numberOfPart ; i++ ) SAVE.force_contributions[k][i]->clear() ;
      forces_t2 = SAVE.force_contributions[k] ;
      INTERACTIONS.type[k]->compute_forces_contribution() ;
    }
    forces_t2 = true_forces_t2 ;
  }
  if( SAVE.MOL_FORCES == 1 ) {
    if( SAVE.MOL_FORCES_SPLIT == 1 ) {
      for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) {
	for( int i=0 ; i<numberOfMol ; i++ ) {
	  SAVE.force_contributions[k][i]->clear() ;
	  for( int j=0 ; j<molecules[i]->Natoms ; j++ ) SAVE.force_contributions[k][i]->position += SAVE.force_contributions[ k-INTERACTIONS.number ][ molecules[i]->atom[j]->id ]->position ;
	}
      }
    } else {
      for( int i=0 ; i<numberOfMol ; i++ ) molecules[i]->fcm.clear() ;
      for( int i=0 ; i<numberOfMol ; i++ ) for( int j=0 ; j<molecules[i]->Natoms ; j++ ) molecules[i]->fcm.position += forces_t1[ molecules[i]->atom[j]->id ]->position ;
    }
  }
  if( SAVE.MOL_VELOCITIES == 1 ) {
    for( int i=0 ; i<numberOfMol ; i++ ) molecules[i]->vcm.clear() ;
    for( int i=0 ; i<numberOfMol ; i++ ) for( int j=0 ; j<molecules[i]->Natoms ; j++ ) molecules[i]->vcm.position += velocities[ molecules[i]->atom[j]->id ]->position ;
  }
  if( SAVE.PER_PART_ENERGY == 1 ) {
    for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_energy[INTERACTIONS.number][i] = 0 ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) INTERACTIONS.type[k]->compute_per_part_energy_contribution( SAVE.per_part_energy[k] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_energy[INTERACTIONS.number][i] += SAVE.per_part_energy[k][i] ;
  }
  if( SAVE.PER_PART_VIRIAL == 1 ) {
    for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_virial[INTERACTIONS.number][i] = 0 ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) INTERACTIONS.type[k]->compute_per_part_virial_contribution( SAVE.per_part_virial[k] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_virial[INTERACTIONS.number][i] += SAVE.per_part_virial[k][i] ;
  }

  // Saving actions :
  if( strcmp( SAVE.configuration_format , "xyz" ) == 0 ) {
    if( SAVE.PARTICLES == 1 || SAVE.VELOCITIES == 1 || SAVE.FORCES == 1 || SAVE.PER_PART_ENERGY == 1 || SAVE.PER_PART_VIRIAL == 1 ) dump_xyz( _particles_file , 12 ) ;

  } else if ( strcmp( SAVE.configuration_format , "ptc" ) == 0 ||
	      strcmp( SAVE.configuration_format , "patch" ) == 0 ) {
    if( SAVE.PARTICLES == 1 ) dump_patch( _particles_file , 12 ) ;

  } else if ( strcmp( SAVE.configuration_format , "sph" ) == 0 ) {
    if( SAVE.PARTICLES == 1 ) dump_sph( _particles_file , 12 ) ;

  } else if ( strcmp( SAVE.configuration_format , "lmp" ) == 0 ) {
    if( SAVE.PARTICLES == 1 ) dump_lmp( _particles_file , 12 ) ;
  }

  if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
    char file_path[300] ;
    strcpy( file_path , _DIRECTORY_NAME ) ;
    if( strcat( file_path , "/unconstrained_CoM.dat" ) == NULL ) cout << " Problems in initialize_storage_file\n" ;
    ofstream _uCoM_file ;
    _uCoM_file.open( file_path , ios_base::app ) ;
    _uCoM_file << setprecision(12) ;
    _uCoM_file << " " << step << " " << global_time << " " << CoM_displacement.position << endl ;
    _uCoM_file.close() ;
  }

  if( SAVE.MOLECULES == 1 || SAVE.MOL_VELOCITIES == 1 || SAVE.MOL_FORCES == 1 ) {
    _molecules_file << setprecision(12) ;
    _molecules_file << "# timestep " << global_time << endl ;
    _molecules_file << "# N " << numberOfMol << endl ;
    _molecules_file << "# box " << box_sides.position << endl ;
    if( SAVE.MOLECULES == 1 && SAVE.MOL_VELOCITIES == 0 && SAVE.MOL_FORCES == 0 ) {
      if( SAVE.MOLECULES_IMG == 1 ) for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << " " << molecules[i]->cm.periodic_box << endl ;
      else for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << endl ;
    } else if( SAVE.MOLECULES == 0 && SAVE.MOL_VELOCITIES == 1 && SAVE.MOL_FORCES == 0 ) for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->vcm.position << endl ;
    else if( SAVE.MOLECULES == 1 && SAVE.MOL_VELOCITIES == 1 && SAVE.MOL_FORCES == 0 ) {
      if( SAVE.MOLECULES_IMG == 1 ) for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << " " << molecules[i]->cm.periodic_box << " " << molecules[i]->vcm.position << endl ;
      else for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << " " << molecules[i]->vcm.position << endl ;
    } else if( SAVE.MOL_FORCES == 1 && SAVE.MOL_FORCES_SPLIT == 1 ) {
      if( SAVE.MOLECULES == 0 && SAVE.MOL_VELOCITIES == 0 ) {
	for( int i=0 ; i<numberOfMol ; i++ ) {
	  for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) _molecules_file << SAVE.force_contributions[k][i]->position << " " ;
	  _molecules_file << endl ;
	}
      } else if( SAVE.MOLECULES == 1 && SAVE.MOL_VELOCITIES == 0 ) {
	if( SAVE.MOLECULES_IMG == 1 ) {
	  for( int i=0 ; i<numberOfMol ; i++ ) {
	    _molecules_file << molecules[i]->cm.position << " " << molecules[i]->cm.periodic_box << " " ;
	    for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) _molecules_file << SAVE.force_contributions[k][i]->position << " " ;
	    _molecules_file << endl ;
	  }
	} else {
	  for( int i=0 ; i<numberOfMol ; i++ ) {
	    _molecules_file << molecules[i]->cm.position << " " ;
	    for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) _molecules_file << SAVE.force_contributions[k][i]->position << " " ;
	    _molecules_file << endl ;
	  }
	}
      } else if( SAVE.MOLECULES == 0 && SAVE.MOL_VELOCITIES == 1 ) {
	for( int i=0 ; i<numberOfMol ; i++ ) {
	  _molecules_file << molecules[i]->vcm.position << " " ;
	  for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) _molecules_file << SAVE.force_contributions[k][i]->position << " " ;
	  _molecules_file << endl ;
	}
      } else if( SAVE.MOLECULES == 1 && SAVE.MOL_VELOCITIES == 1 ) {
	if( SAVE.MOLECULES_IMG == 1 ) {
	  for( int i=0 ; i<numberOfMol ; i++ ) {
	    _molecules_file << molecules[i]->cm.position << " " << molecules[i]->cm.periodic_box << " " << molecules[i]->vcm.position << " " ;
	    for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) _molecules_file << SAVE.force_contributions[k][i]->position << " " ;
	    _molecules_file << endl ;
	  }
	} else {
	  for( int i=0 ; i<numberOfMol ; i++ ) {
	    _molecules_file << molecules[i]->cm.position << " " << molecules[i]->vcm.position << " " ;
	    for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) _molecules_file << SAVE.force_contributions[k][i]->position << " " ;
	    _molecules_file << endl ;
	  }
	}
      }
    } else if( SAVE.MOL_FORCES == 1 ) {
      if( SAVE.MOLECULES == 0 && SAVE.MOL_VELOCITIES == 0 ) for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->fcm.position << endl ;
      else if( SAVE.MOLECULES == 1 && SAVE.MOL_VELOCITIES == 0 ) {
	if( SAVE.MOLECULES_IMG == 1 ) for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << " " << molecules[i]->cm.periodic_box << " " << molecules[i]->fcm.position << endl ;
	else for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << " " << molecules[i]->fcm.position << endl ;
      } else if( SAVE.MOLECULES == 0 && SAVE.MOL_VELOCITIES == 1 ) for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->vcm.position << " " << molecules[i]->fcm.position << endl ;
      else if( SAVE.MOLECULES == 1 && SAVE.MOL_VELOCITIES == 1 ) {
	if( SAVE.MOLECULES_IMG == 1 ) for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << " " << molecules[i]->cm.periodic_box << " " << molecules[i]->vcm.position << " " << molecules[i]->fcm.position << endl ;
	else for( int i=0 ; i<numberOfMol ; i++ ) _molecules_file << molecules[i]->cm.position << " " << molecules[i]->vcm.position << " " << molecules[i]->fcm.position << endl ;
      }
    }
    _molecules_file << endl ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::save_data( int step ) {
  if( SAVE.POTENTIAL_ENERGY == 1 ||
      SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ||
      SAVE.KINETIC_ENERGY == 1 ||
      SAVE.SECONDARY_POTENTIAL == 1 ||
      SAVE.TEMPERATURE == 1 || SAVE.GYR == 1 || SAVE.VOLUME == 1 || SAVE.VIRIAL == 1 ) {

    _thermoinfo_file << step << ' ' << global_time ;
    _thermoinfo_file << setprecision(12) ;
    if( SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ) {
      for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
	if( INTERACTIONS.type[k]->partial_energy.time != global_time ) {
	  INTERACTIONS.type[k]->partial_energy.val = INTERACTIONS.type[k]->compute_energy_contribution() ;
	  INTERACTIONS.type[k]->partial_energy.time = global_time ;
	}
	_thermoinfo_file << ' ' << INTERACTIONS.type[k]->partial_energy.val ;
      }
      if( POTENTIAL_ENERGY.time != global_time ) {
	POTENTIAL_ENERGY.val = 0 ;
	for( int k=0 ; k<INTERACTIONS.number ; k++ ) POTENTIAL_ENERGY.val += INTERACTIONS.type[k]->partial_energy.val ;
	POTENTIAL_ENERGY.time = global_time ;
      }
    }
    if( SAVE.POTENTIAL_ENERGY == 1 ) {
      _thermoinfo_file << ' ' << POTENTIAL_ENERGY.val ;
      if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
	_thermoinfo_file << ' ' << deltaEnergyOfCoM << ' ' << POTENTIAL_ENERGY.val + deltaEnergyOfCoM << ' ' << compute_lambda_derivative_perpart() ;
      }
    }
    if( SAVE.KINETIC_ENERGY == 1 ) {
      _thermoinfo_file << ' ' << compute_kinetic_energy() ;
    }
    if( SAVE.SECONDARY_POTENTIAL == 1 ) {
      _thermoinfo_file << ' ' << calculate_secondary_potential() ;
    }
    if( SAVE.TEMPERATURE == 1 ) {
      _thermoinfo_file << ' ' << ( 2. * compute_kinetic_energy() / DIMENSION / numberOfPart ) ;
    }
    if( SAVE.GYR == 1 ) {
      calculate_molecules_CoM() ;
      _thermoinfo_file << ' ' << calculate_molecules_gyration() ;
    }
    if( SAVE.VOLUME == 1 ) {
      _thermoinfo_file << ' ' << box_sides.position << ' ' << box_sides.volume() ;
    }
    if( SAVE.VIRIAL == 1 ) {
      _thermoinfo_file << ' ' << compute_virial_pressure() ;
    }
    _thermoinfo_file << endl ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::save_runtime_calculations( int step ) {
  if( SAVE.RUNTIME_CALCULATIONS == 1 ) {
    cout << "No runtime calculations have been defined ." << endl ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::clear_runtime_calculations( void ) {
  if( SAVE.RUNTIME_CALCULATIONS == 1 ) {
    cout << "No runtime calculations have been defined ." << endl ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::close_storage_file(void) {
  if( SAVE.PARTICLES == 1 || SAVE.VELOCITIES == 1 || SAVE.FORCES == 1 || SAVE.PER_PART_ENERGY == 1 || SAVE.PER_PART_VIRIAL == 1 ) {
    _particles_file.flush() ;
    _particles_file.close() ;
  }
  if( SAVE.MOLECULES == 1 || SAVE.MOL_VELOCITIES == 1 || SAVE.MOL_FORCES == 1 ) {
    _molecules_file.flush() ;
    _molecules_file.close() ;
  }
  if( SAVE.POTENTIAL_ENERGY == 1 ||
      SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ||
      SAVE.KINETIC_ENERGY == 1 ||
      SAVE.SECONDARY_POTENTIAL == 1 ||
      SAVE.TEMPERATURE == 1 || SAVE.VOLUME == 1 || SAVE.VIRIAL == 1 ) {
    _thermoinfo_file.flush() ;
    _thermoinfo_file.close() ;
  }
  if( SAVE.MSD == 1 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) {
      MSD[i]->_MSD_file.flush() ;
      MSD[i]->_MSD_file.close() ;
      MSD[i]->_saving_steps_file.close() ;
    }
  }
  if( ( SAVE.FORCES == 1 && SAVE.FORCES_SPLIT == 1 ) || ( SAVE.MOL_FORCES == 1 && SAVE.MOL_FORCES_SPLIT == 1 ) ) {
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      for( int i=0 ; i<numberOfPart ; i++ ) delete SAVE.force_contributions[k][i] ;
      delete [] SAVE.force_contributions[k] ;
    }
    if( SAVE.MOL_FORCES == 1 && SAVE.MOL_FORCES_SPLIT == 1 ) for( int k=INTERACTIONS.number ; k<2*INTERACTIONS.number ; k++ ) {
      for( int i=0 ; i<numberOfMol ; i++ ) delete SAVE.force_contributions[k][i] ;
      delete [] SAVE.force_contributions[k] ;
    }
    delete [] SAVE.force_contributions ;
  }
  if( SAVE.PER_PART_ENERGY == 1 ) {
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) delete [] SAVE.per_part_energy[k] ;
    delete [] SAVE.per_part_energy ;
  }
  if( SAVE.PER_PART_VIRIAL == 1 ) {
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) delete [] SAVE.per_part_virial[k] ;
    delete [] SAVE.per_part_virial ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::set_saving_options( const char *MODE ) {
  static save_options OPTIONS ;

  if( strcmp( MODE, "OFF" ) == 0 ) {
    OPTIONS.copy_options( SAVE ) ;
    SAVE.clear() ;
  }else if( strcmp( MODE, "ON" ) == 0 ) {
    SAVE.copy_options( OPTIONS ) ;
  }else if( strcmp( MODE, "EQUILIBRATION" ) == 0 ) {
    SAVE.copy_options( EQ_SAVE_OPT ) ;
  }else if( strcmp( MODE, "PRODUCTION" ) == 0 ) {
    SAVE.copy_options( PROD_SAVE_OPT ) ;
  }else{
    cout << "ERROR: incorrect input for set_saving_options()" << endl;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::save_simulation_conditions( const char *file_name ) {
  int path_length = 0 ;
  while( _DIRECTORY_NAME[path_length] != '\0' ) path_length++ ;
  while( file_name[path_length] != '\0' ) path_length++ ;
  path_length += 200 ;
  char *command = (char *)calloc( path_length , sizeof(char) ) ;
  strcpy( command , _DIRECTORY_NAME ) ;
  strcat( command , "/" ) ;
  strcat( command , file_name ) ;
  ofstream conditions_file ;
  conditions_file.open( command ) ;
  free( command ) ;

  conditions_file << "SIMULATION_MODE " << simulation_mode << endl ;
  conditions_file << "ENSEMBLE " << simulation_ensemble ;
  if( strcmp( simulation_ensemble, "LANGEVIN" ) == 0 || strcmp( simulation_ensemble, "BROWNIAN" ) == 0 ) conditions_file << " , friction: " << langevin_xi ;
  conditions_file << endl ;
  if( strcmp( simulation_ensemble, "NVE" ) == 0 ) {
    conditions_file << "ENERGY_PER_PARTICLE " << energyOfSystem/numberOfPart << endl ;
  } else {
    conditions_file << "KbT " << k_T << endl ;
  }
  conditions_file << "SEED " << seed << endl ;
  conditions_file << "SIGMA " << INTERACTIONS.SIGMA << endl ;
  conditions_file << "EPSILON " << INTERACTIONS.EPSILON << endl ;
  if( particle_mass > 0 ) conditions_file << "PARTICLE_MASS " << particle_mass << endl ;
  else conditions_file << "PARTICLE_MASS either 1 or in the initial configuration file" << endl ;
  conditions_file << "SIDES " << box_sides.position << endl ;
  conditions_file << "INTEGRATION_TIME_STEP " << time_step << endl ;
  conditions_file << "EQUILIBRATION_STEPS " << EQUILIBRATION_STEPS << endl ;
  conditions_file << "PRODUCTION_STEPS " << PRODUCTION_STEPS << endl ;
  conditions_file << "NUMBER_OF_PARTICLES " << numberOfPart << endl ;
  conditions_file << "SAVING_INTERVALS " << SAVE.global_interval << ' ' << SAVE.configuration_interval << ' ' << SAVE.backup_interval << endl ;
  conditions_file << "RUNTIME_CALCULATIONS " << EQ_SAVE_OPT.RUNTIME_CALCULATIONS << ' ' << PROD_SAVE_OPT.RUNTIME_CALCULATIONS << ' ' << SAVE.runtime_interval << endl ;
  if( INITIAL_CONF_BY_FILE == 1 ) conditions_file << "INITIAL_CONFIGURATION_BY_FILE " << configuration_file_name << endl ;
  if( INITIAL_VELOCITIES_BY_FILE == 1 ) conditions_file << "INITIAL_VELOCITIES_BY_FILE " << velocities_file_name << endl ;
  conditions_file << "____ POTENTIALS CONSTANTS " << endl ;
  conditions_file << "SIGMA : " << INTERACTIONS.SIGMA << endl ;
  conditions_file << "EPSILON : " << INTERACTIONS.EPSILON << endl ;
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    conditions_file << "MCS_LENGTH " << MCS_LENGTH << endl ;
    conditions_file << "translation probability, max displacement, max rotation  " << MC_translation_probability << " " << MC_maxdisplacement << " " << MC_maxrotation << endl ;
  }
  conditions_file << INTERACTIONS << endl ;

  conditions_file.close() ;
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::write_backup_data( int current_step ) {
  char file_path[300] ;
  // saving of positions
  strcpy( file_path , _DIRECTORY_NAME ) ;
  strcat( file_path , "/positions_backup.dat" ) ;
  _positions_backup.open( file_path ) ;

  _positions_backup << setprecision(numeric_limits<double>::digits10 + 1) ;
  _positions_backup.seekp( 0 ) ;
  _positions_backup << "# step " << current_step << endl ;
  _positions_backup << "# N " << numberOfPart << endl ;
  _positions_backup << "# box " << box_sides.position << endl ;
  if ( strcmp( SAVE.configuration_format , "ptc" ) == 0 ||
       strcmp( SAVE.configuration_format , "patch" ) == 0 ) {
    _positions_backup << "# format patch" << endl ;
    _positions_backup << "& " << numberOfPart << endl ;
    _positions_backup << setprecision(numeric_limits<double>::digits10 + 1) << box_sides.position << endl ;
    double *patch_parameters = NULL ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) if( strcmp( INTERACTIONS.type[k]->name , "kernfrenkel" ) == 0 ) {
	INTERACTIONS.type[k]->build_bonds_network() ;
	patch_parameters = INTERACTIONS.type[k]->dump_parameters() ;
      }
    if( patch_parameters == NULL ) {
      cout << endl << "*** ERROR: cannot dump a patch backup configuration without a patchy potential !" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    for( int i=0 ; i<numberOfPart ; i++ ) _positions_backup << particles[i]->dump_position_variables( numeric_limits<double>::digits10 + 1 , "patch" , patch_parameters ) << endl ;
    delete [] patch_parameters ;

  } else if ( strcmp( SAVE.configuration_format , "sph" ) == 0 ) {
    _positions_backup << "# format sph" << endl ;
    _positions_backup << numberOfPart << " " << global_time << endl ;
    _positions_backup << setprecision(numeric_limits<double>::digits10 + 1) << box_sides.position << endl ;
    if( DIMENSION ==3 ) for( int i=0 ; i<numberOfPart ; i++ ) _positions_backup << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( numeric_limits<double>::digits10 + 1 , "sph" , NULL ) << " " << particles[i]->radius << endl ;
    else if( DIMENSION == 2 ) for( int i=0 ; i<numberOfPart ; i++ ) _positions_backup << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( numeric_limits<double>::digits10 + 1 , "sph" , NULL ) << " 0.0 " << particles[i]->radius << endl ;

  } else {
    for(int i=0; i<numberOfPart; i++) {
      _positions_backup << particles[i]->dump_position_variables( numeric_limits<double>::digits10 + 1 , "all" ) << " " << particles[i]->periodic_box << endl ;
    }
  }

  // saving of velocities
  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    if( velocities != NULL ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      strcat( file_path , "/velocities_backup.dat" ) ;
      _velocities_backup.open( file_path ) ;

      _velocities_backup << setprecision(numeric_limits<double>::digits10 + 1) ;
      _velocities_backup.seekp( 0 ) ;
      _velocities_backup << "# N " << numberOfPart << endl ;
      _velocities_backup << "# step " << current_step << endl ;
      for(int i=0; i<numberOfPart; i++) {
	_velocities_backup << velocities[i]->position << endl ;
      }
      _velocities_backup << endl ;
      _velocities_backup.close() ;
    } else {
      cout << " *** ERROR: The velocities vector is not allocated, but simulation mode is MD !" << endl ;
      exit( EXIT_FAILURE ) ;
    }
  }

  // saving of information on bonds list and charges
  for( int i=0; i<numberOfPart; i++ ) {
    _positions_backup << i+1 << ' ' << particles[i]->valence << ' ' << particles[i]->charge << endl ;
    if( particles[i]->valence > 0 ) {
      for( int j=0 ; j<particles[i]->valence-1 ; j++ )  _positions_backup << particles[i]->bonded_monomer[j] + 1 << ' ' ;
      _positions_backup << particles[i]->bonded_monomer[ particles[i]->valence-1 ] + 1 << endl ;
    }
  }
  _positions_backup << endl ;

  _positions_backup.close() ;
}
/****************************************************************************************/
