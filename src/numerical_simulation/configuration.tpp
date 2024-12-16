#include "interactions.h"
#include "msd.h"

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
template <typename particle>
configuration<particle>::interactions_array::interactions_array() {
  number = 0 ;
  type = NULL ;
  SIGMA = 0.0 , EPSILON = 0.0 ;
}
/********************************************************************************************************************/

template <typename particle>
configuration<particle>::interactions_array::~interactions_array() {
  for( int i=0 ; i<number ; i++ ) delete type[i] ;
  free( type ) ;
  type = NULL ;
  number = 0 ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::interactions_array::operator=( const interactions_array &other ) {
  number = other.number ;
  type = (interaction<particle> **)calloc( number , sizeof(interaction<particle> *) ) ;
  for(int i=0; i<number; i++) {
    type[i] = other.type[i]->copy() ;
    type[i]->cutoff = other.type[i]->cutoff ;
    type[i]->SHIFT = other.type[i]->SHIFT ;
    type[i]->verlet_list = new verlet<particle>( *( other.type[i]->verlet_list ) ) ;     // copy of verlet list
  }
}
/********************************************************************************************************************/

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
template <typename particle>
configuration<particle>::configuration() {
  strcpy( simulation_ensemble, "0" ) ;
  strcpy( simulation_mode , "00" ) ;
  seed = -1 ;
  DIMENSION = -1 ;
  k_T = 1.0 ;
  NPT_pressure = 1.0 ;
  Lambda3 = 1.0 ;
  fugacity = 1.0 ;
  energyOfSystem = 0.0 ;
  time_step = 1 ;   // The time-step is initialized to 1
  global_time = 0 ;
  EQUILIBRATION_STEPS = 0 ;
  PRODUCTION_STEPS = 0 ;
  MCS_LENGTH = 0 ;
  MC_translation_probability = 0.5 ;
  box_side = 0 ;
  box_sides = particle() ;
  midside = particle() ;
  MC_maxdisplacement = 0.3 ;
  MC_maxboxdisp = 0.3 ;
  MC_maxrotation = 0.2 ;
  MC_rotationjump = 2 * M_PI ;
  // default options saves only configurations and energy
  particle_mass = 0.0 ;
  SAVE.clear() ;
  SAVE.PARTICLES = 1 ;
  SAVE.PARTICLES_IMG = 1 ;
  SAVE.POTENTIAL_ENERGY = 1 ;
  SAVE.configuration_format = new char [4] ;
  sprintf( SAVE.configuration_format , "%s" , "xyz" ) ;
  NHC_thermostat = NULL ;
  vcm_reset_interval = 1000 ;
  numberOfPart = 0 ;
  particles_length = 0 ;
  particles = NULL ;
  velocities = NULL ;
  forces_t1 = NULL ;
  forces_t2 = NULL ;
  dr_gauss = NULL ;
  dv_gauss = NULL ;
  numberOfMol = 0 ;
  molecules_length = 0 ;
  molecules = NULL ;
  compute_MOLcm = 0 ;
  compute_MOLgyr = 0 ;
  monomers_number = 0 ;
  charges_number = 0 ;
  strcpy( _DIRECTORY_NAME , "" ) ;
  inertia_matrix.rebuild( 3 , 3 ) ;
  N_MSDs = 0 ;
  MSD = NULL ;
}
/********************************************************************************************************************/

template <typename particle>
configuration<particle>::configuration( const configuration &other_config ) {
  SAVE.copy_options( other_config.SAVE ) ;
  SAVE.global_interval = other_config.SAVE.global_interval ;
  SAVE.configuration_interval = other_config.SAVE.configuration_interval ;
  SAVE.backup_interval = other_config.SAVE.backup_interval ;
  SAVE.runtime_interval = other_config.SAVE.runtime_interval ;
  sprintf( SAVE.configuration_format , "%s" , other_config.SAVE.configuration_format ) ;
  INTERACTIONS.SIGMA = other_config.INTERACTIONS.SIGMA ;
  INTERACTIONS.EPSILON = other_config.INTERACTIONS.EPSILON ;
  box_side = other_config.box_side ;
  box_sides = other_config.box_sides ;
  midside = other_config.midside ;
  fugacity = other_config.fugacity ;
  k_T = other_config.k_T ;
  energyOfSystem = other_config.energyOfSystem ;
  POTENTIAL_ENERGY = other_config.POTENTIAL_ENERGY ;
  KINETIC_ENERGY = other_config.KINETIC_ENERGY ;
  Lambda3 = other_config.Lambda3 ;
  MC_maxdisplacement = other_config.MC_maxdisplacement ;
  MC_maxboxdisp = other_config.MC_maxboxdisp ;
  MC_maxrotation = other_config.MC_maxrotation ;
  MC_rotationjump = other_config.MC_rotationjump ;
  time_step = other_config.time_step ;
  global_time = other_config.global_time ;
  particle_mass = other_config.particle_mass ;
  EQUILIBRATION_STEPS = other_config.EQUILIBRATION_STEPS ;
  PRODUCTION_STEPS = other_config.PRODUCTION_STEPS ;
  numberOfPart = other_config.numberOfPart ;
  monomers_number = other_config.monomers_number ;
  charges_number = other_config.charges_number ;
  NPT_pressure = other_config.NPT_pressure ;
  virial_press = other_config.virial_press ;
  INTERACTIONS = other_config.INTERACTIONS ; // assignment by copy
  N_MSDs = 0 ;
  MSD = NULL ;

  particles_length = numberOfPart ;
  particles = (particle **)calloc( particles_length , sizeof(particle *) ) ;
  if( particles == NULL ) {
    cout << "\n Allocation error 1 in configuration constructor" ;
    exit(EXIT_FAILURE) ;
  }
  forces_t1 = (particle **)calloc( particles_length , sizeof(particle *) ) ;
  if( forces_t1 == NULL ) {
    cout << "\n Allocation error 1.4 in configuration constructor" ;
    exit( EXIT_FAILURE ) ;
  }
  forces_t2 = (particle **)calloc( particles_length , sizeof(particle *) ) ;
  if( forces_t2 == NULL ) {
    cout << "\n Allocation error 1.5 in configuration constructor" ;
    exit( EXIT_FAILURE ) ;
  }
  // copy of the list of molecules
  numberOfMol = other_config.numberOfMol ;
  molecules_length = numberOfMol ;
  molecules = (molecule<particle> **)calloc( molecules_length , sizeof(molecule<particle> *) ) ;
  if( molecules == NULL ) {
    cout << "\n Allocation error 1.7 in configuration constructor" ;
    exit(EXIT_FAILURE) ;
  }
  for( int i=0 ; i<numberOfMol ; i++ ) {
    molecules[i] = new molecule<particle>( *other_config.molecules[i] ) ;
    for( int j=0 ; j<molecules[i]->Natoms ; j++ )  molecules[i]->atom[j] = particles[ other_config.molecules[i]->atom[j]->id ] ;
  }
  // copy of positions and forces
  for( int i=0 ; i<numberOfPart ; i++ ) {
    particles[i] = new particle ;
    *( particles[i] ) = *( other_config.particles[i] ) ;
    if( other_config.forces_t1[i] != NULL ) {
      forces_t1[i] = new particle ;
      *( forces_t1[i] ) = *( other_config.forces_t1[i] ) ;
    }
    if( other_config.forces_t2[i] != NULL ) {
      forces_t2[i] = new particle ;
      *( forces_t2[i] ) = *( other_config.forces_t2[i] ) ;
    }
  }
  // copy of velocities
  if( other_config.velocities != NULL ) {
    velocities = (particle **)calloc( particles_length , sizeof(particle *) ) ;
    if( velocities == NULL ) {
      cout << "\n Allocation error 1.3 in configuration constructor" ;
      exit( EXIT_FAILURE ) ;
    }
    for( int i=0 ; i<numberOfPart ; i++ ) {
      velocities[i] = new particle ;
      *( velocities[i] ) = *( other_config.velocities[i] ) ;
    }
  }

  // copy of Nose Hoover thermostat
  if( other_config.NHC_thermostat != NULL ) {
    NHC_thermostat = new NoseHooverThermostat( *(other_config.NHC_thermostat) ) ;
  }else{
    NHC_thermostat = NULL ;
  }

}
/********************************************************************************************************************/

template <typename particle>
configuration< particle >::~configuration() {
  // deletion of the array of interactions and verlet lists is automatically done within the destructor of INTERACTIONS

          // deletion of the arrays of particles, velocities and forces
  for( int i=0; i<particles_length; i++ ) {
    delete particles[i] ;
    particles[i] = NULL ;
  }
  free( particles ) ;
  particles = NULL ;
  if( velocities != NULL ) {
    for( int i=0 ; i<numberOfPart ; i++ ) {
      delete velocities[i] ;
      velocities[i] = NULL ;
      delete forces_t1[i] ;
      forces_t1[i] = NULL ;
      delete forces_t2[i] ;
      forces_t2[i] = NULL ;
    }
    free( velocities ) ;
    velocities = NULL ;
    free( forces_t1 ) ;
    forces_t1 = NULL ;
    free( forces_t2 ) ;
    forces_t2 = NULL ;
  }
  if( dr_gauss != NULL ) {
    for( int i=0 ; i<particles_length ; i++ ) delete dr_gauss[i] ;
    delete [] dr_gauss ;
  }
  if( dv_gauss != NULL ) {
    for( int i=0 ; i<particles_length ; i++ ) delete dv_gauss[i] ;
    delete [] dv_gauss ;
  }

  particles_length = 0 ;
  numberOfPart = 0 ;
  monomers_number = 0 ;
  charges_number = 0 ;

          // deletion of molecules, if present
  if( numberOfMol > 0 ) {
    for( int i=0; i<numberOfMol; i++ ) {
      delete molecules[i] ;
      molecules[i] = NULL ;
    }
    free( molecules ) ;
    molecules = NULL ;
    molecules_length = 0 ;
    numberOfMol = 0 ;
  }

          // deletion of the MSD structures
  if( N_MSDs != 0 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) delete MSD[i] ;
    delete MSD ;
    MSD = NULL ;
    N_MSDs = 0 ;
  }

          // deletion of Nose Hoover thermostat
  delete NHC_thermostat ;
}
/****************************************************************************************/
/****************************************************************************************/

/***********************    BLOCK 1   ***************************************************/
/****************************************************************************************/
template <typename particle>
void configuration<particle>::initialize_configuration( const char *input_file_name ) {
  ifstream config_file ;
  if( input_file_name == NULL ) {
    config_file.open( "input_file.dat" ) ;
  }else{
    config_file.open( input_file_name ) ;
  }
  if( config_file.fail() ) {
    cout << "\n  MD_config_file non trovato, apertura fallita" << endl ;
    exit( EXIT_FAILURE ) ;

  } else {
    interpret_input_script( config_file ) ;
    config_file.close() ;
  }

  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    if( INITIAL_VELOCITIES_BY_FILE == 1 ) generate_velocities_by_file() ;
    // generation of the array containing forces
    forces_t1 = new particle* [ particles_length ] ;
    check_memalloc( forces_t1 , "\n Allocation error for forces_t1 in configuration constructor" ) ;
    forces_t2 = new particle* [ particles_length ] ;
    check_memalloc( forces_t2 , "\n Allocation error for forces_t2 in configuration constructor" ) ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      forces_t1[i] = new particle ;
      forces_t2[i] = new particle ;
    }

    if( strcmp( simulation_ensemble , "LANGEVIN" ) == 0 ) {
      // Allocation of the arrays containing the gaussian vectors needed to evolve the system with a Langevin integrator
      dv_gauss = new particle* [ particles_length ] ;
      check_memalloc( dv_gauss , "\n Allocation error for dv_gauss in configuration constructor" ) ;
      dr_gauss = new particle* [ particles_length ] ;
      check_memalloc( dr_gauss , "\n Allocation error for dr_gauss in configuration constructor" ) ;
      for( int i=0 ; i<numberOfPart ; i++ ) {
	dv_gauss[i] = new particle ;
	dr_gauss[i] = new particle ;
      }
    }
  }

  if( compute_MOLcm == 1 ) calculate_molecules_CoM() ;  // This is needed by molecular interactions, if involved (such as Hertzian potential)
  compute_potential_energy() ;           // initial potential energy calculation

  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    calculate_forces_t2() ;                       // initial forces are computed and stored in forces_t1 to start velocity verlet
    support_pp2particle = forces_t1 ;
    forces_t1 = forces_t2 ;
    forces_t2 = support_pp2particle ;
  }

  // creation of the folder containing all saved data
  if( generate_new_directory_for_data() == NULL ) cout << " Problems in configuration constructor\n" ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::interpret_input_script( ifstream &config_file ) {
  double charges_fraction = 0, charges_distribution_sigma = 0 , density = -1.0 ;
  int N_patches = 0 ;
  char charge_distribution[20] ;
  record *file_rec = NULL ;
  file_rec = new record ;
  INITIAL_VELOCITIES_BY_FILE = 0 ;

  while( file_rec->getrecord( config_file ) ) {
    file_rec->split() ;
    if( file_rec->words_number > 0 ) {
      if( file_rec->word[0][0] != '#' ) {

	if( strcmp( file_rec->word[0] , "DIMENSION" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &DIMENSION ) ;
	else if( strcmp( file_rec->word[0] , "PARTICLE_CLASS") == 0 ) {
	  if( strcmp( file_rec->word[1] , "patchy_2D") == 0 )  sscanf( file_rec->word[2] , "%d" , &N_patches ) ;
	}
	else if( strcmp( file_rec->word[0] , "ENSEMBLE" ) == 0 ) {
	  strcpy( simulation_ensemble , file_rec->word[1] ) ;
	  if( strcmp( simulation_ensemble , "NVE" ) != 0 && strcmp( simulation_ensemble , "NVT_NHC" ) != 0 &&
	      strcmp( simulation_ensemble , "NVT_NHE" ) != 0 && strcmp( simulation_ensemble , "LANGEVIN" ) != 0 &&
	      strcmp( simulation_ensemble , "BROWNIAN" ) != 0 && strcmp( simulation_ensemble , "MC_SWAP_CHARGE" ) != 0 &&
	      strcmp( simulation_ensemble , "MC_NVT" ) != 0 && strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) != 0 &&
	      strcmp( simulation_ensemble , "MC_NPT" ) != 0 ) {
	    cout << " *** ERROR : The simulation ensemble " << simulation_ensemble << " does not exist (yet) !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }
	  strcpy( simulation_mode , "MD" ) ;
	  if( strcmp( simulation_ensemble , "BROWNIAN" ) == 0 || strcmp( simulation_ensemble , "LANGEVIN" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%lf" , &langevin_xi ) ;
	  } else if( strcmp( simulation_ensemble , "MC_SWAP_CHARGE" ) == 0  ||
		     strcmp( simulation_ensemble , "MC_NVT" ) == 0  ||
		     strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0  ||
		     strcmp( simulation_ensemble , "MC_NPT" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &MCS_LENGTH ) ;
	    strcpy( simulation_mode , "MC" ) ;
	  }
	}

	else if( strcmp( file_rec->word[0] , "MC_PARAMETER" ) == 0 ) {
	  if( strcmp( file_rec->word[1] , "translational_probability" ) == 0 ) sscanf( file_rec->word[2] , "%lf" , &MC_translation_probability ) ;
	  else if( strcmp( file_rec->word[1] , "max_translation" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%lf" , &MC_maxdisplacement ) ;
	    MC_PARAMS_MAX.MC_maxdisplacement = MC_maxdisplacement ;
	    MC_PARAMS_MIN.MC_maxdisplacement = MC_maxdisplacement ;
	    if( file_rec->words_number > 3 ) {
	      if( strcmp( file_rec->word[3] , "AUTO_ADJUST" ) == 0 ) {
		sscanf( file_rec->word[4] , "%lf" , &(MC_PARAMS_MIN.MC_maxdisplacement) ) ;
		sscanf( file_rec->word[5] , "%lf" , &(MC_PARAMS_MAX.MC_maxdisplacement) ) ;
	      }
	    }
	  } else if( strcmp( file_rec->word[1] , "max_rotation" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%lf" , &MC_maxrotation ) ;
	    MC_PARAMS_MAX.MC_maxrotation = MC_maxrotation ;
	    MC_PARAMS_MIN.MC_maxrotation = MC_maxrotation ;
	    int symfactor = 1 ;
	    for( int w=3; w<file_rec->words_number; w++ ) {
	      if( strcmp( file_rec->word[w] , "AUTO_ADJUST" ) == 0 ) {
		sscanf( file_rec->word[w+1] , "%lf" , &(MC_PARAMS_MIN.MC_maxrotation) ) ;
		sscanf( file_rec->word[w+2] , "%lf" , &(MC_PARAMS_MAX.MC_maxrotation) ) ;
	      } else if( strcmp( file_rec->word[w] , "JUMP" ) == 0 ) sscanf( file_rec->word[w+1] , "%d" , &(symfactor) ) ;
	    }
	    MC_rotationjump = 2 * M_PI / symfactor ;
	  } else if( strcmp( file_rec->word[1] , "max_boxside_change" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%lf" , &MC_maxboxdisp ) ;
	    MC_PARAMS_MAX.MC_maxboxdisp = MC_maxboxdisp ;
	    MC_PARAMS_MIN.MC_maxboxdisp = MC_maxboxdisp ;
	    if( file_rec->words_number > 3 ) {
	      if( strcmp( file_rec->word[3] , "AUTO_ADJUST" ) == 0 ) {
		sscanf( file_rec->word[4] , "%lf" , &(MC_PARAMS_MIN.MC_maxboxdisp) ) ;
		sscanf( file_rec->word[5] , "%lf" , &(MC_PARAMS_MAX.MC_maxboxdisp) ) ;
	      }
	    }
	  }
	}

	else if( strcmp( file_rec->word[0] , "SIGMA" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &INTERACTIONS.SIGMA ) ;
	else if( strcmp( file_rec->word[0] , "EPSILON" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &INTERACTIONS.EPSILON ) ;
	else if( strcmp( file_rec->word[0] , "PARTICLE_MASS" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &particle_mass ) ;
	else if( strcmp( file_rec->word[0] , "SIDE") == 0 ) {
	  sscanf( file_rec->word[1] , "%lf" , &box_side ) ;
	  box_sides.set_equal_comp( box_side ) ;
	  midside.position = (box_sides.position * 0.5) ;
	}
	else if( strcmp( file_rec->word[0], "SIDES") == 0 ) {
	  box_sides.read_by_string( file_rec->word+1 ) ;
	  midside.position = (box_sides.position * 0.5) ;
	  if( box_sides.min_comp() == box_sides.max_comp() ) {
	    box_side = box_sides.min_comp() ;
	  }else{
	    box_side = nanl("") ;
	  }
	}
	else if( strcmp( file_rec->word[0], "KbT" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &k_T ) ;
	else if( strcmp( file_rec->word[0], "PRESSURE" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &NPT_pressure ) ;
	else if( strcmp( file_rec->word[0], "DENSITY" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &density ) ;

	else if( strcmp( file_rec->word[0], "INITIAL_CONFIGURATION" ) == 0 ) {
	  if( strcmp( file_rec->word[1], "file" ) == 0 ) {
	    sscanf( file_rec->word[3] , "%s" , configuration_file_name ) ;
	              // positions at time 0 are generated, giving the file name (configuration file) and the file format
	    generate_config_by_file( configuration_file_name , file_rec->word[2] ) ;
	    if( strcmp( file_rec->word[2] , "xyz" ) == 0 ||
		strcmp( file_rec->word[2] , "nico_bonds" ) == 0 ||
		strcmp( file_rec->word[2] , "nico_bonds_pbc" ) == 0 ) {
	      if( file_rec->words_number > 4 ) {
		sscanf( file_rec->word[4] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
	      }else{
		cout << "ERROR: INITIAL_CONFIGURATION needs a value for the initial energy per particle !" << endl ;
		exit( EXIT_FAILURE ) ;
	      }
	    }
	    if( strcmp( file_rec->word[2] , "nico_rings" ) == 0 ||
		strcmp( file_rec->word[2] , "nico_rings_pbc" ) == 0 ) {
	      if( file_rec->words_number > 4 ) {
	              // molecules are generated, giving the respective file name (polydispersity file)
		sscanf( file_rec->word[4] , "%s" , configuration_file_name ) ;
		generate_molecules_by_file( configuration_file_name ) ;
		if( file_rec->words_number > 5 ) {// here I check if the total energy per particle has been inputed
		  sscanf( file_rec->word[5] , "%lf" , &energyOfSystem ) ;
		  energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
		}else{
		  cout << "ERROR: INITIAL_CONFIGURATION needs a value for the initial energy per particle !" << endl ;
		  exit( EXIT_FAILURE ) ;
		}
	      }else{
		cout << "ERROR: INITIAL_CONFIGURATION needs a polydispersity file !" << endl ;
		exit( EXIT_FAILURE ) ;
	      }
	    }
	  } else {
	    if( file_rec->words_number > 3 ) {
	      // generation of positions' array
	      sscanf( file_rec->word[2] , "%d" , &numberOfPart ) ;
	      particles_length = numberOfPart ;
	      particles = (particle **)calloc( particles_length , sizeof(particle *) ) ;
	      check_memalloc( particles , "\n Allocation error 1 during the generation of the initial configuration" ) ;
	      if( box_side == 0 && box_sides.volume() == 0 ) {
		if( density > 0 ) {
		  box_side = pow( (double)numberOfPart/density , 1./DIMENSION ) ;
		  box_sides.set_equal_comp( box_side ) ;
		  midside.position = (box_sides.position * 0.5) ;
		} else {
		  cout << "*** ERROR (initial configuration): Please, define the box size or the particles density!" << endl ;
		  exit( EXIT_FAILURE ) ;
		}
	      }
	      if( strcmp( file_rec->word[1] , "random" ) == 0 ) {
		sscanf( file_rec->word[3] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
		generate_ControlledEnergy_randomConfig() ;
	      } else if( strcmp( file_rec->word[1] , "fcc" ) == 0 || strcmp( file_rec->word[1] , "bcc" ) == 0 ) {
		sscanf( file_rec->word[3] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
		generate_ideal_lattice_3D( file_rec->word[1] ) ;
	      } else if( strcmp( file_rec->word[1] , "hardspheres" ) == 0 ) {
		if( file_rec->words_number < 5 ) { cout << "*** ERROR : invalid number of options for INITIAL_CONFIGURATION" << endl ; exit(EXIT_FAILURE) ; }
		sscanf( file_rec->word[3] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
		double radius ;
		sscanf( file_rec->word[4] , "%lf" , &radius ) ;
		generate_HS_randomConfig( radius ) ;
	      } else if( strcmp( file_rec->word[1] , "lattice" ) == 0 ) {
		generate_ideal_lattice( file_rec->word[3] ) ;
	      } else {
		cout << endl << " You need to specify the configuration arrangement in the initialization in function main!" << endl;
	      }
	    } else {
	      cout << "ERROR: INITIAL_CONFIGURATION needs the values for the initial energy per particle and number of particles !" << endl ;
	      exit( EXIT_FAILURE ) ;
	    }
	  }
	}

	else if( strcmp( file_rec->word[0], "LENNARD_JONES" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new lennard_jones<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "SOFT" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new soft<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "TEST_POTENTIAL" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new test_potential<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "HARMONIC_SPRING_FIELD" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  } else {
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new harmonic_spring_field<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "HARMONIC_SPRING_FIELD_WIGNERSEITZ" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  } else {
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new harmonic_spring_field_wignerseitz<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "ANGULAR_HARMONIC_SPRING_FIELD" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  } else {
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new angular_h_spring_field<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "ANGULAR_COSINE_SPRING_FIELD" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  } else {
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new angular_c_spring_field<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "HERTZIAN" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << numberOfPart << " " << box_side << " " << box_sides.position << endl ;
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    double u = 0 , sh = 0 ;
	    sscanf( file_rec->word[1] , "%lf" , &u ) ;
	    sscanf( file_rec->word[2] , "%lf" , &sh ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1] = new hertzian_molecule<particle>( u , sh ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , NULL , this , 0 ) ; // In this case the list of particles is not created
	    // it has been put in the hertzian initialization	      compute_MOLcm = 1 ;  // to enable the computation of the molecules center of mass at every time step
	  }
	}

	else if( strcmp( file_rec->word[0], "CONSERVATIVE_HERTZIAN" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << numberOfPart << " " << box_side << " " << box_sides.position << endl ;
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition
	    double u = 0 , sh = 0 ;
	    sscanf( file_rec->word[1] , "%lf" , &u ) ;
	    sscanf( file_rec->word[2] , "%lf" , &sh ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1] = new conservative_hertzian<particle>( u , sh ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , NULL , this , 0 ) ; // In this case the list of particles is not created
	      // it has been put in the hertzian initialization	      compute_MOLcm = 1 ;  // to enable the computation of the molecules center of mass at every time step
	  }
	}

	else if( strcmp( file_rec->word[0], "LENNARD_JONES_POLY" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new lennard_jones_poly<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "KERN_FRENKEL" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new kern_frenkel<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "HARD_SPHERE" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new hard_sphere<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "HS_SQUARE_WELL" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new hs_square_well<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "FENE_POLY" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    double kf , ext_factor ;
	    sscanf( file_rec->word[1] , "%lf" , &kf ) ;
	    sscanf( file_rec->word[2] , "%lf" , &ext_factor ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new fene_poly<particle>( kf , ext_factor ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , monomers_number ) ;
	              // I store the pointer to the list in an array to make consistency checks at the end of simulations
	    CHECKS.bonds_lists_number ++ ;
	    CHECKS.bonds_verlet_link = (verlet<particle> **)realloc( CHECKS.bonds_verlet_link , CHECKS.bonds_lists_number*sizeof(verlet<particle> *) ) ;
	    CHECKS.bonds_verlet_link[ CHECKS.bonds_lists_number-1 ] = INTERACTIONS.type[INTERACTIONS.number-1]->verlet_list ;     // add list the checks
	  }
	}

	else if( strcmp( file_rec->word[0], "WCA_ALPHA" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	              // pair potential definition
	    double alpha_value ;
	    sscanf( file_rec->word[2] , "%lf" , &alpha_value ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1] = new wca_alpha<particle>( alpha_value ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "FENE" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    double kf , r0 ;
	    sscanf( file_rec->word[1] , "%lf" , &kf ) ;
	    sscanf( file_rec->word[2] , "%lf" , &r0 ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new fene<particle>( kf , r0 ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	         // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , monomers_number ) ;
	              // I store the pointer to the list in an array to make consistency checks at the end of simulations
	    CHECKS.bonds_lists_number ++ ;
	    CHECKS.bonds_verlet_link = (verlet<particle> **)realloc( CHECKS.bonds_verlet_link , CHECKS.bonds_lists_number*sizeof(verlet<particle> *) ) ;
	    CHECKS.bonds_verlet_link[ CHECKS.bonds_lists_number-1 ] = INTERACTIONS.type[INTERACTIONS.number-1]->verlet_list ;     // add list the checks
	  }
	}

	else if( strcmp( file_rec->word[0], "CHARGES" ) == 0 ) {
	  strcpy( charge_distribution, file_rec->word[1] ) ;
	  if( strcmp( charge_distribution, "random" ) == 0 ) sscanf( file_rec->word[2] , "%lf" , &charges_fraction ) ;
	  else if( strcmp( charge_distribution, "gauss"  ) == 0  ||  strcmp( charge_distribution, "surface"  ) == 0  ||  strcmp( charge_distribution, "radial"  ) == 0 ) {
	    sscanf( file_rec->word[2] , "%lf" , &charges_fraction ) ;
	    sscanf( file_rec->word[3] , "%lf" , &charges_distribution_sigma ) ;
	  }
	  generate_charge_distribution( charge_distribution, charges_fraction , charges_distribution_sigma ) ;
	}

	else if( strcmp( file_rec->word[0], "DEBYE_HUCKEL" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    if( charges_number == 0 ) cout << " Warning :  You are using the Debye-Huckel potential with no charges in the system." << endl ;
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    double bl , dl ;
	    sscanf( file_rec->word[1] , "%lf" , &bl ) ;
	    sscanf( file_rec->word[2] , "%lf" , &dl ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new debye_huckel<particle>( k_T , bl , dl ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , charges_number ) ;
	    initializeVlist_runtime_calculations( INTERACTIONS.type[INTERACTIONS.number-1]->verlet_list ) ;  // to allow computation of charge profiles
	  }
	}

	else if( strcmp( file_rec->word[0], "DEBYE_HUCKEL_POLY" ) == 0 ) {
	  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
	    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }else{
	    if( charges_number == 0 ) cout << " Warning :  You are using the Debye-Huckel potential with no charges in the system." << endl ;
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    double bl , dl ;
	    sscanf( file_rec->word[1] , "%lf" , &bl ) ;
	    sscanf( file_rec->word[2] , "%lf" , &dl ) ;
	              // pair potential definition
	    INTERACTIONS.type[INTERACTIONS.number-1] = new debye_huckel_poly<particle>( k_T , bl , dl ) ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , charges_number ) ;
	    initializeVlist_runtime_calculations( INTERACTIONS.type[INTERACTIONS.number-1]->verlet_list ) ;  // to allow computation of charge profiles
	  }
	}

	else if( strcmp( file_rec->word[0], "SECONDARY_POTENTIAL" ) == 0 ) {
	  cout << "ERROR: You still have to write the part of code about the computation of the secondary potential . . ." << endl ;
	  exit( EXIT_FAILURE ) ;
	    //	    sscanf( file_rec->word[1] , "%lf" , &SECONDARY_POTENTIAL_SHIFT ) ;
	}
	else if( strcmp( file_rec->word[0], "INTEGRATION_TIME_STEP" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &time_step ) ;
	else if( strcmp( file_rec->word[0], "EQUILIBRATION_STEPS" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &EQUILIBRATION_STEPS ) ;
	else if( strcmp( file_rec->word[0], "PRODUCTION_STEPS" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &PRODUCTION_STEPS ) ;
	else if( strcmp( file_rec->word[0], "CM_RESET_INTERVAL" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &vcm_reset_interval ) ;
	else if( strcmp( file_rec->word[0], "SAVING_INTERVALS" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &SAVE.global_interval ) ;
	  sscanf( file_rec->word[2] , "%d" , &SAVE.configuration_interval ) ;
	  sscanf( file_rec->word[3] , "%d" , &SAVE.backup_interval ) ;
	  SAVE.runtime_interval = INT_MAX ; // Questa cosa è da modificare, la classe runtime deve stare in conf e non in mgel...
	}

	else if( strcmp( file_rec->word[0], "SAVE_FORMAT" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%s" , SAVE.configuration_format ) ;
	  if( strcmp( SAVE.configuration_format , "xyz" ) != 0 &&
	      strcmp( SAVE.configuration_format , "ptc" ) != 0 &&
	      strcmp( SAVE.configuration_format , "patch" ) != 0 &&
	      strcmp( SAVE.configuration_format , "sph" ) != 0 ) {
	    cout << "*** ERROR : Unrecognized output format !" << endl ;
	    exit(EXIT_FAILURE) ;
	  }
	}
	else if( strcmp( file_rec->word[0], "SAVE_ATTRIBUTE" ) == 0 ) {
	  if( strcmp( file_rec->word[1], "particles" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES ) ;
	  } else if( strcmp( file_rec->word[1], "part_img" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES_IMG ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES_IMG ) ;
	  } else if( strcmp( file_rec->word[1], "velocities" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.VELOCITIES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.VELOCITIES ) ;
	  } else if( strcmp( file_rec->word[1], "forces" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.FORCES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.FORCES ) ;
	  } else if( strcmp( file_rec->word[1], "forces_splitted" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.FORCES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.FORCES ) ;
	    EQ_SAVE_OPT.FORCES_SPLIT = EQ_SAVE_OPT.FORCES ;
	    PROD_SAVE_OPT.FORCES_SPLIT = PROD_SAVE_OPT.FORCES ;
	  } else if( strcmp( file_rec->word[1], "mols_com" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.MOLECULES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.MOLECULES ) ;
	  } else if( strcmp( file_rec->word[1], "mols_com_img" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.MOLECULES_IMG ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.MOLECULES_IMG ) ;
	  } else if( strcmp( file_rec->word[1], "mols_velocities" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.MOL_VELOCITIES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.MOL_VELOCITIES ) ;
	  } else if( strcmp( file_rec->word[1], "mols_forces" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.MOL_FORCES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.MOL_FORCES ) ;
	  } else if( strcmp( file_rec->word[1], "mols_forces_splitted" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.MOL_FORCES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.MOL_FORCES ) ;
	    EQ_SAVE_OPT.MOL_FORCES_SPLIT = EQ_SAVE_OPT.MOL_FORCES ;
	    PROD_SAVE_OPT.MOL_FORCES_SPLIT = PROD_SAVE_OPT.MOL_FORCES ;
	  } else if( strcmp( file_rec->word[1], "energy" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PER_PART_ENERGY ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PER_PART_ENERGY ) ;
	  } else if( strcmp( file_rec->word[1], "virial" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PER_PART_VIRIAL ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PER_PART_VIRIAL ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "SAVE_POTENTIAL_ENERGY" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.POTENTIAL_ENERGY ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.POTENTIAL_ENERGY ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_PARTIAL_POTENTIAL_ENERGIES" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.PARTIAL_POTENTIAL_ENERGY ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.PARTIAL_POTENTIAL_ENERGY ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_KINETIC_ENERGY" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.KINETIC_ENERGY ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.KINETIC_ENERGY ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_SECONDARY_POTENTIAL" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.SECONDARY_POTENTIAL ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.SECONDARY_POTENTIAL ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_TEMPERATURE" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.TEMPERATURE ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.TEMPERATURE ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_VOLUME" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.VOLUME ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.VOLUME ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_MSD" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.MSD ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.MSD ) ;
	  if( file_rec->words_number == 4 ) add_initialize_MSD( file_rec->word[3] ) ;
	  else add_initialize_MSD( file_rec->word[3] , file_rec->word[4] ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_VIRIAL_PRESSURE" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.VIRIAL ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.VIRIAL ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_GYR" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.GYR ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.GYR ) ;
	}
	else if( strcmp( file_rec->word[0], "INITIAL_VELOCITIES_BY_FILE" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &INITIAL_VELOCITIES_BY_FILE ) ;
	  sscanf( file_rec->word[2] , "%s" , velocities_file_name );
	}
	else if( strcmp( file_rec->word[0], "SAVING_DIRECTORY_NAME" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%s" , _DIRECTORY_NAME );
	}

	else if( strcmp( file_rec->word[0], "RUNTIME_CALCULATIONS" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.RUNTIME_CALCULATIONS ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.RUNTIME_CALCULATIONS ) ;
	  sscanf( file_rec->word[3] , "%d" , &SAVE.runtime_interval ) ;
	  sscanf( file_rec->word[3] , "%d" , &EQ_SAVE_OPT.runtime_interval ) ;
	  sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.runtime_interval ) ;
	}
      }
    }
    delete file_rec ;
    file_rec = new record ;
  }
  delete file_rec ;
  file_rec = NULL ;

  // initialization of some variables taken in input script, printing some global information
  if( particle_mass > 0 ) for( int i=0; i<numberOfPart; i++ ) particles[i]->mass = particle_mass ;
  if( N_patches > 0 ) for( int i=0; i<numberOfPart; i++ ) particles[i]->set_patches( N_patches ) ;
  // Initialization of verlet lists pointers in particles
  for( int i=0; i<numberOfPart; i++ ) {
    particles[i]->list_ptr = (vlist_element<particle> **)calloc( INTERACTIONS.number , sizeof(vlist_element<particle> *) );
    if( check_memalloc( particles[i]->list_ptr , "\n list pointer in particles within initialization function" ) ) exit( EXIT_FAILURE ) ;
  }
  for( int k=0 ; k<INTERACTIONS.number ; k++ )  INTERACTIONS.type[k]->link_list2particles() ;
  if( DIMENSION != 2 && DIMENSION != 3 ) {
    cout << " *** Specify the geometrical dimension of the system!" << endl ;
    exit(EXIT_FAILURE) ;
  }
  cout << "kbT = " << k_T << endl ;
  midside.position = (box_sides.position * 0.5) ; // redundant, but midside is a pain... find the way to get rid of it
}
/****************************************************************************************/

/***********************    BLOCK 2   ***************************************************/
/****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_potential_energy( void ) {
  POTENTIAL_ENERGY.val = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    INTERACTIONS.type[k]->partial_energy.val = INTERACTIONS.type[k]->compute_energy_contribution() ;
    INTERACTIONS.type[k]->partial_energy.time = global_time ;
    POTENTIAL_ENERGY.val += INTERACTIONS.type[k]->partial_energy.val ;
  }
  POTENTIAL_ENERGY.time = global_time ;
  return POTENTIAL_ENERGY.val ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_potential_energy( particle *part ) {   // MC mainly       // This requires verlet lists constructed in MC mode !
  double energy = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) energy += INTERACTIONS.type[k]->compute_energy_contribution( part ) ;
  return energy ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_potential_energy( particle *part1 , particle *part2 ) {  // MC mainly
  double energy = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) energy += INTERACTIONS.type[k]->return_pair_potential( part1 , part2 ) ;
  if( strcmp( simulation_mode , "MC" ) == 0 ) return energy ;
  else {
    cout << " *** ERROR in function energy_compute( &part1 , &part2 ) : Verlet lists are constructed in MD mode, the single particle energy cannot be computed !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::delta_energy_constrained_CoM( particle *displacement ) {
  double energy = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) energy += INTERACTIONS.type[k]->delta_energy_constrained_CoM( displacement ) ;
  return energy ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_lambda_derivative_perpart( void ) {
  double der = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ )  der += INTERACTIONS.type[k]->lambda_derivative_perpart() ;
  return der ;
}
/*****************************************************************************************/

template <typename particle>
void configuration<particle>::calculate_forces_t2() {     // MD    // calculates the vector of forces acting on particles and place it in forces_t2
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    cout << " *** ERROR : function calculate_forces_t2() cannot be used in MC mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
          // reset of forces_t2
  for( int i=0 ; i<numberOfPart ; i++ ) forces_t2[i]->clear() ;
          // computation of forces due to each type of interaction
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) INTERACTIONS.type[k]->compute_forces_contribution() ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_virial_pressure(void) {      // si può ottimizzare
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    cout << " *** ERROR : to compute the virial in Monte Carlo simulations you have to allocate the array of forces ! " << endl
	 << "             Remember that you should re-organize the information on particles, probably embedding forces in some way," << endl
	 << "             but you have also consider what this would imply for integration with velocity verlet algorithm . . ." << endl
	 << "             Maybe you have to store a variable in the functions that tells what force among two refers to the actual step . . ." << endl ;
    exit( EXIT_FAILURE ) ;
  }

  if( virial_press.time != global_time ) {
    virial_press.val = 0 ;
    /*     SUM OVER ALL THE PAIR INTERACTIONS     */
    /*     **********************************     */
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      INTERACTIONS.type[k]->partial_virial.val = INTERACTIONS.type[k]->compute_virial_contribution() ;
      INTERACTIONS.type[k]->partial_virial.time = global_time ;
      virial_press.val += INTERACTIONS.type[k]->partial_virial.val ;
    }
    virial_press.val /= ( DIMENSION * box_sides.volume() ) ;
    virial_press.time = global_time ;
  }
  return virial_press.val ;
}
/****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_kinetic_energy() {           //MD
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    cout << " *** ERROR : function compute_kinetic_energy() cannot be used in MC mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }

  if( KINETIC_ENERGY.time != global_time ) {
    KINETIC_ENERGY.val = 0 ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      KINETIC_ENERGY.val += velocities[i]->position.square_norm() ;
    }
    KINETIC_ENERGY.val /= 2.0 ;
    KINETIC_ENERGY.time = global_time ;
  }
  return KINETIC_ENERGY.val ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_packing_fraction() {
  double packfrac = 0.0 ;
  for( int i=0 ; i<numberOfPart ; i++ ) packfrac += ( particles[i]->radius * particles[i]->radius ) ;
  packfrac *= ( M_PI / box_sides.volume() ) ;
  return packfrac ;
}
/*****************************************************************************************/

template <>
double configuration<particle_3D>::compute_packing_fraction() {
  double packfrac = 0.0 ;
  for( int i=0 ; i<numberOfPart ; i++ ) packfrac += ( particles[i]->radius * particles[i]->radius * particles[i]->radius ) ;
  packfrac *= ( 4. / 3 * M_PI / box_sides.volume() ) ;
  return packfrac ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::calculate_secondary_potential(void) {                   //MC
  // double energy = 0 ;
  // double dist = 0 ;
  // struct neighbour_list *walker = NULL ;
  // particle apparent_near ;

  cout << " Rewrite this function (secondary potential calculation)" << endl ;
  exit( EXIT_FAILURE ) ;

  /* for(int i=0; i<numberOfPart-1; i++) {                   // attenzione al numberOfPart-1 */
  /*   walker = verlet_list->neighbours[i] ; */
  /*   while( walker != NULL ) {                 //for each i-particle I sum to the total the interaction energy with near j-particle if j>i (see verlet list construction) in order to avoid double counting */
  /*     if( (walker->near) > i ) {                //and if Ri-Rj<Rcutoff */
  /* 	apparent_near = right_copy( walker->near , i ) ;      // I take the coordinates of the nearest copy of the neighbour */
  /* 	if( ( dist = particles[i]->distance( apparent_near ) ) < R_cutoff ) energy += secondary_potential( particles[i] , &apparent_near , dist ) ; */
  /*     } */
  /*     walker = walker->next ; */
  /*   } */
  /* } */
  /* return energy ; */
}
/****************************************************************************************/

template <typename particle>
particle configuration<particle>::compute_mean_position(void) {
  particle mean_pos ;
  for( int i=0; i<numberOfPart; i++ ) mean_pos.position += particles[i]->position ;
  mean_pos.position /= (double)numberOfPart ;

  return mean_pos ;
}
/****************************************************************************************/

template <typename particle>
particle configuration<particle>::compute_centre_of_mass(void) {
  centre_of_mass.position *= 0.0 ;
  double mass_tot = 0 ;
  for( int i=0; i<numberOfPart; i++ ) {
    centre_of_mass.position += ( particles[i]->unwrapped( box_sides.position ) * particles[i]->mass ) ;
    mass_tot += particles[i]->mass ;
  }
  centre_of_mass.position /= mass_tot ;

  return centre_of_mass ;
}
/****************************************************************************************/

template <typename particle>
particle configuration<particle>::compute_mean_velocity( void ) {
  particle mean_vel ;
  for( int i=0 ; i<numberOfPart ; i++ ) mean_vel.position += velocities[i]->position ;
  mean_vel.position /= (double)numberOfPart ;

  return mean_vel ;
}
/****************************************************************************************/

template <typename particle>
struct posizione_3D<double> configuration<particle>::compute_CoM_angular_momentum( void ) {    // in case of identical particles of mass = 1
  particle CoM ;
  posizione_3D<double> ang_mom { 0.0 , 0.0 , 0.0 } ;
  CoM = compute_mean_position() ;
  for( int i=0 ; i<numberOfPart ; i++ ) ang_mom -= ( CoM.position - particles[i]->position ).vectorial_product( velocities[i]->position ) ;

  return ang_mom ;
}
/****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_inertia_momentum( posizione_3D<double> axis ) {    // in case of identical particles of mass = 1
  double I = 0 ;
  particle CoM ;
  CoM = compute_mean_position() ;
  axis /= axis.norm() ;
  for( int i=0; i<numberOfPart; i++ ) {
    I += axis.vectorial_product( particles[i]->position - CoM.position ).vectorial_product( particles[i]->position - CoM.position ).norm() ;
  }

  return I ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::compute_inertia_momentum( void ) {    // in case of identical particles of mass = 1
  particle CoM ;
  posizione_3D<double> r ;
  CoM = compute_mean_position() ;
  for( int i=0; i<3; i++ ) {
    for( int j=0; j<3; j++ ) {
      inertia_matrix.element[i][j] = 0 ;
    }
  }
  for( int i=0; i<numberOfPart; i++ ) {
    r = ( particles[i]->position - CoM.position ) ;
    inertia_matrix.element[0][0] += particles[i]->mass * ( r.y * r.y + r.z * r.z ) ;
    inertia_matrix.element[0][1] -= particles[i]->mass * r.x * r.y ;
    inertia_matrix.element[0][2] -= particles[i]->mass * r.x * r.z ;
    inertia_matrix.element[1][0] -= particles[i]->mass * r.y * r.x ;
    inertia_matrix.element[1][1] += particles[i]->mass * ( r.x * r.x + r.z * r.z ) ;
    inertia_matrix.element[1][2] -= particles[i]->mass * r.y * r.z ;
    inertia_matrix.element[2][0] -= particles[i]->mass * r.z * r.x ;
    inertia_matrix.element[2][1] -= particles[i]->mass * r.z * r.y ;
    inertia_matrix.element[2][2] += particles[i]->mass * ( r.y * r.y + r.x * r.x ) ;
  }
}
/****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_max_bond_distance( void ) {
  neighbour_list<particle> *walker=NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  double max_bonding_distance = 0, dist = 0 ;

  for( int l=0 ; l<CHECKS.bonds_lists_number ; l++ ) {
    for( int i=0 ; i<CHECKS.bonds_verlet_link[l]->particles_in_list ; i++ ) {                 // bonding interactions
      walker = CHECKS.bonds_verlet_link[l]->element[i].neighbours ;
      while( walker != NULL ) {
	right_copy( walker->ptr , CHECKS.bonds_verlet_link[l]->element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
	if( ( dist = particles[i]->position.distance( apparent_near->position ) ) > max_bonding_distance ) max_bonding_distance = dist ;
	walker = walker->next ;
      }
    }
  }

  delete apparent_near ;
  return max_bonding_distance ;
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::calculate_molecules_CoM( void ) {
  cout << "ERROR in calculate_molecules_gyration() : No molecules have been defined within this scope!" << endl ;
  exit( EXIT_FAILURE ) ;
}
/*****************************************************************************************/

template <typename particle>
inline double configuration<particle>::calculate_molecules_gyration( void ) {
  cout << "ERROR in calculate_molecules_gyration() : No molecules have been defined within this scope!" << endl ;
  exit( EXIT_FAILURE ) ;
  return 0 ;
}
/*****************************************************************************************/

template <typename particle>
void configuration<particle>::initializeVlist_runtime_calculations( verlet<particle> *charge_vlist ) {
               // attempt to initialize partially the runtime calc structures
  cout << " Warning : No runtime calculations defined in the base class !" << endl ;
}
/*****************************************************************************************/

template <typename particle>
inline bool configuration<particle>::check_verlet_update(void) {
  bool flag = 0 ;
  verlet<particle> *verlet_list = NULL ;
  vlist_element<particle> *vl_element = NULL ;
  particle displacement ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    verlet_list = INTERACTIONS.type[k]->verlet_list ;
    vl_element = verlet_list->element ;
            // In this version I trigger the update of all the verlet lists at the same time, if needed
    if( verlet_list->disp_on  &&  flag == 0 ) {
      double max_disp = verlet_list->delta_verlet * 0.5 ;
      for( int i=0 ; ( i<verlet_list->particles_in_list )&&( flag == 0 ) ; i++ ) {
	    // guadagno realmente cpu-time con questa cosa ?????????????????????????
            // If any particle has moved more than (rverlet-rcut)/2 I recalculate the verlet list
	displacement.position = vl_element[i].ptr->position - vl_element[i].cell_center->position ;
	if( displacement.position.norm() > max_disp ) flag = 1 ;
      }
    }
  }
  if( flag == 1 ) {
    lists_rebuilds ++ ;
    for( int i=0 ; i<numberOfPart ; i++ ) pbc( particles[i] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      verlet_list = INTERACTIONS.type[k]->verlet_list ;
      if( verlet_list->must_be_updated == 1 ) {
	verlet_list->clear() ;
	if( strcmp( simulation_ensemble , "MC_NPT" ) == 0 ) {
	  if( verlet_list->cells == NULL ) verlet_list->initialize_cellist() ;
	  else verlet_list->reinitialize_cellist() ;
	}
	verlet_list->verlet_creation_bycell( simulation_mode ) ;
      }
    }
  }
  return flag ;
}
/****************************************************************************************/

template <typename particle>
inline bool configuration<particle>::check_verlet_update( particle *part ) {
  bool flag = 0 ;
  double max_disp ;
  verlet<particle> *verlet_list = NULL ;
  for( int k=0 ; k<INTERACTIONS.number && flag == 0  ; k++ ) {
    verlet_list = INTERACTIONS.type[k]->verlet_list ;
    if( verlet_list->disp_on  &&  flag == 0 ) {
      max_disp = verlet_list->delta_verlet * 0.5 ;
      if( ( part->position - part->list_ptr[k]->cell_center->position ).norm() > max_disp ) flag = 1 ;
    }
  }

  if( flag == 1 ) {
    lists_rebuilds ++ ;
    for( int i=0 ; i<numberOfPart ; i++ ) pbc( particles[i] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      verlet_list = INTERACTIONS.type[k]->verlet_list ;
      if( verlet_list->must_be_updated == 1 ) {
	verlet_list->clear() ;
	if( strcmp( simulation_ensemble , "MC_NPT" ) == 0 ) {
	  if( verlet_list->cells == NULL ) verlet_list->initialize_cellist() ;
	  else verlet_list->reinitialize_cellist() ;
	}
	verlet_list->verlet_creation_bycell( simulation_mode ) ;
      }
    }
  }
  return flag ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::add_initialize_MSD( const char *mode , const char *steps_path ) {
  if( mode != NULL && strcmp( mode , "all" ) != 0 && strcmp( mode , "mols_cm" ) != 0 ) {
    cout << "ERROR: MSD initialization mode not recognized !!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  mean_squared_displacement<particle> **temp = NULL ;
  N_MSDs ++ ;
  temp = new mean_squared_displacement<particle> * [N_MSDs] ;
  for( int i=0 ; i<N_MSDs-1 ; i++ ) temp[i] = MSD[i] ;
  if( MSD != NULL ) delete [] MSD ;
  MSD = temp ;
  MSD[N_MSDs-1] = new mean_squared_displacement<particle> ;
  strcpy( MSD[N_MSDs-1]->_steps_path , steps_path ) ;

  if( mode == NULL || strcmp( mode , "all" ) == 0 ) {
    strcpy( MSD[N_MSDs-1]->_msd_filename , "msd_particles.dat" ) ;
    MSD[N_MSDs-1]->Nparts = numberOfPart ;
    MSD[N_MSDs-1]->parts = new particle * [ MSD[N_MSDs-1]->Nparts ] ;
    for( int i=0 ; i<MSD[N_MSDs-1]->Nparts ; i++ ) MSD[N_MSDs-1]->parts[i] = particles[i] ;

  } else if( strcmp( mode , "mols_cm" ) == 0 ) {
    strcpy( MSD[N_MSDs-1]->_msd_filename , "msd_molecules.dat" ) ;
    MSD[N_MSDs-1]->Nparts = numberOfMol ;
    MSD[N_MSDs-1]->parts = new particle * [ MSD[N_MSDs-1]->Nparts ] ;
    for( int i=0 ; i<MSD[N_MSDs-1]->Nparts ; i++ ) MSD[N_MSDs-1]->parts[i] = &( molecules[i]->cm ) ;
  }
};
/****************************************************************************************/
