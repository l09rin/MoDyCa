/* methods to generate initial velocities (for MD) */

/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_MBvelocities( double mean_square_velocity ) {   //MD  // mean_square_velocity is intended as single (x) component mean square
  if( strcmp( simulation_mode , "MC" ) == 0 ) cout << " *** WARNING : you are generating a set of velocities in Monte Carlo mode !" << endl ;
  double dev_stand = sqrt(mean_square_velocity);
  if( velocities == NULL ) {
    velocities = (particle **)calloc(numberOfPart, sizeof(particle *));
  }
  particle total_velocity ;
  if( check_memalloc( velocities , "Allocation error 1 in generate_MBvelocities()" ) ) exit( EXIT_FAILURE ) ;

  // generation of initial velocities
  for(int i=0; i<numberOfPart; i++) {
    velocities[i] = new particle;
    velocities[i]->Gaussian_particle();
    velocities[i]->position *= dev_stand;
  }

  particle CoM_velocity ;
  posizione_3D<double> total_angMomentum ;
  CoM_velocity = compute_mean_velocity() ;
  total_angMomentum = compute_CoM_angular_momentum() ;
  cout << "\n    Velocities have been generated according to a Maxwell-Boltzmann distribution" ;
  cout << "\n  Velocity of center of mass (System units): " << CoM_velocity.dump_position_variables() ;
  cout << "\n  Total angular momentum (System units): " << total_angMomentum << endl ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_velocities_by_file( void ) {   //MD
  if( strcmp( simulation_mode , "MC" ) == 0 ) cout << " *** WARNING : you are reading a set of velocities in Monte Carlo mode !" << endl ;
  //  STRUCTURE OF THE FILE :
  //   first lines : commented (#) or blank lines
  //          Then : # N <number_of_particles>
  // . . . . N velocities . . . .
  if( velocities == NULL ) {
    velocities = (particle **)calloc(numberOfPart, sizeof(particle *));
  } else {
    cout << "\n generate_velocities_by_file :  velocities vector has been already allocated!\n" ;
    exit( EXIT_FAILURE ) ;
  }
  if( check_memalloc( velocities , "Allocation error 1 in generate_velocities_by_file()" ) ) exit( EXIT_FAILURE ) ;

  // generation of initial velocities
  ifstream data_file ;

  data_file.open( velocities_file_name ) ;
  if( data_file.fail() ) {
    cout << "\n  Failure in opening data_file in function generate_velocities_by_file" << endl ;
    exit( EXIT_FAILURE ) ;
  }else{
    record *line = NULL ;
    line = new record ;
    line->getrecord( data_file ) ;
    line->split() ;
    int flag = 0 ;
    while( flag == 0 && !data_file.eof() ) {
      if( line->words_number == 0 ) {
	delete line ;
	line = new record ;
	line->getrecord( data_file ) ;     // Ignoring the record
	line->split() ;
      } else if( line->word[0][0] == '#' && line->words_number < 3 ) {
	delete line ;
	line = new record ;
	line->getrecord( data_file ) ;     // Ignoring the record
	line->split() ;
      } else if( line->word[0][0] == '#' && line->word[1][0] != 'N' ) {
	delete line ;
	line = new record ;
	line->getrecord( data_file ) ;     // Ignoring the record
	line->split() ;
      } else {
	flag = 1 ;
      }
    }
    int N_parts = 0 ;     // Reading number of particles
    if( line->words_number < 3 ) {
      cout << " Wrong file format for inputing the velocities !!" << endl << "It must contain a line containing the number of particles formatted as:" << endl ;
      cout << "  \"# N <number_of_particles>\"" << endl ;
      exit( EXIT_FAILURE ) ;
    } else if( line->word[1][0] != 'N' ) {
      cout << " Wrong file format for inputing the velocities !!" << endl << "It must contain a line containing the number of particles formatted as:" << endl ;
      cout << "  \"# N <number_of_particles>\"" << endl ;
      exit( EXIT_FAILURE ) ;
    } else {
      sscanf( line->word[2] , "%d" , &N_parts ) ;
      if( N_parts != numberOfPart ) {
	cout << "\n ERROR: The number of particles in the velocities' file does not match!" << endl ;
	exit( EXIT_FAILURE ) ;
      }
    }

    delete line ;
    line = new record ;
    for( int i=0; i<numberOfPart; i++ ) {         // Reading particles' velocity
      line->getrecord( data_file ) ;
      velocities[i] = new particle ;
      velocities[i]->read_by_string( line->line ) ;
      delete line ;
      line = new record ;
    }
    delete line ;
    line = NULL ;
  }

  compute_kinetic_energy() ;
  particle CoM_velocity ;
  posizione_3D<double> total_angMomentum ;
  CoM_velocity = compute_mean_velocity() ;
  total_angMomentum = compute_CoM_angular_momentum() ;
  cout << "\n    Velocities have been read by file : \"" << velocities_file_name << "\"" ;
  cout << "\n  Velocity of center of mass (System units): " << CoM_velocity.position ;
  cout << "\n  Total angular momentum (System units): " << total_angMomentum << endl ;
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::CoM_velocity_reset( void ) {
  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    particle total_velocity ;

    for( int i=0 ; i<numberOfPart ; i++ ) total_velocity.position += velocities[i]->position ;
    total_velocity.position /= numberOfPart ;
    for( int i=0 ; i<numberOfPart ; i++ ) velocities[i]->position -= total_velocity.position ;
  } else {
    cout << " *** I cannot reset the velocity of center of mass in Monte Carlo mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::total_angularMomentum_reset( const char *protocol ) {
  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    particle mean_pos ;
    posizione_3D<double> total_ang_mom ;
    total_ang_mom = compute_CoM_angular_momentum() ;
    mean_pos = compute_mean_position() ;

    if( strcmp(protocol, "VTINVERSION") == 0 ) {

      particle radial_velocity , transverse_velocity , new_transverse_velocity , CMposition_versor ;
      double avg_boxL = 0 , avg_modulus_vel = 0 ;
      int iterations = 0 ;
      avg_boxL = box_sides.avg_comp() / 2.0 ;
      avg_modulus_vel = compute_kinetic_energy() * 2.0 / (double)numberOfPart ;

      while( total_ang_mom.norm() > avg_boxL*avg_modulus_vel/numberOfPart && iterations < 100 ) {
	iterations ++ ;
	for( int i=0 ; i<numberOfPart ; i++ ) {
	  if( total_ang_mom.norm() > avg_boxL*avg_modulus_vel/numberOfPart ) {
	    CMposition_versor.position = particles[i]->position - mean_pos.position ;
	    // I subtract the angular momentum of the particle i from the total
	    total_ang_mom -= ( CMposition_versor.position.vectorial_product( velocities[i]->position ) ) ;
	    CMposition_versor.position /= CMposition_versor.position.norm() ;
	    // I compute the transverse velocity
	    radial_velocity.position = ( CMposition_versor.position * ( velocities[i]->position * CMposition_versor.position ) ) ;
	    transverse_velocity.position = ( velocities[i]->position - radial_velocity.position ) ;
	    // I change the transverse velocity in order to reduce the angular momentum
	    new_transverse_velocity.position = CMposition_versor.position.vectorial_product( total_ang_mom ) ;
	    if( new_transverse_velocity.position.norm() > transverse_velocity.position.norm() ) {
	      new_transverse_velocity.position /= ( new_transverse_velocity.position.norm() ) ;
	      new_transverse_velocity.position *= ( transverse_velocity.position.norm() ) ;
	    }
	    velocities[i]->position = ( radial_velocity.position + new_transverse_velocity.position ) ;
	    total_ang_mom += ( ( particles[i]->position - mean_pos.position ).vectorial_product( velocities[i]->position ) ) ;
	  }
	}
      }

    } else if( strcmp(protocol, "RIGID") == 0 || strcmp(protocol, "") == 0 ) {

      compute_inertia_momentum() ;
      matrix inertia_inverse = inertia_matrix.compute_inverse() ;
      posizione_3D<double> angular_velocity { 0.0 , 0.0 , 0.0 } ;
      angular_velocity.x = inertia_inverse.element[0][0] * total_ang_mom.x + inertia_inverse.element[0][1] * total_ang_mom.y + inertia_inverse.element[0][2] * total_ang_mom.z ;
      angular_velocity.y = inertia_inverse.element[1][0] * total_ang_mom.x + inertia_inverse.element[1][1] * total_ang_mom.y + inertia_inverse.element[1][2] * total_ang_mom.z ;
      angular_velocity.z = inertia_inverse.element[2][0] * total_ang_mom.x + inertia_inverse.element[2][1] * total_ang_mom.y + inertia_inverse.element[2][2] * total_ang_mom.z ;
      angular_velocity *= (-1.0) ;
      particle velocity_correction ;
      for( int i=0 ; i<numberOfPart ; i++ ) {
	velocity_correction.position = angular_velocity.vectorial_product( particles[i]->position - mean_pos.position ) ;
	velocities[i]->position += velocity_correction.position ;
      }

    } else {
      cout << "ERROR in function total_angularMomentum_reset(): unknown reset method!" << endl ;
      exit( EXIT_FAILURE ) ;
    }

    total_ang_mom = compute_CoM_angular_momentum() ;
    cout << "  Total angular momentum after reset (System units): " << total_ang_mom << endl ;

  } else {
    cout << " *** I cannot reset the total angular momentum in Monte Carlo mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::set_kinetic_energy( double kin_en ) {
  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    double actual_kin_en_per_part = 0 ;
    for(int i=0; i<numberOfPart; i++) actual_kin_en_per_part += velocities[i]->position.square_norm() / (double)numberOfPart ;
    double actual_dev_std = sqrt( actual_kin_en_per_part/3.0 ), dev_std = sqrt( 2.0*kin_en/3.0/(double)numberOfPart ) ;
    for(int i=0; i<numberOfPart; i++) {
      velocities[i]->position *= ( dev_std/actual_dev_std ) ;
    }
  } else {
    cout << " *** I cannot reset the kinetic energy in Monte Carlo mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/*****************************************************************************************/
