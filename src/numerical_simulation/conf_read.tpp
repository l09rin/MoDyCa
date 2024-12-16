/* method to load configurations from file */

template <typename particle>
particle read_lammps_box( ifstream &data_file , int box_lines ) { cout << "*** Unknown space dimension" << endl ; };

template <typename particle>
void read_lammps_atoms( ifstream &data_file , particle box_sides , int numberOfPart , int particles_length , particle **particles , particle ***velocities , record *line ) { cout << "*** Unknown space dimension" << endl ; };

/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::read_configuration( ifstream &data_file , char *format ) {
  numberOfPart = 0 ;

  if ( strcmp( format , "xyz" ) == 0 ) read_xyz( data_file ) ;

  else if ( strcmp( format , "lmp" ) == 0 ) read_lammps_dump( data_file ) ;

  else if ( strcmp( format , "nico_bonds" ) == 0 ||
	    strcmp( format , "nico_bonds_pbc" ) == 0 ||
	    strcmp( format , "nico_rings" ) == 0 ||
	    strcmp( format , "nico_rings_pbc" ) == 0 ) read_nico_init( data_file , format ) ;

  else if ( strcmp( format , "sph" ) == 0 ||
	    strcmp( format , "flf" ) == 0 ) read_sph( data_file ) ;

  else if ( strcmp( format , "ptc" ) == 0 ||
	    strcmp( format , "patch" ) == 0 ) read_ptc( data_file ) ;

  else {
    cout << "*** ERROR: unrecognized configuration file format!" << endl ;
    exit( EXIT_FAILURE ) ;
  }

  for ( int pidx=0 ; pidx<numberOfPart ; pidx++ ) pbc( particles[pidx] ) ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::read_xyz( ifstream &data_file ) {
  int N = 0 ;
  record *line = NULL ;

  // loading of the configuration parameters
  line = new record ;
  line->getrecord( data_file ) ;
  line->split() ;
  int flagg = 0 ;             // This cicle throw away all blank or comment lines ahead each configuration
  while( flagg == 0 && !(data_file.eof()) ) {
    if( line->line == NULL ) flagg = 1 ;
    else if( line->words_number == 0 ) {
      line->getrecord( data_file ) ;
      line->split() ;
    } else if( strpbrk( line->word[0], "#" ) != NULL ) {
      if( strcmp( line->word[1] , "N" ) == 0 )  sscanf( line->word[2] , "%d" , &N ) ;
      else if( strcmp( line->word[1] , "time" ) == 0 || strcmp( line->word[1] , "timestep" ) == 0 )  sscanf( line->word[2] , "%lf" , &global_time ) ;
      else if( strcmp( line->word[1] , "box") == 0 ) {
	box_sides.read_by_string( line->word + 2 ) ;
	midside.position = (box_sides.position * 0.5) ;
      }
      line->getrecord( data_file ) ;
      line->split() ;
    } else flagg = 1 ;
  }
  numberOfPart = 0 ;
  // loading of all particles in the configuration
  while( line->words_number != 0 && !(data_file.eof()) ) {
    numberOfPart ++ ;
    if( particles_length < numberOfPart ) {
      particles_length += 100 ;
      particles = (particle **)realloc( particles , particles_length*sizeof(particle *) ) ;
      check_memalloc( particles , "Allocation error 1 in read_xyz" ) ;
    }
    particles[numberOfPart-1] = new particle ;
    particles[numberOfPart-1]->id = numberOfPart-1 ;
    if( line->words_number >= 2*DIMENSION ) particles[numberOfPart-1]->read_by_string_pbc( line->line ) ;
    else particles[numberOfPart-1]->read_by_string( line->line ) ;
    particles[numberOfPart-1]->type = 1 ;

    line->getrecord( data_file ) ;
    line->split() ;
  }
  particles_length = numberOfPart ;
  particles = (particle **)realloc( particles , particles_length*sizeof(particle *) ) ;
  if ( N != 0 && N != numberOfPart) printf( "+++ Found a mismatch between the number of read particles and the one indicated in the file! [ %d , %d ]" , numberOfPart , N ) ;
  delete line ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::read_lammps_dump( ifstream &data_file ) {
  record *line = new record ;

  for ( int kk=0 ; kk<4 ; kk++ ) {
    line->getrecord( data_file ) ;
    line->split() ;
    if ( !(data_file.eof()) ) {
      if ( line->words_number > 1 ) {

	if ( strcmp( line->word[1] , "NUMBER" ) == 0 ) {
	  line->getrecord( data_file ) ;
	  line->split() ;
	  sscanf( line->word[0] , "%d" , &numberOfPart ) ;
	  if ( numberOfPart > particles_length ) {
	    particles = (particle **)realloc( particles , numberOfPart*sizeof(particle *) ) ;
	    check_memalloc( particles , "Allocation error 1 in read_lammps_dump" ) ;
	    for( int p=particles_length ; p<numberOfPart ; p++ )  particles[p] = new particle ;
	    if (velocities != NULL) {
	      velocities = (particle **)realloc( velocities , numberOfPart*sizeof(particle *) ) ;
	      check_memalloc( velocities , "Allocation error 2 in read_lammps_dump" ) ;
	      for( int p=particles_length ; p<numberOfPart ; p++ ) velocities[p] = new particle ;
	    }
	    particles_length = numberOfPart ;
	  }

	}else if ( strcmp( line->word[1] , "BOX" ) == 0 ) {
	  int box_lines = 2 ;
	  if( line->words_number > 5 ) box_lines = 3 ;
	  box_sides = read_lammps_box<particle>( data_file , box_lines ) ;
	  midside.position = (box_sides.position * 0.5) ;

	}else if ( strcmp( line->word[1] , "TIMESTEP" ) == 0 ) {
	  line->getrecord( data_file ) ;
	  line->split() ;
	  sscanf( line->word[0] , "%lf" , &global_time ) ;

	}else if ( strcmp( line->word[1] , "ATOMS" ) == 0 ) {
	  read_lammps_atoms<particle>( data_file , box_sides , numberOfPart , particles_length , particles , &velocities , line ) ;
	}
      }
    }
  }

  delete line ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::read_nico_init( ifstream &data_file , const char *format ) {
  record *file_rec = new record ;

  file_rec->getrecord( data_file ) ;     // Ignoring first record

  file_rec->getrecord( data_file ) ;     // Reading number of monomers and constructing configuration's array
  file_rec->split() ;
  sscanf( file_rec->word[0] , "%d" , &monomers_number ) ;
  numberOfPart = monomers_number ;
  particles_length = numberOfPart ;
  particles = new particle* [particles_length] ;

  file_rec->getrecord( data_file ) ;     // Reading box's properties
  file_rec->split() ;
  box_sides.read_by_string( file_rec->line ) ;
  if( box_sides.min_comp() == box_sides.max_comp() ) {
    box_side = box_sides.min_comp() ;
    midside.position = (box_sides.position * 0.5) ;
  }else{
    box_side = nanl("") ;
  }

  for( int i=0; i<monomers_number; i++ ) {         // Reading particles' positions
    file_rec->getrecord( data_file ) ;
    particles[i] = new particle ;
    particles[i]->id = i ;
    if( strcmp( format , "nico_bonds_pbc" ) == 0 || strcmp( format , "nico_rings_pbc" ) == 0 ) particles[i]->read_by_string_pbc( file_rec->line ) ;
    else particles[i]->read_by_string( file_rec->line ) ;
    particles[i]->type = 1 ;
  }

  int monomer_index = 0, bonds_number = 0 ;
  file_rec->getrecord( data_file ) ;
  while( file_rec->char_number != 0 ) {                       // Reading monomers' bonding properties
    file_rec->split() ;
    sscanf( file_rec->word[0] , "%d" , &monomer_index ) ;
    sscanf( file_rec->word[1] , "%d" , &bonds_number ) ;
    particles[monomer_index-1]->valence = bonds_number ;
    if( bonds_number == 0 ) {
      particles[monomer_index-1]->bonded_monomer = NULL ;
    } else {
      particles[monomer_index-1]->bonded_monomer = new int [ bonds_number ] ;
      for( int i=0 ; i<bonds_number ; i++ ) particles[monomer_index-1]->bonded_monomer[i] = -1 ;
      file_rec->getrecord( data_file ) ;
      file_rec->split() ;
      for( int i=0; i<bonds_number; i++ ) {        // For each monomer the array of bonded neighbours is filled
	sscanf( file_rec->word[i] , "%d", particles[monomer_index-1]->bonded_monomer + i ) ;
	particles[monomer_index-1]->bonded_monomer[i] -- ;                  // In configuration file monomers' index starts from 1 instead of 0
      }
    }
    file_rec->getrecord( data_file ) ;
  }

  delete file_rec ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::read_sph( ifstream &data_file ) {
  record *line = new record ;
  char tmp ;

  // loading of the configuration
  global_time = -1 ; // In sph format usually time is not indicated
  line->getrecord( data_file ) ;
  line->split() ;
  int flagg = 0 ;             // This cicle throw away all blank or comment lines ahead each configuration
  while( flagg == 0 && !(data_file.eof()) ) {
    if( line->line == NULL ) flagg = 1 ;
    else if( strpbrk( line->line, "#" ) != NULL || line->words_number == 0 ) {
      line->getrecord( data_file ) ;
      line->split() ;
    } else {
      flagg = 1 ;
    }
  }
  numberOfPart = 0 ;
  // loading of all particles in the configuration
  if( !data_file.eof() ) {
    sscanf( line->word[0] , "%d" , &numberOfPart ) ;
    if ( numberOfPart < 1 ) printf( "*** ERROR: I was not able to read the configuration\n" ) ;
    if ( line->words_number > 1 ) sscanf( line->word[1] , "%lf" , &global_time ) ;
    if( global_time == -1 ) global_time = 0 ;
    line->getrecord( data_file ) ;
    line->split() ;
    box_sides.read_by_string( line->line ) ;
    midside.position = (box_sides.position * 0.5) ;
    // sscanf( line->line , "%lf %lf %lf" , &(box_sides.position.x) , &(box_sides.position.y) , &(box_sides.position.z) ) ;
    if( particles_length < numberOfPart ) {
      int new_length = numberOfPart ;
      particles = (particle **)realloc( particles , new_length*sizeof(particle *) ) ;
      check_memalloc( particles , "Allocation error 1 in read_sph" ) ;
      for( int p=particles_length ; p<new_length ; p++ )  particles[p] = new particle ;
      particles_length = new_length ;
    }
    for ( int pidx=0 ; pidx<numberOfPart ; pidx++ ) {
      line->getrecord( data_file ) ;
      line->split() ;
      if ( data_file.eof() ) {
	cout << "Error in reading the configuration" << endl ;
	fflush(0) ;
	break ;
      }
      sscanf( line->word[0] , "%c" , &tmp ) ;
      particles[pidx]->type = tmp - 'a' + 1 ;
      particles[pidx]->read_by_string( line->word + 1 ) ;
      sscanf( line->word[4] , "%lf" , &(particles[pidx]->radius) ) ;
      particles[pidx]->id = pidx ;
      if ( line->words_number != 5 ) {
	printf( "*** ERROR: Read error (particle) %d \n String: %s\n" , line->words_number , line->line ) ;
	exit( EXIT_FAILURE ) ;
      }
    }
  }

  delete line ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::read_ptc( ifstream &data_file ) {
  cout << "*** ERROR: attempted to use read_ptc method with no patchy particles!" << endl ;
  exit( EXIT_FAILURE ) ;
}

template <>
void configuration<patchy_2D>::read_ptc( ifstream &data_file ) {
  record *line = new record ;
  char tmp ;
  double costheta , sintheta ;

  // loading of the configuration
  line->getrecord( data_file ) ;
  line->split() ;
  int flagg = 0 ;             // This cicle throw away all blank or comment lines ahead each configuration
  while( flagg == 0 && !(data_file.eof()) ) {
    if( line->words_number > 0 ) {
      if( line->word[0][0] == '&' ) flagg = 1 ;
      else {
	line->getrecord( data_file ) ;
	line->split() ;
      }
    }
  }
  numberOfPart = 0 ;
  // loading of all particles in the configuration
  if( !data_file.eof() ) {
    sscanf( line->word[0] , "&%d" , &numberOfPart ) ;
    if ( numberOfPart < 1 ) printf( "*** ERROR: I was not able to read the configuration\n" ) ;
    line->getrecord( data_file ) ;
    line->split() ;
    box_sides.read_by_string( line->line ) ;
    midside.position = (box_sides.position * 0.5) ;
    if( particles_length < numberOfPart ) {
      int new_length = numberOfPart ;
      particles = (patchy_2D **)realloc( particles , new_length*sizeof(patchy_2D *) ) ;
      check_memalloc( particles , "Allocation error 1 in read_sph" ) ;
      for( int p=particles_length ; p<new_length ; p++ )  particles[p] = NULL ;
      particles_length = new_length ;
    }
    for ( int pidx=0 ; pidx<numberOfPart ; pidx++ ) {
      line->getrecord( data_file ) ;
      line->split() ;
      if ( data_file.eof() ) {
	cout << "Error in reading the configuration" << endl ;
	fflush(0) ;
	break ;
      }
      if( particles[pidx] == NULL ) particles[pidx] = new patchy_2D( line->words_number-16 , line->words_number-16 ) ;
      sscanf( line->word[0] , "%c" , &tmp ) ;
      particles[pidx]->type = tmp - 'a' + 1 ;
      particles[pidx]->read_by_string( line->word + 1 ) ;
      sscanf( line->word[4] , "%lf" , &(particles[pidx]->radius) ) ;
      sscanf( line->word[7] , "%lf" , &costheta ) ;
      sscanf( line->word[8] , "%lf" , &sintheta ) ;
      if( sintheta > 0 ) particles[pidx]->rotation = acos( costheta ) ;
      else if( sintheta < 0 ) particles[pidx]->rotation = -acos( costheta ) ;
      else {
	if( costheta > 0 ) particles[pidx]->rotation = 0.0 ;
	else particles[pidx]->rotation = M_PI ;
      }
      particles[pidx]->id = pidx ;
      for( int j=16; j<line->words_number; j++ ) sscanf( line->word[j] , "%d" , &(particles[pidx]->bonded_monomer[j-16]) ) ;
    }
  }

  delete line ;
}
/********************************************************************************************************************/

template <>
particle_3D read_lammps_box( ifstream &data_file , int box_lines ) {
  record *line = new record ;
  particle_3D box_sides ;

  line->getrecord( data_file ) ;
  line->split() ;
  double inf = 0 , sup = 0 ;
  sscanf( line->word[0] , "%lf" , &inf ) ;
  sscanf( line->word[1] , "%lf" , &sup ) ;
  box_sides.position.x = sup - inf ;
  line->getrecord( data_file ) ;
  line->split() ;
  sscanf( line->word[0] , "%lf" , &inf ) ;
  sscanf( line->word[1] , "%lf" , &sup ) ;
  box_sides.position.y = sup - inf ;
  line->getrecord( data_file ) ;
  line->split() ;
  sscanf( line->word[0] , "%lf" , &inf ) ;
  sscanf( line->word[1] , "%lf" , &sup ) ;
  box_sides.position.z = sup - inf ;
  return box_sides ;
}

template <>
particle_2D read_lammps_box( ifstream &data_file , int box_lines ) {
  record *line = new record ;
  particle_2D box_sides ;

  line->getrecord( data_file ) ;
  line->split() ;
  double inf = 0 , sup = 0 ;
  sscanf( line->word[0] , "%lf" , &inf ) ;
  sscanf( line->word[1] , "%lf" , &sup ) ;
  box_sides.position.x = sup - inf ;
  line->getrecord( data_file ) ;
  line->split() ;
  sscanf( line->word[0] , "%lf" , &inf ) ;
  sscanf( line->word[1] , "%lf" , &sup ) ;
  box_sides.position.y = sup - inf ;
  if( box_lines == 3 ) {
    line->getrecord( data_file ) ;
    line->split() ;
  }
  return box_sides ;
}

template <>
patchy_2D read_lammps_box( ifstream &data_file , int box_lines ) {
  record *line = new record ;
  patchy_2D box_sides ;

  line->getrecord( data_file ) ;
  line->split() ;
  double inf = 0 , sup = 0 ;
  sscanf( line->word[0] , "%lf" , &inf ) ;
  sscanf( line->word[1] , "%lf" , &sup ) ;
  box_sides.position.x = sup - inf ;
  line->getrecord( data_file ) ;
  line->split() ;
  sscanf( line->word[0] , "%lf" , &inf ) ;
  sscanf( line->word[1] , "%lf" , &sup ) ;
  box_sides.position.y = sup - inf ;
  if( box_lines == 3 ) {
    line->getrecord( data_file ) ;
    line->split() ;
  }
  return box_sides ;
}
/********************************************************************************************************************/

template <>
void read_lammps_atoms( ifstream &data_file , particle_3D box_sides , int numberOfPart , int particles_length , particle_3D **particles , particle_3D ***velocities , record *line ) {
  int type_col = -1 , mol_col = -1 , x_col = -1 , y_col = -1 , z_col = -1 , vx_col = -1 , vy_col = -1 , vz_col = -1 , ix_col = -1 , iy_col = -1 , iz_col = -1 , id_col = -1 , q_col = -1 ;
  for ( int i=2 ; i<line->words_number ; i++ ) {
    if ( strcmp( line->word[i] , "type" ) == 0 ) type_col = i-2 ;
    else if ( strcmp( line->word[i] , "mol" ) == 0 ) mol_col = i-2 ;
    else if ( strcmp( line->word[i] , "x" ) == 0 || strcmp( line->word[i] , "xu" ) == 0 ) x_col = i-2 ;
    else if ( strcmp( line->word[i] , "y" ) == 0 || strcmp( line->word[i] , "yu" ) == 0 ) y_col = i-2 ;
    else if ( strcmp( line->word[i] , "z" ) == 0 || strcmp( line->word[i] , "zu" ) == 0 ) z_col = i-2 ;
    else if ( strcmp( line->word[i] , "vx" ) == 0 ) vx_col = i-2 ;
    else if ( strcmp( line->word[i] , "vy" ) == 0 ) vy_col = i-2 ;
    else if ( strcmp( line->word[i] , "vz" ) == 0 ) vz_col = i-2 ;
    else if ( strcmp( line->word[i] , "ix" ) == 0 ) ix_col = i-2 ;
    else if ( strcmp( line->word[i] , "iy" ) == 0 ) iy_col = i-2 ;
    else if ( strcmp( line->word[i] , "iz" ) == 0 ) iz_col = i-2 ;
    else if ( strcmp( line->word[i] , "id" ) == 0 ) id_col = i-2 ;
    else if ( strcmp( line->word[i] , "q" ) == 0 ) q_col = i-2 ;
  }
  if ( *velocities == NULL && (vx_col != -1 || vy_col != -1 || vz_col != -1) ) {
    *velocities = (particle_3D **)calloc(particles_length, sizeof(particle_3D *));
    for(int i=0; i<particles_length; i++)  (*velocities)[i] = new particle_3D ;
  }
  for ( int i=0 ; i<numberOfPart ; i++ ) {
    line->getrecord( data_file ) ;
    line->split() ;
    if ( id_col != -1 ) sscanf( line->word[id_col] , "%d" , &(particles[i]->id) ) ;
    if ( type_col != -1 ) sscanf( line->word[type_col] , "%d" , &(particles[i]->type) ) ;
    if ( mol_col != -1 ) sscanf( line->word[mol_col] , "%d" , &(particles[i]->mol) ) ;
    if ( x_col != -1 ) sscanf( line->word[x_col] , "%lf" , &(particles[i]->position.x) ) ;
    if ( y_col != -1 ) sscanf( line->word[y_col] , "%lf" , &(particles[i]->position.y) ) ;
    if ( z_col != -1 ) sscanf( line->word[z_col] , "%lf" , &(particles[i]->position.z) ) ;
    if ( q_col != -1 ) sscanf( line->word[q_col] , "%lf" , &(particles[i]->charge) ) ;
    if ( ix_col != -1 ) sscanf( line->word[ix_col] , "%d" , &(particles[i]->periodic_box.x) ) ;
    if ( iy_col != -1 ) sscanf( line->word[iy_col] , "%d" , &(particles[i]->periodic_box.y) ) ;
    if ( iz_col != -1 ) sscanf( line->word[iz_col] , "%d" , &(particles[i]->periodic_box.z) ) ;
    particles[i]->position.x = particles[i]->position.x + ( box_sides.position.x * particles[i]->periodic_box.x ) ;
    particles[i]->position.y = particles[i]->position.y + ( box_sides.position.y * particles[i]->periodic_box.y ) ;
    particles[i]->position.z = particles[i]->position.z + ( box_sides.position.z * particles[i]->periodic_box.z ) ;
    particles[i]->periodic_box.x = 0 ;
    particles[i]->periodic_box.y = 0 ;
    particles[i]->periodic_box.z = 0 ;
    if ( vx_col != -1 ) sscanf( line->word[vx_col] , "%lf" , &((*velocities)[i]->position.x) ) ;
    if ( vy_col != -1 ) sscanf( line->word[vy_col] , "%lf" , &((*velocities)[i]->position.y) ) ;
    if ( vz_col != -1 ) sscanf( line->word[vz_col] , "%lf" , &((*velocities)[i]->position.z) ) ;
  }
}

template <>
void read_lammps_atoms( ifstream &data_file , particle_2D box_sides , int numberOfPart , int particles_length , particle_2D **particles , particle_2D ***velocities , record *line ) {
  int type_col = -1 , mol_col = -1 , x_col = -1 , y_col = -1 , vx_col = -1 , vy_col = -1 , ix_col = -1 , iy_col = -1 , id_col = -1 , q_col = -1 ;
  for ( int i=2 ; i<line->words_number ; i++ ) {
    if ( strcmp( line->word[i] , "type" ) == 0 ) type_col = i-2 ;
    else if ( strcmp( line->word[i] , "mol" ) == 0 ) mol_col = i-2 ;
    else if ( strcmp( line->word[i] , "x" ) == 0 || strcmp( line->word[i] , "xu" ) == 0 ) x_col = i-2 ;
    else if ( strcmp( line->word[i] , "y" ) == 0 || strcmp( line->word[i] , "yu" ) == 0 ) y_col = i-2 ;
    else if ( strcmp( line->word[i] , "vx" ) == 0 ) vx_col = i-2 ;
    else if ( strcmp( line->word[i] , "vy" ) == 0 ) vy_col = i-2 ;
    else if ( strcmp( line->word[i] , "ix" ) == 0 ) ix_col = i-2 ;
    else if ( strcmp( line->word[i] , "iy" ) == 0 ) iy_col = i-2 ;
    else if ( strcmp( line->word[i] , "id" ) == 0 ) id_col = i-2 ;
    else if ( strcmp( line->word[i] , "q" ) == 0 ) q_col = i-2 ;
  }
  if ( *velocities == NULL && (vx_col != -1 || vy_col != -1) ) {
    *velocities = (particle_2D **)calloc(particles_length, sizeof(particle_2D *));
    for(int i=0; i<particles_length; i++)  (*velocities)[i] = new particle_2D ;
  }
  for ( int i=0 ; i<numberOfPart ; i++ ) {
    line->getrecord( data_file ) ;
    line->split() ;
    if ( id_col != -1 ) sscanf( line->word[id_col] , "%d" , &(particles[i]->id) ) ;
    if ( type_col != -1 ) sscanf( line->word[type_col] , "%d" , &(particles[i]->type) ) ;
    if ( mol_col != -1 ) sscanf( line->word[mol_col] , "%d" , &(particles[i]->mol) ) ;
    if ( x_col != -1 ) sscanf( line->word[x_col] , "%lf" , &(particles[i]->position.x) ) ;
    if ( y_col != -1 ) sscanf( line->word[y_col] , "%lf" , &(particles[i]->position.y) ) ;
    if ( q_col != -1 ) sscanf( line->word[q_col] , "%lf" , &(particles[i]->charge) ) ;
    if ( ix_col != -1 ) sscanf( line->word[ix_col] , "%d" , &(particles[i]->periodic_box.x) ) ;
    if ( iy_col != -1 ) sscanf( line->word[iy_col] , "%d" , &(particles[i]->periodic_box.y) ) ;
    particles[i]->position.x = particles[i]->position.x + ( box_sides.position.x * particles[i]->periodic_box.x ) ;
    particles[i]->position.y = particles[i]->position.y + ( box_sides.position.y * particles[i]->periodic_box.y ) ;
    particles[i]->periodic_box.x = 0 ;
    particles[i]->periodic_box.y = 0 ;
    if ( vx_col != -1 ) sscanf( line->word[vx_col] , "%lf" , &((*velocities)[i]->position.x) ) ;
    if ( vy_col != -1 ) sscanf( line->word[vy_col] , "%lf" , &((*velocities)[i]->position.y) ) ;
  }
}
