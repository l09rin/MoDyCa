/* Methods to generate the initial configuration */

/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::generate_ControlledEnergy_randomConfig( void ) {
  cout << "   The function configuration<particle>::generate_ControlledEnergy_randomConfig() generates a random" << endl ;
  cout << "     distribution of particles as if they interacted with a WCA repulsive potential." << endl ;

  double distance = 0, particle_potential_energy = 0, max_potential_energy = energyOfSystem/numberOfPart;   // potential energy must be less than total energy in order to obtain positive kinetic energy
  int too_high_energy = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  lennard_jones<particle> lj_interaction ;
  // generation of initial positions
  cout << endl << "  Generation of random positions via controlled energy insertion" << endl;
  for( int i=0 ; i<numberOfPart ; i++ ) {
    particles[i] = new particle ;
    particles[i]->id = i ;
    do {
      particles[i]->rndPos_generate( box_sides.position ) ;
      particles[i]->type = 1 ;
      too_high_energy = 0 ;
      for( int j=0; j<i; j++) {
	right_copy( particles[j] , particles[i] , apparent_near ) ;                     // PBC
	distance = particles[i]->position.distance( apparent_near->position ) ;
	particle_potential_energy = ( distance > 1.122462 ) ? 0 : ( lj_interaction.return_pair_potential( distance )+1.0 ) ;
	if( particle_potential_energy > max_potential_energy ) too_high_energy = 1 ;
      }
    }while( too_high_energy == 1 ) ;
  }
  delete apparent_near ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_HS_randomConfig( double radius ) {
  cout << "   The function configuration<particle>::generate_HS_randomConfig() generates a random" << endl ;
  cout << "     distribution of particles as if they interacted with a HS infinitely repulsive potential." << endl ;
  if ( double packfrac = numberOfPart / box_sides.volume() * 4./3 * M_PI * radius * radius * radius  >  0.6 ) {
    cout << "The packing fraction is too high ! [" << packfrac << "]" << endl ;
    exit(EXIT_FAILURE) ;
  }

  int keep_trying = 1 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;

  verlet<particle> * verlet_list = new verlet<particle>( this , numberOfPart , 2.0*radius*1.00001 , 2.0*radius*0.00001 ) ;
  verlet_list->initialize_cellist() ;
  struct posizione_3D<int> *particleScell = verlet_list->cells->particleScell ;
  struct clist_element *particle_postIt = verlet_list->cells->particle_postIt ;
  struct cell_array_element ***index = verlet_list->cells->index ;
  struct posizione_3D<int> nmax = verlet_list->cells->nmax , nearCell ;
  struct clist_element *walker = NULL ;

  // generation of initial positions
  cout << endl << "  Generation of random positions" << endl;
  int max_cz = (nmax.z > 1) ? 2 : 1 ;
  int min_cz = (nmax.z > 1) ? -1 : 0 ;
  for( int i=0 ; i<numberOfPart ; i++ ) {
    particles[i] = new particle ;
    particles[i]->id = i ;
    particles[i]->type = 1 ;
    do {
      particles[i]->rndPos_generate( box_sides.position ) ;
      particleScell[i] = ( particles[i]->position / verlet_list->cells->cell_side ).convert2int_floor() ;
      // I start with the exploration of the interactions inside the cell
      keep_trying = 0 ;
      for( int cx=-1; cx<2 && (!keep_trying); cx++ ) {
	for( int cy=-1; cy<2 && (!keep_trying); cy++ ) {
	  for( int cz=min_cz; cz<max_cz && (!keep_trying); cz++ ) {
	    nearCell.x = ( particleScell[i].x + cx + nmax.x ) % nmax.x ; // storing the coordinates of the near cell to analyse
	    nearCell.y = ( particleScell[i].y + cy + nmax.y ) % nmax.y ;
	    nearCell.z = ( particleScell[i].z + cz + nmax.z ) % nmax.z ;
	    walker = index[ nearCell.x ][ nearCell.y ][ nearCell.z ].list ;
	    while( walker != NULL ) {          // this cycle scans the list of a cell to individuate neighbours for i particle
	      right_copy( particles[walker->host] , particles[i] , apparent_near ) ;
	      if( particles[i]->position.distance(apparent_near->position) < 2.0*radius ) keep_trying = 1 ;
	      walker = walker->next ;
	    }
	  }
	}
      }

      if( ! keep_trying ) {
	particle_postIt[i].next = index[ particleScell[i].x ][ particleScell[i].y ][ particleScell[i].z ].list ;
	index[ particleScell[i].x ][ particleScell[i].y ][ particleScell[i].z ].list = ( particle_postIt + i ) ;
      }
    }while( keep_trying == 1 ) ;
  }
  delete apparent_near ;
  delete verlet_list ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_ideal_lattice_3D( const char *ltype ) {               // AGGIORNARE PER IL CASO SCATOLA NON CUBICA
  int inserted_particles = 0, flag = 0, cells_per_side = 0 ;
  double cell_side = 0 ;                  // cell_side = side length of the unit cubic cell

  if( strcmp( ltype , "fcc" ) == 0 ) {
    if( numberOfPart % 4 != 0 ) {
      cout << "*** ERROR: To generate an fcc configuration N needs to be multiple of 4 !" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    cells_per_side = round( pow((double)numberOfPart/4.0, 1.0/3.0) ) ;    // cells_per_side = number of fcc unit cells (4 particles) per side
    cell_side = (double)box_side/cells_per_side ;
    cout << endl << "  Generation of fcc-crystal; cells per side : " << cells_per_side << " , cell's side : " << cell_side << endl;

    flag = 0;
    for(int i=0; i<cells_per_side && flag==0; i++) {
      for(int j=0; j<cells_per_side && flag==0; j++) {
	for(int k=0; k<cells_per_side && flag==0; k++) {
	  if( inserted_particles < numberOfPart ) {            // particle 1
	    particles[inserted_particles] = new particle( cell_side*i , cell_side*j , cell_side*k ) ;
	    particles[inserted_particles]->id = inserted_particles ;
	    particles[inserted_particles]->type = 1 ;
	    inserted_particles++ ;
	  }else{
	    flag=1 ;
	  }
	  if( inserted_particles < numberOfPart ) {            // particle 2
	    particles[inserted_particles] = new particle( cell_side*(i+0.5) , cell_side*(j+0.5) , cell_side*k ) ;
	    particles[inserted_particles]->id = inserted_particles ;
	    particles[inserted_particles]->type = 1 ;
	    inserted_particles++ ;
	  }else{
	    flag=1 ;
	  }
	  if( inserted_particles < numberOfPart ) {            // particle 3
	    particles[inserted_particles] = new particle( cell_side*(i+0.5) , cell_side*j , cell_side*(k+0.5) ) ;
	    particles[inserted_particles]->id = inserted_particles ;
	    particles[inserted_particles]->type = 1 ;
	    inserted_particles++ ;
	  }else{
	    flag=1 ;
	  }
	  if( inserted_particles < numberOfPart ) {            // particle 4
	    particles[inserted_particles] = new particle( cell_side*i , cell_side*(j+0.5) , cell_side*(k+0.5) ) ;
	    particles[inserted_particles]->id = inserted_particles ;
	    particles[inserted_particles]->type = 1 ;
	    inserted_particles++ ;
	  }else{
	    flag=1 ;
	  }
	}
      }
    }

  } else if( strcmp( ltype , "bcc" ) == 0 ) {
    if( numberOfPart % 2 != 0 ) {
      cout << "*** ERROR: To generate a bcc configuration N needs to be multiple of 2 !" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    cells_per_side = round( pow((double)numberOfPart/2.0, 1.0/3.0) ) ;    // cells_per_side = number of bcc unit cells (2 particles) per side
    cell_side = (double)box_side/cells_per_side ;
    cout << endl << "  Generation of bcc-crystal; cells per side : " << cells_per_side << " , cell's side : " << cell_side << endl;

    flag = 0;
    for(int i=0; i<cells_per_side && flag==0; i++) {
      for(int j=0; j<cells_per_side && flag==0; j++) {
	for(int k=0; k<cells_per_side && flag==0; k++) {
	  if( inserted_particles < numberOfPart ) {            // particle 1
	    particles[inserted_particles] = new particle( cell_side*i , cell_side*j , cell_side*k ) ;
	    particles[inserted_particles]->id = inserted_particles ;
	    particles[inserted_particles]->type = 1 ;
	    inserted_particles++ ;
	  }else{
	    flag=1 ;
	  }
	  if( inserted_particles < numberOfPart ) {            // particle 2
	    particles[inserted_particles] = new particle( cell_side*(i+0.5) , cell_side*(j+0.5) , cell_side*(k+0.5) ) ;
	    particles[inserted_particles]->id = inserted_particles ;
	    particles[inserted_particles]->type = 1 ;
	    inserted_particles++ ;
	  }else{
	    flag=1 ;
	  }
	}
      }
    }
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_ideal_lattice( const char *ltype ) {               // AGGIORNARE PER IL CASO SCATOLA NON CUBICA

  if( strcmp( ltype , "square" ) == 0 ) {
    int inserted_particles = 0, flag = 0, cells_per_side = round( pow((double)numberOfPart, 1.0/2.0) ) ; // cells_per_side = number of unit cells per side
    double cell_side=(double)box_side/cells_per_side;                  // cell_side = side length of the unit square cell
    cout << endl << "  Generation of a square lattice; cells per side : " << cells_per_side << " , cell's side : " << cell_side << endl;

    flag = 0;
    for(int i=0; i<cells_per_side && flag==0; i++) {
      for(int j=0; j<cells_per_side && flag==0; j++) {
	if( inserted_particles < numberOfPart ) {            // particle 1
	  particles[inserted_particles] = new particle( cell_side*i , cell_side*j ) ;
	  pbc( particles[inserted_particles] ) ;
	  particles[inserted_particles]->id = inserted_particles ;
	  particles[inserted_particles]->type = 1 ;
	  inserted_particles++ ;
	}else{
	  flag=1 ;
	}
      }
    }

  } else if( strcmp( ltype , "hcp2D" ) == 0 ) {
    int inserted_particles = 0 ;
    double lattice_step = sqrt( box_sides.volume() * 2.0 / sqrt(3) / numberOfPart ) ;
    int ncells_y = 2 * round( sqrt( numberOfPart / sqrt(3) ) ) , ncells_x = 0 ;
    while( numberOfPart % ncells_y != 0 && ncells_y > 2 ) ncells_y -= 2 ;
    ncells_x = numberOfPart / ncells_y ;
    if( ncells_y < 0.8 * 2 * ncells_x / sqrt(3) || numberOfPart % ncells_x != 0 ) {
      cout << "*** ERROR: Please, choose a different particles number" << endl ;
      exit( EXIT_FAILURE ) ;
    } else {
      cout << endl << "  Generation of a hcp 2D lattice; cells per side : " << ncells_x << "," << ncells_y << " ; lattice constant : " << lattice_step << endl ;
    }
    box_sides.position.x = lattice_step * ncells_x ;
    box_sides.position.y = lattice_step * ncells_y * 0.5*sqrt(3.0) ;
    midside.position = (box_sides.position * 0.5) ;

    for(int i=0; i<ncells_x; i++) {
      for(int j=0; j<ncells_y; j++) {
	// particle 1
	particles[inserted_particles] = new particle( lattice_step * ( i + j * 0.5 ) , 0.5*sqrt(3) * j * lattice_step ) ;
	pbc( particles[inserted_particles] ) ;
	particles[inserted_particles]->id = inserted_particles ;
	particles[inserted_particles]->type = 1 ;
	inserted_particles++ ;
      }
    }

  } else {
    cout << "*** ERROR: " << ltype << " is not a recognized lattice type!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_config_by_file( const char *file_name , const char *format ) {
  ifstream data_file ;

  data_file.open( file_name ) ;
  if( data_file.fail() ) {
    cout << "\n  Failure in opening data_file in function generate_config_by_file" << endl ;
    exit( EXIT_FAILURE ) ;

         //  FORMAT :  XYZ
  } else if( strcmp( format , "xyz" ) == 0 ) read_xyz( data_file ) ;
         //  FORMAT :  LAMMPS DUMP
  else if ( strcmp( format , "lmp" ) == 0 ) read_lammps_dump( data_file ) ;
         //  FORMAT :  SPH
  else if ( strcmp( format , "sph" ) == 0 ) read_sph( data_file ) ;
         //  FORMAT :  PTC (3D), PATCH (2D)
  else if ( strcmp( format , "ptc" ) == 0 ||
	    strcmp( format , "patch" ) == 0 ) read_ptc( data_file ) ;
         //  FORMAT :  NICO BONDS or NICO RINGS
  else if( strcmp( format , "nico_bonds" ) == 0 ||
	   strcmp( format , "nico_bonds_pbc" ) == 0 ||
	   strcmp( format , "nico_rings" ) == 0 ||
	   strcmp( format , "nico_rings_pbc" ) == 0 ) read_nico_init( data_file , format ) ;
         //  FORMAT :  UNKNOWN
  else {
    cout << "Configuration loading failed ! You need to specify the file format (xyz|nico_bonds|nico_rings) !" << endl ;
    exit( EXIT_FAILURE ) ;
  }

  if( numberOfPart == 0 ) {
    cout << "Configuration loading failed." << endl ;
    exit( EXIT_FAILURE ) ;
  } else for( int i=0; i<numberOfPart; i++ ) pbc( particles[i] ) ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_molecules_by_file( const char *file_name ) {
  ifstream data_file ;
  int mon_per_ring = 0 ;

  // This first routine on data_file calculate the exact number of molecules
  data_file.open( file_name ) ;
  if( data_file.fail() ) {
    cout << "\n  Failure in opening data_file in function generate_config_by_file" << endl ;
    exit( EXIT_FAILURE ) ;
  } else {

    if( numberOfMol != 0 ) {
      cout << "ERROR: molecules vector already initialized !!" << endl ;
      exit ( EXIT_FAILURE ) ;
    }
    record *line = NULL ;
    line = new record ;
    line->getrecord( data_file ) ;
    while( line->char_number != 0 ) {
      numberOfMol ++ ;
      delete line ;
      line = new record ;
      line->getrecord( data_file ) ;
    }
    mon_per_ring = numberOfPart / numberOfMol ;
    if( mon_per_ring * numberOfMol != numberOfPart ) cout << "warning: rings seem not to have all the same number of monomers !!" << endl ;
    molecules_length = numberOfMol ;
    molecules = new molecule<particle> * [ molecules_length ] ;
  }
  data_file.close() ;
  data_file.clear() ;

  // Now I read the polydispersity of the monomers diameter in each molecule
  data_file.open( file_name ) ;
  if( data_file.fail() ) {
    cout << "\n  Failure in opening data_file in function generate_config_by_file" << endl ;
    exit( EXIT_FAILURE ) ;
  } else {

    int bad_flag = 0 , ordering_flag = 0 , molecule_index = -1 ;
    double monomers_diameter = 0 ;
    record *line = NULL ;
    line = new record ;
    line->getrecord( data_file ) ;
    for( int i=0 ; i<molecules_length ; i++ ) {
      line->split() ;
      if( line->words_number < 2 ) {
	cout << "ERROR: the file containing polydispersities is in a wrong format !!" << endl ;
	exit( EXIT_FAILURE ) ;
      } else if( line->words_number > 2 ) bad_flag = 1 ;
      if( bad_flag == 1 ) cout << "Warning: the file containing polydispersities could be in a wrong format !!" << endl ;
      sscanf( line->word[0] , "%d" , &molecule_index ) ;
      if( i != molecule_index-1 ) ordering_flag = 1 ;
      if( ordering_flag == 1 ) cout << "Warning: the molecules in the file are not ordered !!" << endl ;
      molecules[ molecule_index-1 ] = new molecule<particle> ;
      molecules[ molecule_index-1 ]->Natoms = mon_per_ring ;
      molecules[ molecule_index-1 ]->atom = new particle * [ mon_per_ring ] ;
      sscanf( line->word[1] , "%lf" , &monomers_diameter ) ;
      molecules[ molecule_index-1 ]->monomers_diameter = monomers_diameter ;
      for( int j=0 ; j<molecules[ molecule_index-1 ]->Natoms ; j++ ) {
	molecules[ molecule_index-1 ]->atom[j] = particles[ (molecule_index-1)*mon_per_ring + j ] ;
	molecules[ molecule_index-1 ]->atom[j]->radius = 0.5 * monomers_diameter ;
      }
      delete line ;
      line = new record ;
      line->getrecord( data_file ) ;
    }
    delete line ;
    line = NULL ;

  }
  data_file.close() ;

  if( numberOfMol == 0 ) {
    cout << "Molecules loading failed." << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/
