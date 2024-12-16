cell_list::cell_list() {
  particles_in_list = 0 ;
  index = NULL ;           // tridimensional array of lists of particles contained in each cell
  particleScell = NULL ;   // array containing coordinates of the cell wherein each particle is located
  particle_postIt = NULL ; // array of the list elements associated to each particle
  nmax.x = 1 ;             // number of pieces in which is divided each side of the box
  nmax.y = 1 ;
  nmax.z = 1 ;
  NearCells2explore = NULL ; // array of triplets containing the relative position of the cells to explore during the verlet list construction
  N_of_NearCells2explore = 0 ;
  cell_side.x = 0 ;                    // length of side of each cell
  cell_side.y = 0 ;
  cell_side.z = 0 ;

};
/****************************************************************************************/

//  THIS WILL WORK ONLY FOR PBCs
template< typename T >
cell_list::cell_list( T &box_sides_vec , double r_verlet , int parts_in_list ) {
  nmax = ( box_sides_vec / r_verlet ).convert2int_floor() ;
  if( nmax.x < 4 ) nmax.x = 1 ;
  if( nmax.y < 4 ) nmax.y = 1 ;
  if( nmax.z < 4 ) nmax.z = 1 ;

  if( nmax.x == 1 && nmax.y == 1 && nmax.z == 1 ) {
    particles_in_list = 0 ;
    index = NULL ;
    particleScell = NULL ;
    particle_postIt = NULL ;
    NearCells2explore = NULL ;
    N_of_NearCells2explore = 0 ;
    cell_side.x = 0 ;
    cell_side.y = 0 ;
    cell_side.z = 0 ;

  }else{
    particles_in_list = parts_in_list ;
    cell_side = ( box_sides_vec / nmax ) ;
    set_nearCells2explore() ;

    index = (struct cell_array_element ***)calloc( nmax.x, sizeof(struct cell_array_element **) ) ;        // allocation of the cell index
    if( check_memalloc( index , "error in allocation of cells->index_1 in initialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
    for( int i=0; i<nmax.x; i++ ) {          // initialization of the 3D array of cells
      index[i] = (struct cell_array_element **)calloc( nmax.y, sizeof(struct cell_array_element *) ) ;
      if( check_memalloc( index[i] , "\n error in allocation of cells->index_2 in initialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
      for( int j=0; j<nmax.y; j++ ) {
	index[i][j] = (struct cell_array_element *)calloc( nmax.z, sizeof(struct cell_array_element) ) ;
	if( check_memalloc( index[i][j] , "\n error in allocation of cells->index_3 in initialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
	for( int k=0; k<nmax.z; k++ ) {            // initialization to NULL value of the cell index
	  index[i][j][k].list = NULL ;
	  index[i][j][k].cell_is_explored = 0 ;
	}
      }
    }

    if( particles_in_list == 0 ) {                        // initialization of the vector containing the coordinates of the particle's cell
      particleScell = NULL ;
      particle_postIt = NULL ;
    }else{                                      // initialization of the vector containing the cell in which a particle stay
      particleScell = (struct posizione_3D<int> *)calloc( particles_in_list , sizeof(struct posizione_3D<int>) ) ;
      if( check_memalloc( particleScell , "error in allocation of cell->particleScell in initialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
      particle_postIt = (struct clist_element *)calloc( particles_in_list , sizeof(struct clist_element) ) ;
      if( check_memalloc( particle_postIt , "error in allocation of cell->particle_postIt in initialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
      for( int i=0 ; i<particles_in_list ; i++ ) {
	particle_postIt[i].host = i ;            // initialization of the vector containing the elements of the cell list
	particle_postIt[i].next = NULL ;
      }
    }

  }
}
/****************************************************************************************/

cell_list::cell_list( const cell_list &other ) {
  particles_in_list = 0 ;
  index = NULL ;
  particleScell = NULL ;
  particle_postIt = NULL ;
  nmax.x = other.nmax.x ;
  nmax.y = other.nmax.y ;
  nmax.z = other.nmax.z ;
  NearCells2explore = NULL ;
  N_of_NearCells2explore = 0 ;
  cell_side.x = other.cell_side.x ;
  cell_side.y = other.cell_side.y ;
  cell_side.z = other.cell_side.z ;

  cout << "*** ERROR: Coding of copy constructor for the class cell_list has not been yet completed!!" << endl ;
  exit( EXIT_FAILURE ) ;
};
/****************************************************************************************/

cell_list::~cell_list() {
  clear() ;      // clean all pointers in cell structure
  free_index() ;

  free( particleScell ) ;
  particleScell = NULL ;
  free( particle_postIt ) ;
  particle_postIt = NULL ;
  free( NearCells2explore ) ;
  NearCells2explore = NULL ;
  N_of_NearCells2explore = 0 ;

  cell_side.x = 0 ;
  cell_side.y = 0 ;
  cell_side.z = 0 ;
  nmax.x = 1 ;
  nmax.y = 1 ;
  nmax.z = 1 ;
  particles_in_list = 0 ;
};
/****************************************************************************************/

/****************************************************************************************/

template < typename T >            //  THIS WILL WORK ONLY FOR PBCs
void cell_list::reinitialize_cellist_index( T &box_sides_vec , double r_verlet  ) {
  struct posizione_3D<int> act_nmax = ( box_sides_vec / r_verlet ).convert2int_floor() ;
  if( act_nmax.x < 4 ) act_nmax.x = 1 ;
  if( act_nmax.y < 4 ) act_nmax.y = 1 ;
  if( act_nmax.z < 4 ) act_nmax.z = 1 ;
  if( act_nmax == nmax ) cell_side = ( box_sides_vec / nmax ) ;
  else {  // . . . or if we need to rebuild the array of cells

    if( act_nmax.x == 1 && act_nmax.y == 1 && act_nmax.z == 1 ) {
      free_index() ;
      nmax = act_nmax ;
    } else {
      // if the topology of the space changes, I also change the way to explore cells
      bool TOPOLOGY_CHANGE = ( (act_nmax.x==1) + (act_nmax.y==1) + (act_nmax.z==1) ) != ( (nmax.x==1) + (nmax.y==1) + (nmax.z==1) ) ;

      // reallocation of the cells list, in case their number changes
      free_index() ;
      nmax = act_nmax ;
      cell_side = ( box_sides_vec / nmax ) ;

      if( TOPOLOGY_CHANGE ) {
	free( NearCells2explore ) ;
	NearCells2explore = NULL ;
	set_nearCells2explore() ;
      }

      index = (struct cell_array_element ***)calloc( nmax.x, sizeof(struct cell_array_element **) ) ;        // allocation of the cell index
      if( check_memalloc( index , "error in allocation of cells->index_1 in reinitialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
      for( int i=0; i<nmax.x; i++ ) {          // initialization of the 3D array of cells
	index[i] = (struct cell_array_element **)calloc( nmax.y, sizeof(struct cell_array_element *) ) ;
	if( check_memalloc( index[i] , "\n error in allocation of cells->index_2 in reinitialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
	for( int j=0; j<nmax.y; j++ ) {
	  index[i][j] = (struct cell_array_element *)calloc( nmax.z, sizeof(struct cell_array_element) ) ;
	  if( check_memalloc( index[i][j] , "\n error in allocation of cells->index_3 in reinitialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
	  for( int k=0; k<nmax.z; k++ ) {            // initialization to NULL value of the cell index
	    index[i][j][k].list = NULL ;
	    index[i][j][k].cell_is_explored = 0 ;
	  }
	}
      }
    }

  }
}
/****************************************************************************************/

//  THIS WILL WORK ONLY FOR PBCs
void cell_list::set_nearCells2explore( void ) {
  if( nmax.x == 1 && nmax.y == 1 ) {
    N_of_NearCells2explore = 1 ;
    NearCells2explore = (struct posizione_3D<int> *)calloc( N_of_NearCells2explore, sizeof(struct posizione_3D<int>) ) ;
    if( check_memalloc( NearCells2explore , "error 1 in allocation of NearCells2explore in set_nearCells2explore()" ) ) exit( EXIT_FAILURE ) ;
    NearCells2explore[0].x = 0 , NearCells2explore[0].y = 0 , NearCells2explore[0].z = 1 ;
  }else if( nmax.x == 1 && nmax.z == 1 ) {
    N_of_NearCells2explore = 1 ;
    NearCells2explore = (struct posizione_3D<int> *)calloc( N_of_NearCells2explore, sizeof(struct posizione_3D<int>) ) ;
    if( check_memalloc( NearCells2explore , "error 2 in allocation of NearCells2explore in set_nearCells2explore()" ) ) exit( EXIT_FAILURE ) ;
    NearCells2explore[0].x = 0 , NearCells2explore[0].y = 1 , NearCells2explore[0].z = 0 ;
  }else if( nmax.y == 1 && nmax.z == 1 ) {
    N_of_NearCells2explore = 1 ;
    NearCells2explore = (struct posizione_3D<int> *)calloc( N_of_NearCells2explore, sizeof(struct posizione_3D<int>) ) ;
    if( check_memalloc( NearCells2explore , "error 3 in allocation of NearCells2explore in set_nearCells2explore()" ) ) exit( EXIT_FAILURE ) ;
    NearCells2explore[0].x = 1 , NearCells2explore[0].y = 0 , NearCells2explore[0].z = 0 ;

  }else if( nmax.z == 1 ) {
    N_of_NearCells2explore = 4 ;
    NearCells2explore = (struct posizione_3D<int> *)calloc( N_of_NearCells2explore, sizeof(struct posizione_3D<int>) ) ;
    if( check_memalloc( NearCells2explore , "error 4 in allocation of NearCells2explore in set_nearCells2explore()" ) ) exit( EXIT_FAILURE ) ;
    NearCells2explore[0].x = 0 , NearCells2explore[0].y = 1 , NearCells2explore[0].z = 0 ;
    NearCells2explore[1].x = 1 , NearCells2explore[1].y = 1 , NearCells2explore[1].z = 0 ;
    NearCells2explore[2].x = 1 , NearCells2explore[2].y = 0 , NearCells2explore[2].z = 0 ;
    NearCells2explore[3].x = 1 , NearCells2explore[3].y = -1 , NearCells2explore[3].z = 0 ;
  }else if( nmax.y == 1 ) {
    N_of_NearCells2explore = 4 ;
    NearCells2explore = (struct posizione_3D<int> *)calloc( N_of_NearCells2explore, sizeof(struct posizione_3D<int>) ) ;
    if( check_memalloc( NearCells2explore , "error 5 in allocation of NearCells2explore in set_nearCells2explore()" ) ) exit( EXIT_FAILURE ) ;
    NearCells2explore[0].x = 0 , NearCells2explore[0].y = 0 , NearCells2explore[0].z = 1 ;
    NearCells2explore[1].x = 1 , NearCells2explore[1].y = 0 , NearCells2explore[1].z = 1 ;
    NearCells2explore[2].x = 1 , NearCells2explore[2].y = 0 , NearCells2explore[2].z = 0 ;
    NearCells2explore[3].x = 1 , NearCells2explore[3].y = 0 , NearCells2explore[3].z = -1 ;
  }else if( nmax.x == 1 ) {
    N_of_NearCells2explore = 4 ;
    NearCells2explore = (struct posizione_3D<int> *)calloc( N_of_NearCells2explore, sizeof(struct posizione_3D<int>) ) ;
    if( check_memalloc( NearCells2explore , "error 6 in allocation of NearCells2explore in set_nearCells2explore()" ) ) exit( EXIT_FAILURE ) ;
    NearCells2explore[0].x = 0 , NearCells2explore[0].y = 0 , NearCells2explore[0].z = 1 ;
    NearCells2explore[1].x = 0 , NearCells2explore[1].y = 1 , NearCells2explore[1].z = 1 ;
    NearCells2explore[2].x = 0 , NearCells2explore[2].y = 1 , NearCells2explore[2].z = 0 ;
    NearCells2explore[3].x = 0 , NearCells2explore[3].y = 1 , NearCells2explore[3].z = -1 ;

  }else{
    N_of_NearCells2explore = 13 ;
    NearCells2explore = (struct posizione_3D<int> *)calloc( N_of_NearCells2explore, sizeof(struct posizione_3D<int>) ) ;
    if( check_memalloc( NearCells2explore , "error 7 in allocation of NearCells2explore in set_nearCells2explore()" ) ) exit( EXIT_FAILURE ) ;
    NearCells2explore[0].x = 0 , NearCells2explore[0].y = 0 , NearCells2explore[0].z = 1 ;
    NearCells2explore[1].x = -1 , NearCells2explore[1].y = 1 , NearCells2explore[1].z = 1 ;
    NearCells2explore[2].x = -1 , NearCells2explore[2].y = 1 , NearCells2explore[2].z = 0 ;
    NearCells2explore[3].x = -1 , NearCells2explore[3].y = 1 , NearCells2explore[3].z = -1 ;
    NearCells2explore[4].x = 0 , NearCells2explore[4].y = 1 , NearCells2explore[4].z = 1 ;
    NearCells2explore[5].x = 0 , NearCells2explore[5].y = 1 , NearCells2explore[5].z = 0 ;
    NearCells2explore[6].x = 0 , NearCells2explore[6].y = 1 , NearCells2explore[6].z = -1 ;
    NearCells2explore[7].x = 1 , NearCells2explore[7].y = 1 , NearCells2explore[7].z = 1 ;
    NearCells2explore[8].x = 1 , NearCells2explore[8].y = 1 , NearCells2explore[8].z = 0 ;
    NearCells2explore[9].x = 1 , NearCells2explore[9].y = 1 , NearCells2explore[9].z = -1 ;
    NearCells2explore[10].x = 1 , NearCells2explore[10].y = 0 , NearCells2explore[10].z = 1 ;
    NearCells2explore[11].x = 1 , NearCells2explore[11].y = 0 , NearCells2explore[11].z = 0 ;
    NearCells2explore[12].x = 1 , NearCells2explore[12].y = 0 , NearCells2explore[12].z = -1 ;
  }
}
/****************************************************************************************/

void cell_list::free_index() {
  for( int i=0 ; i<nmax.x ; i++ ) {
    for( int j=0 ; j<nmax.y ; j++ ) {
      free( index[i][j] ) ;
      index[i][j] = NULL ;
    }
    free( index[i] ) ;
    index[i] = NULL ;
  }
  free( index ) ;
  index = NULL ;
}
/****************************************************************************************/

inline void cell_list::clear( void ) {
  for( int i=0 ; i<nmax.x ; i++ ) {
    for( int j=0 ; j<nmax.y ; j++ ) {
      for( int k=0 ; k<nmax.z ; k++ ) {
	index[i][j][k].list = NULL ;
	index[i][j][k].cell_is_explored = 0 ;
      }
    }
  }

  for( int i=0 ; i<particles_in_list ; i++ ) {
    particle_postIt[i].next = NULL ;
  }
}
/****************************************************************************************/
