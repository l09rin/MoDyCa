/* Here the constructor, distructors and main methods of the class verlet for the verlet lists are defined */
#include "molecules.h"

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
template <typename particle>
verlet<particle>::verlet() {
  cout << "Do not use the default constructor for the verlet() list, please." << endl ;
  exit( EXIT_FAILURE ) ;
}
/********************************************************************************************/

template <typename particle>
verlet<particle>::verlet( configuration<particle> *conf_ptr ) {
  DOUBLE_PAIRS = 0 ;
  if( conf_ptr == NULL ) {
    cout << "ERROR in verlet() constructor: The pointer to God cannot be NULL!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  god = conf_ptr ;
  cells = NULL ;
  //  NearCells2explore = NULL ;
  //  N_of_NearCells2explore = 0 ;
  r_verlet = 0 ;
  delta_verlet = 0 ;
  particles_in_list = 0 ;
  arrays_length = 0 ;
  element = NULL ;
  must_be_updated = 0 ;
  disp_on = 0 ;
  mols = NULL ;
  Nmols = 0 ;
  mols_length = 0 ;
}
/********************************************************************************************/

template <typename particle>
verlet<particle>::verlet( configuration<particle> *conf_ptr , int initial_length , double verlet_ray , double verlet_delta , const char *disp_mode ) {
  DOUBLE_PAIRS = 0 ;
  if( conf_ptr == NULL ) {
    cout << "ERROR in verlet() constructor: The pointer to God cannot be NULL!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  god = conf_ptr ;
  cells = NULL ;
  r_verlet = verlet_ray ;
  delta_verlet = verlet_delta ;
  particles_in_list = initial_length ;
  arrays_length = initial_length ;
  if( arrays_length == 0 ) {
    element = NULL ;
    disp_on = 0 ;
  }else{
    element = (vlist_element<particle> *)calloc( arrays_length , sizeof(vlist_element<particle>) ) ;
    if( check_memalloc( element , "Error in verlet constructor" ) ) exit( EXIT_FAILURE ) ;
    for( int i=0 ; i<arrays_length ; i++ ) {
      element[i].ptr = NULL ;
      element[i].cell_center = NULL ;
      element[i].neighbours = NULL ;
    }
    if( disp_mode == NULL ) {
      for( int i=0 ; i<arrays_length ; i++ ) element[i].cell_center = new particle ;
      disp_on = 1 ;
      must_be_updated = 1 ;
    }else if( strcmp( disp_mode, "NO_DISPLACEMENT_VECTOR" ) == 0 ) {
      disp_on = 0 ;
    }else{
      cout << " ERROR in verlet creation function !" << endl ;
      fflush(0) ;
      exit( EXIT_FAILURE ) ;
    }
  }
  mols = NULL ;
  Nmols = 0 ;
  mols_length = 0 ;
}
/********************************************************************************************/

template <typename particle>
verlet<particle>::verlet( const verlet &other_verlet ) {
  DOUBLE_PAIRS = other_verlet->DOUBLE_PAIRS ;
  //  You can copy the list only within the same configuration !!
  god = other_verlet.god ;
  cells = NULL ;
  r_verlet = other_verlet.r_verlet ;
  delta_verlet = other_verlet.delta_verlet ;
  particles_in_list = other_verlet.particles_in_list ;
  arrays_length = other_verlet.arrays_length ;
  must_be_updated = other_verlet.must_be_updated ;
  disp_on = other_verlet.disp_on ;

  element = (vlist_element<particle> *)calloc( arrays_length , sizeof(vlist_element<particle>) ) ;
  if( check_memalloc( element , "Error in verlet constructor" ) ) exit( EXIT_FAILURE ) ;
  for( int i=0 ; i<arrays_length ; i++ ) {
    element[i].ptr = other_verlet.element[i].ptr ;
    element[i].cell_center = NULL ;
    copy_neighbour_list( element+i , other_verlet.element[i].neighbours ) ;
  }
  if( disp_on ) {
    for( int i=0 ; i<arrays_length ; i++ ) {
      element[i].cell_center = new particle ;
      element[i].cell_center = other_verlet.element[i].cell_center ;
    }
  }

  Nmols = other_verlet.Nmols ;
  mols_length = Nmols ;
  mols = new molecule<particle>* [ mols_length ] ;
  for( int i=0 ; i<mols_length ; i++ ) mols[i] = other_verlet.mols[i] ;
}
/********************************************************************************************/

template <typename particle>
verlet<particle>::~verlet() {
  clear() ;
  for( int i=0 ; i<arrays_length ; i++ ) {
    element[i].ptr = NULL ;
    if( element[i].cell_center != NULL ) {
      delete element[i].cell_center ;
      element[i].cell_center = NULL ;
    }
  }
  delete [] mols ;
  mols = NULL ;
  Nmols = 0 ;
  mols_length = 0 ;
  if( element != NULL ) free( element ) ;
  element = NULL ;
  if( cells != NULL ) delete cells ;
  cells = NULL ;
  arrays_length = 0 ;
  god = NULL ;
}
/********************************************************************************************/
/****************************************************************************************/


template <typename particle>
inline void verlet<particle>::initialize_particles_list( const char *kind ) {
  strcpy( INITIALIZATION_MODE , kind ) ;
  if( kind == NULL || strcmp(kind, "ORDERED") == 0 ) {
    for( int i=0 ; i<particles_in_list ; i++ ) element[i].ptr = god->particles[i] ;
  }else if( strcmp(kind, "BONDS") == 0 ) {
    int j = 0 ;
    for( int i=0 ; i<god->numberOfPart ; i++ ) {
      if( god->particles[i]->valence > 0 ) {
	j ++ ;
	if( j > particles_in_list ) {
	  cout << "ERROR: Mismatching in the number of bonded atoms (in function initialize_particles_list()) !" << endl ;
	  exit( EXIT_FAILURE ) ;
	}
	element[j-1].ptr = god->particles[i] ;
      }
    }
  }else if( strcmp(kind, "CHARGES") == 0 ) {
    if( particles_in_list == 0 ) {
      cout << "  ERROR: The number of charges cannot be 0 if you want to simulate a charged microgel!!" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    int j = 0 ;
    for( int i=0 ; i<god->numberOfPart ; i++ ) {
      if( god->particles[i]->charge != 0 ) {
	j ++ ;
	if( j > particles_in_list ) {
	  cout << "ERROR: Mismatching in the number of charged atoms (in function initialize_particles_list()) !" << endl ;
	  exit( EXIT_FAILURE ) ;
	}
	element[j-1].ptr = god->particles[i] ;
      }
    }
  }else{
    cout << "ERROR: Unrecognized type of listing in function initialize_particles_list() !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void verlet<particle>::initialize_molecules_list( const char *kind ) {
  DOUBLE_PAIRS = 0 ;
  if( kind == NULL ) {
    mols_length = god->molecules_length ;
    Nmols = god->numberOfMol ;
    mols = new molecule<particle> * [ mols_length ] ;
    for( int i=0 ; i<Nmols ; i++ ) mols[i] = god->molecules[i] ;
    for( int i=Nmols ; i<mols_length ; i++ ) mols[i] = NULL ;

  } else {
    ifstream data_file ;
    data_file.open( kind ) ;
    if( data_file.fail() ) {
      cout << "\n  Failure in opening molecules file in function verlet::initialize_molecules_list() !!" << endl ;
      exit( EXIT_FAILURE ) ;
    } else {
      int flag = 0 , molecule_index = -1 ;
      molecule<particle> **service_pointer = NULL ;
      record *line = NULL ;
      line = new record ;
      line->getrecord( data_file ) ;
      while( line->char_number != 0 ) {                       // Reading monomers' bonding properties
	line->split() ;
	if( line->words_number > 1 ) flag = 1 ;      // control on the format of the file
	if( flag == 1 ) cout << "warning: the file could contain wrongly formatted information !!" << endl ;
	sscanf( line->word[0] , "%d" , &molecule_index ) ;    // reading the molecule index
	Nmols ++ ;
	if( Nmols >= mols_length ) {      // enlargement of the molecule pointers list
	  mols_length += 100 ;
	  service_pointer = new molecule<particle> * [ mols_length ] ;
	  if( check_memalloc( service_pointer , "Allocation error 1 in function verlet::initialize_molecules_list()" ) ) exit( EXIT_FAILURE ) ;
	  for( int i=0 ; i<Nmols-1 ; i++ ) service_pointer[i] = mols[i] ;
	  delete [] mols ;
	  mols = service_pointer ;
	  service_pointer = NULL ;
	}
	mols[Nmols-1] = god->molecules[molecule_index-1] ;
	delete line ;
	line = new record ;
	line->getrecord( data_file ) ;
      }
      delete line ;
      line = NULL ;
      mols_length = Nmols ;
      service_pointer = new molecule<particle> * [ mols_length ] ;
      if( check_memalloc( service_pointer , "Allocation error 2 in function verlet::initialize_molecules_list()" ) ) exit( EXIT_FAILURE ) ;
      for( int i=0 ; i<Nmols ; i++ ) service_pointer[i] = mols[i] ;
      delete [] mols ;
      mols = service_pointer ;
      service_pointer = NULL ;

    }
  }
}
/****************************************************************************************/


// This version of the function works only if particles coordinates are always in the box
template <typename particle>
inline void verlet<particle>::add_connection_between( int i , int j ) {
  neighbour_list<particle> *newNear = NULL ;
                     //  I insert j in the i-list ONLY
  newNear = (neighbour_list<particle> *)calloc( 1 , sizeof(neighbour_list<particle>) ) ;
  if( check_memalloc( newNear , "Allocation failed! newNear not allocated in function add_connection_between()" ) ) exit( EXIT_FAILURE ) ;
  newNear->ptr = element[j].ptr ;
  newNear->next = element[i].neighbours ;
  element[i].neighbours = newNear ;
  newNear = NULL ;
}
/****************************************************************************************/

template <typename particle>
inline void verlet<particle>::clear( void ) {
  if( ! disp_on ) {
    for( int i=0 ; i<particles_in_list ; i++ ) {
      free_neigh_list( element[i].neighbours ) ;
      element[i].neighbours = NULL ;
    }
  } else {
    for( int i=0 ; i<particles_in_list ; i++ ) {
      // element[i].cell_center->clear() ;
      free_neigh_list( element[i].neighbours ) ;
      element[i].neighbours = NULL ;
    }
  }
}
/****************************************************************************************/

template <typename particle>
inline void verlet<particle>::free_neigh_list( neighbour_list<particle> *lista ) {
  if( lista != NULL ) {
    neighbour_list<particle> *walker ;
    while( lista != NULL ) {
      walker = lista->next ;
      free( lista ) ;
      lista = walker ;
    }
  }
}
/****************************************************************************************/

template <typename particle>
void verlet<particle>::copy_neighbour_list( vlist_element<particle> *gancio , neighbour_list<particle> *lista ) {
          // return the pointer to a copy of lista
  if( lista == NULL ) {
    gancio->neighbours = NULL ;
  }else{
    neighbour_list<particle> *walker = lista , *new_element = NULL ;
    while( walker != NULL ) {
      new_element = (neighbour_list<particle> *)calloc( 1 , sizeof(neighbour_list<particle>) ) ;
      new_element->ptr = walker->ptr ;
      new_element->next = gancio->neighbours ;
      gancio->neighbours = new_element ;
      new_element = NULL ;
      walker = walker->next ;
    }
  }
}
/****************************************************************************************/

template <typename particle>
void verlet<particle>::bonding_verlet_creation( const char *simulation_kind ) {
  if( strcmp( simulation_kind, "MD" ) == 0 ) {
    DOUBLE_PAIRS = 0 ;

    int near_id = -1 ;
    for( int i=0 ; i<particles_in_list ; i++ ) {
      for( int j=0 ; j<element[i].ptr->valence ; j++ ) {
	near_id = element[i].ptr->bonded_monomer[j] ;
	for( int k=i+1 ; k<particles_in_list ; k++ ) {
	  if( element[k].ptr->id == near_id ) add_connection_between( i , k ) ;
	}
      }
    }

  }else if( strcmp( simulation_kind, "MC" ) == 0 ) {
    DOUBLE_PAIRS = 1 ;

    int near_id = -1 ;
    for( int i=0 ; i<particles_in_list ; i++ ) {
      for( int j=0 ; j<element[i].ptr->valence ; j++ ) {
	near_id = element[i].ptr->bonded_monomer[j] ;
	for( int k=i+1 ; k<particles_in_list ; k++ ) {
	  if( element[k].ptr->id == near_id ) {
	    add_connection_between( i , k ) ;
	    add_connection_between( k , i ) ;
	  }
	}
      }
    }

  }else{
    cout << "ERROR: You must specify the kind of simulation (MC|MD) when using function bonding_verlet_creation() !!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/
