/* this library file contains the definitions of all the functions related to the construction of the verlet list from a cell list */
#include "cells.h"

/****************************************************************************************/

template <typename particle>
inline void verlet<particle>::verlet_creation_bycell( const char *simulation_kind ) {
  if( disp_on )  for( int i=0 ; i<particles_in_list ; i++ )  element[i].cell_center->position = element[i].ptr->position ;

  if( strcmp( simulation_kind, "MD" ) == 0 ) {
    DOUBLE_PAIRS = 0 ;

    particle *j_apparent_particle = NULL ;
    j_apparent_particle = new particle ;

    if( cells == NULL ) {
      for( int i=0 ; i<particles_in_list - 1 ; i++ ) {
	for( int j=i+1 ; j<particles_in_list ; j++ ) {
	  god->right_copy( element[j].ptr , element[i].ptr , j_apparent_particle ) ;
	  if( element[i].ptr->position.distance(j_apparent_particle->position) < r_verlet ) {  // if i-j distance is less than verlet radius I insert the connection
	    add_connection_between( i , j ) ;
	  }
	}
      }
    }else{
      cells->clear() ;
      create_cellist() ;
      struct posizione_3D<int> nmax = cells->nmax , nearCell ;
      struct posizione_3D<int> *currentCell = cells->particleScell ;
      struct cell_array_element ***cells_index = cells->index ;
      struct clist_element *first_walker = NULL , *second_walker = NULL ;
      for( int i=0 ; i<particles_in_list ; i++ ) {       // this cycle scan the particle to insert in the verlet list
	if( cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].cell_is_explored == 0 ) {
	  // I start with the exploration of the interactions inside the cell
	  first_walker = cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].list ;
	  while( first_walker != NULL ) {          // this cycle scans the list of a cell to individuate neighbours for i particle
	    second_walker = first_walker->next ;       // for avoiding double counting
	    while( second_walker != NULL ) {
	      god->right_copy( element[second_walker->host].ptr , element[first_walker->host].ptr , j_apparent_particle ) ;      // I put in point2 the scaled coordinates of the nearest copy of the neighbour
	      if( element[first_walker->host].ptr->position.distance(j_apparent_particle->position) < r_verlet ) {        // if i-j distance is less than verlet radius I insert the connection
		add_connection_between( first_walker->host , second_walker->host ) ;
	      }
	      second_walker = second_walker->next ;
	    }
	    first_walker = first_walker->next ;
	  }
	  // Now I look for the interactions among particles in the cell and particles in the neighbouring ones
	  for( int j=0 ; j<cells->N_of_NearCells2explore ; j++ ) {       // for avoiding double counting
	    nearCell.x = ( currentCell[i].x + cells->NearCells2explore[j].x + nmax.x ) % nmax.x ; // storing the coordinates of the near cell to analyse
	    nearCell.y = ( currentCell[i].y + cells->NearCells2explore[j].y + nmax.y ) % nmax.y ;
	    nearCell.z = ( currentCell[i].z + cells->NearCells2explore[j].z + nmax.z ) % nmax.z ;
	    first_walker = cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].list ;
	    while( first_walker != NULL ) {          // this cycle scans the list of a cell to individuate neighbours for i particle
	      second_walker = cells_index[ nearCell.x ][ nearCell.y ][ nearCell.z ].list ;
	      while( second_walker != NULL ) {          // this cycle scans the list of a cell to individuate neighbours for i particle
		god->right_copy( element[second_walker->host].ptr , element[first_walker->host].ptr , j_apparent_particle ) ;      // I put in point2 the scaled coordinates of the nearest copy of the neighbour
		if( element[first_walker->host].ptr->position.distance(j_apparent_particle->position) < r_verlet ) {        // if i-j distance is less than verlet radius I insert the connection
		  add_connection_between( first_walker->host , second_walker->host ) ;
		}
		second_walker = second_walker->next ;
	      }
	      first_walker = first_walker->next ;
	    }
	  }
	  cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].cell_is_explored = 1 ;
	}
      }
    }

    delete j_apparent_particle ;

  }else if( strcmp( simulation_kind, "MC" ) == 0 ) {
    DOUBLE_PAIRS = 1 ;

    particle *j_apparent_particle = NULL ;
    j_apparent_particle = new particle ;

    if( cells == NULL ) {
      for( int i=0 ; i<particles_in_list - 1 ; i++ ) {
	for( int j=i+1 ; j<particles_in_list ; j++ ) {
	  god->right_copy( element[j].ptr , element[i].ptr , j_apparent_particle ) ;
	  if( element[i].ptr->position.distance(j_apparent_particle->position) < r_verlet ) {  // if i-j distance is less than verlet radius I insert the connection
	    add_connection_between( i , j ) ;
	    add_connection_between( j , i ) ;
	  }
	}
      }

    } else {
      cells->clear() ;
      create_cellist() ;
      struct posizione_3D<int> nmax = cells->nmax, nearCell ;
      struct posizione_3D<int> *currentCell = cells->particleScell ;
      struct cell_array_element ***cells_index = cells->index ;
      struct clist_element *first_walker = NULL , *second_walker = NULL ;
      for( int i=0 ; i<particles_in_list ; i++ ) {       // this cycle scan the particle to insert in the verlet list
	if( cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].cell_is_explored == 0 ) {
	  // I start with the exploration of the interactions inside the cell
	  first_walker = cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].list ;
	  while( first_walker != NULL ) {          // this cycle scans the list of a cell to individuate neighbours for i particle
	    second_walker = first_walker->next ;       // for avoiding double counting
	    while( second_walker != NULL ) {
	      god->right_copy( element[second_walker->host].ptr , element[first_walker->host].ptr , j_apparent_particle ) ;      // I put in point2 the scaled coordinates of the nearest copy of the neighbour
	      if( element[first_walker->host].ptr->position.distance(j_apparent_particle->position) < r_verlet ) {        // if i-j distance is less than verlet radius I insert the connection
		add_connection_between( first_walker->host , second_walker->host ) ;
		add_connection_between( second_walker->host , first_walker->host ) ;
	      }
	      second_walker = second_walker->next ;
	    }
	    first_walker = first_walker->next ;
	  }
	  // Now I look for the interactions among particles in the cell and particles in the neighbouring ones
	  for( int j=0 ; j<cells->N_of_NearCells2explore ; j++ ) {
	    nearCell.x = ( currentCell[i].x + cells->NearCells2explore[j].x + nmax.x ) % nmax.x ; // storing the coordinates of the near cell to analyse
	    nearCell.y = ( currentCell[i].y + cells->NearCells2explore[j].y + nmax.y ) % nmax.y ;
	    nearCell.z = ( currentCell[i].z + cells->NearCells2explore[j].z + nmax.z ) % nmax.z ;
	    first_walker = cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].list ;
	    while( first_walker != NULL ) {          // this cycle scans the list of a cell to individuate neighbours for i particle
	      second_walker = cells_index[ nearCell.x ][ nearCell.y ][ nearCell.z ].list ;
	      while( second_walker != NULL ) {          // this cycle scans the list of a cell to individuate neighbours for i particle
		god->right_copy( element[second_walker->host].ptr , element[first_walker->host].ptr , j_apparent_particle ) ;      // I put in point2 the scaled coordinates of the nearest copy of the neighbour
		if( element[first_walker->host].ptr->position.distance(j_apparent_particle->position) < r_verlet ) {        // if i-j distance is less than verlet radius I insert the connection
		  add_connection_between( first_walker->host , second_walker->host ) ;
		  add_connection_between( second_walker->host , first_walker->host ) ;
		}
		second_walker = second_walker->next ;
	      }
	      first_walker = first_walker->next ;
	    }
	  }
	  cells_index[ currentCell[i].x ][ currentCell[i].y ][ currentCell[i].z ].cell_is_explored = 1 ;
	}
      }
    }

    delete j_apparent_particle ;
  }else{
    cout << "ERROR: You must specify the kind of simulation (MC|MD) when using function verlet_creation_bycell() !!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/

template <typename particle>
inline int verlet<particle>::get_neighbour_indices( particle *point , int *neigh_array_length , int **neigh_array ) {
  int Nneighs = 0 ;
  int dx = (cells->nmax.x > 1) , dy = (cells->nmax.y > 1) , dz = (cells->nmax.z > 1) ;
  struct posizione_3D<int> nearCell , point_cell = ( point->position / cells->cell_side ).convert2int_floor() ;
  struct clist_element *walker = NULL ;
  particle apparent_particle ;

  for( int i=point_cell.x-dx ; i<point_cell.x+dx+1 ; i++ ) {
    nearCell.x = ( i + cells->nmax.x ) % cells->nmax.x ;
    for( int j=point_cell.y-dy ; j<point_cell.y+dy+1 ; j++ ) {
      nearCell.y = ( j + cells->nmax.y ) % cells->nmax.y ;
      for( int k=point_cell.z-dz ; k<point_cell.z+dz+1 ; k++ ) {
	nearCell.z = ( k + cells->nmax.z ) % cells->nmax.z ;
	walker = cells->index[ nearCell.x ][ nearCell.y ][ nearCell.z ].list ;
	while( walker != NULL ) {
	  god->right_copy( element[walker->host].ptr , point , &apparent_particle ) ;
	  if( point->position.distance(apparent_particle.position) < r_verlet ) {
	    if( Nneighs >= *neigh_array_length ) {
	      (*neigh_array_length)++ ;
	      *neigh_array = (int *)realloc( *neigh_array , *neigh_array_length * sizeof(int*) ) ;
	    }
	    (*neigh_array)[Nneighs] = walker->host ;
	    Nneighs++ ;
	  }
	  walker = walker->next ;
	}
      }
    }
  }
  return Nneighs ;
}
/****************************************************************************************/

template <typename particle>
inline void verlet<particle>::initialize_cellist( void ) {
  cells = new cell_list( god->box_sides.position , r_verlet , particles_in_list ) ;
  if( check_memalloc( cells , "error in allocation of cell list in initialize_cellist()" ) ) exit( EXIT_FAILURE ) ;
  if( cells->nmax.x == 1 && cells->nmax.y == 1 && cells->nmax.z == 1 ) {
    delete cells ;
    cells = NULL ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void verlet<particle>::reinitialize_cellist( void ) {
  if( check_memalloc( cells , "\n*** ERROR: cells not allocated!" ) ) exit( EXIT_FAILURE ) ;
  cells->reinitialize_cellist_index( god->box_sides.position , r_verlet ) ;
  if( cells->nmax.x == 1 && cells->nmax.y == 1 && cells->nmax.z == 1 ) {
    delete cells ;
    cells = NULL ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void verlet<particle>::create_cellist( void ) {
  for( int i=0 ; i<particles_in_list ; i++ ) {    // calculating the coordinates of the cell where the i-th particle is in
    cells->particleScell[i] = ( element[i].ptr->position / cells->cell_side ).convert2int_floor() ;
    cells->particle_postIt[i].next = cells->index[ cells->particleScell[i].x ][ cells->particleScell[i].y ][ cells->particleScell[i].z ].list ;    // connecting i-th particle's list element to the index of her cell
    cells->index[ cells->particleScell[i].x ][ cells->particleScell[i].y ][ cells->particleScell[i].z ].list = ( cells->particle_postIt + i ) ;
  }
}
/****************************************************************************************/
