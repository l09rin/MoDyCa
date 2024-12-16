/* methods to print system and simulation information (useful for debugging purposes) */

template <typename particle>
void configuration<particle>::print_all_interactions( void ) {
  cout << "******************************************************************************\n";
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) cout << INTERACTIONS.type[k]->ostr().str() << "   " << INTERACTIONS.type[k]->compute_energy_contribution() << endl ;
  cout << "******************************************************************************\n";
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::print_all_information( void ) {
  neighbour_list<particle> *walker=NULL;
  cout << "******************************************************************************\n";
  if( strcmp( simulation_ensemble , "MD" ) == 0 ) {
    for( int i=0 ; i<numberOfPart ; i++ ) {
      cout << "particle n. " << i << " :" << endl ;
      cout << "           position:  " << *(particles[i]) ;
      cout << "           velocity:  " << *(velocities[i]) ;
      cout << "        total force:  " << *(forces_t1[i]) ;
      cout << "   prev. step force:  " << *(forces_t2[i]) ;
    }
  } else {
    for( int i=0 ; i<numberOfPart ; i++ ) {
      cout << "particle n. " << i << " :" << endl ;
      cout << "           position:  " << *(particles[i]) ;
    }
  }
  cout << endl ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    int particles_in_list = INTERACTIONS.type[k]->verlet_list->particles_in_list ;
    verlet<particle> *verlet_list = INTERACTIONS.type[k]->verlet_list ;
    vlist_element<particle> *element = INTERACTIONS.type[k]->verlet_list->element ;
    cout << "Verlet list n. " << k << "  " << verlet_list->r_verlet << "  " << verlet_list->delta_verlet << endl ;
    for( int l=0 ; l<CHECKS.bonds_lists_number ; l++ ) {
      if( verlet_list == CHECKS.bonds_verlet_link[l] ) cout << "Bonding verlet list" << endl ;
    }
    //    if( INTERACTIONS.type[k]->verlet_list == RUNTIME_CALCULATIONS.charge_vlist ) cout << "Charge verlet list" << endl ;
    for( int i=0 ; i<particles_in_list ; i++ ) {
      cout << "particle " << i << ": ";
      walker = element[i].neighbours ;
      while( walker != NULL ) {
	cout << " " << walker->ptr->id << "," ;
	walker = walker->next ;
      }
      cout << endl ;
    }
  }
  cout << "******************************************************************************\n";
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::print_displacements( void ) {
  cout << "******************************************************************************\n";
  verlet<particle> *verlet_list = NULL ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    cout << "  Verlet displacements, pot: " << k << endl ;
    verlet_list = INTERACTIONS.type[k]->verlet_list ;
    if( verlet_list->disp_on == 1 ) {
      cout << "    max displacement: " << verlet_list->delta_verlet * 0.5 << endl ;
      for( int i=0 ; i<verlet_list->particles_in_list ; i++ ) cout << i << "\t" << verlet_list->element[i].ptr->position - verlet_list->element[i].cell_center->position << endl ;
    }
    cout << endl << endl ;
  }
  cout << "******************************************************************************\n\n" ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::print_max_velocity( void ) {
  cout << "max velocity :  " << max_component_of_vector( particles , numberOfPart ) << "  (sqrt(EPSILON/M))" ;
}
/****************************************************************************************/

template <typename particle>
double configuration<particle>::max_component_of_vector( particle **vec , int length ) {
  double max = 0 , max_component = 0 ;
  for( int i=0 ; i<length ; i++ ) {
    max_component = vec[i]->norm() ;
    max = (max_component>max) ? max_component : max ;
  }
  return max ;
}
/****************************************************************************************/

template <typename particle>
double *configuration<particle>::min_part_distance( void ) {
  double min = box_sides.norm() , dist = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  int id1 = -1 , id2 = -1 ;

  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    int particles_in_list = INTERACTIONS.type[k]->verlet_list->particles_in_list ;
    vlist_element<particle> *element = INTERACTIONS.type[k]->verlet_list->element ;
    for( int i=0 ; i<particles_in_list ; i++ ) {
      walker = element[i].neighbours ;
      while( walker != NULL ) {
	right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
	dist = element[i].ptr->distance( apparent_near ) ;
	if ( dist < min ) {
	  min = dist ;
	  id1 = element[i].ptr->id ;
	  id2 = walker->ptr->id ;
	}
	walker = walker->next ;
      }
    }
  }
  delete apparent_near ;
  double *a = new double[3] ;
  a[0] = min ;
  a[1] = id1 ;
  a[2] = id2 ;
  return a ;
}
/****************************************************************************************/

template <typename particle>
double *configuration<particle>::max_part_distance( void ) {
  double max = 0 , dist = 0 ;
  neighbour_list<particle> *walker = NULL ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  int id1 = -1 , id2 = -1 ;

  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    int particles_in_list = INTERACTIONS.type[k]->verlet_list->particles_in_list ;
    vlist_element<particle> *element = INTERACTIONS.type[k]->verlet_list->element ;
    for( int i=0 ; i<particles_in_list ; i++ ) {
      walker = element[i].neighbours ;
      while( walker != NULL ) {
	right_copy( walker->ptr , element[i].ptr , apparent_near ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
	dist = element[i].ptr->distance( apparent_near ) ;
	if ( dist > max ) {
	  max = dist ;
	  id1 = element[i].ptr->id ;
	  id2 = walker->ptr->id ;
	}
	walker = walker->next ;
      }
    }
  }
  delete apparent_near ;
  double *a = new double[3] ;
  a[0] = max ;
  a[1] = id1 ;
  a[2] = id2 ;
  return a ;
}
/****************************************************************************************/

template <typename particle>
void verlet<particle>::print_cell_info( void ) {
  cout << endl << "nmax : " << cells->nmax << endl ;
  cout << "cells sides : " << cells->cell_side << endl ;
  cout << "index pointer : " << cells->index << endl ;
  cout << "near cells to explore : " << cells->N_of_NearCells2explore << endl ;
  for( int i=0 ; i<cells->N_of_NearCells2explore ; i++ ) cout << "   " << cells->NearCells2explore[i] << endl ;
  cout << endl ;
  for( int i=0 ; i<particles_in_list ; i++ ) {
    cout << "particle " << cells->particle_postIt[i].host << " ; next " << cells->particle_postIt[i].next << " ; position ( " << element[i].ptr->position.x << " , " << element[i].ptr->position.y << " , " << element[i].ptr->position.z << " ) ; cell ( " << cells->particleScell[i].x << " , " << cells->particleScell[i].y << " , " << cells->particleScell[i].z << " )" << endl ;
  }
  cout << endl << endl ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::interactions_array::print_on_file( double delta , int steps ) {
  // This function writes the pair potential and pair force values as a function of the distance on a file
  particle zero , target ;
  char filename[500] ;
  double dist = 0 ;
  for( int j=0 ; j<number ; j++ ) {
    ofstream pot_out , force_out ;
    sprintf( filename , "potential_%d.dat" , j ) ;
    pot_out.open( filename ) ;
    sprintf( filename , "force_%d.dat" , j ) ;
    force_out.open( filename ) ;
    for( int i=1 ; i<=steps ; i++ ) {
      target.position.x = i * delta ;
      dist = zero.position.distance( target.position ) ;
      if( i*delta < type[j]->cutoff ) {
	pot_out << i*delta << " " << type[j]->return_pair_potential( dist ) << endl ;
	// This should be modified in case of external fields, for which pair force gives 0 everywhere
	force_out << i*delta << " " << type[j]->return_pair_force( &zero , &target ).position.x << endl ;
      } else {
	pot_out << i*delta << " " << "0" << endl ;
	force_out << i*delta << " " << zero << endl ;
      }
    }
  }
}
/****************************************************************************************/
