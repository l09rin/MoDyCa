#ifndef CELLS_H
#define CELLS_H
/* Here the structs to deal with cell lists are defined */
/*    MAYBE MODIFY TO MAKE IT MORE EFFICIENT, SUCH AS SINGLE ARRAY FOR INDEX . . . */
#include "vector3D.h"

struct clist_element{
  struct clist_element *next;
  int host;                            // index of the particle
};

struct cell_array_element{     // array of the cells
  struct clist_element *list ;
  int cell_is_explored = 0 ;
};

class cell_list {
 public :
  int particles_in_list = 0 ;
  struct cell_array_element ***index ;      // tridimensional array of lists of particles contained in each cell
  struct posizione_3D<int> *particleScell ;                 // array containing coordinates of the cell wherein each particle is located
  struct clist_element *particle_postIt ;        // array of the list elements associated to each particle
  struct posizione_3D<int> nmax ;                            // number of pieces in which is divided each side of the box
  struct posizione_3D<int> *NearCells2explore ; // array of triplets containing the relative position of the cells to explore during the verlet list construction
  int N_of_NearCells2explore ;
  posizione_3D<double> cell_side ;                    // length of side of each cell

  template <typename T>
  void reinitialize_cellist_index( T &box_sides_vec , double r_verlet ) ;  // calculate the number of cells and initialize the list
  void set_nearCells2explore( void ) ;  // set the array of neighboring cells to explore for each one
  void free_index( void ) ;       // free the memory of the entire cell list structure
  inline void clear( void ) ;      // assign NULL to all elements of cell-index and clist-next

  cell_list() ;
  template <typename T>
  cell_list( T &box_sides_vec , double r_verlet , int parts_in_list ) ;
  cell_list( const cell_list &other ) ;
  ~cell_list() ;
};


#include "cells.cpp"

#endif
