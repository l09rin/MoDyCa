#ifndef VERLET_H
#define VERLET_H

template <typename particle> class configuration ;
template <typename particle> class molecule ;
class cell_list ;

template <typename particle>
struct neighbour_list {
  particle *ptr = NULL ;
  //    signed int ID_in_list ;
  struct neighbour_list<particle> *next ;
};

template <typename particle>
struct vlist_element {
  particle *ptr = NULL ;
  particle *cell_center = NULL ;
  struct neighbour_list<particle> *neighbours = NULL ;
};


template <typename particle>
class verlet {
public :
  struct vlist_element<particle> *element = NULL ;
  double r_verlet , delta_verlet ;
  int particles_in_list ;
  int arrays_length ;                  // length of the array element[]
  bool must_be_updated , disp_on ;
  bool DOUBLE_PAIRS ;     // This will tell if the list has been constructed in MC (=1) or MD (=0, default) mode
  char INITIALIZATION_MODE[10] ;

  configuration<particle> *god = NULL ;
  cell_list *cells = NULL ;

  molecule<particle> **mols = NULL ;          // this is useful in case of molecular potentials are employed (say rings...)
  int Nmols = 0 , mols_length = 0 ;

  verlet( void ) ;
  verlet( configuration<particle> *conf_ptr ) ;
  verlet( configuration<particle> *conf_ptr , int initial_length , double verlet_ray , double verlet_delta , const char *disp_mode = NULL ) ;
  verlet( const verlet &other_verlet ) ;
  ~verlet() ;

  void initialize_particles_list( const char *kind = NULL ) ;  // calculate the number of cells and initialize the list
  void initialize_molecules_list( const char *kind = NULL ) ;  // put to 0 the particles list and initialize the molecules list
  inline void verlet_creation_bycell( const char *simulation_kind = NULL ) ;   // creation of verlet list using the cell list method
  void bonding_verlet_creation( const char *simulation_kind = NULL );   // creation of verlet list using the cell list method (bonding and non bonding lists are created)
  inline int get_neighbour_indices( particle *point , int *neigh_array_length , int **neigh_array ) ;  // returns the number of neighboring particle inices to the point point, stored in the array neigh_array
  void copy_neighbour_list( struct vlist_element<particle> *gancio , struct neighbour_list<particle> *lista ) ; // return the pointer to a copy of lista
  inline void clear( void ) ;     // this member clears the verlet list

  // PRINT INFORMATION
  void print_cell_info( void ) ;

  inline void initialize_cellist( void ) ;      // try to allocate and initialize a cells list
  inline void reinitialize_cellist( void ) ;      // try to re-initialize a cells list
private :
  inline void add_connection_between( int i , int j ) ;     // Put particle j in the verlet list of particle i
  void free_neigh_list( struct neighbour_list<particle> *lista ) ;     // free a dynamically allocated list of concatenated elements of struct neighbour_list kind
  inline void create_cellist( void ) ;      // rebuild the cells list

};


#include "verlet.tpp"
#include "verlet_cells.tpp"

#endif
