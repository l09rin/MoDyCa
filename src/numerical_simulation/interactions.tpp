#include "configuration.h"
#include "verlet.h"

template <typename particle>
interaction<particle>::~interaction() {
  delete verlet_list ;
  verlet_list = NULL ;
};

template <typename particle>
void interaction<particle>::link_list2particles( void ) {
  if( verlet_list ) {
    for( int i=0; i<verlet_list->particles_in_list; i++ ) verlet_list->element[i].ptr->list_ptr[idx] = verlet_list->element + i ;
  }
};
