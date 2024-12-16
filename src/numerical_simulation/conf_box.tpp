/* This library file contains the functions to deal with periodic boundary conditions in the simulation box */

/****************************************************************************************/
template <typename particle>
inline void configuration<particle>::right_copy( particle_3D *jptr , particle_3D *iptr , particle_3D *copy ) {
  double dx = jptr->position.x - iptr->position.x ;
  double dy = jptr->position.y - iptr->position.y ;
  double dz = jptr->position.z - iptr->position.z ;

  copy->position = jptr->position ;
  if( fabs(dx) > midside.position.x ) {   // Do i-particle interact with j-particle rather than one of its copies?  x axis...
    if( dx > 0 )  copy->position.x -= box_sides.position.x ;
    else  copy->position.x += box_sides.position.x ;
  }
  if( fabs(dy) > midside.position.y ) {   // y axis...
    if( dy > 0 )  copy->position.y -= box_sides.position.y ;
    else  copy->position.y += box_sides.position.y ;
  }
  if( fabs(dz) > midside.position.z ) {    // z axis...
    if( dz > 0 )  copy->position.z -= box_sides.position.z ;
    else  copy->position.z += box_sides.position.z ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::right_copy( particle_2D *jptr , particle_2D *iptr , particle_2D *copy ) {
  double dx = jptr->position.x - iptr->position.x ;
  double dy = jptr->position.y - iptr->position.y ;

  copy->position = jptr->position ;
  if( fabs(dx) > midside.position.x ) {   // Do i-particle interact with j-particle rather than one of its copies?  x axis...
    if( dx > 0 )  copy->position.x -= box_sides.position.x ;
    else  copy->position.x += box_sides.position.x ;
  }
  if( fabs(dy) > midside.position.y ) {   // y axis...
    if( dy > 0 )  copy->position.y -= box_sides.position.y ;
    else  copy->position.y += box_sides.position.y ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::pbc( particle_3D *part ) {
  if ( part->position.x < 0 ) {
    part->position.x += box_sides.position.x ;
    (part->periodic_box.x) -- ;
  } else if ( part->position.x >= box_sides.position.x  ) {
    part->position.x -= box_sides.position.x ;
    (part->periodic_box.x) ++ ;
  }
  if ( part->position.y < 0 ) {
    part->position.y += box_sides.position.y ;
    (part->periodic_box.y) -- ;
  } else if ( part->position.y >= box_sides.position.y  ) {
    part->position.y -= box_sides.position.y ;
    (part->periodic_box.y) ++ ;
  }
  if ( part->position.z < 0 ) {
    part->position.z += box_sides.position.z ;
    (part->periodic_box.z) -- ;
  } else if ( part->position.z >= box_sides.position.z  ) {
    part->position.z -= box_sides.position.z ;
    (part->periodic_box.z) ++ ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::pbc( particle_2D *part ) {
  if ( part->position.x < 0 ) {
    part->position.x += box_sides.position.x ;
    (part->periodic_box.x) -- ;
  } else if ( part->position.x >= box_sides.position.x  ) {
    part->position.x -= box_sides.position.x ;
    (part->periodic_box.x) ++ ;
  }
  if ( part->position.y < 0 ) {
    part->position.y += box_sides.position.y ;
    (part->periodic_box.y) -- ;
  } else if ( part->position.y >= box_sides.position.y  ) {
    part->position.y -= box_sides.position.y ;
    (part->periodic_box.y) ++ ;
  }
}
/****************************************************************************************/
