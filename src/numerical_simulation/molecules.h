#ifndef MOLECULES_H
#define MOLECULES_H

/* This library defines the class molecule, to deal with systems with intra-molecular potentials */

template <typename particle>
class molecule {
 public:
  particle cm , vcm , fcm , omega ;
  double gyration = 0 , gyration_square = 0 , monomers_diameter = 0 ;
  particle **atom = NULL ;
  int Natoms = 0 ;

  // Variables for the hertzian potential calculation
  double hertzian_radius = 0 ; // distance and force units respective to the molecule's parameters

  inline void compute_CoM( particle box_sides ) ;
  inline void compute_gyration_square( particle box_sides ) ;
  inline void compute_gyration( void ) ;

  molecule() ;
  molecule( const molecule &other ) ;
  ~molecule() ;
};


#include "molecules.tpp"

#endif
