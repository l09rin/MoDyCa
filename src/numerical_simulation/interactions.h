#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include <iostream>
#include <string>
#include "../record.h"
#include "../system_comm.h"

template <typename particle> class verlet ;
template <typename particle> class configuration ;

template <typename particle>
class interaction {
 public:
  int idx = -1 ;
  char name[15] = "" ;
  double cutoff = 0 ;
  double SHIFT = 0 ;
  verlet<particle> *verlet_list = NULL ;

  // in generate_by_record function the verlet list is also initialized, using conf_ptr and parts_number, considering simulation_kind (MD|MC)
  virtual void generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) = 0 ;
  virtual double return_pair_potential( double distance ) = 0 ;
  virtual double return_pair_potential( particle *part1 , particle *part2 ) = 0 ; // This version accounts for the cutoff, same molecule's belonging and PBCs
  virtual particle return_pair_force( particle *i_part , particle *j_part ) = 0 ;
 public:
  time_dep_variable<double> partial_energy , partial_virial ;
  virtual double compute_energy_contribution( void ) = 0 ;
  virtual double compute_energy_contribution( particle *part ) = 0 ;
  virtual double delta_energy_constrained_CoM( particle *disp ) = 0 ;  // non trivial only for external fields, used in MC_NVT_CMfixed ensemble
  virtual double lambda_derivative_perpart( void ) = 0 ;  // non trivial only for external fields, used in MC_NVT_CMfixed ensemble, lambda derivative of the potential
  virtual void compute_forces_contribution( void ) = 0 ;
  virtual double compute_virial_contribution( void ) = 0 ;
  virtual void compute_per_part_energy_contribution( double *envec ) = 0 ;
  virtual void compute_per_part_virial_contribution( double *virvec ) = 0 ;
  virtual string ostr( void ) = 0 ;
  virtual double *dump_parameters( void ) = 0 ;
  virtual void build_bonds_network( double bonds_cutoff = 0.0 ) = 0 ;
  void link_list2particles( void ) ;
  virtual ~interaction() ;
};


#include "interactions.tpp"

#include "Interactions/lennard_jones.tpp"
#include "Interactions/lennard_jones_poly.tpp"
#include "Interactions/soft.tpp"
#include "Interactions/fene.tpp"
#include "Interactions/fene_poly.tpp"
#include "Interactions/solvophobic.tpp"
#include "Interactions/wca_alpha.tpp"
#include "Interactions/debye_huckel.tpp"
#include "Interactions/debye_huckel_poly.tpp"
#include "Interactions/hertzian_molecule.tpp"
#include "Interactions/conservative_hertzian.tpp"
#include "Interactions/kern_frenkel.tpp"
#include "Interactions/hard_sphere.tpp"
#include "Interactions/hs_square_well.tpp"
#include "Interactions/harmonic_spring_field.tpp"
#include "Interactions/harmonic_spring_field_wignerseitz.tpp"
#include "Interactions/angular_harmonic_spring_field.tpp"
#include "Interactions/angular_cosine_spring_field.tpp"
#include "Interactions/test_potential.tpp"

#endif
