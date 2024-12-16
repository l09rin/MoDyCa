/*****************************************************************************************/

template <typename particle>
class conservative_hertzian : public interaction<particle> {  // MD  // energy in units of EPSILON
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::idx ;

 public:
  double U , S_H ;  // amplitude of the hertzian and scaling factor to define the molecule diameter with respect to the monomer's one (sigma_H = S_H * sigma_monomer)
  conservative_hertzian( double amplitude , double SH ) ;

 private:
  particle force ;
  double multiplier = 0 , sqrt_factor = 0 ;
  inline double pair_potential( particle *i_part , particle *j_part , double distance , double mol_radius ) ;
  inline particle pair_force( particle *i_part , particle *j_part , double distance , double mol_radius ) ;

 public:
  double return_pair_potential( double distance ) ;
  double return_pair_potential( double distance , double mol_radius ) ;
  double return_pair_potential( particle *part1 , particle *part2 ) ;
  particle return_pair_force( particle *i_part , particle *j_part ) ;
  particle return_pair_force( particle *i_part , particle *j_part , double mol_radius ) ;
  void generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) ;

  double compute_energy_contribution( void ) ;
  inline double compute_energy_contribution( particle *part ) ;
  double delta_energy_constrained_CoM( particle *disp ) ;
  double lambda_derivative_perpart( void ) ;
  void compute_forces_contribution( void ) ;
  double compute_virial_contribution( void ) ;
  void compute_per_part_energy_contribution( double *envec ) ;
  void compute_per_part_virial_contribution( double *virvec ) ;

  string ostr( void ) ;
  void build_bonds_network( double bonds_cutoff = 0.0 ) ;
  double *dump_parameters( void ) ;
  conservative_hertzian *copy( void ) ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

template <typename particle>
conservative_hertzian<particle>::conservative_hertzian( double amplitude , double SH ) {
  U = amplitude ;
  S_H = SH ;
};
/********************************************************************************************/

template <typename particle>
void conservative_hertzian<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
	              // CUTOFF
  cutoff = 0 ;
	              // SHIFT
  SHIFT = 0 ;
	              // generation of the verlet list
  verlet_list = new verlet<particle>( conf_ptr ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 0 ;
	              // initialization of the verlet list
  if( strcmp( file_rec->word[3] , "all" ) == 0 ) {
    verlet_list->initialize_molecules_list( NULL ) ;
  } else if( strcmp( file_rec->word[3] , "file" ) == 0 ) {
    verlet_list->initialize_molecules_list( file_rec->word[4] ) ;
  } else {
    cout << "\n Molecules list format not recognized !!" << endl ;
    exit( EXIT_FAILURE ) ;
  }       // Here I tell to the integrator that I need to update the molecules positions when updating particles ones
  verlet_list->god->compute_MOLcm = 1 ;

  for( int i=0 ; i<verlet_list->Nmols ; i++ ) verlet_list->mols[i]->hertzian_radius = verlet_list->mols[i]->monomers_diameter * S_H * 0.5 ;
};
/********************************************************************************************/

template <typename particle>
inline double conservative_hertzian<particle>::pair_potential( particle *i_part , particle *j_part , double distance , double mol_radius ) {
    // energy in units of EPSILON
    //        U * (1-dist)**(5/2)
    /***********   HERTZIAN POTENTIAL          **********/
  multiplier = distance / mol_radius ;
  sqrt_factor = sqrt( 1.0 - multiplier ) ;
  return U * ( 1.0 - multiplier ) * ( 1.0 - multiplier ) * sqrt_factor ;
};
/********************************************************************************************/

template <typename particle>
inline particle conservative_hertzian<particle>::pair_force( particle *i_part , particle *j_part , double distance , double mol_radius ) {
    // returns the force j_part applies on i_part, in units of EPSILON
    //        5/2 * U / sigma_H * (1-dist)**(3/2)
    /***********   HERTZIAN FORCE          **********/
  multiplier = distance / mol_radius ;
  sqrt_factor = sqrt( 1.0 - multiplier ) ;
  multiplier = 2.5 * U * ( 1.0 - multiplier ) * sqrt_factor / mol_radius / distance ;
  force.position = i_part->position - j_part->position ;
  force.position *= multiplier ;

  return force ;
};
/********************************************************************************************/

template <typename particle>
double conservative_hertzian<particle>::return_pair_potential( double distance ) {
  cout << "To return the value of the hertzian potential the molecule's diameter is needed !!" << endl ;
  return nanl("") ;
};
/********************************************************************************************/

template <typename particle>
double conservative_hertzian<particle>::return_pair_potential( double distance , double mol_radius ) {
  return pair_potential( NULL , NULL , distance , mol_radius ) + SHIFT ;
};
/********************************************************************************************/

template <typename particle>
double conservative_hertzian<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  double energy = 0 ;
  double dist = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  particle *mol_cm = NULL , **atom = NULL ;
  int molecules_in_list = verlet_list->Nmols , atoms_in_mol = 0 ;
  molecule<particle> **mols = verlet_list->mols ;
  double mol_radius = 0 ;
  bool keep_running = 1 ;

  // Here should put info on mol idx in particles, to make it faster...
  for( int i=0 ; i<molecules_in_list && keep_running ; i++ ) {
    atoms_in_mol = mols[i]->Natoms ;
    atom = mols[i]->atom ;
    for( int j=0 ; j<atoms_in_mol ; j++ ) {
      if( atom[j]->id == part1->id ) keep_running = 0 ;
    }
    if( ! keep_running ) {
      for( int j=0 ; j<atoms_in_mol ; j++ ) {
	if( atom[j]->id == part2->id ) {
	  mol_cm = &(mols[i]->cm) ;
	  mol_radius = mols[i]->hertzian_radius ;
	  for( int k=0 ; k<atoms_in_mol ; k++ ) {
	    verlet_list->god->right_copy( atom[k] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
	    if( ( dist = mol_cm->position.distance( apparent_near->position ) ) < mol_radius ) {
	      energy += ( pair_potential( apparent_near , mol_cm , dist , mol_radius ) + SHIFT ) ;
	    }
	  }
	}
      }
    }
  }

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
particle conservative_hertzian<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  cout << "To return the value of the hertzian force the molecule's diameter is needed !!" << endl ;
  particle force ;
  force.position = force.position * nanl("") ;
  return force ;
};
/********************************************************************************************/

template <typename particle>
particle conservative_hertzian<particle>::return_pair_force( particle *i_part , particle *j_part , double mol_radius ) {
  double dist = i_part->position.distance( j_part->position ) ;
  return pair_force( i_part , j_part , dist , mol_radius ) ;
};
/********************************************************************************************/

template <typename particle>
double conservative_hertzian<particle>::compute_energy_contribution( void ) {
  double energy = 0 ;
  double dist = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  particle *mol_cm = NULL , **atom = NULL ;
  int molecules_in_list = verlet_list->Nmols , atoms_in_mol = 0 ;
  molecule<particle> **mols = verlet_list->mols ;
  double mol_radius = 0 ;

  for( int i=0 ; i<molecules_in_list ; i++ ) {
    atoms_in_mol = mols[i]->Natoms ;
    atom = mols[i]->atom ;
    mol_cm = &(mols[i]->cm) ;
    mol_radius = mols[i]->hertzian_radius ;
    for( int j=0 ; j<atoms_in_mol ; j++ ) {
      verlet_list->god->right_copy( atom[j] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
      if( ( dist = mol_cm->position.distance( apparent_near->position ) ) < mol_radius ) {
	energy += ( pair_potential( apparent_near , mol_cm , dist , mol_radius ) + SHIFT ) ;
      }
    }
  }

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
inline double conservative_hertzian<particle>::compute_energy_contribution( particle *part ) {
  // This function computes the energy of ONE particle only
  double energy = 0 ;
  double dist = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  particle *mol_cm = NULL , **atom = NULL ;
  int molecules_in_list = verlet_list->Nmols , atoms_in_mol = 0 ;
  molecule<particle> **mols = verlet_list->mols ;
  double mol_radius = 0 ;
  bool keep_running = 1 ;

  for( int i=0 ; i<molecules_in_list && keep_running ; i++ ) {
    atoms_in_mol = mols[i]->Natoms ;
    atom = mols[i]->atom ;
    for( int j=0 ; j<atoms_in_mol ; j++ ) {
      if( atom[j] == part ) keep_running = 0 ;
    }
    if( ! keep_running ) {
      mol_cm = &(mols[i]->cm) ;
      mol_radius = mols[i]->hertzian_radius ;
      for( int j=0 ; j<atoms_in_mol ; j++ ) {
	verlet_list->god->right_copy( atom[j] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
	if( ( dist = mol_cm->position.distance( apparent_near->position ) ) < mol_radius ) {
	  energy += ( pair_potential( apparent_near , mol_cm , dist , mol_radius ) + SHIFT ) ;
	}
      }
    }
  }

  delete apparent_near ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
double conservative_hertzian<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double conservative_hertzian<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from ideal lattice sites with the external harmonic spring potential
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
void conservative_hertzian<particle>::compute_forces_contribution(void) {
  double dist = 0 ;
  particle *apparent_near = NULL , deltaF , **forces_t2 = verlet_list->god->forces_t2 ;
  apparent_near = new particle ;
  particle *mol_cm = NULL , **atom = NULL ;
  int molecules_in_list = verlet_list->Nmols , atoms_in_mol = 0 , j = 0 , k = 0 ;
  molecule<particle> **mols = verlet_list->mols ;
  double mol_radius = 0 ;

  for( int i=0 ; i<molecules_in_list ; i++ ) {
    atoms_in_mol = mols[i]->Natoms ;
    atom = mols[i]->atom ;
    mol_cm = &(mols[i]->cm) ;
    mol_radius = mols[i]->hertzian_radius ;
    for( j=0 ; j<atoms_in_mol ; j++ ) {
      verlet_list->god->right_copy( atom[j] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
      if( ( dist = mol_cm->position.distance( apparent_near->position ) ) < mol_radius ) {
	deltaF = pair_force( apparent_near , mol_cm , dist , mol_radius ) ;
	forces_t2[ atom[j]->id ]->position += deltaF.position ;         // force applied on the atom j
	deltaF.position /= atoms_in_mol ;
	for( k=0 ; k<atoms_in_mol ; k++ ) forces_t2[ atom[k]->id ]->position -= deltaF.position ;         // force redistributed over all the molecule's atoms
      }
    }
  }
  delete apparent_near ;
};
/********************************************************************************************/

template <typename particle>
double conservative_hertzian<particle>::compute_virial_contribution( void ) {
  double dist = 0 ;
  particle *apparent_near = NULL , deltaF ;
  apparent_near = new particle ;
  double virial = 0 ;
  particle *mol_cm = NULL , **atom = NULL ;
  int molecules_in_list = verlet_list->Nmols , atoms_in_mol = 0 , j = 0 , k = 0 ;
  molecule<particle> **mols = verlet_list->mols ;
  double mol_radius = 0 ;

  for( int i=0 ; i<molecules_in_list ; i++ ) {
    atoms_in_mol = mols[i]->Natoms ;
    atom = mols[i]->atom ;
    mol_cm = &(mols[i]->cm) ;
    mol_radius = mols[i]->hertzian_radius ;
    for( j=0 ; j<atoms_in_mol ; j++ ) {
      verlet_list->god->right_copy( atom[j] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
      if( ( dist = mol_cm->position.distance( apparent_near->position ) ) < mol_radius ) {
	deltaF = pair_force( apparent_near , mol_cm , dist , mol_radius ) ;
	virial += ( deltaF.position * ( apparent_near->position - mol_cm->position ) ) ;
	deltaF.position /= atoms_in_mol ;
	for( k=0 ; k<atoms_in_mol ; k++ ) {
	  verlet_list->god->right_copy( atom[k] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
	  virial -= ( deltaF.position * ( apparent_near->position - mol_cm->position ) ) ;         // force redistributed over all the molecule's atoms
	}
      }
    }
  }
  delete apparent_near ;
  return virial ;
};
/********************************************************************************************/

template <typename particle>
string conservative_hertzian<particle>::ostr( void ) {
  ostringstream st ;
  st <<  "CONSERVATIVE HERTZIAN RINGS parameters: ( U=" << U << " ; sigma_H/sigma_m=" << S_H << " )" << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
conservative_hertzian<particle> *conservative_hertzian<particle>::copy( void ) {
  conservative_hertzian *copy = new conservative_hertzian() ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void conservative_hertzian<particle>::compute_per_part_energy_contribution( double *envec ) {
  double energy = 0 ;
  double dist = 0 ;
  particle *apparent_near = NULL ;
  apparent_near = new particle ;
  particle *mol_cm = NULL , **atom = NULL ;
  int molecules_in_list = verlet_list->Nmols , atoms_in_mol = 0 ;
  molecule<particle> **mols = verlet_list->mols ;
  double mol_radius = 0 ;

  for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) envec[i] = 0 ;
  for( int i=0 ; i<molecules_in_list ; i++ ) {
    atoms_in_mol = mols[i]->Natoms ;
    atom = mols[i]->atom ;
    mol_cm = &(mols[i]->cm) ;
    mol_radius = mols[i]->hertzian_radius ;
    energy = 0 ;
    for( int j=0 ; j<atoms_in_mol ; j++ ) {
      verlet_list->god->right_copy( atom[j] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
      if( ( dist = mol_cm->position.distance( apparent_near->position ) ) < mol_radius ) {
	energy += ( pair_potential( apparent_near , mol_cm , dist , mol_radius ) + SHIFT ) ;
      }
    }
    energy /= atoms_in_mol ;
    for( int j=0 ; j<atoms_in_mol ; j++ ) envec[ atom[j]->id ] += energy ;
  }

  delete apparent_near ;
};
/********************************************************************************************/

template <typename particle>
void conservative_hertzian<particle>::compute_per_part_virial_contribution( double *virvec ) {
  double dist = 0 ;
  particle *apparent_near = NULL , deltaF ;
  apparent_near = new particle ;
  particle *mol_cm = NULL , **atom = NULL ;
  int molecules_in_list = verlet_list->Nmols , atoms_in_mol = 0 , j = 0 , k = 0 ;
  molecule<particle> **mols = verlet_list->mols ;
  double mol_radius = 0 ;

  for( int i=0 ; i<verlet_list->god->numberOfPart ; i++ ) virvec[i] = 0 ;
  for( int i=0 ; i<molecules_in_list ; i++ ) {
    atoms_in_mol = mols[i]->Natoms ;
    atom = mols[i]->atom ;
    mol_cm = &(mols[i]->cm) ;
    mol_radius = mols[i]->hertzian_radius ;
    for( j=0 ; j<atoms_in_mol ; j++ ) {
      verlet_list->god->right_copy( atom[j] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
      if( ( dist = mol_cm->position.distance( apparent_near->position ) ) < mol_radius ) {
	deltaF = pair_force( apparent_near , mol_cm , dist , mol_radius ) ;
	virvec[ atom[j]->id ] += ( deltaF.position * ( apparent_near->position - mol_cm->position ) ) ;
	deltaF.position /= atoms_in_mol ;
	for( k=0 ; k<atoms_in_mol ; k++ ) {
	  verlet_list->god->right_copy( atom[k] , mol_cm , apparent_near ) ;      // I take the coordinates of the nearest copy to the mol center of mass.
	  virvec[ atom[k]->id ] -= ( deltaF.position * ( apparent_near->position - mol_cm->position ) ) ;         // force redistributed over all the molecule's atoms
	}
      }
    }
  }
  delete apparent_near ;
};
/********************************************************************************************/

template <typename particle>
void conservative_hertzian<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * conservative_hertzian<particle>::dump_parameters( void ) {
  double *params = new double [2] ;
  params[0] = U ;
  params[1] = S_H ;
  return params ;
};
/********************************************************************************************/

/********************************************************************************************/
