#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <climits>
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "../matrix.h"
#include "../record.h"
#include "particles2D.h"
#include "particles3D.h"
#include "patchy2D.h"
#include "molecules.h"
#include "NoseHoover_thermostat.h"
#include "cells.h"
#include "verlet.h"

#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

template <typename particle> class interaction ;
template <typename particle> class mean_squared_displacement ;

// DEFINITION OF THE CLASS CONFIGURATION, CONTAINING ALL THE INFORMATION ABOUT THE SYSTEM AND THE SIMULATION METHODS


template <typename particle>
class configuration {

 public :

  int seed ;                    // It can be used to initialize srand48()

  /*******************************************   SIMULATION BOX AND ENSEMBLES CONSTANTS   *******************************************/

  int DIMENSION ;
  double k_T ;
  double energyOfSystem ;               // total value, sum over all particles
  double box_side ;
  particle box_sides ; // If the box is not cubic    // AGGIORNARE LE FUNZIONI IN CUI RIENTRA
  particle midside ;
  inline void pbc( particle_3D *part ) ;   // imposes the periodic boundary conditions
  inline void pbc( particle_2D *part ) ;   // imposes the periodic boundary conditions
  char simulation_ensemble[15] ;
  char simulation_mode[3] ;  // It can be MC or MD and determines how the verlet lists are constructed
         // MC-ONLY VARIABLES
  double NPT_pressure ;  // NPT ensemble
  double Lambda3 ;       // NPT, GC ensembles
  double fugacity ;       // NPT, GC ensembles, chemical potential
  double MC_maxdisplacement ;     // in MonteCarlo metropolis position's changes are drawn randomly in a cube of side MC_maxdisplacement
  // in MonteCarlo metropolis orientation's changes are drawn randomly within an angle opening of MC_maxrotation, rotationjump performs angle jumps
  //       according to the symmetry of the molecules to rightly sample angles in low temp / highly constrained crystals. . .
  double MC_maxrotation , MC_rotationjump ;
  double MC_maxboxdisp ;         // in MonteCarlo metropolis NPT each box sides coordinate change is drawn randomly within an interval [-MC_maxboxdisp;MC_maxboxdisp]
  struct MC_PARAMS {
    double MC_maxdisplacement ;
    double MC_maxrotation ;
    double MC_maxboxdisp ;

    MC_PARAMS( double val_disp = 1.0 , double val_rot = M_PI , double val_boxdisp = 0.1 ) {
      MC_maxdisplacement = val_disp , MC_maxrotation = val_rot , MC_maxboxdisp = val_boxdisp ;
    };
  } MC_PARAMS_MAX , MC_PARAMS_MIN ;
         // MD-ONLY VARIABLES
  double time_step , global_time ;           // global_time is the global time of the simulation, time_step is the integration step for MD
  double particle_mass ;

  /*******************************************   THERMODYNAMIC QUANTITIES (COMPUTATION NEEDED)   *******************************************/

  time_dep_variable<double> virial_press ;
  time_dep_variable<double> POTENTIAL_ENERGY , KINETIC_ENERGY ;               // total values, sum over all particles


  /*******************************************   MEAN SQUARED DISPLACEMENT (msd.h)   *******************************************/
  mean_squared_displacement<particle> **MSD = NULL ;
  int N_MSDs = 0 ;
  void add_initialize_MSD( const char *mode , const char *steps_path = "msd_saving_steps.dat" ) ;

                                               // READ CONFIGURATIONS BY FILES ( read_data.cpp / dump_data.cpp )
  void read_configuration( ifstream &data_file , char *format ) ;
  void read_xyz( ifstream &data_file ) ;
  void read_lammps_dump( ifstream &data_file ) ;
  void read_nico_init( ifstream &data_file , const char *format ) ;
  void read_sph( ifstream &data_file ) ;
  void read_ptc( ifstream &data_file ) ;
  void dump_configuration( ofstream &out_file , char *format , int Ndigits = 8 ) ;
  void dump_molecules( ofstream *out_file , char *format , int yes_mol , int no_type ) ;
  void dump_patch( ofstream &out_file , int Ndigits = 8 ) ;
  void dump_sph( ofstream &out_file , int Ndigits = 8 ) ;
  void dump_xyz( ofstream &out_file , int Ndigits = 8 ) ;
  void dump_lmp( ofstream &out_file , int Ndigits = 8 ) ;

  /*******************************************   INTERACTIONS ( folder Interactions )   *******************************************/

  class interactions_array {
  public:
    int number ;
    double SIGMA, EPSILON ;
    interaction<particle> **type ;

    void print_on_file( double delta , int steps ) ;
    //friend ostream & operator<< <particle>( ostream &out , const interactions_array &int_array ) ;
    friend ostream & operator<<( ostream &out , const interactions_array &int_array ) {
      for( int i=0 ; i<int_array.number ; i++ ) {
	out << int_array.type[i]->ostr() ;
      }
      return out ;
    };
    void operator=( const interactions_array &other ) ;
    interactions_array() ;
    ~interactions_array() ;
  } INTERACTIONS , SECONDARY_INTERACTIONS ;

  /*******************************************   CHECKS   *******************************************/

  struct check {
    int bonds_lists_number = 0 ;
    verlet<particle> **bonds_verlet_link = NULL ;
  } CHECKS ;
  int lists_rebuilds = 0 ;


  /*******************************************   MONTE CARLO STATISTICS   *******************************/

  struct MC_stats {  // keeps track of the moves attempted/performed to change the system's :
    long int NoP_attempted = 0 , NoP_accepted = 0 ;  // number of particles
    long int pos_attempted = 0 , pos_accepted = 0 ;  // particle's position
    long int vol_attempted = 0 , vol_accepted = 0 ;  // box sides
    long int rot_attempted = 0 , rot_accepted = 0 ;  // particle's orientation
    long int chr_swap_attempted = 0 , chr_swap_accepted = 0 ;  // particle's charge

    void clear( void ) {
      NoP_attempted = 0 ; NoP_accepted = 0 ; pos_attempted = 0 ; pos_accepted = 0 ; vol_attempted = 0 ; vol_accepted = 0 ; rot_attempted = 0 ; rot_accepted = 0 ; chr_swap_attempted = 0 ; chr_swap_accepted = 0 ;
    };
  } MonteCarlo_stats ;

  /*******************************************   SYSTEM VECTORS   *******************************************/

  int numberOfPart , monomers_number , charges_number , particles_length ;   // number of all particles, monomers fene-bonded, charged monomers in the system
  particle **particles = NULL ;                       // It is a vector containing pointers to particle positions
  int numberOfMol , molecules_length , compute_MOLcm = 0 , compute_MOLgyr = 0 ;
  molecule<particle> **molecules = NULL ;
  matrix inertia_matrix ;
  particle **velocities = NULL ;                      // It is a vector containing pointers to particle velocities
  particle **forces_t1 = NULL, **forces_t2 = NULL ;     // It is a vector containing total forces acting on each particle in the system
  particle **support_pp2particle = NULL ;
      // intended as calculated with unwrapped coordinates, used in ensemble MC_NVT_CMfixed, for the moment
      // CoM_displacement is the displacement of the centre of mass since the beginning of the simulation
  particle centre_of_mass , CoM_displacement ;
  double deltaEnergyOfCoM = 0.0 ;   // energy difference among the system with and without a constrained centre of mass in MC NVT ensemble
          // INPUT - OUTPUT
  int INITIAL_CONF_BY_FILE ;             // if I want to start from a configuration of stored positions
  int INITIAL_VELOCITIES_BY_FILE ;             // if I want to start from a configuration of stored velocities
  char configuration_file_name[300], velocities_file_name[300] ;

  struct save_options {      // I use this variables in the main() to set if storing this data in files produced
    int PARTICLES , PARTICLES_IMG , VELOCITIES , FORCES , FORCES_SPLIT , PER_PART_ENERGY , PER_PART_VIRIAL , MOL_IDX , IDX , TYPE , CHARGE ;
    int MOLECULES , MOLECULES_IMG , MOL_VELOCITIES , MOL_FORCES , MOL_FORCES_SPLIT ;
    int POTENTIAL_ENERGY , PARTIAL_POTENTIAL_ENERGY , KINETIC_ENERGY , SECONDARY_POTENTIAL , MSD ;
    int VOLUME , VIRIAL , TEMPERATURE , GYR , RUNTIME_CALCULATIONS ;      // Data are saved if correspondent values are set to 1, 0 by default
    particle ***force_contributions = NULL ;   // This is needed in case of forces_splitted are used, to split the several contribution of the forces
    double **per_part_energy = NULL , **per_part_virial = NULL ;  // This is needed if you want to dump the values of energy and/or virial for each particle
    int global_interval ;               // number of steps between successive savings of global quantities
    int configuration_interval ;        // number of steps between successive savings of positions and velocities
    int backup_interval ;               // number of steps between storage of positions and velocities for a restart of the simulation
    int runtime_interval ;              // number of steps between successive savings of runtime calculations
    char *configuration_format ;
    void copy_options( const save_options &other ) ;
    void clear( void ) ;
    save_options( void ) ;
  } SAVE, EQ_SAVE_OPT, PROD_SAVE_OPT ;

  configuration( void ) ;
  configuration( const configuration &other_config ) ;   // copy constructor
  ~configuration() ;


  /*******************************************   BLOCK 1 - DEFINITIONS   *******************************************/

  void initialize_configuration( const char *input_file_name = NULL ) ;
  void interpret_input_script( ifstream &config_file ) ;


  /*******************************************   INITIAL POSITIONS SETTINGS ( initial_configuration.cpp )   *******************************************/

  void generate_ControlledEnergy_randomConfig( void ) ;     // this function generate a random configuration through controlled energy insertion for MD simulations
  void generate_HS_randomConfig( double radius = 0.5 ) ;  // this function generate a random configuration as if the particles were HS of determined radius
  void generate_ideal_lattice_3D( const char *ltype ) ;
  void generate_ideal_lattice( const char *ltype ) ;
  void generate_config_by_file( const char *file_name , const char *format = NULL ) ;
  void generate_molecules_by_file( const char *file_name ) ;


  /*******************************************   INITIAL VELOCITY SETTINGS ( velocities.h )   *******************************************/
  void generate_MBvelocities( double mean_square_velocity ) ; // this function generates a vector of velocity according to a
                                                            // Maxwell-Boltzmann distribution with mean value = compute_mean_velocity and mean square value = mean_square_velocity
  void generate_velocities_by_file( void ) ;
  inline void CoM_velocity_reset( void ) ;
  void total_angularMomentum_reset( const char protocol[] = "" ) ;  // It cancels the total angular momentum in the CM (protocol=RIGID|VTINVERSION|NULL)
                        // VTINVERSION = inverting the transverse component of velocities contributing to L ; RIGID = as if it were a rigid body : NULL = same as RIGID
  void set_kinetic_energy( double kin_en ) ;   // fix to kin_en the overall value of kinetic energy of system


  /*******************************************   INITIAL CHARGE ASSIGNMENT ( charges.h )   *******************************************/

              // generates the distribution of charges on the microgel network
  virtual void generate_charge_distribution( char *distribution, double charges_fraction = 0 , double charges_distribution_sigma = 0 ) ;


  /*******************************************   BLOCK 2 - COMPUTATIONS   *******************************************/

  double compute_potential_energy( void ) ;            // calculates the total energy of the configuration
  double compute_potential_energy( particle *part ) ;      // calculates the energy of a single particle (all the amount of pair and molecular contributions are assigned)
  double compute_potential_energy( particle *part1 , particle *part2 ) ;   // calculates the energy of a single pair of particles (molecular contributions are assigned only if the pair belongs to the same molecule)
  double delta_energy_constrained_CoM( particle *displacement ) ;            // calculates the energy difference among the constrained and unconstrained CoM system after the displacement move of a particle in the MC NVT ensemble with fixed CoM integration scheme
  double compute_lambda_derivative_perpart( void ) ;            // calculates the lambda derivative of the potentials in the spring lattice fields with centre of mass constraint
  void calculate_forces_t2() ;              // calculates the vector of forces acting on particles and place it in forces_t2
  double compute_virial_pressure( void ) ;
  double compute_kinetic_energy() ;            // calculates the kinetic energy of the configuration
  double compute_packing_fraction() ;
  double calculate_secondary_potential( void ) ;
  particle compute_centre_of_mass( void ) ;
  particle compute_mean_position( void ) ;                    // center of mass in case of identical particles
  particle compute_mean_velocity( void ) ;                    // velocity of center of mass in case of identical particles
  posizione_3D<double> compute_CoM_angular_momentum( void ) ;                    // in case of identical particles of mass 1
  double compute_inertia_momentum( posizione_3D<double> axis ) ;                    // in case of identical particles of mass 1
  void compute_inertia_momentum( void ) ;
  double compute_max_bond_distance( void ) ;         // calculates the maximume stretching of bonds in the configuration

  inline virtual double calculate_molecules_gyration( void ) ;
  inline virtual void calculate_molecules_CoM( void ) ;
  virtual void initializeVlist_runtime_calculations( verlet<particle> *charge_vlist = NULL ) ;  // attempt to initialize partially the runtime calc structures


  /*******************************************   SIMULATION ENSEMBLES ( Integrators/ )   *******************************************/

  int EQUILIBRATION_STEPS ;
  int PRODUCTION_STEPS ;
  void MonteCarlo( const char *questions ) ;                      // MC_SWAP_CHARGE, NVT , NPT ensembles
  void MD( const char *questions ) ;    // Molecular Dynamics simulation of NVE, NVT_NHC, NVT_NHE, BROWNIAN ensembles. IF Questions == NO_STOP equilibration and production steps cannot be modified

 protected:

  /*******************************************   OUTPUT   *******************************************/

  char _DIRECTORY_NAME[300] ;
  ofstream _thermoinfo_file ;
  ofstream _particles_file ;
  ofstream _molecules_file ;
  ofstream _positions_backup ;
  ofstream _velocities_backup ;

  /*******************************************   EVOLUTION VARIABLES   *******************************************/
  /* They are not defined as static in their reference functions because they would be slower (more checks needed at run-time . . .) */
  /* Indeed local declaration for standard types could be even better, I do not know . . . */

          // MC variables
  int MCS_LENGTH ;   // number of changing attempts for each MonteCarlo step
  double MC_translation_probability ;     // variable defining the probability that each MonteCarlo attempt would be a translational move
  double *bias_NoP_potential ;  // for GC US, not used until now
  int bias_Nmin , bias_Nmax ;
          // MD variables
  NoseHooverThermostat *NHC_thermostat ;                 // This array gathers all the indipendent NH chains coupled to the system. For NVT_NHchains propagation method.
  double NoseHooverFriction , NoseHooverMass ;       // These variables are used as thermostat in NVT_NHequation method
  int NoseHoover_degreesOFfreedom , vcm_reset_interval ;  // vcm_reset_interval is the number of steps between successive reset of vcm and center of mass during NH equation integration
  double langevin_xi ;    // Diffusion coefficient and friction used to integrate the system with smoluchovski and langevin dynamics, respectively
  double sqrt_2kT_over_xidt ;         // Supporting variables for Brownian (Smoluchovski) dynamics integrator
  particle **dr_gauss = NULL , **dv_gauss = NULL , brown_dv_gauss ;           // Supporting variables for Langevin dynamics integrator
  double langevin_C0 , langevin_C1 , langevin_C2 , sigma_dr , sigma_dv , correlation_dv_dr , uncorrelation_dv_dr ; // Supporting variables for Langevin dynamics integrator


  /*******************************************   CELL - VERLET METHOD   *******************************************/

public:
  inline bool check_verlet_update();  // This member checks if the verlet list has to be updated, triggering it if needed
  inline bool check_verlet_update( particle *part );  // same as previous method, only checks one particle
                                                  // VERLET METHOD
  inline void right_copy( particle_3D *jptr , particle_3D *iptr , particle_3D *copy ) ;  // In a box with periodic boundary conditions this member returns the j-th part copy that is nearest to i-th part
  inline void right_copy( particle_2D *jptr , particle_2D *iptr , particle_2D *copy ) ;  // In a box with periodic boundary conditions this member returns the j-th part copy that is nearest to i-th part


  /*******************************************   BLOCK 2 - DEFINITIONS   *******************************************/
protected:
  inline void check_action( int current_step , int starting_time , int screen_print_interval ) ;   // This function contains all the routines deciding if saving at each integration cycle


  /*******************************************   INTEGRATORS ( Integrators/ )   *******************************************/

                                                  // BLOCK 3 - MD NVE
  inline void VelocityVerletIntegrationStep( void ) ;
                                                 // BLOCK 4 - MD NVT
  void initialize_NHC_thermostat( void ) ;         // in this function thermostats are created and basis for the integration in Yoshida scheme are defined
  inline void NHC_integration_step( double time_delta ) ;
  inline void Factorized_NHC_evolutionOperator( double time_delta ) ;            // this algorithm follows factorization of NHC operator given in Tuckerman's book (p. 194)
  inline void NoseHooverEquation_direct_integration( void ) ;              // this is a version of velocity verlet algorithm for Nosè-Hoover equation integration
  inline virtual void NHE_cm_reset( void ) ;   // This function resets the velocity of the center of mass for the correct NHE integration
                                                 // BLOCK 5 - MD BROWNIAN DYNAMICS
  inline void BROWN_integration_step( double time_delta ) ;
                                                 // BLOCK 6 - MD LANGEVIN DYNAMICS
  inline void LANGEVIN_integration_step( double time_delta ) ;

                                                 // MC dynamics
  inline void MonteCarloSwapChargeStep( void ) ;
  inline void MonteCarlo_NVT( particle old_position ) ;       // translational moves only in the NVT ensemble (also rotational if particle is patchy)
  inline void MonteCarlo_NVT_CMfixed( particle old_position ) ;       // MC NVT with the constraint of fixed centre of mass
  inline void MonteCarlo_NPT( particle old_position , double *old_coordinates ) ;       // translational and volume moves in the NPT ensemble (also rotational if particle is patchy)
  inline void volume_move( int direction , double *old_coordinates ) ;
  inline void volume_move_isotropic( void ) ;


  /*******************************************   DATA SAVING ( data_saving.h )   *******************************************/

  void plot_equilibration_data( void ) ;       // produce plots of energy and number of particles versus MC-time
  char *generate_new_directory_for_data( void ) ;              // this function generate a new directory for saving data
  char *generate_equilibration_subdirectory( char *root ) ;        // this function generate a subdirectory in which are generated storage files for equilibration data. Returns Root-directory name
  void initialize_storage_file( int mcsteps, double seed ) ;  // opens the storage file and prints the headlines
  inline void save_configuration( int step ) ;                        // saves the configuration in positions' and velocities' files
  inline void save_data( int step ) ;                        // saves global quantities in file thermoinfo_file
  inline virtual void save_runtime_calculations( int step ) ;        // attempt to save quantities calculated and averaged on the fly
  virtual void clear_runtime_calculations( void ) ;        // attempt to clear the quantities calculated and averaged on the fly
  void set_saving_options( const char *MODE ) ;                     // This function allows to stop and resume the saving data process
  void close_storage_file( void ) ;
  void save_simulation_conditions( const char *file_name ) ;  // this function save in the simulation folder a copy of configuration file with a new name
  inline void write_backup_data( int current_step ) ;                    // This fuction update two files in which the program stores velocities and positions every (backup_interval) steps


  /*******************************************   PRINT DEBUG INFORMATION ( print_info.h )   *******************************************/

  void print_all_interactions( void ) ;
  void print_all_information( void ) ;
  void print_displacements( void ) ;
  void print_max_velocity( void ) ;
  double max_component_of_vector( particle **vec , int length ) ;
  double *min_part_distance( void ) ;
  double *max_part_distance( void ) ;

};


#include "configuration.tpp"

#include "Integrators/molecular_dynamics.tpp"
#include "Integrators/montecarlo.tpp"
#include "conf_box.tpp"
#include "conf_init.tpp"
#include "conf_velocities.tpp"
#include "conf_charges.tpp"
#include "conf_save.tpp"
#include "conf_print_info.tpp"
#include "conf_read.tpp"
#include "conf_dump.tpp"


#endif
