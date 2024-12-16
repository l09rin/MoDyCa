#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
using namespace std;
#include <iomanip>
#include "system_comm.h"
#include "numerical_simulation/configuration.h"

#define PARTICLE_3D 0
#define PARTICLE_2D 1
#define PATCHY_2D 2

// This function read the input script to get the space dimension and the particle class that have to be used to build the configuration class
void preprocess_inputscript( const char* input_script_name , int *D , int *PARTICLE_CLASS , int *SEED ) {
  char part_class_name[16] = "" ;

  read_parameter<int>( "DIMENSION" , input_script_name , D ) ;
  read_parameter<int>( "SEED" , input_script_name , SEED ) ;
  read_parameter<char>( "PARTICLE_CLASS" , input_script_name , part_class_name ) ;
  if( strcmp( part_class_name , "particle_3D") == 0 )  *PARTICLE_CLASS = PARTICLE_3D ;
  else if( strcmp( part_class_name , "particle_2D") == 0 )  *PARTICLE_CLASS = PARTICLE_2D ;
  else if( strcmp( part_class_name , "patchy_2D") == 0 )  *PARTICLE_CLASS = PATCHY_2D ;
  cout << " Particle class : " << *PARTICLE_CLASS << endl ;
};



template <typename particle>
void run_simulation( int D , int seed , char *input_script_file ) {
  configuration<particle> configOFsystem ;
  configOFsystem.seed = seed ;
  srand48( configOFsystem.seed );                      // initialization of the random numbers generator

  configOFsystem.initialize_configuration( input_script_file );    // initialization of the positions' vector and entire configuration of the system

  // generation of the velocities, if needed
  if( strcmp( configOFsystem.simulation_mode , "MD" ) == 0 &&
      configOFsystem.INITIAL_VELOCITIES_BY_FILE == 0 ) {
    if( strcmp( configOFsystem.simulation_ensemble, "NVE" ) == 0 ) {
      configOFsystem.KINETIC_ENERGY.val = configOFsystem.energyOfSystem - configOFsystem.POTENTIAL_ENERGY.val ;
      if( configOFsystem.KINETIC_ENERGY.val < 0 ) {
	cout << endl << "***ERROR: This fixed energy initial configuration has a negative kinetic energy" << endl;
	exit( EXIT_FAILURE ) ;
      }
      // mean square velocity has been settled to fit the total energy initial value
      configOFsystem.generate_MBvelocities( 2.*configOFsystem.KINETIC_ENERGY.val/3/configOFsystem.numberOfPart ) ;
      // kinetic energy has been fixed exactly rescaling velocities
      configOFsystem.set_kinetic_energy( configOFsystem.KINETIC_ENERGY.val ) ;

    } else {
      // mean square velocity has been set to fit the total energy initial value
      configOFsystem.generate_MBvelocities( configOFsystem.k_T ) ;
    }
    configOFsystem.CoM_velocity_reset() ;
    //    configOFsystem.total_angularMomentum_reset() ;
    cout << "  Velocity of center of mass after reset (System units): " << configOFsystem.compute_mean_velocity() << endl << endl ;
  }

  // printing some info
  if( strcmp( configOFsystem.simulation_ensemble, "NVE" ) == 0 ) {
    cout << " * NVE ensemble." << endl ;
    cout << " Total energy per particle = " << configOFsystem.energyOfSystem/configOFsystem.numberOfPart << " (epsilon)" << endl ;
  } else {
    cout << " * NVT ensemble." << endl ;
    cout << " Simulation temperature = " << configOFsystem.k_T << " (epsilon/Kb)" << endl ;
  }
  printf(" Total potential energy    = %lf (epsilon)\n"
	 " Numerical density = %lf (sigma)^(-%d)\n"
	 " packing fraction = %lf \n"
	 " Number of particles in the box = %d\n" , configOFsystem.POTENTIAL_ENERGY.val/configOFsystem.numberOfPart , configOFsystem.numberOfPart/configOFsystem.box_sides.volume() , D , configOFsystem.compute_packing_fraction() , configOFsystem.numberOfPart ) ;
  cout << " Box sides : " << configOFsystem.box_sides.position << " (sigma)\n" ;
  cout << configOFsystem.INTERACTIONS ;

  // configOFsystem.INTERACTIONS.print_on_file( 0.01 , 1000 ) ;

  // running the simulation
  if( strcmp( configOFsystem.simulation_mode, "MC" ) == 0 ) configOFsystem.MonteCarlo( "NO_STOP" ) ;
  else configOFsystem.MD( "NO_STOP" ) ;
};


int main( int argc, char **argv ) {
  int D = 3 , seed = -1 , seed_from_file = -1 , PARTICLE_CLASS = 0 ;
  char input_script_file[300] ;
  strcpy( input_script_file, "MD_input_file.dat" ) ;

  //  reading of the external variables
  for( int i=1; i<argc; i++ ) {
    if( strcmp( argv[i], "-seed" ) == 0 ) {
      sscanf( argv[i+1], "%d", &seed ) ;
      i++ ;
    }else if( strcmp( argv[i], "-in" ) == 0 ) {
      strcpy( input_script_file, argv[i+1] ) ;
      i++ ;
    }
  }

  printf( "=============================================================\n"
	  " Molecular Dynamics and MonteCarlo simulation tool.\n"
	  " Periodic boundary conditions are used.\n" ) ;

  preprocess_inputscript( input_script_file , &D , &PARTICLE_CLASS , &seed_from_file ) ;
  if( seed == -1 ) {
    if( seed_from_file != -1 )  seed = seed_from_file ;
    else  seed = time(0) ;
  }
  if( !(D==2 && (PARTICLE_CLASS==PARTICLE_2D || PARTICLE_CLASS==PATCHY_2D)) &&
      !(D==3 && PARTICLE_CLASS==PARTICLE_3D) ) {
    cout << "*** ERROR : Inconsistent space dimension" << endl ;
    exit(EXIT_FAILURE) ;
  }

  
  if ( PARTICLE_CLASS == PARTICLE_3D )  run_simulation<particle_3D>( D , seed , input_script_file ) ;
  else if ( PARTICLE_CLASS == PARTICLE_2D )  run_simulation<particle_2D>( D , seed , input_script_file ) ;
  else if ( PARTICLE_CLASS == PATCHY_2D )  run_simulation<patchy_2D>( D , seed , input_script_file ) ;

  printf("=============================================================\n");
  exit(EXIT_SUCCESS);
}
