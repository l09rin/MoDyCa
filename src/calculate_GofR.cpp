#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <thread>
#include <mutex>
#include <chrono>
using namespace std;
#include "system_comm.h"
#include "numerical_simulation/configuration.h"

class filelist {
public:
  char *address = NULL ;
  char format[5] ;       // lmp or xyz
  int NofConfigurations ;
};


mutex mtx ;



template <typename particle>
void program_body( int argc, char *argv[] ) ;

template <typename particle>
double * gR_calculation( int points , double delta_R , double max_distance , configuration<particle> &config , int NO_PBC , int CROSS_CALC , int type1 = 0 , int type2 = 0 ) ;
template <typename particle>
inline double shell_volume( int index , double delta_R ) { cout << "*** Unrecognized particle class !" << endl ; exit(EXIT_FAILURE) ;};

template <typename particle>
double ** gRtheta_calculation( int r_pts , int theta_pts , double delta_R , double delta_theta , double max_distance , configuration<particle> &config , int NO_PBC , int CROSS_CALC , int type1 = 0 , int type2 = 0 ) { cout << "*** This function has not yet implemented !" << endl ; exit(EXIT_FAILURE) ;};



template <typename particle>
void calculate_singleconf_gr( int &conf_counter , filelist *configuration_files , int &current_file , int nofFiles , ifstream &data_file , int &analysed_confs , int &bad_confs , int &totalNconfs , int YES_THETA , double *** &grtheta_of_conf , double ** &gr_of_conf , int r_pts , int theta_pts , double delta_R , double delta_theta , double max_distance , int NO_PBC , int CROSS_CALC , int type1 , int type2 , configuration<particle> *configs , int Nthread ) {
  int conf_number = 0 ;
  char file_path[600] ;
  bool STOP_FLAG = 0 ;
  mtx.lock() ;
  cout << "Thread n. " << Nthread << " started." << endl ;
  mtx.unlock() ;
  while ( STOP_FLAG == 0 ) {
    mtx.lock() ;
    if ( data_file.eof() ) {
      if ( current_file < nofFiles - 1 ) {
	current_file ++ ;
	data_file.close() ;
	data_file.open( configuration_files[current_file].address ) ;
	if ( ! data_file.is_open() ) { printf("*** ERROR: The specified file does not exist! [%s]" , configuration_files[current_file].address) ; STOP_FLAG = 1 ; }
	conf_counter = 0 ;
      } else STOP_FLAG = 1 ;
    } else {
      conf_counter ++ ;
      // if configuration_files[i].NofConfigurations == 0 the program reads all the configurations stored in that file
      if( configuration_files[current_file].NofConfigurations > 0 && conf_counter > configuration_files[current_file].NofConfigurations ) {
	if ( current_file < nofFiles - 1 ) {
	  current_file ++ ;
	  data_file.close() ;
	  data_file.open( configuration_files[current_file].address ) ;
	  if ( ! data_file.is_open() ) { printf("*** ERROR: The specified file does not exist! [%s]" , configuration_files[current_file].address) ; STOP_FLAG = 1 ; }
	  conf_counter = 0 ;
	} else STOP_FLAG = 1 ;
      }
    }

    if ( STOP_FLAG == 0 ) {
      analysed_confs ++ ;
      conf_number = analysed_confs ;
      if ( analysed_confs > totalNconfs ) {
	totalNconfs += 100 ;
	if( YES_THETA == 1 ) {
	  grtheta_of_conf = ( double ***)realloc( grtheta_of_conf , totalNconfs*sizeof(double **) ) ;
	  for( int j=totalNconfs-100 ; j<totalNconfs ; j++ ) grtheta_of_conf[j] = NULL ;
	} else {
	  gr_of_conf = ( double **)realloc( gr_of_conf , totalNconfs*sizeof(double *) ) ;
	  for( int j=totalNconfs-100 ; j<totalNconfs ; j++ ) gr_of_conf[j] = NULL ;
	}
      }

      // loading of the configuration
      configs[Nthread].read_configuration( data_file , configuration_files[current_file].format ) ;
      if ( data_file.eof() ) {
	if ( current_file < nofFiles - 1 ) {
	  current_file ++ ;
	  data_file.close() ;
	  data_file.open( configuration_files[current_file].address ) ;
	  if ( ! data_file.is_open() ) { printf("*** ERROR: The specified file does not exist! [%s]" , configuration_files[current_file].address) ; STOP_FLAG = 1 ; }
	  conf_counter = 0 ;
	} else STOP_FLAG = 1 ;
      }
      if( configs[Nthread].numberOfPart <= 0 ) analysed_confs -- ;
      if( STOP_FLAG == 0 && configs[Nthread].numberOfPart <= 0 ) bad_confs ++ ;
    }
    mtx.unlock() ;
 
    // computation of the pair distribution function for the configuration loaded
    if( configs[Nthread].numberOfPart > 0 ) {
      mtx.lock() ;
      printf( "I'm scanning configuration n. %d, containing %d particles\n" , conf_number , configs[Nthread].numberOfPart ) ; fflush(0) ;
      mtx.unlock() ;
      sprintf( file_path , ".grtmp/conf_%d.dat" , conf_number ) ;   // Temporary files where g(r)s are stored for each configuration
      ofstream tmp_file ;
      tmp_file.open( file_path ) ;
      if ( YES_THETA == 1 ) {
	double **tmpp = gRtheta_calculation<particle>( r_pts , theta_pts , delta_R , delta_theta , max_distance , configs[Nthread] , NO_PBC , CROSS_CALC , type1 , type2 ) ;
	mtx.lock() ;
	grtheta_of_conf[conf_number-1] = tmpp ;
	mtx.unlock() ;
	tmp_file << "# r\t theta\t g(r,theta)\n" ;
	for( int r_i=0 ; r_i<r_pts ; r_i++ ) {
	  for( int theta_j=0 ; theta_j<theta_pts ; theta_j++ ) tmp_file << (r_i+0.5)*delta_R << '\t' << (theta_j+0.5)*delta_theta << '\t' << tmpp[r_i][theta_j] << '\n' ;
	}
      } else {
	double *tmpp = gR_calculation<particle>( r_pts , delta_R , max_distance , configs[Nthread] , NO_PBC , CROSS_CALC , type1 , type2 ) ;
	mtx.lock() ;
	gr_of_conf[conf_number-1] = tmpp ;
	mtx.unlock() ;
	tmp_file << "# r\t g(r)\n" ;
	for( int index=0 ; index<r_pts ; index++ ) tmp_file << (index+0.5)*delta_R << '\t' << tmpp[index] << '\n' ;
      }
      tmp_file.close() ;
      configs[Nthread].numberOfPart = 0 ;
      mtx.lock() ;
      cout << "Configuration " << conf_number << " done." << endl ;
      fflush( 0 ) ;
      mtx.unlock() ;
    }

  }
}



int main( int argc, char *argv[] ) {

  cout << "========================================================================================" << endl ;
  cout << "This program calculate the pair distribution function on a set of configurations stored in " << endl ;
  cout << "one or more files in which each configuration:" << endl ;
  cout << "      @-  starts with a row containing the character '#' or with a blank line" << endl ;
  cout << "      @-  each row contains a single particle's coordinates x,y,z in columns 0,1,2" << endl ;
  cout << "      @-  ends with a row containing the character '#' or with a blank line" << endl ;
  cout << endl ;
  cout << " Usage tips:" << endl ;
  cout << "      ./calculate_GofR.exe   -box_sides <X_length> <Y_length> {<Z_length>} | -box_side <L_in_SIGMA_units>" << endl
       << "                           { -D 2|3    [space dimension, default = 3] }" << endl
       << "                           { -max <L_in_SIGMA_units>}" << endl
       << "                           [ -delta <deltaR_in_SIGMA_units> | -N <average_number_of_particles> | -points <number_of_points> ]" << endl
       << "                           { -cross <part_type_1> <part_type_2> }" << endl
       << "                           { -nopbc }" << endl
       << "                           { -theta_D <delta_theta> | -theta_D <delta_theta> }" << endl
       << "                           { -out <out_file_name> }" << endl
       << "                           { -parallel <n_threads> }" << endl
       << "                           { -format xyz|lmp|sph }" << endl
       <<"                              -files configuration_file_1 { -NoC <number_of_configurations_in_file> | -all } configuration_file_2 ......." << endl << endl ;
  cout << "    last option must be -files ;" << endl ;
  cout << "    the option -nopbc can be used to avoid considering nearest neighbours to compute pair distances ." << endl ;
  cout << "    the options -theta_D or -theta_N can be used to compute the pair distribution function as a function of both angular and radial pair distance ." << endl << endl ;

  int number_of_arguments = argc , D = 3 ;

  for( int i=1; i<number_of_arguments; i++ ) {
    if( strcmp(argv[i], "-D") == 0 ) {
      sscanf( argv[i+1], "%d", &D ) ;
      cout << "Changed dimension of the space : " << D << endl ;
      i++ ;
    }
  }

  if( D == 3 ) program_body<particle_3D>( argc , argv ) ;
  else if( D == 2 ) program_body<particle_2D>( argc , argv ) ;
  else cout << "*** " << D << " is not a valid dimension !!" << endl ;
}


template <typename particle>
void program_body( int argc, char *argv[] ) {

  int number_of_arguments = argc ;
  int command_output = 0 , averageN = 0 , r_pts = 0 , theta_pts = 0 ;
  double delta_R = -1 , max_distance = -1 , delta_theta = -1 ;
  configuration<particle> *configs = NULL ;
  particle box_sides ;
  box_sides.set_equal_comp( -1.0 ) ;
  filelist *configuration_files = NULL ;
  configuration_files = new filelist [number_of_arguments] ;
  char FORMAT[4] = "xyz" ;
  int nofFiles = 0 ;
  int NO_PBC = 0 , YES_THETA = 0 ;
  int CROSS_CALC = 0 , type1 = 0 , type2 = 0 ;
  int PARALLEL_THREADS = 1 ;
  thread **multi_thread = NULL ;
  char out_name[600] = "pair_dist_function.dat" ;

  for( int i=1; i<number_of_arguments; i++ ) {

    if( strcmp(argv[i], "-box_sides") == 0 ) {
      box_sides.read_by_string( argv+i+1 ) ;
      cout << "Inserted box edges length: " << box_sides.position << endl ;
      i++ ;

    } else if( strcmp(argv[i], "-box_side") == 0 ) {
      sscanf( argv[i+1] , "%lf" , &(box_sides.position.x) ) ;
      box_sides.set_equal_comp( box_sides.position.x ) ;
      cout << "Inserted box edges length: " << box_sides.position << endl ;
      i++ ;

    } else if( strcmp(argv[i], "-out") == 0 ) {
      sscanf( argv[i+1], "%s", out_name ) ;
      i++ ;

    } else if( strcmp(argv[i], "-nopbc") == 0 ) {
      NO_PBC = 1 ;

    } else if( strcmp(argv[i], "-max") == 0 ) {
      sscanf( argv[i+1], "%lf", &max_distance ) ;
      cout << "Inserted maximum distance: " << max_distance << endl ;
      i++ ;

    } else if( strcmp(argv[i], "-delta") == 0 ) {
      sscanf( argv[i+1], "%lf", &delta_R ) ;
      cout << "Inserted bin length: " << delta_R << endl ;
      i++ ;

    } else if( strcmp(argv[i], "-N") == 0 ) {
      sscanf( argv[i+1], "%d", &averageN ) ;
      cout << "Inserted average number of particles: " << averageN << endl ;
      i++ ;

    } else if( strcmp(argv[i], "-points") == 0 ) {
      sscanf( argv[i+1], "%d", &r_pts ) ;
      cout << "Number of points for which g(r) will be computed: " << r_pts << endl ;
      i++ ;

    } else if( strcmp(argv[i], "-parallel") == 0 ){
      sscanf( argv[i+1], "%d", &PARALLEL_THREADS ) ;
      i++ ;

    } else if( strcmp(argv[i], "-files") == 0 ) {
      i ++ ;
      while( i < number_of_arguments ) {
	configuration_files[nofFiles].address = argv[i] ;
	sprintf( configuration_files[nofFiles].format , "%s", FORMAT ) ;
	if( i+1 < number_of_arguments ) {
	  if( strcmp(argv[i+1], "-NoC") == 0 ) {
	    if( i+2 < number_of_arguments ) {
	      sscanf( argv[i+2], "%d", &(configuration_files[nofFiles].NofConfigurations) ) ;
	      i += 2 ;
	    } else {
	      cout << "Invalid format in the specification of configuration files" << endl ;
	      exit(EXIT_FAILURE) ;
	    }
	  } else if( strcmp(argv[i+1], "-all") == 0 ) {
	    configuration_files[nofFiles].NofConfigurations = 0 ;
	    i += 1 ;
	  } else {
	    configuration_files[nofFiles].NofConfigurations = 1 ;
	  }
	} else {
	  configuration_files[nofFiles].NofConfigurations = 1 ;
	}
	nofFiles ++ ;
	i += 1 ;
	cout << "Opening: " << configuration_files[nofFiles-1].address << '\n' ;
	fflush(0) ;
      }

    } else if( strcmp(argv[i], "-format") == 0 ) {
      if( strcmp(argv[i+1], "lmp") == 0 ) {
	sprintf( FORMAT , "lmp" ) ;
	// box sides will be read from the configurations file
	box_sides.set_equal_comp( 0.0 ) ;
	cout << "Box edges will be read by file" << endl ;
      } else if( strcmp(argv[i+1], "xyz") == 0 ) {
	sprintf( FORMAT , "xyz" ) ;
      } else if( strcmp(argv[i+1], "sph") == 0 ) {
	sprintf( FORMAT , "sph" ) ;
      } else {
	printf("+++ Format not recognized, default option: xyz\n");
      }
      i ++ ;

    } else if( strcmp(argv[i], "-cross") == 0 ) {
      CROSS_CALC = 1 ;
      sscanf( argv[i+1], "%d", &type1 ) ;
      sscanf( argv[i+2], "%d", &type2 ) ;
      cout << "The cross radial distribution function will be calculated among pairs of type " << type1 << " and " << type2 << endl ;
      i+=2 ;

    } else if( strcmp(argv[i], "-theta_D") == 0 ) {
      YES_THETA = 1 ;
      sscanf( argv[i+1], "%lf", &delta_theta ) ;
      i++ ;

    } else if( strcmp(argv[i], "-theta_N") == 0 ) {
      YES_THETA = 1 ;
      sscanf( argv[i+1], "%d", &theta_pts ) ;
      i++ ;

    } else if( strcmp(argv[i], "-help") == 0  ||  strcmp(argv[i], "--help") == 0  ) {
      exit( EXIT_SUCCESS ) ;
    }
  }
  fflush(0) ;

  if( box_sides.min_comp() < 0 && box_sides.max_comp() < 0 ) {
    cout << "ERROR:  You have to insert the values of the Box sides to properly account for PBC !" << endl ;
    exit(EXIT_FAILURE) ;
  }
  if( max_distance < 0 ) {
    if ( box_sides.min_comp() == 0 ) {
      cout << "ERROR:  You have to insert the value of the maximum distance to compute the g(r) !" << endl ;
      exit(EXIT_FAILURE) ;
    } else max_distance = box_sides.min_comp() / 2.0 ;
  }
  if( averageN > 0 ) {
    delta_R = max_distance / averageN ;
    r_pts = averageN ;
  }else if( r_pts > 0 ) {
    delta_R = max_distance / r_pts ;
  }else if( delta_R > 0 ) {
    r_pts = (int)ceil(max_distance/delta_R) ;
    max_distance = delta_R * r_pts ;
  }else if( delta_R <= 0 ) {
    cout << "ERROR:  You must insert a value for the bin-length or a value for the average number of particles !" << endl ;
    exit(EXIT_FAILURE) ;
  }
  if( YES_THETA == 1 ) {
    if ( delta_theta > 0 ) {
      if ( delta_theta < M_PI ) {
	theta_pts = (int) ceil( M_PI / delta_theta ) ;
	delta_theta = M_PI / theta_pts ;
	cout << "The angular distribution function will be calculated with " << theta_pts << " angular bins of width " << delta_theta << endl ;
      } else {
	cout << " *** ERROR: the angular bin must be positive and smaller than pi !" << endl ;
	exit( EXIT_FAILURE ) ;
      }
    } else if ( theta_pts > 0 ) {
      delta_theta = M_PI / theta_pts ;
      cout << "The angular distribution function will be calculated with " << theta_pts << " angular bins of width " << delta_theta << endl ;
    } else {
      cout << " *** ERROR: to compute the angular pair distribution function insert a valid value for the number or width of bins !" << endl ;
      exit( EXIT_FAILURE ) ;
    }
  }
  // Preparation of the parallel calculus
  multi_thread = new thread* [PARALLEL_THREADS] ;
  configs = new configuration<particle> [PARALLEL_THREADS] ;
  for ( int i=0 ; i<PARALLEL_THREADS ; i++ ) configs[i].box_sides = box_sides ;

  // Creation of the array containing the pair distribution functions for each configuration
  // calculation of the form factor g(r)
  int totalNconfs = 0 ;
  double **gr_of_conf = NULL , ***grtheta_of_conf = NULL ;

  // Now, for each file and each configuration I compute the pair distribution function and save it in the array gr_of_conf[]
  int analysed_confs = 0 , bad_confs = 0 ;
  // I create a hidden directory where I save partial data, to recover in case of bad termination
  command_output = system( "mkdir .grtmp" ) ;
  if (command_output < 0) cout << "*** ERROR: mkdir .grtmp " << endl ;
  int current_file=0 ;
  ifstream data_file ;
  data_file.open( configuration_files[current_file].address ) ;
  if ( ! data_file.is_open() ) { printf("*** ERROR: The specified file does not exist! [%s]" , configuration_files[current_file].address) ; exit(EXIT_FAILURE); }

  int conf_counter = 0 ;
  for ( int l=0 ; l<PARALLEL_THREADS ; l++ ) {
    multi_thread[l] = new thread( calculate_singleconf_gr<particle> , ref(conf_counter) , configuration_files , ref(current_file) , nofFiles , ref(data_file) , ref(analysed_confs) , ref(bad_confs) , ref(totalNconfs) , YES_THETA , ref(grtheta_of_conf) , ref(gr_of_conf) , r_pts , theta_pts , delta_R , delta_theta , max_distance , NO_PBC , CROSS_CALC , type1 , type2 , configs , l ) ;
    this_thread::sleep_for( chrono::milliseconds(100) ) ;	
  }
  for ( int l=0 ; l<PARALLEL_THREADS ; l++ ) {
    multi_thread[l]->join() ;	
  }

  data_file.close() ;
  cout << endl ;


  // calculation of averaged g(r)
  cout << endl << "I'm calculating the averaged g(r) over " << analysed_confs << " configurations; I found " << bad_confs << " bad configurations" << endl ;
  fflush(0) ;
  double *total_gr = NULL , **total_grtheta = NULL ;
  if ( YES_THETA == 1 ) {
    total_grtheta = new double* [r_pts] ;
    for( int i=0; i<r_pts; i++ ) {
      total_grtheta[i] = new double [theta_pts] ;
      for( int j=0; j<theta_pts; j++ ) {
	total_grtheta[i][j] = 0 ;
	for( int l=0 ; l<analysed_confs ; l++ )  total_grtheta[i][j] += grtheta_of_conf[l][i][j] ;
	total_grtheta[i][j] /= analysed_confs ;
      }
    }
  } else {
    total_gr = new double [r_pts] ;
    for( int i=0; i<r_pts; i++ ) {
      total_gr[i] = 0 ;
      for( int l=0 ; l<analysed_confs ; l++ )  total_gr[i] += gr_of_conf[l][i] ;
      total_gr[i] /= analysed_confs ;
    }
  }

  // saving data on average g(r)
  cout << endl << "I'm saving final data . . . ." << endl ;
  fflush(0) ;
  ofstream out_file ;
  out_file.open( out_name ) ;
  if ( YES_THETA == 1 ) {
    out_file << "# r\t theta\t g(r,theta)\n" ;
    for( int i=0; i<r_pts; i++ ) {
      for( int j=0; j<theta_pts; j++ )  out_file << ( (i+0.5)*delta_R ) << '\t' << ( (j+0.5)*delta_theta ) << '\t' << total_grtheta[i][j] << '\n' ;
    }
  } else {
    out_file << "#r\t g(r)\n" ;
    for( int i=0; i<r_pts; i++ )  out_file << ( (i+0.5)*delta_R ) << '\t' << total_gr[i] << '\n' ;
  }
  out_file.close() ;


  // Now I free all the occupied memory and...
  delete[] configuration_files ;
  delete [] configs ;
  if ( YES_THETA == 1 ) {
    for( int l=0 ; l<analysed_confs ; l++ ) {
      for( int i=0 ; i<r_pts ; i++ ) delete[] grtheta_of_conf[l][i] ;
      delete[] grtheta_of_conf[l] ;
    }
    free( grtheta_of_conf ) ;
    for( int i=0 ; i<r_pts ; i++ ) delete[] total_grtheta[i] ;
    delete[] total_grtheta ;
  } else {
    for( int j=0 ; j<analysed_confs ; j++ ) delete[] gr_of_conf[j] ;
    free( gr_of_conf ) ;
    delete[] total_gr ;
  }
  // ... then I finally quit
  command_output = system( "rm -r .grtmp" ) ;   // If nothing has gone wrong I can eliminate the temporary directory containing single-configuration g(r)s
  if (command_output < 0) cout << "*** ERROR: rm -r .grtmp " << endl ;

  exit( EXIT_SUCCESS ) ;
}


/***************************************************************************************************************/


template <typename particle>
double * gR_calculation( int points , double delta_R , double max_distance , configuration<particle> &config , int NO_PBC , int CROSS_CALC , int type1 , int type2 ) {
  long int type1_N = 0 , type2_N = 0 , N_distinct_pairs = 0 ;
  double *gr_array = NULL ;
  gr_array = new double [points] ;
  for( int e=0 ; e<points ; e++ ) gr_array[e] = 0 ;

  // calculation of g(r) for a single configuration
  double dist = 0. ;
  int index = 0 ;

  if( NO_PBC == 1 && CROSS_CALC == 0 ) {
    for( int Ipart=0 ; Ipart<config.numberOfPart-1 ; Ipart++ ) {
      for( int Jpart=Ipart+1 ; Jpart<config.numberOfPart ; Jpart++ ) {
	dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->position ) ;
	if( dist < max_distance ) {
	  index = (int) ceil( dist/delta_R ) ;
	  gr_array[index-1] += 1.0 ;
	}
      }
    }
    N_distinct_pairs = (long int)config.numberOfPart*(config.numberOfPart-1)/2 ;

  } else if( NO_PBC == 0 && CROSS_CALC == 0 ) {
    long int yes=0,no=0;
    for( int Ipart=0 ; Ipart<config.numberOfPart-1 ; Ipart++ ) {
      for( int Jpart=Ipart+1 ; Jpart<config.numberOfPart ; Jpart++ ) {
	dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->PBC_NearNeigh2( *config.particles[Ipart] , config.box_sides ).position ) ;
	if( dist < max_distance ) {
	  index = (int) ceil( dist/delta_R ) ;
	  gr_array[index-1] += 1.0 ;
	  yes++;
	}else no++;
      }
    }
    long int tot=0;
  for( index=0 ; index<points ; index++ ) {
    tot+=gr_array[index];
  }
    N_distinct_pairs = (long int)config.numberOfPart*(config.numberOfPart-1)/2 ;

  } else if( CROSS_CALC == 1 ) {
    for( int Ipart=0 ; Ipart<config.numberOfPart ; Ipart++ ) {
      if ( config.particles[Ipart]->type == type1 ) {
	type1_N ++ ;
	for( int Jpart=0 ; Jpart<config.numberOfPart ; Jpart++ ) {
	  if ( config.particles[Jpart]->type == type2 ) {
	    if ( Ipart != Jpart ) {
	      if( NO_PBC == 0 ) dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->PBC_NearNeigh2( *config.particles[Ipart] , config.box_sides ).position ) ;
	      else dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->position ) ;
	      if( dist < max_distance ) {
		index = (int) ceil( dist/delta_R ) ;
		gr_array[index-1] += 1.0 ;
	      }
	    }
	  }
	}
      } else if ( config.particles[Ipart]->type == type2 ) type2_N ++ ;
    }
    if ( type1 != type2 ) N_distinct_pairs = type1_N * type2_N ;
    else N_distinct_pairs = type1_N*(type1_N-1) ;

  }
  if ( N_distinct_pairs < 0 ) { printf("*** ERROR: Overflow of long int variable!\n") ; fflush(0) ; exit(EXIT_FAILURE) ; }

  for( index=0 ; index<points ; index++ ) {
    gr_array[index] /= ( shell_volume<particle>(index , delta_R) * N_distinct_pairs / config.box_sides.volume() ) ;
  }

  return gr_array ;
}
/***************************************************************************************************************/

template<>
inline double shell_volume<particle_3D>( int index , double delta_R ) {
  return 4.0/3.0*M_PI*( pow(index+1,3.0)-pow(index,3.0) )*pow(delta_R,3.0) ;
}

template<>
inline double shell_volume<particle_2D>( int index , double delta_R ) {
  return M_PI*( pow(index+1,2.0)-pow(index,2.0) )*pow(delta_R,2.0) ;
}
/***************************************************************************************************************/

template <>
double ** gRtheta_calculation<particle_3D>( int r_pts , int theta_pts , double delta_R , double delta_theta , double max_distance , configuration<particle_3D> &config , int NO_PBC , int CROSS_CALC , int type1 , int type2 ) {
  long int type1_N = 0 , type2_N = 0 , N_distinct_pairs = 0 ;
  double **grtheta_array = NULL ;
  grtheta_array = new double* [r_pts] ;
  for( int i=0 ; i<r_pts ; i++ ) {
    grtheta_array[i] = new double [theta_pts] ;
    for( int j=0 ; j<theta_pts ; j++ ) grtheta_array[i][j] = 0 ;
  }

  // calculation of g(r) for a single configuration
  double dist = 0. ;
  int r_i = 0 , theta_j = 0 ;

  if( NO_PBC == 1 && CROSS_CALC == 0 ) {
    for( int Ipart=0 ; Ipart<config.numberOfPart-1 ; Ipart++ ) {
      for( int Jpart=Ipart+1 ; Jpart<config.numberOfPart ; Jpart++ ) {
	dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->position ) ;
	if( dist < max_distance ) {
	  r_i = (int) ceil( dist/delta_R ) ;
	  theta_j = (int) ceil( ( config.particles[Jpart]->position - config.particles[Ipart]->position ).zenit() / delta_theta ) ;
	  grtheta_array[r_i-1][theta_j-1] += 1.0 ;
	}
      }
    }
    N_distinct_pairs = config.numberOfPart*(config.numberOfPart-1)/2 ;

  } else if( NO_PBC == 0 && CROSS_CALC == 0 ) {
    for( int Ipart=0 ; Ipart<config.numberOfPart-1 ; Ipart++ ) {
      for( int Jpart=Ipart+1 ; Jpart<config.numberOfPart ; Jpart++ ) {
	dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->PBC_NearNeigh2( *config.particles[Ipart] , config.box_sides ).position ) ;
	if( dist < max_distance ) {
	  r_i = (int) ceil( dist/delta_R ) ;
	  theta_j = (int) ceil( ( config.particles[Jpart]->PBC_NearNeigh2( *config.particles[Ipart] , config.box_sides ).position - config.particles[Ipart]->position ).zenit() / delta_theta ) ;
	  grtheta_array[r_i-1][theta_j-1] += 1.0 ;
	}
      }
    }
    N_distinct_pairs = config.numberOfPart*(config.numberOfPart-1)/2 ;

  } else if( CROSS_CALC == 1 ) {
    for( int Ipart=0 ; Ipart<config.numberOfPart ; Ipart++ ) {
      if ( config.particles[Ipart]->type == type1 ) {
	type1_N ++ ;
	for( int Jpart=0 ; Jpart<config.numberOfPart ; Jpart++ ) {
	  if ( config.particles[Jpart]->type == type2 ) {
	    if ( Ipart != Jpart ) {
	      if( NO_PBC == 0 ) {
		dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->PBC_NearNeigh2( *config.particles[Ipart] , config.box_sides ).position ) ;
		theta_j = (int) ceil( ( config.particles[Jpart]->PBC_NearNeigh2( *config.particles[Ipart] , config.box_sides ).position - config.particles[Ipart]->position ).zenit() / delta_theta ) ;
	      } else {
		dist = config.particles[Ipart]->position.distance( config.particles[Jpart]->position ) ;
		theta_j = (int) ceil( ( config.particles[Jpart]->position - config.particles[Ipart]->position ).zenit() / delta_theta ) ;
	      }
	      if( dist < max_distance ) {
		r_i = (int) ceil( dist/delta_R ) ;
		grtheta_array[r_i-1][theta_j-1] += 1.0 ;
	      }
	    }
	  }
	}
      } else if ( config.particles[Ipart]->type == type2 ) type2_N ++ ;
    }
    if ( type1 != type2 ) N_distinct_pairs = type1_N * type2_N ;
    else N_distinct_pairs = type1_N*(type1_N-1) ;

  }
  if ( N_distinct_pairs < 0 ) { printf("*** ERROR: Overflow of long int variable!\n") ; fflush(0) ; exit(EXIT_FAILURE) ; }

  for( r_i=0 ; r_i<r_pts ; r_i++ ) {
    for( theta_j=0 ; theta_j<theta_pts ; theta_j++ ) {
      grtheta_array[r_i][theta_j] /= ( 2.0/3.0*M_PI*( cos(theta_j*delta_theta)-cos((theta_j+1)*delta_theta) )*( pow(r_i+1,3.0)-pow(r_i,3.0) )*pow(delta_R,3.0)*N_distinct_pairs/config.box_sides.volume() ) ;
    }
  }

  return grtheta_array ;
}
/***************************************************************************************************************/
