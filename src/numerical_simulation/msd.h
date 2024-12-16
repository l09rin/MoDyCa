#ifndef MSD_H
#define MSD_H
/* methods to compute the mean squared displacement on the fly */

template <typename particle>
class mean_squared_displacement {
public:
  int Nparts = 0 , saving_step ;
  particle *box_sides = NULL , **parts = NULL ;
  ofstream _MSD_file ;
  char _msd_filename[600] ;
  ifstream _saving_steps_file ;
  char _steps_path[600] ;

  void save_msd_conf( double tempo , particle boxSides ) ;          // mode = "all" for all the particles, "mols_cm" for molecules centers of mass
  void initialize_MSDfiles( const char *directory_path ) ;
  ~mean_squared_displacement() {
    if( parts != NULL ) delete parts ;
    parts = NULL ;
    Nparts = 0 ;
  };
};
/****************************************************************************************/

template <typename particle>
void mean_squared_displacement<particle>::save_msd_conf( double tempo , particle boxSides ) {
  _MSD_file << "# N " << Nparts << endl ;
  _MSD_file << "# timestep " << tempo << endl ;
  for( int i=0 ; i<Nparts ; i++ ) _MSD_file << parts[i]->unwrapped( boxSides.position ) << endl ;
  _MSD_file << endl ;
  if( _saving_steps_file.eof() ) {
    saving_step = -1 ;
  } else {
    _saving_steps_file >> saving_step ;
  }
};
/****************************************************************************************/

template <typename particle>
void mean_squared_displacement<particle>::initialize_MSDfiles( const char *directory_path ) {
  char file_path[300] ;
  strcpy( file_path , directory_path ) ;
  strcat( file_path , "/" ) ;
  if( strcat( file_path , _msd_filename ) == NULL) cout << " Problems in initialize() of mean_squared_displacement object\n" ;
  _MSD_file.open( file_path ) ;
  if( _MSD_file.fail() ) {
    cout << "Failure in the opening of the MSD configurations file !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  _MSD_file << setprecision(12) ;
  _saving_steps_file.open( _steps_path ) ;
  if( _saving_steps_file.fail() ) {
    cout << "Failure in the opening of the MSD saving-steps file ! File " << _steps_path << " not found !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  _saving_steps_file >> saving_step ;
};
/****************************************************************************************/

#endif
