/* method to set charges to particles */

/****************************************************************************************/

template <typename particle>
void configuration<particle>::generate_charge_distribution( char *distribution, double charges_fraction , double charges_distribution_sigma ) {
  if( strcmp( distribution, "all" ) == 0 ) {

    if( monomers_number == 0 ) charges_number = numberOfPart ;
    else charges_number = monomers_number ;
    for( int i=0; i<charges_number; i++ ) particles[i]->charge = 1 ;

  }else if( strcmp( distribution, "all_noXl" ) == 0 ) {
    int xlinkers = 0 ;
    for( int i=0; i<numberOfPart; i++ ) {
      if( particles[i]->valence > 2 ) xlinkers ++ ;
    }
    charges_number = numberOfPart - xlinkers ;
    for( int i=0; i<numberOfPart; i++ ) {
      if( particles[i]->valence < 3 ) {
	particles[i]->charge = 1 ;
      }
    }

  }else if( strcmp( distribution, "random" ) == 0 ) {
    charges_number = floor( charges_fraction*numberOfPart ) ;
    int xlinkers = 0 ;
    for( int i=0; i<numberOfPart; i++ ) {
      if( particles[i]->valence > 2 ) xlinkers ++ ;
    }
    if( charges_number > (numberOfPart-xlinkers) ) charges_number = numberOfPart - xlinkers ;
    for( int i=0; i<charges_number; i++ ) {
      int randn = floor( drand48()*numberOfPart ) ;
      if( ( particles[randn]->charge == 0 ) && ( particles[randn]->valence < 3 ) ) {
	particles[randn]->charge = 1 ;
      }else{
	i -- ;
      }
    }

  } else if( strcmp( distribution, "file" ) == 0 ) {

    int monomer_lines = 0 ;
    ifstream data_file ;

    data_file.open( configuration_file_name ) ;
    if( data_file.fail() ) {
      cout << "\n  Failure in opening data_file in function generate_config_by_file" << endl ;
      exit( EXIT_FAILURE ) ;

    }else{
      record *line = NULL ;
      line = new record ;

      line->getrecord( data_file ) ;     // Ignoring first record
      delete line ;
      line = new record ;
      line->getrecord( data_file ) ;     // Reading number of monomers
      line->split() ;
      sscanf( line->word[0] , "%d" , &monomer_lines ) ;
      delete line ;
      line = new record ;
      line->getrecord( data_file ) ;     // Discarding box's properties
      line->split() ;
      delete line ;
      line = new record ;
      for( int i=0; i<monomer_lines; i++ ) {         // Discarding particles' positions
	line->getrecord( data_file ) ;
	delete line ;
	line = new record ;
      }

      int monomer_index = 0, bonds_number = 0 ;
      charges_number = 0 ;
      for( int i=0; i<monomer_lines; i++ ) {         // Reading charges
	line->getrecord( data_file ) ;
	line->split() ;
	if( line->words_number < 3 ) {
	  cout << "ERROR: The file format is incorrect, it does not contain information on charges !" << endl ;
	  exit( EXIT_FAILURE ) ;
	}
	sscanf( line->word[0] , "%d" , &monomer_index ) ;
	sscanf( line->word[1] , "%d" , &bonds_number ) ;
	sscanf( line->word[2] , "%lf" , &( particles[monomer_index-1]->charge ) ) ;
	delete line ;
	line = new record ;
	if( particles[monomer_index-1]->charge != 0 ) {
	  charges_number ++ ;
	}
	if( bonds_number > 0 ) {
	  line->getrecord( data_file ) ;
	  delete line ;
	  line = new record ;
	}
      }
      delete line ;
      line = NULL ;
    }

    data_file.close() ;

  }else{
    cout << "The selected distribution charge is not known" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/****************************************************************************************/
