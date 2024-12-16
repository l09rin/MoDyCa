#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;
#include <iomanip>
#include "record.h"




int main( int argc, char **argv ) {

  char *infile = NULL , *outfile = NULL , *binfile = NULL ;
  int col = -1 , Nbins = 0 , N = 0 ;
  double binwidth = 0 , min = 0 , max = 0 , value = 0 ;
  ifstream data_file , bin_file ;
  ofstream out_file ;

  printf( "  This script computes the probability distribution of a series of values contained\n" ) ;
  printf( "     in a file starting with an arbitrary number of commented (# ) or blank lines.\n\n" ) ;
  printf( "       Usage :\n" ) ;
  printf( "                hystogram.sh -in <file_name> -col <column_number> {-binN <number_of_bins> | -binW <bin_width>} {-binfile <file_w_bin_extremes>} -interval <min>:<max> -out <output_file_name>\n" ) ;
  printf( "         if character 'x' is given as an argument for -binN, the number of bins is set automatically to sqrt(N), with N the number of data in the input file\n" ) ;
  printf( "         if character 'x' is given as <min> or <max> argument for -interval, the max/min value over all data in the input file is automatically calculated\n" ) ;

  for( int i=1; i<argc; i++ ) {

    if( strcmp( argv[i], "-col" ) == 0 ) {
      sscanf( argv[i+1], "%d", &col ) ;
      i++ ;

    }else if( strcmp( argv[i], "-in" ) == 0 ) {
      infile = argv[i+1] ;
      i++ ;

    }else if( strcmp( argv[i], "-out" ) == 0 ) {
      outfile = argv[i+1] ;
      i++ ;

    }else if( strcmp( argv[i], "-binfile" ) == 0 ) {
      binfile = argv[i+1] ;
      i++ ;

    }else if( strcmp( argv[i], "-binN" ) == 0 ) {
      if( argv[i+1][0] == 'x' ) {
	data_file.open( infile ) ;
	if( data_file.fail() ) {
	  cout << "\n  Failure in opening the input file" << endl ;
	  exit( EXIT_FAILURE ) ;
	}
	record *line = NULL ;
	line = new record ;
	line->getrecord( data_file ) ;
	while( !data_file.eof() ) {
	  line->split() ;
	  if( line->words_number >= col ) if( strcmp( line->word[0] , "#" ) != 0 ) N ++ ;
	  delete line ;
	  line = new record ;
	  line->getrecord( data_file ) ;
	}
	data_file.close() ;
	data_file.clear() ;
	Nbins = sqrt( (double)N ) ;
      } else sscanf( argv[i+1], "%d", &Nbins ) ;
      i++ ;

    }else if( strcmp( argv[i], "-binW" ) == 0 ) {
      sscanf( argv[i+1], "%lf", &binwidth ) ;
      i++ ;

    }else if( strcmp( argv[i], "-interval" ) == 0 ) {
      record *interv = NULL ;
      interv = new record( argv[i+1] ) ;
      interv->split( ':' ) ;
      sscanf( interv->word[0], "%lf", &min ) ;
      sscanf( interv->word[1], "%lf", &max ) ;
      if( interv->word[0][0] == 'x' || interv->word[1][0] == 'x' ) {
	data_file.open( infile ) ;
	if( data_file.fail() ) {
	  cout << "\n  Failure in opening the input file" << endl ;
	  exit( EXIT_FAILURE ) ;
	}
	record *line = NULL ;
	line = new record ;
	line->getrecord( data_file ) ;
	int firstline = 1 ;
	double tmp = 0 ;
	while( !data_file.eof() ) {
	  line->split() ;
	  if( line->words_number >= col ) {
	    if( strcmp( line->word[0] , "#" ) != 0 ) {
	      if( firstline == 1 ) {
		if( interv->word[0][0] == 'x' ) sscanf( line->word[col-1] , "%lf" , &min ) ;
		if( interv->word[1][0] == 'x' ) sscanf( line->word[col-1] , "%lf" , &max ) ;
		firstline = 0 ;
	      }
	      if( firstline == 0 ) {
		sscanf( line->word[col-1] , "%lf" , &tmp ) ;
		if( tmp < min && interv->word[0][0] == 'x' ) min = tmp ;
		if( tmp > max && interv->word[1][0] == 'x' ) max = tmp ;
	      }
	    }
	  }
	  delete line ;
	  line = new record ;
	  line->getrecord( data_file ) ;
	}
	data_file.close() ;
	data_file.clear() ;
      }
      i++ ;

    }
  }

  char default_name[50] = "hystogram.dat" ;
  if( outfile == NULL ) outfile = default_name ;

  if( binfile != NULL ) {

    Nbins = 10 ;
    int i = 0 ;
    double *bin_extremes = NULL ;
    bin_extremes = (double *) realloc( bin_extremes , Nbins * sizeof(double) ) ;
    bin_file.open( binfile ) ;
    if( bin_file.fail() ) {
      cout << "\n  Failure in opening the bin file" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    record *line = NULL ;
    line = new record ;
    line->getrecord( bin_file ) ;
    while( !bin_file.eof() ) {
      line->split() ;
      if( strcmp( line->word[0] , "#" ) != 0 ) {
	if( i == Nbins ) {
	  Nbins += 10 ;
	  bin_extremes = (double *) realloc( bin_extremes , Nbins * sizeof(double) ) ;
	}
	sscanf( line->word[0] , "%lf" , bin_extremes+i ) ;
	i ++ ;
      }
      delete line ;
      line = new record ;
      line->getrecord( bin_file ) ;
    }
    bin_file.close() ;
    Nbins = i ;
    bin_extremes = (double *) realloc( bin_extremes , Nbins * sizeof(double) ) ;

    double *hyst = new double[Nbins-1] ;
    int j = 0 ;
    N = 0 ;
    for( int i=0 ; i<Nbins-1 ; i++ ) hyst[i] = 0 ;
    data_file.open( infile ) ;
    if( data_file.fail() ) {
      cout << "\n  Failure in opening the input file" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    line = NULL ;
    line = new record ;
    line->getrecord( data_file ) ;
    while( !data_file.eof() ) {
      line->split() ;
      if( line->words_number >= col ) {
	if( strcmp( line->word[0] , "#" ) != 0 ) {
	  sscanf( line->word[col-1] , "%lf" , &value ) ;
	  if( value >= bin_extremes[0] && value <= bin_extremes[Nbins-1] ) {
	    j = Nbins-1 ;
	    while( bin_extremes[j] >= value ) j-- ;
	    if( j < 0 ) j = 0 ;
	    hyst[j] ++ ;
	    N ++ ;
	  }
	}
      }
      delete line ;
      line = new record ;
      line->getrecord( data_file ) ;
    }
    data_file.close() ;

    out_file.open( outfile ) ;
    if( out_file.fail() ) {
      cout << "\n  Failure in opening the output file" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    out_file << setprecision(18) ;
    out_file << "# bin custom" << endl ;
    for( int i=0 ; i<Nbins-1 ; i++ ) out_file << 0.5*( bin_extremes[i] + bin_extremes[i+1] ) << " " << hyst[i]/N/( bin_extremes[i+1] - bin_extremes[i] ) << " " << hyst[i] << endl ;
    out_file.close() ;

  } else {

    if( Nbins > 0 ) binwidth = ( max - min ) / Nbins ;
    else if( binwidth > 0 ) {
      Nbins = (int)floor( ( max - min ) / binwidth ) ;
      binwidth = ( max - min ) / Nbins ;
    }

    printf( "\nInput file\t%s\nOutput file\t%s\nNumber of bins\t%d\nbin width\t%lf\nfile column\t%d\nmin val\t%lf\nmax val\t%lf\n\n" , infile , outfile , Nbins , binwidth , col , min , max ) ;

    double *hyst = new double[Nbins] ;
    int j = 0 ;
    N = 0 ;
    for( int i=0 ; i<Nbins ; i++ ) hyst[i] = 0 ;
    data_file.open( infile ) ;
    if( data_file.fail() ) {
      cout << "\n  Failure in opening the input file" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    record *line = NULL ;
    line = new record ;
    line->getrecord( data_file ) ;
    while( !data_file.eof() ) {
      line->split() ;
      if( line->words_number >= col ) {
	if( strcmp( line->word[0] , "#" ) != 0 ) {
	  sscanf( line->word[col-1] , "%lf" , &value ) ;
	  if( value >= min && value <= max ) {
	    j = floor( ( value - min ) / binwidth ) ;
	    if( j == Nbins ) j-- ;
	    hyst[j] ++ ;
	    N ++ ;
	  }
	}
      }
      delete line ;
      line = new record ;
      line->getrecord( data_file ) ;
    }
    data_file.close() ;

    out_file.open( outfile ) ;
    if( out_file.fail() ) {
      cout << "\n  Failure in opening the output file" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    out_file << setprecision(12) ;
    out_file << "# bin " << binwidth << endl ;
    for( int i=0 ; i<Nbins ; i++ ) out_file << min+(i+0.5)*binwidth << " " << hyst[i]/N/binwidth << " " << hyst[i] << endl ;
    out_file.close() ;

  }

  printf("=============================================================\n");
  exit(EXIT_SUCCESS);
}
