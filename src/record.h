#ifndef RECORD_H
#define RECORD_H
/* This library contains the class record, which allows to read files line by line and to split them in words */
using namespace std;

class record{
 public:
  char *line = NULL ;
  char **word = NULL ;
  int char_number ;
  int words_number ;

  record( int numberOfChar = 0 ) ;
  record( const char *stringa ) ;
  record( const record &other_rec ) ;
  ~record() ;

  bool getrecord( ifstream &infile ) ;
  void split( void ) ;
  void split( char separation ) ;
};

#endif
