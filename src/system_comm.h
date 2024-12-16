#ifndef SYSTEM_COMM_H
#define SYSTEM_COMM_H
/* some useful function to interact with the user and the system */
#include <iostream>
using namespace std;

template < class T >
inline int check_memalloc( T *root , const char *point ) ;      // Checks if memory is correctly allocated at the pointer root

char *rename(const char *old_name);    // changes the name/path of the file old_name, taking the new one from keyboard
char request(const char *question, const char *options, char separator);      /* Takes as arguments a question, possible answers stored in the string options,
										 separated by the char separator. */
template <typename T>
void read_parameter( const char *param_name , const char *input_script_name , T *param_val ) ;
template <typename T>
void get_param_from_string( const char *param_val , T *val ) ;

template<class T>
const string ClassName( const T *ptr ) ;
class particle_3D ;
class particle_2D ;
class patchy_2D ;


template < typename T >
struct time_dep_variable {
  T val = 0 ;
  double time = -1 ;

  struct time_dep_variable operator=( struct time_dep_variable other ) {
    val = other.val ;
    time = other.time ;
  };
};


#include "system_comm.tpp"

#endif
