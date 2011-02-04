/************************** SVN Revision Information **************************
 **    $Id: salloc.h 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/* 
 Smart ALLOCation 

 Written by Filipe Ribeiro 2006

 salloc.h: 
  
  Includes macro definitions and function prototypes 
  to deal and keep track of dynamically allocated memory.

  You should NOT use any of the prototyped functions directly. 
  Please use ONLY the following macros:

    my_malloc(...)
    my_calloc(...)
    my_malloc_init(...)
    my_free(...)
    my_free_all()
    my_alloc_report(...)
 
  Also, avoid using malloc(), calloc(), and free() directly.
  Example:

#include "salloc.h"
int main()
{
  double *a;
  int *b;
  my_calloc( a, 100, double );
  my_malloc_init( b, 1000, int );
  . 
  .
  .
  my_free(a);
  my_alloc_report("test");
}
 
This program allocates a 100 double block and sets 'a' to point to it.
Then it allocates a 1000 int block, initializes it to zero, and sets 'b' to point to it.
At the end it frees the memory block to which 'a' points.
Finally it outputs a report to the file 'test.alloc.log' about the allocated memory 
that was not deallocated (in this case the memory block pointed to by 'b') 
plus some statistics
 
*/

/* 
 The following preprocessor statement make sure that 
 this header file is only include once 
*/

#ifndef SALLOC_H
#  define SALLOC_H 1





/*
 ALLOC_M is used to tell salloc() to use malloc()
 ALLOC_C is used to tell salloc() to use calloc()
*/
#  define ALLOC_M  0
#  define ALLOC_C  1


#  define ALLOC_NO_INIT 0
#  define ALLOC_INIT    1






/*
 This is a wrapper macro for all the allocation macros
 the allocation command is actually a parameter of the macro! 

 If nobj == 0 then nothing happens and the code just keeps running.

 Do NOT use this macro directly!
*/
#  define __alloc_wrapper__( _ptr_, _nobj_, _type_, _flag_, _init_ ) \
            ( ((_nobj_) != 0 && \
               (_ptr_ = (_type_*) salloc(_nobj_, sizeof(_type_), #_ptr_, #_type_, \
                                        __FILE__, __LINE__, _flag_, _init_ ) ) == NULL) ? \
               printf("attempting to allocate %d objects of type '%s' of size %d bytes\n",\
                      (int) _nobj_, #_type_, (int) sizeof(_type_)),\
               error_handler("can't allocate memory") : 0 )







/* 
 Use ONLY the following macros when allocating and deallocating memory.
*/ 

#  define my_malloc( _ptr_, _nobj_, _type_ ) \
      __alloc_wrapper__( _ptr_, _nobj_, _type_, ALLOC_M, ALLOC_NO_INIT )

#  define my_calloc( _ptr_, _nobj_, _type_ ) \
      __alloc_wrapper__( _ptr_, _nobj_, _type_, ALLOC_C, ALLOC_NO_INIT )


#  define my_malloc_init( _ptr_, _nobj_, _type_ ) \
      __alloc_wrapper__( _ptr_, _nobj_, _type_, ALLOC_M, ALLOC_INIT )



#  define my_free( _ptr_ )  \
                    ( _ptr_ == NULL ? \
                      error_handler("can't free NULL pointer") : \
                      (sfree(_ptr_,#_ptr_,__FILE__,__LINE__), _ptr_ = NULL,0) )

#  define my_free_all() sfree_all(__FILE__,__LINE__)

#  define my_alloc_report(_tag_) salloc_report(__FILE__,__LINE__, pct.thispe, _tag_)






/* 
 Because SALLOC_C is defined in salloc.c before the #include "salloc.h" statement
 the following preprocessor statement makes sure that the remainder is only included 
 in source files other than salloc.c.
*/

#  ifndef SALLOC_C




/* 
 The following three macros redefine "malloc", "calloc", and "free" 
 so that if you use them you'll get a warning message.
*/

#    define malloc(_s_) \
         (printf("**** please don't use 'malloc' (@ %s,%d) use 'my_malloc' instead! ****\n", __FILE__, __LINE__),\
          malloc(_s_))

#    define calloc(_n_,_s_) \
         (printf("**** please don't use 'calloc' (@ %s,%d) use 'my_calloc' instead! ****\n", __FILE__, __LINE__),\
          calloc(_n_,_s_))

#    define free(_s_) \
         (printf("**** please don't use 'free' (@ %s,%d) use 'my_free' instead! ****\n", __FILE__, __LINE__),\
          free(_s_))






/* 
 We must declare the prototypes of the functions directly used 
 by the macros my_malloc(), my_calloc(), my_free(), etc. defined above.

 Do NOT use these functions directly. Use ONLY the macros.
*/

extern void salloc_debug_level( int d );
extern void *salloc( size_t nobj, size_t size1, 
                     char *vn, char *type, char *fn, int ln, 
                     int flag, int init );
extern void sfree( void *ptr, char *vn, char *fn, int ln );
extern void sfree_all( char *fn, int ln );
extern void salloc_report( char *fn, int ln, int ipe, char *tag );



#  endif


#endif
