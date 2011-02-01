/************************** SVN Revision Information **************************
 **    $Id: salloc.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/* 
 Smart ALLOCation 

 Written by Filipe Ribeiro 2006

 salloc.c: 
  
  Includes functions to deal and keep track of 
  dynamically allocated memory.

  You should NOT use any of these functions directly. 
  Please use ONLY the macros defined in salloc.h:

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

#define SALLOC_C 1

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>


/* 
 Include md.h because of the macro variables
 NEGF1, ORDER_N, ON2, REAL_SPACE, RS2, SMP
*/
#include "md.h"


/* 
 USE_SALLOC is a switch to turn on/off the use of Smart ALLOCation.
 Let's set the default to 1.
 Add USE_SALLOC=0 to the Makefile is you want to switch SALLOC off 
 (although there's really no good reason not to use SALLOC ;-).
*/

#ifndef USE_SALLOC
#  define USE_SALLOC 1
#endif





/* 
 Definition of the node type salloc_node_t.
 Each node points to the previous and next node, and stores information
 on the pointer, size, variable name, filename and line number where the 
 variable was allocated.
*/

typedef struct salloc_node_st *salloc_nodep_t;


typedef struct salloc_node_st 
{ 
    char *fn;          /* filename where salloc was called */
    int ln;            /* line number of call to salloc */
    char *vn;          /* string with the pointer variable name */
    char *type;        /* string with the type name */
    size_t nobj;       /* number of objs */
    size_t size1;      /* size in bytes of 1 obj */
    size_t size;       /* size in bytes */
    void *ptr;
    salloc_nodep_t prev;        /* pointer to the previous node */
    salloc_nodep_t next;        /* pointer to the next node */
    double time;
} salloc_node_t;





/* 
 Declarations of various variables. 
*/


static salloc_nodep_t lastnode = NULL;
static size_t nnodes = 0;
static size_t maxnnodes = 0;
static size_t totalmem = 0;
static size_t maxtotalmem = 0;
/*static size_t totalalloc = 0;*/

static size_t salloc_ncalls = 0;
static double salloc_time = 0.0;
static size_t sfree_ncalls = 0;
static double sfree_time = 0.0;


/* 
 Setting the debug level to 1 makes the code 
 print out a message every time something is allocated or deallocated.
*/
static int debug = 0;
void salloc_debug_level( int d ) { debug = d; }





/* 
 The following routine initializes the allocated data to 0.0 
 for selected types.
*/

static void initialize( void *ptr, size_t n, 
                        char *vn, char *type, char *fn, int ln  )
{
    size_t i;
    if ( strcmp( type, "int" ) == 0 ) 
    {
        int *p = (int *) ptr;
        for ( i = 0; i < n; i ++ ) 
            p[i] = 0;
    }
    else if ( strcmp( type, "double" ) == 0 ) 
    {
        double *p = (double *) ptr;
        for ( i = 0; i < n; i ++ ) 
            p[i] = 0.0;
    }
#if ON2 || RS2
    else if ( strcmp( type, "real_t" ) == 0 ) 
    {
        real_t *p = (real_t *) ptr;
        for ( i = 0; i < n; i ++ ) 
            p[i] = 0.0;
    }
#else
    else if ( strcmp( type, "REAL" ) == 0 ) 
    {
        REAL *p = (REAL *) ptr;
        for ( i = 0; i < n; i ++ ) 
            p[i] = 0.0;
    }
#endif
#if NEGF1 || ORDER_N
    else if ( strcmp( type, "fftw_complex" ) == 0 ) 
    {
        fftw_complex *p = (fftw_complex *) ptr;
        for ( i = 0; i < n; i ++ ) 
            p[i].re = p[i].im = 0.0;
    }
#endif
#if NEGF1  
    else if ( strcmp( type, "doublecomplex" ) == 0 ) 
    {
        doublecomplex *p = (doublecomplex *) ptr;
        for ( i = 0; i < n; i ++ ) 
            p[i].r = p[i].i = 0.0;
    }
#endif

    else 
    {
        printf( "!!!! warning: requested initialization at '%s:%d' of data pointed to by '%s' not done.\n",
                fn, ln, vn );
        printf( "              please inplement initialization procedure for type '%s' in '%s:%d'.\n",
                type, __FILE__, __LINE__ );
    }

}









/* 
 The following are suppose to be a SMP-safe wrappers for 
 malloc(), calloc(), and free() routines (not actually tested with SMP).
*/ 

static void *smp_alloc( size_t nobj, size_t size1, 
                        char *vn, char *type, char *fn, int ln, 
                        int flag, int init )
{
    void *ptr = NULL;
    
    if (nobj > 0) 
    {

#       if SMP
        pthread_mutex_lock(&malloc_mutex);
#       endif        

        if      (flag == ALLOC_M)  ptr = malloc( nobj * size1 );
        else if (flag == ALLOC_C)  ptr = calloc( nobj,  size1 );
        else                       ptr = NULL;
        
#       if SMP
        pthread_mutex_unlock(&malloc_mutex);
#       endif

        if (ptr && init)    
            initialize( ptr, nobj, vn, type, fn, ln );
    }

    return ptr;
}


static void smp_free( void *ptr )
{
#   if SMP
    pthread_mutex_lock(&malloc_mutex);
#   endif        

    free( ptr );

#   if SMP
    pthread_mutex_unlock(&malloc_mutex);
#   endif
}









/*
 Main routine for memory allocation.
 Creates a new list node, takes care of statistics, and calls allocation routine.
*/

void *salloc( size_t nobj, size_t size1, 
              char *vn, char *type, char *fn, int ln, 
              int flag, int init )
{

#if ! USE_SALLOC

    return smp_alloc( nobj, size1, vn, type, fn, ln, flag, init );

#else
    

    double time = -my_crtc();
    salloc_nodep_t node;

    node = (salloc_nodep_t) smp_alloc( 1, sizeof(salloc_node_t), 
                                       "node", "salloc_node_t", __FILE__, __LINE__,
                                       ALLOC_M, ALLOC_NO_INIT );

    if (node == NULL)
        return NULL;
  
    node->prev = node->next = NULL;
  
    if (nnodes > 0) 
    {
        node->prev = lastnode;
        node->prev->next = node;
    }

    lastnode = node;
    nnodes++;

    if (nnodes > maxnnodes) 
        maxnnodes = nnodes;
  
    node->size1 = size1;
    node->nobj = nobj;
    node->size = size1 * nobj;
  
    /* allocate memory */
    node->ptr = smp_alloc( nobj, size1, vn, type, fn, ln, flag, init );

    totalmem += node->size + sizeof(salloc_node_t); /* size in bytes */
  
    if (totalmem > maxtotalmem) 
        maxtotalmem = totalmem;

    /* strcpy() is not necessary because the file name, variable name, and the type */
    /* are always constant strings created by the preprocessor */
    /* with the macros my_malloc() or my_calloc() in salloc.h */
    node->vn = vn;  
    node->type = type;  
    node->fn = fn; 
    node->ln = ln;
    
    node->time = (time += my_crtc());
    if (debug) 
        printf( "salloc: allocated: %p %10d b = %10d x %10d b (%s* %s @ %s:%d %.6fs)\n", 
                node->ptr, (int) node->size, (int) nobj, (int) size1, 
                node->type, node->vn, node->fn, node->ln,
                node->time);
    
    salloc_ncalls ++;
    salloc_time += time;

#  if NEGF1    
    rmg_timings( ALLOC_TIME, time);
#  elif REAL_SPACE || ORDER_N
    rmg_timings( ALLOC_TIME, time, 0);
#  endif     

    return node->ptr;

#endif

}



/*
 Main routine for memory deallocation.
 Searches the list of nodes for the node associated with the specified pointer.
 Destroys one node in the list, takes care of statistics, and calls free routine.
*/


void sfree( void *ptr, char *vn, char *fn, int ln )
{

#if ! USE_SALLOC

    smp_free( ptr );
   
#else

    double time = -my_crtc();
    if (ptr != NULL) 
    {
        salloc_nodep_t node = lastnode;
        int done = 0;
  
        do 
        {
            /*    printf ("    %p == %p ? %s\n", ptr, node->ptr, node->label);*/
            
            if (ptr == node->ptr) 
            {
                if (lastnode == node) 
                    lastnode = node->prev;

                if (node->prev) 
                    node->prev->next = node->next;

                if (node->next) 
                    node->next->prev = node->prev;
      
                smp_free( ptr ); 
                totalmem -= node->size;
      
                if (debug) 
                {
                    printf( "sfree: deallocated: %p %10d b = %10d x %10d b (%s* %s @ %s:%d)\n", 
                            node->ptr, (int) node->size, (int) node->nobj, (int) node->size1, 
                            node->type, vn, fn, ln );
                }

                smp_free( node );
                nnodes--;
                done = 1;
            }
            else
                node = node->prev;

        } while ( node && !done );

        if (!done) 
            printf("sfree: can't find pointer '%p' to free (%s @ %s:%d)\n", 
                   ptr, vn, fn, ln );
        
    }
    time += my_crtc();
    
    sfree_ncalls ++;
    sfree_time += time;

#  if NEGF1    
    rmg_timings( ALLOC_TIME, time);
#  elif REAL_SPACE || ORDER_N
    rmg_timings( ALLOC_TIME, time, 0);
#  endif     

#endif
}




/* 
 Function to deallocate all the memory space allocated using salloc.
 Should only be used at the end of the program.
 In fact, I would argue that this should not be used at all. 
 You should sfree all the variables allocated with salloc individually.
 Use at your own risk.
*/

void sfree_all( char *fn, int ln )
{

#if USE_SALLOC

    if (nnodes > 0) 
    {
        salloc_nodep_t prev;


        if (debug) 
        { 
            printf( "sfree_all: (@ %s:%d)\n", fn, ln );
            printf( "     total of %10d bytes (%d pointers)\n",
                    (int) totalmem, (int) nnodes );
        }

        while (lastnode) 
        {
            smp_free( lastnode->ptr );
            totalmem -= lastnode->size;

            if (debug) 
                printf( "              %10d b = %10d x %10d b (%s* %s @ %s:%d)\n", 
                        (int) lastnode->size, (int) lastnode->nobj, (int) lastnode->size1, 
                        lastnode->type, lastnode->vn, lastnode->fn, lastnode->ln );
      
            prev = lastnode->prev;

            smp_free( lastnode ); 
            nnodes--;

            lastnode = prev;

        }


        if ( totalmem != 0 ) 
            printf( "sfree_all: totalmem = %d b != 0 (@ %s:%d)\n", 
                    (int) totalmem, fn, ln );

        if ( nnodes != 0 ) 
            printf( "sfree_all: nnodes = %d != 0 (@ %s:%d)\n", 
                    (int) nnodes, fn, ln );
    }

#endif
}







/* 
 This function outputs a report of all the memory space that is currently allocated 
 to the file "<tag>.alloc.log".
 Very useful at the end of the code to figure out what memory you forgot to free 
 or if you have a memory leak.
*/


void salloc_report( char *fn, int ln, int ipe, char *tag )
{
#if USE_SALLOC
    
    salloc_nodep_t node;
    FILE *repf;
    char repfn[200];
    int jpe, npe;

    MPI_Comm_size(MPI_COMM_WORLD, &npe);


    sprintf( repfn, "%s.alloc.log", tag );

    for ( jpe = 0; jpe < npe; jpe ++ ) 
    {
        if (jpe == ipe) 
        {
            if (ipe == 0) 
                my_fopen( repf, repfn, "w" );
            else 
                my_fopen( repf, repfn, "a" );
            
  
            fprintf( repf, "\n" );
            fprintf( repf, "\n" );
            fprintf( repf, "salloc_report: (@ %s:%d, pe=%3d)\n", fn, ln, ipe );
      
            for ( node = lastnode; node; node = node->prev ) 
            {
                fprintf( repf, "pe=%3d %p %10d b = %10d x %10d b (%s* %s @ %s:%d)\n",
                         ipe, node->ptr, 
                         (int) node->size, (int) node->nobj, (int) node->size1,
                         node->type, node->vn, node->fn, node->ln);
            }
            fprintf( repf, "\n" );
            fprintf( repf, "pe=%3d total of %10d b not deallocated\n", ipe, (int) totalmem );
            fprintf( repf, "pe=%3d total of %10d   blocks not deallocated\n", ipe, (int) nnodes );
            fprintf( repf, "\n" );
            fprintf( repf, "pe=%3d total of %10d b peak memory\n", ipe, (int) maxtotalmem );
            fprintf( repf, "pe=%3d total of %10d   peak allocated blocks\n", ipe, (int) maxnnodes );
            fprintf( repf, "\n" );
            fprintf( repf, "pe=%3d total of %10.4g s allocating memory\n", ipe, salloc_time );
            fprintf( repf, "pe=%3d total of %10d   memory allocations\n", ipe, (int) salloc_ncalls );
            fprintf( repf, "\n" );
            fprintf( repf, "pe=%3d total of %10.4g s deallocating memory\n", ipe, sfree_time );
            fprintf( repf, "pe=%3d total of %10d   memory deallocations\n", ipe, (int) sfree_ncalls );
            fprintf( repf, "\n" );
            fprintf( repf, "\n" );

            fclose( repf );
        }
        my_barrier();
    }
#endif
}




