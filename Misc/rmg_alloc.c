

#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "rmg_error.h"
#include "fftw3.h"


static void rmg_alloc_initialize (void *ptr, size_t n, char *type);

void *rmg_malloc(int n, size_t size )
{
    void *ptr;
    if(NULL == (ptr = malloc(n * size))) {
           printf("\n memory size required %d x %lu", n, size);
           rmg_error_handler("can't allocate memory malloc ");
    }
    return ptr;
}

void *rmg_malloc_init(int n, size_t size, char *type )
{
    void *ptr;
    if(NULL == (ptr = malloc(n * size))) {
            printf("\n memory size required %d x %lu", n, size);
           rmg_error_handler("can't allocate memory malloc_init " );
    }
    rmg_alloc_initialize (ptr, n, type);
    return ptr;
}

void *rmg_calloc(int n, size_t size )
{
    void *ptr;
    if(NULL == (ptr = calloc(n, size))) {
            printf("\n memory size required %d x %lu", n, size);
           rmg_error_handler("can't allocate memory calloc ");
    }
    return ptr;
}

void rmg_free( void *ptr )
{

    if(ptr != NULL)
    {
         free( ptr );
    }
    else
    {
           rmg_error_handler("don't try to free null pointer");
    }

}


/* 
 The following routine initializes the allocated data to 0.0 
 for selected types.
*/

static void rmg_alloc_initialize (void *ptr, size_t n, char *type)
{
    size_t i;
    if (strcmp (type, "int") == 0)
    {
        int *p = (int *) ptr;
        for (i = 0; i < n; i++)
            p[i] = 0;
    }
    else if (strcmp (type, "double") == 0)
    {
        double *p = (double *) ptr;
        for (i = 0; i < n; i++)
            p[i] = 0.0;
    }
    else if (strcmp (type, "double") == 0)
    {
        double *p = (double *) ptr;
        for (i = 0; i < n; i ++ )
            p[i] = 0.0;
    }
    else if (strcmp (type, "fftw_complex") == 0)
    {
        fftw_complex *p = (fftw_complex *) ptr;
        for (i = 0; i < n; i++) {
            double *p1 = (double *)p;
            p1[0] = p1[1] = 0.0;
            p++;
        }
    }
    else if ( strcmp( type, "complex double" ) == 0 )
    {
        complex double *p = (complex double *) ptr;
        for ( i = 0; i < n; i ++ )
            p[i] = 0.0;
    }
    else
    {
        printf ("!!!! warning: requested initialization of data pointed to by '%p' not done.\n", ptr);
        printf ("              please inplement initialization procedure for type '%s' in '%s:%d'.\n",
             type, __FILE__, __LINE__);

    }

}

