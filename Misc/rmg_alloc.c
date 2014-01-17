

#include <stdlib.h>
#include <string.h>
#include "rmgtypes.h"
#include "rmg_error.h"
#include "fftw.h"


static void rmg_alloc_initialize (void *ptr, size_t n, char *type);

void rmg_malloc(void *ptr, int n, size_t size )
{
    if(NULL == (ptr = malloc(n * size))) {
           rmg_error_handler("can't allocate memory");
    }
}

void rmg_malloc_init(void *ptr, int n, size_t size, char *type )
{
    if(NULL == (ptr = malloc(n * size))) {
           rmg_error_handler("can't allocate memory");
    }
    rmg_alloc_initialize (ptr, n, type);

}

void rmg_calloc(void *ptr, int n, size_t size )
{
    if(NULL == (ptr = calloc(n, size))) {
           rmg_error_handler("can't allocate memory");
    }
}

void rmg_free( void *ptr )
{

    free( ptr );

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
    else if (strcmp (type, "rmg_double_t") == 0)
    {
        rmg_double_t *p = (rmg_double_t *) ptr;
        for (i = 0; i < n; i ++ )
            p[i] = 0.0;
    }
    else if (strcmp (type, "fftw_complex") == 0)
    {
        fftw_complex *p = (fftw_complex *) ptr;
        for (i = 0; i < n; i++)
            p[i].re = p[i].im = 0.0;
    }
#if 0
    else if ( strcmp( type, "_Complex double" ) == 0 )
    {
        complex double *p = (complex double *) ptr;
        for ( i = 0; i < n; i ++ )
            p[i] = 0.0;
    }
#endif
    else
    {
        printf ("!!!! warning: requested initialization of data pointed to by '%s' not done.\n", ptr);
    }

}

