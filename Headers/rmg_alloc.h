
#ifndef RMG_ALLOC_H
#define RMG_ALLOC_H 1

void rmg_malloc(void *ptr, int n, size_t size );
void rmg_malloc_init(void *ptr, int n, size_t size, char *type );
void rmg_free( void *ptr );
void my_free( void *ptr );

#define my_free rmg_free
#define my_malloc( _ptr_, _nobj_, _type_ ) \
        rmg_malloc(_ptr_, (int) _nobj_, sizeof(_type_))

#define my_calloc( _ptr_, _nobj_, _type_ ) \
        rmg_calloc(_ptr_, (int) _nobj_, sizeof(_type_))

#define my_malloc_init( _ptr_, _nobj_, _type_ ) \
        rmg_malloc_init(_ptr_, (int) _nobj_, sizeof(_type_), "_type_")

#endif
