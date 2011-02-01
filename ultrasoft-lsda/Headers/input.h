/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#ifndef NODELIB_H_INCLUDED
#define NODELIB_H_INCLUDED

#define require( getdata ) if( ! getdata ) error_handler( "No input and no default in "#getdata)

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "const.h"


typedef enum {
	INIT = 1,
	TAGS = 2,
	INFO = 4,
	END = 8,
	ITEM = 16, /* member in a LIFO stack of items */
	LIST = 32, 
	INT = 64, /* a fundamental type, currently an int */
	BOOL = 128, /* a fundamental type of true/false (stdbool.h) */
	DBL = 256, /* a fundamental type, currently a double */
	STR = 512, /* a fundamental type, currently a char* */
	OPT = TAGS + STR, /*alias for a validated string*/
    LINE = ITEM + STR, /*alias for a line item string*/
	LINES = LIST + STR, /* alias for new-line delimited list */
	VEC = LIST + DBL, /*alias for a 3 item list of doubles*/
} flags_t;

typedef struct item_t {
	flags_t as;
	union {
		struct node_t *list;
		bool boolean;
		int integer;
		double rational;
		char *string;
	} the;
	struct item_t *next;
} item_t;

typedef struct node_t {
	char *name;
	item_t *is; /* Content *is* data or index pointer */
	struct node_t *last, *next;
} node_t;

/* Needed in verify as pointer, see Input/input.c or Input/verify.c */
extern const bool SET;
extern const bool USET;

extern node_t Root;
#define this Root.is->the.list

/* Node/Item functionality with default behavior for "this" */
item_t	*newItem( flags_t type, void *value );
bool	 pushItem( item_t *new, node_t *here );
#define	 push( DATA ) pushItem( DATA, this )
bool	 appendItem( item_t *new, node_t *here );
#define	 append( DATA ) appendItem( DATA, this )
item_t	*popItem( node_t *here );
#define	 pop() popItem ( this )
int		 itemize( char delimiter );
int		 countItems( void );
void	 fcatItem( FILE *stream, item_t *item );
#define  catItem( ITEM ) fcatItem ( ct.logfile , ITEM )
#define  DcatItem( ITEM ) DEBUG ? fcatItem ( stderr, ITEM ), fflush(NULL): 0

node_t	*newNode( char *name, item_t *item );
node_t	*linkNode( node_t *new );
node_t	*unlinkNode( void );
bool	 killNode( node_t *node );
#define  kill() killNode( unlinkNode() )
bool	 findNode( char *name );
void	 fcatNode( FILE *stream, node_t *here );
#define	 catNode( NODE ) fcatNode( ct.logfile, NODE );
void	 catNodes( void );
int		 tagsload( void );			
int		 tagstrip( void );
bool	 validate( char *optlist );
bool	 verify( char *tagname, const void *optvalue );
bool	 get_data(char *meta, void *dest, flags_t flags, char *data);
bool	 cast( flags_t flags );
bool	 assign( flags_t flaga, void *dest );

/* Input processing that is node/item independent */
char	*filetostr(char *fname);
void	 read_atoms ( void );
void	 read_pdb ( void );
void	 pdb_clean( void );

#endif /* NODELIB_H_INCLUDED */
