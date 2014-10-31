/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

struct BOX {

	int 	x1, x2, y1, y2, z1, z2; 
}; 
typedef	struct 	BOX 	BOX; 

struct COMPASS {

	int 	type; 
	BOX  	box1, box2, box3, box4;  
}; 
typedef	struct	COMPASS	COMPASS;


extern COMPASS 	chargeDensityCompass; 
extern COMPASS 	potentialCompass; 
