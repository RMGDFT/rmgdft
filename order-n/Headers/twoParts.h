/************************** SVN Revision Information **************************
 **    $Id: twoParts.h 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
struct BOX
{

    int x1, x2, y1, y2, z1, z2;
};
typedef struct BOX BOX;

struct COMPASS
{

    int type;
    BOX box1, box2, box3, box4;
};
typedef struct COMPASS COMPASS;


COMPASS chargeDensityCompass;
COMPASS potentialCompass;
