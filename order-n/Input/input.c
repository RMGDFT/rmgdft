/************************** SVN Revision Information **************************
 *  **    $Id: input.c 1041 2009-04-07 01:19:24Z froze $    **
 *  ******************************************************************************/


#include "input.h"

const bool SET = true;
const bool UNSET = false;

item_t RootIndex = { LIST, {NULL}, NULL };
node_t Root = { ":", &RootIndex, &Root, &Root };
