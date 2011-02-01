/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "input.h"

const bool SET = true;
const bool UNSET = false;

item_t RootIndex = { LIST, {NULL}, NULL };
node_t Root = { ":", &RootIndex, &Root, &Root };
