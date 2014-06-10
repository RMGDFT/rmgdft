/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"


void read_init(char *meta)
{
    int tmp, is;
    char *tbuf, *tptr, *dstr;
    float ftmp;
    int status,size, num_image;
    struct stat buffer;

    my_malloc (tptr, MAX_PATH, char);
    if((status = stat(meta, &buffer)) == -1)
    {
        dprintf("\n using default path and input"); 
        strncpy(pct.image_path[0], "./", MAX_PATH);
        strncpy(pct.image_input[0], "input", MAX_PATH);
        pct.images = 1;
        pct.image_npes[0] = 1;
        ct.images_per_node = 1;
        ct.spin_flag = 0;
        pct.pe_kpoint = 1;
        return;
    }

    if (pct.worldrank == 0)
    {
        newNode (meta, newItem (STR, filetostr (meta)));
        size = tagstrip ();
        dstr = this->is->the.string;
    }
    MPI_Bcast (&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (pct.worldrank != 0)
    {
        my_malloc (dstr, size, char);
    }
    MPI_Bcast (dstr, size, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (pct.worldrank != 0)
    {
        newNode (meta, newItem (STR, dstr));
    }
    size = tagsload();

    /* do spin polarized calculation? */
    get_data ("spin_polarization", &ct.spin_flag, BOOL, "false");
    dprintf("\n ct.spin_flag  %d", ct.spin_flag); 

    get_data ("num_images", &pct.images, INT, "1");
    get_data ("image_per_node", &ct.images_per_node, INT, "1");
    get_data ("pegroup_for_kpoint", &pct.pe_kpoint, INT, "1");

    dprintf("\n pct.images %d", pct.images);
    num_image = 0;

    get_data ("image_infos", &tmp, INIT | LIST, NULL);
    tbuf = tptr;
    int tot_pe = 0;
    while (get_data ("image_infos", tbuf, ITEM | STR, NULL))
    {
        if(sscanf (tbuf, "%s %s %d",  
                    pct.image_path[num_image],
                    pct.image_input[num_image], 
                    &pct.image_npes[num_image]) != 3)
        {
            printf("\n imgae info wrong");
            exit(0);
        }

        tot_pe += pct.image_npes[num_image];
        num_image++;
    }

    if(tot_pe != pct.total_npes) 
    {
        dprintf("\n require %d npes != %d from job", tot_pe, pct.total_npes);
        exit(0);
    }
    if(num_image != pct.images) 
    {
        dprintf("\n number of image and image info not consistent");
        exit(0);
    }
    if(pct.images > MAX_IMGS)
    {
        dprintf("\n numbef of image %d > %d MAX_IMGS", pct.images, MAX_IMGS);
        exit(0);
    }

    if(pct.images % ct.images_per_node != 0)
    {
        dprintf("you are crazy to use the multiple image per node feature"); 
        exit(0);
    }
    if(ct.images_per_node >1 )
    {
        for(tmp = 1; tmp < pct.images; tmp++) 
            if(pct.image_npes[tmp] != pct.image_npes[0])
            {
                dprintf("\n image %d has different NPES %d", pct.image_npes[tmp], pct.image_npes[0]);
                exit(0);
            }
    }
}




