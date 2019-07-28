#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

read_orbital.c

This routine reads the orbitals of all probes and the central part. 
Reading data for the left probe, right probe and the central part is 
(from ON cal) straightforward. However, for the down and up probes
the operation orbitals_NEGF(x,y,z) = orbitals_ON(y,x,z) is needed.

*/



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"


void read_orbital (STATE * states)
{
    int fhand;
    int st, st1, st_new;
    long nbytes;
    char newname[MAX_PATH + 200];
    char msg[200];

    int idx, idx0, idx2, subsystem;
    int ixmin, ixmax, iymin, iymax;
    int incx, incy, ix, iy, iz;
    double *array_tmp;

 
    /* Wait until everybody gets here */
    MPI_Barrier(pct.img_comm);

/*
    printf ("state_begin, state_end %d %d \n", ct.state_begin, ct.state_end);
*/
    

    for (st = ct.state_begin; st < ct.state_end; st++)
    {

        st_new = 0;
        for (subsystem = 0; subsystem < cei.num_subsystem; subsystem++)
        {
            idx0 = cei.subsystem_idx[subsystem];

            for (st1 = lcr[idx0].state_begin; st1 < lcr[idx0].state_end; st1++)
            {
                if(st_new == st)
                {
                    sprintf (newname, "%s%s%d", lcr[idx0].name, ".orbit_", st1);
                    my_open( fhand, newname, O_RDWR, S_IREAD | S_IWRITE );
                    if (fhand < 0)
                    {
                        printf ("\n %s, st1 = %d %d", newname, st1, st);
                        error_handler (" Unable to open file ");
                    }

                    idx = states[st].size * (int) sizeof (double);

                    /* ====================== Reading orbitals ======================= */                    

                    if(idx0 <= 2 | idx0 > cei.num_probe) /* Satisfies left probe, central parta & right probe */
                    {
 
                        nbytes = read (fhand, states[st].psiR, idx);
                        if (nbytes != idx)
                        {
                            printf ("\n read %d is different from %d for state %d", (int) nbytes, idx, st);
                            error_handler ("Unexpected end of file orbit");
                        }

                        nbytes = read (fhand, &ixmin, sizeof (int));
                        nbytes = read (fhand, &ixmax, sizeof (int));
                        nbytes = read (fhand, &iymin, sizeof (int));
                        nbytes = read (fhand, &iymax, sizeof (int));
                        if (nbytes != sizeof (int))
                        {
                            printf ("\n read %d is different from %d for state %d", (int) nbytes,
                                    (int) sizeof (int), st);
                            error_handler ("Unexpected end of file orbit");
                        }

                        states[st].ixmin_old = ixmin;
                        states[st].ixmax_old = ixmax;
                        states[st].iymin_old = iymin;
                        states[st].iymax_old = iymax;
/*                   
                        if(idx0 == 1) 
                        printf (" state_pos1  %d %d %d %d %d %d \n", st, idx0, states[st].ixmin_old, 
                        states[st].ixmax_old, states[st].iymin_old, states[st].iymax_old); 

                        if(idx0 == 0) 
                        printf (" state_pos0  %d %d %d %d %d %d \n", st, idx0, states[st].ixmin_old, 
                        states[st].ixmax_old, states[st].iymin_old, states[st].iymax_old); 
*/                   
                    }
                    else /* Satisfies down and up probes */ 
                    {
/*
                    printf ("reading orbitals for 3rd/4th probe, idx0 =  %d %s \n", idx0, newname);
*/
                        my_malloc_init(array_tmp, idx, double);

                        nbytes = read (fhand, array_tmp, idx);
                        if (nbytes != idx)
                        {
                            printf ("\n read %d is different from %d for state %d", (int) nbytes, idx, st);
                            error_handler ("Unexpected end of file orbit");
                        }

                        incx = states[st].orbit_nz * states[st].orbit_ny;
                        incy = states[st].orbit_nz * states[st].orbit_nx;
                        for(ix = 0; ix < states[st].orbit_nx; ix++)
                        {
                            for(iy = 0; iy < states[st].orbit_ny; iy++)
                            {
                                for(iz = 0; iz < states[st].orbit_nz; iz++)
                                {
                                    idx = iz + iy * states[st].orbit_nz + ix * incx;
                                    idx2= iz + ix * states[st].orbit_nz + iy * incy; /* Check */ 
                                    states[st].psiR[idx] = array_tmp[idx2];
/*
                       if(pct.gridpe ==0) printf (" urgent  %d %d %d %f \n", ix, iy, iz, 
                                  states[st].psiR[idx]); 
                       printf (" urgent  %d %d %d %f \n", ix, iy, iz, states[st].psiR[idx]); 
*/
                                }
                            }
                        }
                        my_free(array_tmp);

                        nbytes = read (fhand, &ixmin, sizeof (int));
                        nbytes = read (fhand, &ixmax, sizeof (int));
                        nbytes = read (fhand, &iymin, sizeof (int));
                        nbytes = read (fhand, &iymax, sizeof (int));
                        if (nbytes != sizeof (int))
                        {
                            printf ("\n read %d is different from %d for state %d", (int) nbytes,
                                    (int) sizeof (int), st);
                            error_handler ("Unexpected end of file orbit");
                        }

                        states[st].ixmin_old = iymin;
                        states[st].ixmax_old = iymax;
                        states[st].iymin_old = ixmin;
                        states[st].iymax_old = ixmax;
/*                        
                        printf (" state_pos3  %d %d %d %d %d %d \n", st, idx0, states[st].ixmin_old, 
                        states[st].ixmax_old, states[st].iymin_old, states[st].iymax_old); 
*/
                    }   /* if statement ends */
                    /* ==================== Reading finished ==================== */                    

#if 	DEBUG
                    sprintf (msg, "State %d sum ", st);
                    print_sum_idx (states[st].size, states[st].psiR, msg);
#endif

                    close (fhand);


                }   /* if statement ends */
                st_new++;
            }   /* st1 loop ends */
        }   /* subsystem loop ends */
    }   /* st loop ends */


    MPI_Barrier(pct.img_comm);

    fflush (NULL);

}                               /* end read_data */
