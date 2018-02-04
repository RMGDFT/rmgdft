#include <vector>
#include <cstdlib>
#include <cmath>
#include "const.h"
#include "params.h"
#include "Atomic.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "transition.h"
#include "RmgSumAll.h"
#include "GlobalSums.h"
#include "boundary_conditions.h"


/*we can consider all ions (expensive) or just those nearby (may not be always sufficient)
 * Using local ions only, error is expected to be very small */

using namespace std;

void Vdd(double * rho)
{

    double hxgrid, hygrid, hzgrid, vel, *atomic_rho, distmin, del_vdd, check, dist;
    double xc, yc, zc, *xtal_x, *xtal_y, *xtal_z;
    int i, j, FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET, FP0_BASIS, dimx, dimy, dimz, min_ion=0, num_ions;
    int *rel_index, rel_index_size, *rel_index_size_all, *rel_index_recv=NULL, recv_buff_size, offset;
    int ion_index, count, ion_multiply, *xtal_index, ix, iy, iz, ionct;
    int loopx_max, loopx_min, loopy_max, loopy_min, loopz_max, loopz_min;
    ION *iptr;
    double *loc_array_recv=NULL, x[3], cartesian[3];
    MPI_Request *recv_request=NULL, send_request1=NULL, send_request2=NULL;
    double a00, a01, a02;
    double a10, a11, a12;
    double a20, a21, a22;

    //double timex = my_crtc ();

    hxgrid = get_hxxgrid(); //The fine grid spacing in x
    hygrid = get_hyygrid(); //The fine grid spacing in y
    hzgrid = get_hzzgrid(); //The fine grid spacing in z

    vel = get_vel_f(); //the volume element of the grid

    dimx = get_FPX0_GRID(); //The number of fine grid points in the x direction for this processor
    dimy = get_FPY0_GRID(); //The number of fine grid points in the y direction for this processor
    dimz = get_FPZ0_GRID(); //The number of fine grid points in the z direction for this processor

    FPX_OFFSET = get_FPX_OFFSET(); //The processor grid offset from pe 0 of this processor in the x direction
    FPY_OFFSET = get_FPY_OFFSET(); //The processor grid offset from pe 0 of this processor in the y direction
    FPZ_OFFSET = get_FPZ_OFFSET(); //The processor grid offset from pe 0 of this processor in the z direction

    FP0_BASIS = get_FP0_BASIS(); //The number of grid points per processor

    /*Lattice parameters to avoid calling to_cartesian repeatedly*/
    a00 = Rmg_L.get_a0(0);
    a01 = Rmg_L.get_a0(1);
    a02 = Rmg_L.get_a0(2);

    a10 = Rmg_L.get_a1(0);
    a11 = Rmg_L.get_a1(1);
    a12 = Rmg_L.get_a1(2);

    a20 = Rmg_L.get_a2(0);
    a21 = Rmg_L.get_a2(1);
    a22 = Rmg_L.get_a2(2);

    /*setup ions, we do it here, so that we avoid extra work when looping over all gridpoints*/
    /*Only loop over local ions, i.e.ions whose local potentials have overlap with current domain*/
    num_ions = pct.num_loc_ions;


    /*For non-cluster (periodic) boundary conditions we have to conside periodic images of ions*/
    ion_multiply = 27;
    loopx_min = -1;
    loopx_max = 1;

    loopy_min = -1;
    loopy_max = 1;

    loopz_min = -1;
    loopz_max = 1;

    if (ct.boundaryflag == CLUSTER)
    {
	ion_multiply = 1;

	loopx_min = 0;
	loopx_max = 0;

	loopy_min = 0;
	loopy_max = 0;

	loopz_min = 0;
	loopz_max = 0;
    }

    if (ct.boundaryflag == SURFACE)
    {
	ion_multiply = 9;

	loopz_min = 0;
	loopz_max = 0;
    }

    /*Crystal coordinates of ions and their indices*/
    xtal_x = new double[ion_multiply * num_ions];
    xtal_y = new double[ion_multiply * num_ions];
    xtal_z = new double[ion_multiply * num_ions];
    xtal_index = new int [ion_multiply * num_ions];

    //rmg_printf("\n Vdd setup took %f seconds\n", my_crtc () - timex);


    //timex = my_crtc ();
    /*Store xtal coordinates of ions (including periodic images) into a local array*/
    count = 0;
    for (ionct = 0; ionct < num_ions; ionct++) //for each ion number
    {
	/*Absolute ionic index*/
	ion_index = pct.loc_ions_list[ionct];

	iptr = &ct.ions[ion_index];

	/*For non-cluster (periodic) boundary conditions we have to conside periodic images of ions*/
	for (ix = loopx_min; ix<= loopx_max; ix++)
	{ 
	    for (iy = loopy_min; iy<= loopy_max; iy++)
	    { 
		for (iz = loopz_min; iz<= loopz_max; iz++)
		{

		    xtal_x[count] = iptr->xtal[0]  + (double) ix;
		    xtal_y[count] = iptr->xtal[1]  + (double) iy;
		    xtal_z[count] = iptr->xtal[2]  + (double) iz;

		    xtal_index[count] = ion_index;

		    count++;
		}
	    }
	}
    }

    //rmg_printf("\n Vdd loop 1 took %f seconds\n", my_crtc () - timex);




    /*Get atomic charge density*/
    atomic_rho = new double[FP0_BASIS];
    LcaoGetAtomicRho(atomic_rho);


    /* Array for storing indices of ions*/
    rel_index = new int[num_ions];
    rel_index_size = 0;


    /* Array for storing partial integration in a voronoi domain*/
    double * loc_array;
    loc_array = new double[num_ions]; 

    for (int i = 0; i < num_ions; i++)
	loc_array[i] = 0.0;


    /* Loop over gridpoints*/
    /* Get the offset of the x coordinate of this processor*/

    xc = FPX_OFFSET * hxgrid;
    //Go through each each grid point by looping through each cartesian direction
    j = 0; //index for use in all three dimensions

    //timex = my_crtc ();
    for (int ix = 0; ix < dimx; ix++)
    {

	yc = FPY_OFFSET * hygrid; //Get the offset of the y coordinate of this processor
	for (int iy = 0; iy < dimy; iy++)
	{

	    zc = FPZ_OFFSET * hzgrid; //Get the offset of the z coordinate of this processor
	    for (int iz = 0; iz < dimz; iz++)
	    {

		del_vdd = atomic_rho[j] - rho[j]; //voronoi density deformation integral at this gridpoint 

		/*Loop over ions (including periodic images, find distances from current grid point
		 * and find which ion is the closest*/
		distmin = 1000000.0;
		for (ionct = 0; ionct < ion_multiply * num_ions; ionct++) //for each ion number
		{
		    /*Vector between current grid point and ion in crystal coordinates*/
		    x[0] = xc - xtal_x[ionct];
		    x[1] = yc - xtal_y[ionct];
		    x[2] = zc - xtal_z[ionct];


		    //to_cartesian (x, cartesian);
		    //Do not call this function, overhead is too expensive this deep
		    //in a 4 time nested loop, calculate things explicitly instead
		    cartesian[0] = x[0] * a00 + x[1] * a10 + x[2] * a20;  
		    cartesian[1] = x[0] * a01 + x[1] * a11 + x[2] * a21;  
		    cartesian[2] = x[0] * a02 + x[1] * a12 + x[2] * a22;  



		    dist  = cartesian[0] * cartesian[0];
		    dist += cartesian[1] * cartesian[1];
		    dist += cartesian[2] * cartesian[2];

		    if (dist < distmin)
		    {
			distmin = dist;
			min_ion = xtal_index[ionct];
		    }


		}


		/*Check if index is already recorded*/
		for (i = 0; i < rel_index_size; i++) 
		{
		    if (rel_index[i] == min_ion)
		    {
			loc_array[i] += del_vdd;
			break;
		    }
		}

		/*If not, record*/
		if (i == rel_index_size)
		{
		    rel_index_size += 1;
		    rel_index[i] = min_ion; 
		    loc_array[i] += del_vdd;
		}


		zc += hzgrid; //advance to the next gridpoint in the z direction
		j++; //incrememnt 3-d index


	    }
	    yc += hygrid; //advance to the next gridpoint in the y direction
	}
	xc += hxgrid; //advance to the next gridpoint in the x direction
    }

    //rmg_printf("\n Vdd loop 2 took %f seconds\n", my_crtc () - timex);


    /*Array to store number of voronoi domains on each processor*/
    if(pct.gridpe == 0)
	rel_index_size_all = new int[pct.grid_npes];
    else
	rel_index_size_all = NULL;


    /*Processor 0 gets sizes of local arrays, i.e. how many atoms does each 
     * processor handle, store it as an array in rel_index_size */
    MPI_Gather(&rel_index_size, 1, MPI_INT, rel_index_size_all, 1, MPI_INT, 0, pct.grid_comm);



    if(pct.gridpe == 0)
    {
	//timex = my_crtc ();

	/*Get total size of the receive buffer (not counting local data)*/
	recv_buff_size = 0;

	for (j=1; j<pct.grid_npes; j++)
	    recv_buff_size += rel_index_size_all[j];

	/*Memory for receive buffers*/
	rel_index_recv = new int[recv_buff_size];
	loc_array_recv = new double[recv_buff_size];

	recv_request = new MPI_Request[2*pct.grid_npes];


	offset = 0;

	for (j=1; j<pct.grid_npes; j++)
	{
	    /*First communication sends indices of ions that a given processor has info about, INT only*/
	    MPI_Irecv (&rel_index_recv[offset], rel_index_size_all[j], MPI_INT, j, 111, pct.grid_comm, &recv_request[j]);
	    
	    /*Second call sends values for each ion*/
	    MPI_Irecv (&loc_array_recv[offset], rel_index_size_all[j], MPI_DOUBLE, j, 111, pct.grid_comm, &recv_request[j + pct.grid_npes]);

	    offset += rel_index_size_all[j];
	}





	for (i=0; i<ct.num_ions; i++)
	    ct.ions[i].partial_charge = 0.0;

	
	/*PE 0 own contribution, dealing with this separately*/
	for (i=0; i<rel_index_size; i++)
	{

	    ion_index = rel_index[i];

	    ct.ions[ion_index].partial_charge +=  loc_array[i];
	}

	/*Wait until communication is finished, these are non-blocking calls*/
	/*Integer communication*/
	MPI_Waitall(pct.grid_npes -1,  &recv_request[1], MPI_STATUSES_IGNORE);
	/*Floating point array communication*/
	MPI_Waitall(pct.grid_npes -1,  &recv_request[pct.grid_npes + 1], MPI_STATUSES_IGNORE);

	/*Add received data*/
	for (j=0; j< recv_buff_size; j++)
	{
	    ion_index = rel_index_recv[j];

	    ct.ions[ion_index].partial_charge += loc_array_recv[j];

	}
	//rmg_printf("\n Vdd p2p communication took %f seconds\n", my_crtc () - timex);





	check = 0.0;
	for (i = 0; i < ct.num_ions; i++)
	{

	    iptr = &ct.ions[i];
	    iptr->partial_charge *= vel;

	    check +=  iptr->partial_charge;
	}

	rmg_printf("\n\nVDD: Summation of partial charges is %e (report if substantially different from 0)", check);


    }

    else
    {

	//timex = my_crtc ();
	MPI_Isend (rel_index, rel_index_size, MPI_INT,    0, 111, pct.grid_comm, &send_request1);
	MPI_Isend (loc_array, rel_index_size, MPI_DOUBLE, 0, 111, pct.grid_comm, &send_request2);

	/*Wait until communication is finished*/
	/*Integer communication*/
	MPI_Waitall(1,  &send_request1, MPI_STATUSES_IGNORE);
	/*Floating point array communication*/
	MPI_Waitall(1,  &send_request2, MPI_STATUSES_IGNORE);

	//printf("\n Vdd communication took %f seconds\n", my_crtc () - timex);

    }


    delete [] atomic_rho;
    delete [] loc_array;
    delete [] rel_index;
    delete [] xtal_x;
    delete [] xtal_y;
    delete [] xtal_z;
    delete [] xtal_index;

    if(pct.gridpe == 0)
    {
	delete [] rel_index_recv;
	delete [] loc_array_recv;
	delete [] recv_request;
	delete [] rel_index_size_all;
    }

}
