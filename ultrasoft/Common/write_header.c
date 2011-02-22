/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/write_header.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void write_header(void)
 *   Writes out header information 
 * INPUTS
 *   nothing
 * OUTPUT
 *   print out header information
 * PARENTS
 *   main.c 
 * CHILDREN
 *   write_pos.c
 * SOURCE
 */


#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"
#include "svnrev.h"

static void init_write_pos (void);


char *lattice_type[] = {
    "",
    "Cubic_primitive",
    "Cubic_FC",
    "Cubic_BC",
    "Hexagonal",
    "Trigonal_primitive",
    "Tetragonal_primitive",
    "Tetragonal_BC",
    "Orthorhombic_primitive",
    "Orthorhombic_base_centred",
    "Orthorhombic_BC",
    "Orthorhombic_FC",
    "Monoclinic_primitive",
    "Monoclinic_base_centred",
    "Triclinic_primitive"
};


/* Writes out header information */
void write_header (void)
{

    int kpt, idx;
    time_t tt;
    REAL t1;


    char *timeptr;
    time (&tt);
    timeptr = ctime (&tt);


    printf ("\n");
    printf (" Real-space finite difference molecular dynamics\n\n");
    switch (ct.boundaryflag)
    {

    case PERIODIC:
        printf (" Using periodic boundary conditions\n");
        break;

    case CLUSTER:
        printf (" Using cluster boundary conditions\n");
        break;

    case SURFACE:
        printf (" Using surface boundary conditions\n");
        break;

    }                           /* end switch */

    printf ("\n");
    /*printf (" Code %s\n\n", QMD_VERSION); */
    printf ("    Code Revision %d, Last change on %s", SVN_REV, SVN_REVDATE);
    printf ("    %s\n", ct.description);
    printf ("    Run started at %s\n", timeptr);

    printf ("\n\n");
    printf ("    Files used\n");
    printf ("       Control input file       %s\n", ct.cfile);
    printf ("       Data input file          %s\n", ct.infile);
    printf ("       Data output file         %s\n", ct.outfile);
    /*  printf("       Pseudopotential file     %s\n", ct.pspfile); */

    printf ("\n#######Control input file parameters#####\n\n");
    findNode ("description");
    get_data (NULL, NULL, END, NULL);
    printf ("\n\n#######End input file parameters#########\n");

    printf ("\n\n");

    printf ("    Maximum number of md steps = %d\n", ct.max_md_steps);
    printf ("    Maximum number of scf steps = %d\n", ct.max_scf_steps);
    printf ("    Output every %d steps\n", ct.outcount);
    printf ("    Checkpoint every %d steps\n", ct.checkpoint);
    printf ("    RMS criterion %14.10f\n", ct.thr_rms);

    printf ("\n");
    if (ct.rmvmovie == 0)
    {
        printf ("    Rotmovie not outputted\n");
    }
    else
    {
        printf ("    Rotmovie outputted every %d steps\n", ct.rmvmovie);
    }
    if (ct.chmovie == 0)
    {
        printf ("    DX charge movie not outputted\n");
    }
    else
    {
        printf ("    DX charge movie outputted every %d steps\n", ct.chmovie);
    }


    printf ("\n");
    if (ct.runflag == 1)
        printf ("    Restart from %s\n", ct.infile);

    else
        printf ("    Initial run\n");



    if (ct.initdiag)
        printf ("    Initial subspace diagonalization\n");



    printf ("\n");
    printf ("    Grid discretization:\n");
    printf ("        Hx  = %12.6f a0\n", ct.hxgrid * ct.xside);
    printf ("        Hy  = %12.6f a0\n", ct.hygrid * ct.yside);
    printf ("        Hz  = %12.6f a0\n", ct.hzgrid * ct.zside);
    printf ("        NX  = %d\n", ct.psi_nxgrid);
    printf ("        NY  = %d\n", ct.psi_nygrid);
    printf ("        NZ  = %d\n", ct.psi_nzgrid);

    printf ("\n");
    printf ("    Bravais lattice type is %s\n", lattice_type[ct.ibrav]);
    printf ("    Cell volume      = %12.6f a0^3\n", ct.vel * ct.psi_nbasis);
    printf ("    Grid anisotropy  = %12.6f\n", ct.anisotropy);

    printf ("\n");
    printf ("    Processor topology:  total PE's = %d\n", (PE_X * PE_Y * PE_Z));
    printf ("       PE_X  = %d\n", PE_X);
    printf ("       PE_Y  = %d\n", PE_Y);
    printf ("       PE_Z  = %d\n", PE_Z);

    printf ("\n");
    printf ("    Fine grid (for charge density):\n");
    printf ("       FG_NX / CG_NX  = %d\n", FG_NX);
    printf ("       FG_NY / CG_NY  = %d\n", FG_NY);
    printf ("       FG_NZ / CG_NZ  = %d\n", FG_NZ);
    printf ("    Interpolation type:\n");
    if (!ct.interp_flag)
        printf ("       Cubic interpolation\n");
    else
        printf ("       B-spline interpolation of %d order, using trade_image%d\n",
                ct.interp_order, ct.interp_trade);


    printf ("    Double grid (for non-local projectors):\n");
    printf ("       nxfgrid / CG_NX  = %d\n", ct.nxfgrid);
    printf ("       nyfgrid / CG_NY  = %d\n", ct.nyfgrid);
    printf ("       nzfgrid / CG_NZ  = %d\n", ct.nzfgrid);


    printf ("\n");
    printf ("    Energy cutoff  parameter  %12.6f  %12.6f  %12.6f\n", ct.cparm, ct.betacparm,
            ct.qcparm);

    /* We compute the equivalent energy cutoff using the density of grid
     * points in the cell with a correction for the grid anisotropy.
     */
    t1 = pow (ct.vel, 0.333333333333);
    t1 = PI / (t1 * ct.anisotropy);
    t1 = t1 * t1 / 2.0;
    printf ("       Equivalent energy cutoff  %12.6f Ry\n", t1);

    printf ("\n");
    printf ("    Basis vectors:\n");
    printf ("       A1  = %15.6f  %15.6f  %15.6f a0\n", ct.a0[0], ct.a0[1], ct.a0[2]);
    printf ("       A2  = %15.6f  %15.6f  %15.6f a0\n", ct.a1[0], ct.a1[1], ct.a1[2]);
    printf ("       A3  = %15.6f  %15.6f  %15.6f a0\n", ct.a2[0], ct.a2[1], ct.a2[2]);

#if 0
    /* check minimum imaging assumption */
    printf ("\n");
    printf ("\n    TEST OF MINIMUM IMAGING ASSUMPTION");
    printf ("\n          SPEC-PAIR  ENERGY CONTRIBUTION");
    t1 = (ct.xside < ct.yside) ? ct.xside : ct.yside;
    t1 = (t1 < ct.zside) ? t1 : ct.zside;
    t1 = 0.5 * t1;
    for (i = 0; i < ct.num_species; i++)
    {
        for (j = i; j < ct.num_species; j++)
        {
            printf ("\n            %2d %2d    %13.7e", i, j,
                    ct.sp[i].zvalence * ct.sp[i].zvalence *
                    erfc (t1 / sqrt (ct.sp[i].rc * ct.sp[i].rc + ct.sp[j].rc * ct.sp[i].rc)) / t1);
        }
    }
#endif


    printf ("\n");

    if (pct.spin_flag)
    {
    	printf ("    This is a spin polarized calculation \n");
    	printf ("    Number of spin up states   = %d\n", ct.num_states);
    	printf ("    Number of spin down states   = %d\n", ct.num_states_oppo);
    }
    else
    {
    	printf ("    This is NOT a spin polarized calculation \n");
    printf ("    Number of states   = %d\n", ct.num_states);
    }	
    printf ("    Number of species  = %d\n", ct.num_species);
    printf ("    Number of ions     = %d\n", ct.num_ions);
    printf ("    Density mixing     = %12.6f\n", ct.mix);
    printf ("    Projector mixing   = %12.6f\n", ct.prjmix);
    printf ("    Comp charge        = %12.6f\n", ct.crho);
    if (ct.background_charge)
    {
        printf ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        printf ("    BACKGROUND CHARGE NOT ZERO!!!\n");
        printf ("    background charge  = %12.6f\n", ct.background_charge);
        printf ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }


    switch (ct.xctype)
    {
    case LDA_PZ81:             /* LDA Perdew Zunger 81 */
        printf ("    XC using LDA with Perdew-Zunger 81\n");
        break;

    case GGA_BLYP:
        printf ("    XC using GGA with BLYP\n");
        break;

    case GGA_XB_CP:            /* GGA X-Becke C-Perdew */
        printf ("    XC using GGA with X-Becke and C-Perdew\n");
        break;

    case GGA_XP_CP:            /* GGA X-Perdew C-Perdew */
        printf ("    XC using GGA with X-Perdew and C-Perdew\n");
        break;

    case GGA_PBE:
        printf ("    XC using GGA with PBE\n");
        break;

    default:
        error_handler ("Unknown exchange-correlation functional");

    }                           /* end switch */


    printf ("\n");
    printf ("    Poisson global step        = %12.6f\n", ct.poi_parm.gl_step);
    printf ("    Poisson pre                = %d\n", ct.poi_parm.gl_pre);
    printf ("    Poisson post               = %d\n", ct.poi_parm.gl_pst);
    printf ("    Poisson multigrid levels   = %d\n", ct.poi_parm.levels);

    printf ("\n");
    printf ("    Psi global step            = %12.6f\n", ct.eig_parm.gl_step);
    printf ("    Psi pre                    = %d\n", ct.eig_parm.gl_pre);
    printf ("    Psi post                   = %d\n", ct.eig_parm.gl_pst);
    printf ("    Psi multigrid levels       = %d\n", ct.eig_parm.levels);

    printf ("\n");
    printf ("    Brillouin Zone sampling at K-points with weight:\n");
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        printf ("    %12.6f   %12.6f   %12.6f   %12.6f\n",
                ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);

    }

    printf ("\n");
    switch (ct.forceflag)
    {
    case MD_QUENCH:
        printf ("    Quenching electrons only. Ions are not allowed to move.\n");
        break;

    case MD_FASTRLX:
        printf ("    Molecular dynamics using fast relax.\n");
        break;

    case BAND_STRUCTURE:
        printf ("    Band structure calculation.\n");
        break;

    case MD_CVE:
        printf ("    Molecular dynamics - CVE\n");
        break;

    case MD_CVT:
        printf ("    Molecular dynamics - CVT\n");
        break;

    case MD_CPT:
        printf ("    Molecular dynamics - CPT\n");
        break;

    case PLOT:
        printf ("    Plotting density in DX form.\n");
        break;

    case PSIPLOT:
        printf ("    Plotting Psi^2 in DX form.\n");
        break;

    case NEB_RELAX:
        printf ("    Molecular dynamics using Nudged Elastic Band.\n");
        break;

    default:
        error_handler ("Unknown molecular dynamics method.");
    }

    /* Forces are updated under normalized constraint field */
    if (verify ("atom_constraints", NULL))
    {
        printf ("    Constrained per atom dynamics vector field.\n");
        for (idx = 0; idx < ct.num_ions; idx++)
        {
            printf ("       % 10f % 10f % 10f % 10f\n",
					ct.ions[idx].constraint.setA_coord[0],
					ct.ions[idx].constraint.setA_coord[1],
					ct.ions[idx].constraint.setA_coord[2],
                    ct.ions[idx].constraint.setA_weight);
        }
    }

    if (ct.forceflag != MD_QUENCH)
    {
        printf ("    Timestep for molecular dynamics = %12.8f\n", ct.iondt);
        printf ("    SCF count per MD step           = %d\n", ct.max_scf_steps);
    }



    printf ("\n");
    printf ("    Atomic Species Information:\n");
    printf ("      Maximum number of non-local projectors = %d\n", ct.max_nl);
    for (idx = 0; idx < ct.num_species; idx++)
    {

        SPECIES *sp;
        sp = &ct.sp[idx];

        printf ("\n");
        printf ("      Species %d '%s'\n", idx + 1, sp->pseudo_symbol);
        printf ("      %s\n", sp->description);
        printf ("        Atomic symbol   = %s\n", sp->atomic_symbol);
        printf ("        Atomic number   = %d\n", sp->atomic_number);
        printf ("        Atomic mass     = %12.6f mu\n", sp->atomic_mass);
        printf ("        Zvalence        = %12.6f\n", sp->zvalence);
        printf ("        Core rc         = %12.6f\n", sp->rc);
        printf ("        L potentials    = %d\n", sp->num_potentials);
        printf ("        L local         = %d\n", sp->local);
        printf ("        lrcut           = %12.6f\n", sp->lrcut);
        printf ("        nlrcut          = %12.6f  %12.6f  %12.6f\n", sp->nlrcut[0], sp->nlrcut[1],
                sp->nlrcut[2]);
        printf ("        rwidth          = %12.6f\n", sp->rwidth);
        printf ("        gwidth          = %12.6f\n", sp->gwidth);


    }                           /* end for(idx = 0;idx < ct.num_species;idx++) */


    /* Write out the ionic positions and displacements */
    init_write_pos ();


    if ((pct.imgpe == 0) && (verify ("pdb_atoms", NULL)))
        write_pdb ();




}                               /* end write_header */


/********************************/









static void init_write_pos (void)
{
    int ion;
    ION *iptr;


    printf ("\n\n\n  INITIAL IONIC POSITIONS AND DISPLACEMENTS:\n");


    printf ("\nSpecies   X           Y           Z           dX          dY          dZ");

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];

        printf ("\n  %d   %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f",
                iptr->species + 1,
                iptr->crds[0], iptr->crds[1], iptr->crds[2],
                iptr->crds[0] - iptr->icrds[0],
                iptr->crds[1] - iptr->icrds[1], iptr->crds[2] - iptr->icrds[2]);

    }                           /* end for */

    printf ("\n");

}                               /* end write_pos */

/******/
