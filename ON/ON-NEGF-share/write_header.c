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
 *   md.c 
 * CHILDREN
 *   write_pos.c
 * SOURCE
 */


#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"
#include "grid.h"
#include "prototypes_on.h"


char *lattice_type[] = {
    "",
    "CUBIC_PRIMITIVE",
    "CUBIC_FC",
    "CUBIC_BC",
    "HEXAGONAL",
    "TRIGONAL_PRIMITIVE",
    "TETRAGONAL_PRIMITIVE",
    "TETRAGONAL_BC",
    "ORTHORHOMBIC_PRIMITIVE",
    "ORTHORHOMBIC_BASE_CENTRED",
    "ORTHORHOMBIC_BC",
    "ORTHORHOMBIC_FC",
    "MONOCLINIC_PRIMITIVE",
    "MONOCLINIC_BASE_CENTRED",
    "TRICLINIC_PRIMITIVE"
};



/* Writes out header information */
void write_header(void)
{

    int kpt, idx;
    time_t tt;
    double t1;

    time(&tt);


    switch (ct.boundaryflag)
    {

    case PERIODIC:
        printf("\n USING PERIODIC BOUNDARY CONDITIONS\n");
        break;

    case CLUSTER:
        printf("\n USING CLUSTER BOUNDARY CONDITIONS\n");
        break;

    case SURFACE:
        printf("\n USING SURFACE BOUNDARY CONDITIONS\n");
        break;

    }                           /* end switch */

    printf ("    Code Revision %s, Last change on %s", SVN_REV, SVN_REVDATE);
    printf("\n\n    %s", ct.description);

    printf("\n\n    FILES USED");
    printf("\n       Control input file       %s", ct.cfile);
    printf("\n       Data input file          %s", ct.infile);
    printf("\n       Data output file         %s", ct.outfile);


//    printf ("\n#######Control input file parameters#####\n\n");
//    findNode ("description");
//    get_data (NULL, NULL, END, NULL);

    printf("\n\n    Number of steps = %d", ct.max_md_steps);
    printf("\n    Output every %d steps", ct.outcount);
    printf("\n    Checkpoint every %d steps", ct.checkpoint);
    printf("\n    rms criterion %e", ct.thr_rms);


    printf("\n\n    GRID DISCRETIZATION:");
    printf("\n        Hx  = %12.6f  bohr", get_hxgrid() * get_xside());
    printf("\n        Hy  = %12.6f  bohr", get_hygrid() * get_yside());
    printf("\n        Hz  = %12.6f  bohr", get_hzgrid() * get_zside());
    printf("\n        NX  = %d", get_NX_GRID());
    printf("\n        NY  = %d", get_NY_GRID());
    printf("\n        NZ  = %d\n", get_NZ_GRID());
    printf("\n        FG_RATIO  = %d", get_FG_RATIO());

    printf("\n    BRAVAIS LATTICE TYPE IS %s", lattice_type[get_ibrav_type()]);
    printf("\n    Cell volume      = %16.8f", get_vel() * ct.psi_nbasis);
    printf("\n    Grid anisotropy  = %16.8f", get_anisotropy());

    printf("\n\n    PROCESSOR TOPOLOGY:  Total PE's = %d", NPES);
    printf("\n       PE_KPOINT  = %d", pct.pe_kpoint);
    printf("\n       pct.pe_x  = %d", pct.pe_x);
    printf("\n       pct.pe_y  = %d", pct.pe_y);
    printf("\n       pct.pe_z  = %d\n", pct.pe_z);


    printf("\n\n    ENERGY CUTOFF  PARAMETER  %12.6f", ct.cparm);

    /* We compute the equivalent energy cutoff using the density of grid
     * points in the cell with a correction for the grid anisotropy.
     */
    t1 = pow(get_vel(), 0.333333333333);
    t1 = PI / (t1 * get_anisotropy());
    t1 = t1 * t1 / 2.0;
    printf("\n    EQUIVALENT ENERGY CUTOFF  %12.6f Rydbergs", t1);


    printf("\n\n    Basis vectors:");

//    printf("\n       A1  = %15.6f,  %15.6f,  %15.6f", ct.a0[0], ct.a0[1], ct.a0[2]);
//    printf("\n       A2  = %15.6f,  %15.6f,  %15.6f", ct.a1[0], ct.a1[1], ct.a1[2]);
//    printf("\n       A3  = %15.6f,  %15.6f,  %15.6f\n\n", ct.a2[0], ct.a2[1], ct.a2[2]);

#if 0
    /* check minimum imaging assumption */
    printf("\n");
    printf("\n    TEST OF MINIMUM IMAGING ASSUMPTION");
    printf("\n          SPEC-PAIR  ENERGY CONTRIBUTION");
    t1 = (get_xside() < get_yside()) ? get_xside() : get_yside();
    t1 = (t1 < get_zside()) ? t1 : get_zside();
    t1 = 0.5 * t1;
    for (i = 0; i < ct.num_species; i++)
    {
        for (j = i; j < ct.num_species; j++)
        {
            printf("\n            %2d %2d    %13.7e", i, j,
                    ct.sp[i].zvalence * ct.sp[i].zvalence *
                    erfc(t1 / sqrt(ct.sp[i].rc * ct.sp[i].rc + ct.sp[j].rc * ct.sp[i].rc)) / t1);
        }                       /* j */
    }                           /* i */
#endif
    printf("\n");


    printf("\n    Number of states   = %d", ct.num_states);
    printf("\n    Number of species  = %d", ct.num_species);
    printf("\n    Number of ions     = %d", ct.num_ions);
    printf("\n    Density mixing     = %12.6f", ct.mix);
    printf("\n    Projector mixing   = %12.6f", ct.prjmix);
    printf("\n    Comp charge        = %12.6f", ct.crho);
    if (ct.background_charge)
    {
        printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        printf("\n    BACKGROUND CHARGE NOT ZERO!!!       ");
        printf("\n    background charge  = %12.6f", ct.background_charge);
        printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    }                           /* end if */


    switch (ct.xctype)
    {
        case LDA_PZ81:             /* LDA Perdew Zunger 81 */
            printf("\n    XC USING LDA WITH PERDEW-ZUNGER 81");
            break;

        case GGA_BLYP:
            printf("\n    XC USING GGA WITH BLYP");
            break;

        case GGA_XB_CP:            /* GGA X-Becke C-Perdew */
            printf("\n    XC USING GGA WITH X-BECKE AND C-PERDEW");
            break;

        case GGA_XP_CP:            /* GGA X-Perdew C-Perdew */
            printf("\n    XC USING GGA WITH X-PERDEW AND C-PERDEW");
            break;

        case GGA_PBE:
            printf("\n    XC USING GGA WITH PBE");
            break;

        default:
            error_handler("Unknown exchange-correlation functional");

    }                           /* end switch */


    printf("\n    POISSON GLOBAL STEP        = %12.6f", ct.poi_parm.gl_step);
    printf("\n    POISSON PRE                = %d", ct.poi_parm.gl_pre);
    printf("\n    POISSON POST               = %d", ct.poi_parm.gl_pst);
    printf("\n    POISSON MULTIGRID LEVELS   = %d", ct.poi_parm.levels);

    printf("\n    PSI GLOBAL STEP            = %12.6f", ct.eig_parm.gl_step);
    printf("\n    PSI PRE                    = %d", ct.eig_parm.gl_pre);
    printf("\n    PSI POST                   = %d", ct.eig_parm.gl_pst);
    printf("\n    PSI MULTIGRID LEVELS   = %d", ct.eig_parm.levels);

    printf("\n\n    BRILLOUIN ZONE SAMPLING AT K-POINTS WITH WEIGHT:");
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        printf("\n    %12.6f   %12.6f   %12.6f   %12.6f",
                ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);

    }                           /* end for */

    if (ct.forceflag)
    {

        switch (ct.forceflag)
        {

            case MD_FASTRLX:
                printf("\n\n    MOLECULAR DYNAMICS USING FAST RELAX");
                break;


            case BAND_STRUCTURE:
                printf("\n\n    BAND STRUCTURE CALCULATION");
                break;

            case MD_CVE:
                printf("\n\n    MOLECULAR DYNAMICS - CVE");
                break;

            case MD_CVT:
                printf("\n\n    MOLECULAR DYNAMICS - CVT");
                break;

            case MD_CPT:
                printf("\n\n    MOLECULAR DYNAMICS - CPT");
                break;

            case PLOT:
                printf("\n\n    PLOTTING DENSITY IN DX FORM");
                break;

            case PSIPLOT:
                printf("\n\n    PLOTTING PSI^2 IN DX FORM");
                break;
            case TDDFT:
                printf("\n\n    TDDFT");
                break;

            default:
                error_handler("UNKNOWN MOLECULAR DYNAMICS METHOD");

        }                       /* end switch */

        printf("\n\n    TIMESTEP FOR MOLECULAR DYNAMICS = %12.8f", ct.iondt);
        printf("\n\n    SCF COUNT PER MD STEP           = %d", ct.max_scf_steps);


    }
    else
    {

        printf("\n\n    IONS ARE NOT ALLOWED TO MOVE");

    }                           /* end if */



    printf("\n\n\n    ATOMIC SPECIES INFORMATION:");
    printf("\n\n\n      Maximum number of non-local projectors = %d", ct.max_nl);
    for (idx = 0; idx < ct.num_species; idx++)
    {

        SPECIES *sp;
        sp = &ct.sp[idx];

        printf("\n\n      Species %d", idx + 1);
        printf("\n      %s", sp->description);
        printf("\n        ATOMIC NUMBER   = %d", sp->atomic_number);
        printf("\n        ATOMIC MASS     = %12.6f", sp->atomic_mass);
        printf("\n        ZVALENCE        = %12.6f", sp->zvalence);
        printf("\n        CORE RC         = %12.6f", sp->rc);
        printf("\n        L POTENTIALS    = %d", sp->num_potentials);
        printf("\n        LOCAL           = %d", sp->local);
        printf("\n        RWIDTH          = %12.6f", sp->rwidth);
        printf("\n        GWIDTH          = %12.6f", sp->gwidth);





    }                           /* idx */


    /* Write out the ionic postions and displacements */
    write_pos();

    printf("\n        Number of TF ions is %d", ct.num_tfions);

    if (ct.num_tfions > 0)
    {
	write_tfions();
    }





}                               /* end write_header */

/******/
