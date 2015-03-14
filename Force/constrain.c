/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/***** Common-MGDFT/constrain.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *                      Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void constraint(ION *)
 *   Enforces per atom constraints as set in the ct.ions[ion].constraint 3-vector.
 *   Atomic forces are projected onto the constraint vector in one of three ways
 *   based on the ct.ions[ion].contraint_type integer.
 *   constraint_type == 0::do nothing,
 *   constraint_type == 1::restrict to plane as defined by the constraint vector as normal,
 *   constraint_type == 2::restrict to direction as defined by the constraint vector.
 * INPUTS
 *   pointer to the ion that needs to have its constraint enforced.
 * OUTPUT
 *   no output
 * PARENTS
 *   rmg_fastrlx.c
 * CHILDREN
 *   
 * SOURCE
 */

#define X 0
#define Y 1
#define Z 2


#include "const.h"
#include "rmgtypes.h"
#include "rmgtypedefs.h"
#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include <math.h>

void constrain (void)
{
    int ion;
    ION *iptr;
    printf("Entering constrained forces for image %d", pct.thisimg+1);
    double *Tau, *Img_L, *Img_R;
    my_malloc(Tau, 3*ct.num_ions, double);
    my_malloc(Img_L, 3*ct.num_ions, double);
    my_malloc(Img_R, 3*ct.num_ions, double);

    switch (ct.constrainforces)
    {
        case 5:   /* NEB tangent to higher energy adjacent image with climbing/descending extrema */
            {
                double Mag_T = 0.0;
                double Mag_L = 0.0;
                double Mag_R = 0.0;
                double FdotT = 0.0;

                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = get_ion(ion);

                    /*Calculate displacement vectors from self to left and right image coords */
                    Img_L[3*ion+X] = iptr->crds[X] - iptr->constraint.setA_coord[X];
                    Img_L[3*ion+Y] = iptr->crds[Y] - iptr->constraint.setA_coord[Y];
                    Img_L[3*ion+Z] = iptr->crds[Z] - iptr->constraint.setA_coord[Z];

                    Img_R[3*ion+X] = iptr->constraint.setB_coord[X] - iptr->crds[X];
                    Img_R[3*ion+Y] = iptr->constraint.setB_coord[Y] - iptr->crds[Y];
                    Img_R[3*ion+Z] = iptr->constraint.setB_coord[Z] - iptr->crds[Z];

                    /*Calculate vector magnitudes (squared) */
                    Mag_L += Img_L[3*ion+X]*Img_L[3*ion+X]
                        + Img_L[3*ion+Y]*Img_L[3*ion+Y]
                        + Img_L[3*ion+Z]*Img_L[3*ion+Z];

                    Mag_R += Img_R[3*ion+X]*Img_R[3*ion+X]
                        + Img_R[3*ion+Y]*Img_R[3*ion+Y]
                        + Img_R[3*ion+Z]*Img_R[3*ion+Z];

                }

                if( ( iptr->constraint.setA_weight >= ct.TOTAL && ct.TOTAL <= iptr->constraint.setB_weight) || \
                        ( iptr->constraint.setA_weight < ct.TOTAL && ct.TOTAL > iptr->constraint.setB_weight) )
                {   /* this image energy is an extrema along the band */
                    /* Calculate tangent vector Tau */
                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] = Img_L[3*ion+X]/Mag_L + Img_R[3*ion+X]/Mag_R;
                        Tau[3*ion+Y] = Img_L[3*ion+Y]/Mag_L + Img_R[3*ion+Y]/Mag_R;
                        Tau[3*ion+Z] = Img_L[3*ion+Z]/Mag_L + Img_R[3*ion+Z]/Mag_R;

                        Mag_T += Tau[3*ion+X]*Tau[3*ion+X]
                            + Tau[3*ion+Y]*Tau[3*ion+Y]
                            + Tau[3*ion+Z]*Tau[3*ion+Z];
                    }

                }
                else if ( iptr->constraint.setA_weight > iptr->constraint.setB_weight)
                {   /* this image energy is in a decreasing to the right section */
                    /* Calculate tangent vector Tau */
                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] = Img_L[3*ion+X];
                        Tau[3*ion+Y] = Img_L[3*ion+Y];
                        Tau[3*ion+Z] = Img_L[3*ion+Z];
                    }

                    Mag_T = Mag_L;
                }
                else if ( iptr->constraint.setA_weight < iptr->constraint.setB_weight)
                {   /* this image energy is in a decreasing to the left section */
                    /* Calculate tangent vector Tau */
                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] = Img_R[3*ion+X];
                        Tau[3*ion+Y] = Img_R[3*ion+Y];
                        Tau[3*ion+Z] = Img_R[3*ion+Z];
                    }

                    Mag_T = Mag_R;
                }


                /* check tau tangent vector size */
                if ( Mag_T == 0 ) {
                    error_handler("Image collision this(%d) image.", pct.thisimg);
                } else {
                    Mag_T = sqrt(Mag_T);
                }

                /* Normalize tangent vector Tau */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    Tau[3*ion+X] /= Mag_T;
                    Tau[3*ion+Y] /= Mag_T;
                    Tau[3*ion+Z] /= Mag_T;
                }

                /* Find the amount of force along the Tau vector */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    FdotT += iptr->force[ct.fpt[0]][X] * Tau[3*ion+X]
                        + iptr->force[ct.fpt[0]][Y] * Tau[3*ion+Y]
                        + iptr->force[ct.fpt[0]][Z] * Tau[3*ion+Z];
                }

                /* determine image motion along band */
                if  ( iptr->constraint.setA_weight < ct.TOTAL && ct.TOTAL > iptr->constraint.setB_weight)
                {   /* this image energy is a maxima - climb image */
                    Mag_T = -FdotT;
                }
                else if ( iptr->constraint.setA_weight > ct.TOTAL && ct.TOTAL < iptr->constraint.setB_weight)
                { /* this image energy is a minima - descend image */
                    Mag_T = 2*FdotT;
                }
                else
                { /* not extrema - keep evenly spaced (RMSD) images */
                    /*Calculate vector norms */
                    Mag_L =  sqrt(Mag_L) ;
                    Mag_R =  sqrt(Mag_R) ;
                    Mag_T = ct.neb_spring_constant*(Mag_R - Mag_L);
                }

                /* Remove physical force along Tau, replace it with the restoring force */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr->constraint.forcemask[X] = (Mag_T - FdotT) * Tau[3*ion+X];
                    iptr->constraint.forcemask[Y] = (Mag_T - FdotT) * Tau[3*ion+Y];
                    iptr->constraint.forcemask[Z] = (Mag_T - FdotT) * Tau[3*ion+Z];
                }
            }
            break;
        case 4:                    /* NEB tangent to higher energy adjacent image */
            {
                double Mag_T = 0.0;
                double Mag_L = 0.0;
                double Mag_R = 0.0;
                double FdotT = 0.0;

                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = &ct.ions[ion];

                    /*Calculate displacement vectors from self to left and right image coords */
                    Img_L[3*ion+X] = iptr->crds[X] - iptr->constraint.setA_coord[X];
                    Img_L[3*ion+Y] = iptr->crds[Y] - iptr->constraint.setA_coord[Y];
                    Img_L[3*ion+Z] = iptr->crds[Z] - iptr->constraint.setA_coord[Z];

                    Img_R[3*ion+X] = iptr->constraint.setB_coord[X] - iptr->crds[X];
                    Img_R[3*ion+Y] = iptr->constraint.setB_coord[Y] - iptr->crds[Y];
                    Img_R[3*ion+Z] = iptr->constraint.setB_coord[Z] - iptr->crds[Z];

                    /*Calculate vector magnitudes (squared) */
                    Mag_L += Img_L[3*ion+X]*Img_L[3*ion+X]
                        + Img_L[3*ion+Y]*Img_L[3*ion+Y]
                        + Img_L[3*ion+Z]*Img_L[3*ion+Z];

                    Mag_R += Img_R[3*ion+X]*Img_R[3*ion+X]
                        + Img_R[3*ion+Y]*Img_R[3*ion+Y]
                        + Img_R[3*ion+Z]*Img_R[3*ion+Z];

                }
                /*Calculate vector norms, protect against zeros */
                Mag_L = Mag_L > 0.0 ? sqrt(Mag_L) : error_handler("Left image(%d) collision in NEB", pct.thisimg);
                Mag_R = Mag_R > 0.0 ? sqrt(Mag_R) : error_handler("Right image(%d) collision in NEB", pct.thisimg);

                if( ( iptr->constraint.setA_weight >= ct.TOTAL && ct.TOTAL <= iptr->constraint.setB_weight) || \
                        ( iptr->constraint.setA_weight <= ct.TOTAL && ct.TOTAL >= iptr->constraint.setB_weight) )
                {   /* this image energy is an extrema along the band */
                    /* Calculate tangent vector Tau */
                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] = Img_L[3*ion+X]/Mag_L + Img_R[3*ion+X]/Mag_R;
                        Tau[3*ion+Y] = Img_L[3*ion+Y]/Mag_L + Img_R[3*ion+Y]/Mag_R;
                        Tau[3*ion+Z] = Img_L[3*ion+Z]/Mag_L + Img_R[3*ion+Z]/Mag_R;

                        Mag_T += Tau[3*ion+X]*Tau[3*ion+X]
                            + Tau[3*ion+Y]*Tau[3*ion+Y]
                            + Tau[3*ion+Z]*Tau[3*ion+Z];
                    }

                    Mag_T = Mag_T > 0.0 ? sqrt(Mag_T) : error_handler("Left/Right image collision in NEB", pct.thisimg);
                }
                else if ( iptr->constraint.setA_weight > ct.TOTAL && ct.TOTAL > iptr->constraint.setB_weight)
                {   /* this image energy is in a decreasing to the right section */
                    /* Calculate tangent vector Tau */
                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] = Img_L[3*ion+X];
                        Tau[3*ion+Y] = Img_L[3*ion+Y];
                        Tau[3*ion+Z] = Img_L[3*ion+Z];
                    }

                    Mag_T = Mag_L;
                }
                else if ( iptr->constraint.setA_weight < ct.TOTAL && ct.TOTAL < iptr->constraint.setB_weight)
                {   /* this image energy is in a decreasing to the left section */
                    /* Calculate tangent vector Tau */
                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] = Img_R[3*ion+X];
                        Tau[3*ion+Y] = Img_R[3*ion+Y];
                        Tau[3*ion+Z] = Img_R[3*ion+Z];
                    }

                    Mag_T = Mag_R;
                }


                /* Normalize tangent vector Tau */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    Tau[3*ion+X] /= Mag_T;
                    Tau[3*ion+Y] /= Mag_T;
                    Tau[3*ion+Z] /= Mag_T;
                }

                /* Find the amount of force along the Tau vector */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    FdotT += iptr->force[ct.fpt[0]][X] * Tau[3*ion+X]
                        + iptr->force[ct.fpt[0]][Y] * Tau[3*ion+Y]
                        + iptr->force[ct.fpt[0]][Z] * Tau[3*ion+Z];
                }

                /* Calculate the necessary restoring force to keep images from instersecting */
                Mag_T = ct.neb_spring_constant*(Mag_R - Mag_L);

                /* Remove physical force along Tau, replace it with the restoring force */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr->constraint.forcemask[X] = (Mag_T - FdotT) * Tau[3*ion+X];
                    iptr->constraint.forcemask[Y] = (Mag_T - FdotT) * Tau[3*ion+Y];
                    iptr->constraint.forcemask[Z] = (Mag_T - FdotT) * Tau[3*ion+Z];
                }
            }
            break;

        case 3:                    /* NEB tangent to normalized adjacent images */
            {
                double Mag_T = 0.0;
                double Mag_L = 0.0;
                double Mag_R = 0.0;
                double FdotT = 0.0;

                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = &ct.ions[ion];

                    /*Calculate displacement vectors from self to left and right image coords */
                    Img_L[3*ion+X] = iptr->crds[X] - iptr->constraint.setA_coord[X];
                    Img_L[3*ion+Y] = iptr->crds[Y] - iptr->constraint.setA_coord[Y];
                    Img_L[3*ion+Z] = iptr->crds[Z] - iptr->constraint.setA_coord[Z];

                    Img_R[3*ion+X] = iptr->constraint.setB_coord[X] - iptr->crds[X];
                    Img_R[3*ion+Y] = iptr->constraint.setB_coord[Y] - iptr->crds[Y];
                    Img_R[3*ion+Z] = iptr->constraint.setB_coord[Z] - iptr->crds[Z];

                    /*Calculate vector magnitudes (squared) */
                    Mag_L += Img_L[3*ion+X]*Img_L[3*ion+X]
                        + Img_L[3*ion+Y]*Img_L[3*ion+Y]
                        + Img_L[3*ion+Z]*Img_L[3*ion+Z];

                    Mag_R += Img_R[3*ion+X]*Img_R[3*ion+X]
                        + Img_R[3*ion+Y]*Img_R[3*ion+Y]
                        + Img_R[3*ion+Z]*Img_R[3*ion+Z];

                }
                /*Calculate vector norms, protect against zeros */
                Mag_L = Mag_L > 0.0 ? sqrt(Mag_L) : 1.0;
                Mag_R = Mag_R > 0.0 ? sqrt(Mag_R) : 1.0;

                /* Calculate tangent vector Tau */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    Tau[3*ion+X] = Img_L[3*ion+X]/Mag_L + Img_R[3*ion+X]/Mag_R;
                    Tau[3*ion+Y] = Img_L[3*ion+Y]/Mag_L + Img_R[3*ion+Y]/Mag_R;
                    Tau[3*ion+Z] = Img_L[3*ion+Z]/Mag_L + Img_R[3*ion+Z]/Mag_R;

                    Mag_T += Tau[3*ion+X]*Tau[3*ion+X]
                        + Tau[3*ion+Y]*Tau[3*ion+Y]
                        + Tau[3*ion+Z]*Tau[3*ion+Z];

                }

                /* check tau tangent vector size */
                if ( Mag_T == 0 ) {
                    error_handler("Image collision in both left and right images.");
                } else {
                    Mag_T = sqrt(Mag_T);
                }

                /* Normalize tangent vector Tau */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    Tau[3*ion+X] /= Mag_T;
                    Tau[3*ion+Y] /= Mag_T;
                    Tau[3*ion+Z] /= Mag_T;
                }

                /* Find the amount of force along the Tau vector */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    FdotT += iptr->force[ct.fpt[0]][X] * Tau[3*ion+X]
                        + iptr->force[ct.fpt[0]][Y] * Tau[3*ion+Y]
                        + iptr->force[ct.fpt[0]][Z] * Tau[3*ion+Z];
                }

                /* Calculate the necessary restoring force to keep images from instersecting */
                Mag_T = ct.neb_spring_constant*(Mag_R - Mag_L);

                /* Remove physical force along Tau, replace it with the restoring force */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr->force[ct.fpt[0]][X] += (Mag_T - FdotT) * Tau[3*ion+X];
                    iptr->force[ct.fpt[0]][Y] += (Mag_T - FdotT) * Tau[3*ion+Y];
                    iptr->force[ct.fpt[0]][Z] += (Mag_T - FdotT) * Tau[3*ion+Z];
                }
            }
            break;

        case 2:                    /* NEB tangent to adjacent images */
            {
                double Mag_T = 0.0;
                double Mag_L = 0.0;
                double Mag_R = 0.0;
                double FdotT = 0.0;

                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = &ct.ions[ion];

                    /*Calculate displacement vectors from self to left and right image coords */
                    Img_L[3*ion+X] = iptr->crds[X] - iptr->constraint.setA_coord[X];
                    Img_L[3*ion+Y] = iptr->crds[Y] - iptr->constraint.setA_coord[Y];
                    Img_L[3*ion+Z] = iptr->crds[Z] - iptr->constraint.setA_coord[Z];

                    Img_R[3*ion+X] = iptr->constraint.setB_coord[X] - iptr->crds[X];
                    Img_R[3*ion+Y] = iptr->constraint.setB_coord[Y] - iptr->crds[Y];
                    Img_R[3*ion+Z] = iptr->constraint.setB_coord[Z] - iptr->crds[Z];

                    /*Calculate displacement vector between left and right image coords */
                    Tau[3*ion+X] = Img_L[3*ion+X] + Img_R[3*ion+X];
                    Tau[3*ion+Y] = Img_L[3*ion+Y] + Img_R[3*ion+Y];
                    Tau[3*ion+Z] = Img_L[3*ion+Z] + Img_R[3*ion+Z];

                    /*Calculate vector magnitudes (squared) */
                    Mag_L += Img_L[3*ion+X]*Img_L[3*ion+X]
                        + Img_L[3*ion+Y]*Img_L[3*ion+Y]
                        + Img_L[3*ion+Z]*Img_L[3*ion+Z];

                    Mag_R += Img_R[3*ion+X]*Img_R[3*ion+X]
                        + Img_R[3*ion+Y]*Img_R[3*ion+Y]
                        + Img_R[3*ion+Z]*Img_R[3*ion+Z];

                    Mag_T += Tau[3*ion+X]*Tau[3*ion+X]
                        + Tau[3*ion+Y]*Tau[3*ion+Y]
                        + Tau[3*ion+Z]*Tau[3*ion+Z];
                }
                /*Calculate vector magnitudes */
                Mag_L = sqrt(Mag_L);
                Mag_R = sqrt(Mag_R);

                /* Normalize tangent vector Tau */
                if (Mag_T != 0.0)
                {
                    Mag_T = sqrt(Mag_T);
                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] /= Mag_T;
                        Tau[3*ion+Y] /= Mag_T;
                        Tau[3*ion+Z] /= Mag_T;
                    }
                }

                /* Find the amount of force along the Tau vector */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    FdotT += iptr->force[ct.fpt[0]][X] * Tau[3*ion+X]
                        + iptr->force[ct.fpt[0]][Y] * Tau[3*ion+Y]
                        + iptr->force[ct.fpt[0]][Z] * Tau[3*ion+Z];
                }

                /* Calculate the necessary restoring force to keep images from instersecting */
                Mag_T = ct.neb_spring_constant*(Mag_R - Mag_L);

                /* Remove physical force along Tau, replace it with the restoring force */
                printf("Applying NEB force constraints\n");

                for (ion=0; ion < ct.num_ions; ion++)
                {
                    //    printf("ION %3d F_orig: %g, %g, %g\n",ion, iptr->force[ct.fpt[0]][X], iptr->force[ct.fpt[0]][Y] ,iptr->force[ct.fpt[0]][Z]);
                    //    printf("ION %3d FdotT : %g, %g, %g\n",ion, FdotT*iptr->force[ct.fpt[0]][X], FdotT*iptr->force[ct.fpt[0]][Y] ,FdotT*iptr->force[ct.fpt[0]][Z]);
                    //    printf("ION %3d F_rest: %g, %g, %g\n\n",ion, Mag_T*iptr->force[ct.fpt[0]][X], Mag_T*iptr->force[ct.fpt[0]][Y] ,Mag_T*iptr->force[ct.fpt[0]][Z]);


                    iptr->constraint.forcemask[X] = (Mag_T - FdotT) * Tau[3*ion+X];
                    iptr->constraint.forcemask[Y] = (Mag_T - FdotT) * Tau[3*ion+Y];
                    iptr->constraint.forcemask[Z] = (Mag_T - FdotT) * Tau[3*ion+Z];
                }
            }
            break;

        case 1:                    /* In plane, 2-D restriction (typical). */
            {
                double FdotC;
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = &ct.ions[ion];

                    FdotC = iptr->force[ct.fpt[0]][X] * iptr->constraint.setA_coord[X]
                        + iptr->force[ct.fpt[0]][Y] * iptr->constraint.setA_coord[Y]
                        + iptr->force[ct.fpt[0]][Z] * iptr->constraint.setA_coord[Z];

                    iptr->force[ct.fpt[0]][X] -= FdotC * iptr->constraint.setA_coord[X];
                    iptr->force[ct.fpt[0]][Y] -= FdotC * iptr->constraint.setA_coord[Y];
                    iptr->force[ct.fpt[0]][Z] -= FdotC * iptr->constraint.setA_coord[Z];
                }
            }
            break;
        case 0:                    /* Constraint disabled, should not get here */
            break;
    }

    my_free(Img_R);
    my_free(Img_L);
    my_free(Tau);
    return;
}


