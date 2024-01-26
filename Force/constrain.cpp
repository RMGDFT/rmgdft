/***** Common-MGDFT/constrain.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *                      Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void constraint(ION *)
 *   Enforces per atom constraints as set in the Atoms[ion].constraint 3-vector.
 *   Atomic forces are projected onto the constraint vector in one of three ways
 *   based on the Atoms[ion].contraint_type integer.
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


#include <float.h>
#include <math.h>
#include "const.h"
#include "rmgtypedefs.h"
#include "common_prototypes.h"
#include "main.h"

#define X 0
#define Y 1
#define Z 2

void constrain (void)
{

    int ion;
    ION *iptr=NULL;
    if(pct.imgpe==0) fprintf(ct.logfile, "Entering constrained forces for image %d", pct.thisimg+1);

    double *Tau, *Img_L, *Img_R;
    Tau = new double[3*ct.num_ions];
    Img_L = new double[3*ct.num_ions];
    Img_R = new double[3*ct.num_ions];

    switch (ct.constrainforces)
    {
        case 5:   /* NEB tangent to higher energy adjacent image with climbing/descending extrema */
            {
                double Mag_T = 0.0;
                double Mag_L = 0.0;
                double Mag_R = 0.0;
                double FdotT = 0.0;

                double Eleft = Atoms[0].constraint.setA_weight;
                double Eright = Atoms[0].constraint.setB_weight;
                double Eself = ct.TOTAL;
                   
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = &Atoms[ion];

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

                if( ( Eleft >= Eself && Eself <= Eright) || ( Eleft < Eself && Eself > Eright) )
                {   /* this image energy is an extrema along the band */
                    /* Calculate tangent vector Tau */
                    double abs_vl = std::abs(Eleft - Eself);
                    double abs_vr = std::abs(Eright - Eself);
                    double vmax = std::max(abs_vl, abs_vr);
                    double vmin = std::min(abs_vl, abs_vr);
                    double delta_left, delta_right;

                    if (Eleft > Eright) 
                    {
                        delta_left = vmax;
                        delta_right = vmin;
                    }
                    else
                    {
                        delta_left = vmin;
                        delta_right = vmin;
                    }

                    for (ion=0; ion < ct.num_ions; ion++)
                    {
                        Tau[3*ion+X] = Img_L[3*ion+X] * delta_left + Img_R[3*ion+X] * delta_right;
                        Tau[3*ion+Y] = Img_L[3*ion+Y] * delta_left + Img_R[3*ion+Y] * delta_right;
                        Tau[3*ion+Z] = Img_L[3*ion+Z] * delta_left + Img_R[3*ion+Z] * delta_right;

                        Mag_T += Tau[3*ion+X]*Tau[3*ion+X]
                            + Tau[3*ion+Y]*Tau[3*ion+Y]
                            + Tau[3*ion+Z]*Tau[3*ion+Z];
                    }

                }
                else if ( Eleft> Eright)
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
                else if ( Eleft < Eright )
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
                    std::string errmsg = "Image collision this(" + std::to_string(pct.thisimg) + ") image.";
                    rmg_error_handler(__FILE__, __LINE__, errmsg.c_str());
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
                    iptr = &Atoms[ion];
                    FdotT += iptr->force[ct.fpt[0]][X] * Tau[3*ion+X]
                        + iptr->force[ct.fpt[0]][Y] * Tau[3*ion+Y]
                        + iptr->force[ct.fpt[0]][Z] * Tau[3*ion+Z];
                }

                /* determine image motion along band */
                if(  Eleft < Eself && Eself > Eright) 
                {   /* this image energy is a maxima - climb image */
                    Mag_T = -FdotT;
                }
                if(  Eleft > Eself && Eself < Eright) 
                { /* this image energy is a minima - descend image */
                    Mag_T = 2*FdotT;
                }
                else
                { /* not extrema - keep evenly spaced (RMSD) images */
                    /*Calculate vector norms */
                    Mag_L =  sqrt(Mag_L) ;
                    Mag_R =  sqrt(Mag_R) ;
                    Mag_T = ct.neb_spring_constant*(Mag_R -Mag_L);
                }

                /* Remove physical force along Tau, replace it with the restoring force */
                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = &Atoms[ion];
                    iptr->constraint.forcemask[X] = (Mag_T - FdotT) * Tau[3*ion+X];
                    iptr->constraint.forcemask[Y] = (Mag_T - FdotT) * Tau[3*ion+Y];
                    iptr->constraint.forcemask[Z] = (Mag_T - FdotT) * Tau[3*ion+Z];
                    iptr->force[ct.fpt[0]][X] +=iptr->constraint.forcemask[X];
                    iptr->force[ct.fpt[0]][Y] +=iptr->constraint.forcemask[Y]; 
                    iptr->force[ct.fpt[0]][Z] +=iptr->constraint.forcemask[Z];
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
                    iptr = &Atoms[ion];

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
                if(Mag_L > 0.0)
                {
                    Mag_L = sqrt(Mag_L);
                }
                else
                {
                    if(pct.gridpe==0) printf("Left image(%d) collision in NEB", pct.thisimg);
                    rmg_error_handler(__FILE__, __LINE__, "Terminating");
                }
                if(Mag_R > 0.0)
                {
                    Mag_R = sqrt(Mag_R);
                }
                else
                {
                    if(pct.gridpe==0) printf("Right image(%d) collision in NEB", pct.thisimg);
                    rmg_error_handler(__FILE__, __LINE__, "Terminating");
                }

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

                    if(Mag_T > 0.0)
                    {
                        Mag_T = sqrt(Mag_T);
                    }
                    else
                    {
                        rmg_error_handler(__FILE__, __LINE__, "Left/Right image collision in NEB");
                    }
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
                    iptr = &Atoms[ion];

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
                    rmg_error_handler(__FILE__, __LINE__, "Image collision in both left and right images.");
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
                    iptr = &Atoms[ion];

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
                if(pct.imgpe==0) fprintf(ct.logfile, "Applying NEB force constraints\n");

                for (ion=0; ion < ct.num_ions; ion++)
                {
                    iptr = &Atoms[ion];

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
                    iptr = &Atoms[ion];

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

    delete [] Img_R;
    delete [] Img_L;
    delete [] Tau;
    return;

}


