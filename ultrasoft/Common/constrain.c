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

#include <float.h>
#include <math.h>
#include "main.h"

void constrain (void)
{

	switch (ct.constrainforces)
	{
		case 3:                    /* NEB tangent to normalized adjacent images */
			{
				int ion;
				ION *iptr;
				REAL Mag_T = 0.0;
				REAL Mag_L = 0.0;
				REAL Mag_R = 0.0;
				REAL FdotT = 0.0;
				REAL Tau[3*ct.num_ions];
				REAL Img_L[3*ct.num_ions];
				REAL Img_R[3*ct.num_ions];

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
				int ion;
				ION *iptr;
				REAL Mag_T = 0.0;
				REAL Mag_L = 0.0;
				REAL Mag_R = 0.0;
				REAL FdotT = 0.0;
				REAL Tau[3*ct.num_ions];
				REAL Img_L[3*ct.num_ions];
				REAL Img_R[3*ct.num_ions];

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
				for (ion=0; ion < ct.num_ions; ion++)
				{
					iptr->force[ct.fpt[0]][X] += (Mag_T - FdotT) * Tau[3*ion+X];
					iptr->force[ct.fpt[0]][Y] += (Mag_T - FdotT) * Tau[3*ion+Y];
					iptr->force[ct.fpt[0]][Z] += (Mag_T - FdotT) * Tau[3*ion+Z];
				}
			}
			break;
		case 1:                    /* In plane, 2-D restriction (typical). */
			{
				int ion;
				REAL FdotC;
				ION *iptr;
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

	return;
}

#if 0	

	switch (iptr->constraint.setA_type)
	{
		case 4:    /* NEB (simple) tangent to adjacent images */
			{   
				REAL MagT, MagL, MagR, FdotT, Tau[3], PosL[3], PosR[3];

				PosL[0] = iptr->crds[0] - iptr->constraint.setA_coord[0];
				PosL[1] = iptr->crds[1] - iptr->constraint.setA_coord[1];
				PosL[2] = iptr->crds[2] - iptr->constraint.setA_coord[2];

				PosR[0] = iptr->constraint.setB_coord[0] - iptr->crds[0];
				PosR[1] = iptr->constraint.setB_coord[1] - iptr->crds[1];
				PosR[2] = iptr->constraint.setB_coord[2] - iptr->crds[2];

				MagL = sqrt(PosL[0]*PosL[0] + PosL[1]*PosL[1] + PosL[2]*PosL[2]);
				MagR = sqrt(PosR[0]*PosR[0] + PosR[1]*PosR[1] + PosR[2]*PosR[2]);

				Tau[0] = PosL[0] + PosR[0];
				Tau[1] = PosL[1] + PosR[1];
				Tau[2] = PosL[2] + PosR[2];

				MagT = sqrt(Tau[0]*Tau[0] + Tau[1]*Tau[1] + Tau[2]*Tau[2]);

				if (MagT != 0.0)
				{
					Tau[0] /= MagT;
					Tau[1] /= MagT;
					Tau[2] /= MagT;
				}

				FdotT = iptr->force[fpt][0] * Tau[0]
					+ iptr->force[fpt][1] * Tau[1]
					+ iptr->force[fpt][2] * Tau[2];

				MagT = ct.neb_spring_constant*(MagL - MagR);

				iptr->force[fpt][0] += (MagT - FdotT) * Tau[0];
				iptr->force[fpt][1] += (MagT - FdotT) * Tau[1];
				iptr->force[fpt][2] += (MagT - FdotT) * Tau[2];

			}
			break;
		case 3:    /* NEB (enhanced) tangent to adjacent images */
			{   
				REAL MagT, MagL, MagR, FdotT, Tau[3], PosL[3], PosR[3];

				PosL[0] = iptr->crds[0] - iptr->constraint.setA_coord[0];
				PosL[1] = iptr->crds[1] - iptr->constraint.setA_coord[1];
				PosL[2] = iptr->crds[2] - iptr->constraint.setA_coord[2];

				PosR[0] = iptr->constraint.setB_coord[0] - iptr->crds[0];
				PosR[1] = iptr->constraint.setB_coord[1] - iptr->crds[1];
				PosR[2] = iptr->constraint.setB_coord[2] - iptr->crds[2];

				MagL = sqrt(PosL[0]*PosL[0] + PosL[1]*PosL[1] + PosL[2]*PosL[2]);
				MagR = sqrt(PosR[0]*PosR[0] + PosR[1]*PosR[1] + PosR[2]*PosR[2]);

				if (MagL == 0.0) MagL = 1.0;
				if (MagR == 0.0) MagR = 1.0;

				Tau[0] = PosL[0]/MagL + PosR[0]/MagR;
				Tau[1] = PosL[1]/MagL + PosR[1]/MagR;
				Tau[2] = PosL[2]/MagL + PosR[2]/MagR;

				MagT = sqrt(Tau[0]*Tau[0] + Tau[1]*Tau[1] + Tau[2]*Tau[2]);

				if (MagT != 0.0)
				{
					Tau[0] /= MagT;
					Tau[1] /= MagT;
					Tau[2] /= MagT;
				}

				FdotT = iptr->force[fpt][0] * Tau[0]
					+ iptr->force[fpt][1] * Tau[1]
					+ iptr->force[fpt][2] * Tau[2];

				MagT = ct.neb_spring_constant*(MagL - MagR);
				printf("\n\tNat force            is :%f\t%f\t%f", iptr->force[fpt][0], iptr->force[fpt][1], iptr->force[fpt][2]);
				printf("\n\tNEB force correction is :%f\t%f\t%f", (MagT - FdotT) * Tau[0], (MagT - FdotT) * Tau[1], (MagT - FdotT) * Tau[2]);

				iptr->force[fpt][0] += (MagT - FdotT) * Tau[0];
				iptr->force[fpt][1] += (MagT - FdotT) * Tau[1];
				iptr->force[fpt][2] += (MagT - FdotT) * Tau[2];

			}
			break;

		case 2:                    /* NEB tangent to adjacent images */
			{
				REAL Mag, MagL, MagR, FdotT;
				REAL Tau[3], TauL[3], TauR[3];
				REAL Left_E = iptr->constraint.setA_weight;
				REAL This_E = ct.TOTAL;
				REAL Right_E = iptr->constraint.setB_weight;


				/*Calculate displacement vectors to left and right image coords */
				TauL[0] = iptr->crds[0] - iptr->constraint.setA_coord[0];
				TauL[1] = iptr->crds[1] - iptr->constraint.setA_coord[1];
				TauL[2] = iptr->crds[2] - iptr->constraint.setA_coord[2];

				TauR[0] = iptr->constraint.setB_coord[0] - iptr->crds[0];
				TauR[1] = iptr->constraint.setB_coord[1] - iptr->crds[1];
				TauR[2] = iptr->constraint.setB_coord[2] - iptr->crds[2];

				/*Calculate vector magnitudes to left and right image coords */
				MagL = sqrt(TauL[0]*TauL[0] + TauL[1]*TauL[1] + TauL[2]*TauL[2]);
				MagR = sqrt(TauR[0]*TauR[0] + TauR[1]*TauR[1] + TauR[2]*TauR[2]);

				/* Calculate energy weighted tangent vector */
				/* Is this image an extrema ? */
				if ( (This_E > Left_E && This_E > Right_E)||(This_E < Left_E && This_E < Right_E) )	
				{
					REAL dEmax = max( fabs(This_E - Left_E), fabs(This_E - Right_E) );
					REAL dEmin = min( fabs(This_E - Left_E), fabs(This_E - Right_E) );

					if ( Left_E > Right_E )
					{
						Tau[0] = dEmax*TauL[0] + dEmin*TauR[0];
						Tau[1] = dEmax*TauL[1] + dEmin*TauR[1];
						Tau[2] = dEmax*TauL[2] + dEmin*TauR[2];
					}
					else
					{
						Tau[0] = dEmax*TauR[0] + dEmin*TauL[0];
						Tau[1] = dEmax*TauR[1] + dEmin*TauL[1];
						Tau[2] = dEmax*TauR[2] + dEmin*TauL[2];
					}
				}
				else /* Image is in a monotonic region */
				{
					if ( Left_E > Right_E )
					{
						Tau[0] = TauL[0];
						Tau[1] = TauL[1];
						Tau[2] = TauL[2];
					}
					else
					{
						Tau[0] = TauR[0];
						Tau[1] = TauR[1];
						Tau[2] = TauR[2];
					}
				}

				/* Normalize Tau */
				Mag = sqrt(Tau[0]*Tau[0] + Tau[1]*Tau[1] + Tau[2]*Tau[2]);
				if (Mag != 0.0)
				{
					Tau[0] /= Mag;
					Tau[1] /= Mag;
					Tau[2] /= Mag;
				}

				/* Find amount of force along tangent */
				FdotT = iptr->force[fpt][0] * Tau[0] +
					iptr->force[fpt][1] * Tau[1] +
					iptr->force[fpt][2] * Tau[2];

				/* Find restoring force necessary to maintain image separation */
				Mag = -( MagR - MagL )*(ct.neb_spring_constant);

				/* Remove Tau parallel add restoring component of force */
				iptr->force[fpt][0] += (Mag - FdotT) * Tau[0];
				iptr->force[fpt][1] += (Mag - FdotT) * Tau[1];
				iptr->force[fpt][2] += (Mag - FdotT) * Tau[2];

				printf("\tNEB force correction is :%f\t%f\t%f\n", (Mag - FdotT) * Tau[0], (Mag - FdotT) * Tau[1], (Mag - FdotT) * Tau[2]);

			}
			break;
		case 1:                    /* In plane, 2-D restriction (typical). */
			{
				REAL fdotc = iptr->force[fpt][0] * iptr->constraint.setA_coord[0] +
					iptr->force[fpt][1] * iptr->constraint.setA_coord[1] +
					iptr->force[fpt][2] * iptr->constraint.setA_coord[2];

				iptr->force[fpt][0] -= fdotc * iptr->constraint.setA_coord[0];
				iptr->force[fpt][1] -= fdotc * iptr->constraint.setA_coord[1];
				iptr->force[fpt][2] -= fdotc * iptr->constraint.setA_coord[2];
			}
			break;
		case 0:                    /* Disable constraint */
			break;
	}
}
else
{
	/*constrained to position of ion number -(setA_type) */
	;
}
}
#endif

