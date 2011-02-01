/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/pack_doubles.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (c) 1999. Emil Briggs
 *     This code is free software. You are free to modify it or distribute
 *     it in any way you see fit.
 * FUNCTION
 *   int pack_doubles_int(double *inb, unsigned *outb, int count, int precision)
 *   void unpack_doubles_int(double *outb, unsigned *inb, int count, int precision)
 *
 * Description:
 *   Functions to convert a buffer of 8 byte double precision floating
 *   point numbers to and from a packed 5 or 6 byte format.
 *   
 * Usage:
 *   When doing parallel numerical calculations it is sometimes the 
 *   case that single precision arithmetic does not provide enough
 *   accuracy but double precision arithmetic provides much more
 *   than is needed. If the calculations are run on a cluster 
 *   architecture then the extra bandwidth required for sending
 *   a full 8 bytes per floating point number can be expensive.
 *   These routines convert a buffer of 64 bit double precision
 *   floating point numbers in IEEE format to a packed 5 or 6 byte
 *   format. The conversion causes some loss of precision but can
 *   improve execution speed considerably.
 *
 *   A normal IEEE double precision number consists of
 *     mantissa		53 bits
 *     exponent		11 bits
 *   The packed 5 byte format uses
 *     mantissa         29 bits
 *     exponent         11 bits
 *   The packed 6 byte format uses
 *     mantissa         37 bits
 *     exponent         11 bits
 *
 *   
 * A typical calling sequence would be
 *
 *    bidx = pack_doublesc(nmat1, (unsigned *)nmat3, idx, 6);
 *    MPI_Sendrecv(nmat3, bidx, MPI_BYTE, nb_ids[NB_U], 1, nmat4, bidx,
 *                 MPI_BYTE, nb_ids[NB_D], 1, pct.thisgrp_comm, &mstatus);
 *    unpack_doublesc(nmat2, (unsigned *)nmat4, idx, 6);                           
 *
 *
 * Please email any bug reports to:
 *   briggs@tick.physics.ncsu.edu
 *
 *
 * Function prototypes (C calling convention)
 *   pack_doublesc(double *inb, unsigned *outb, int count, int precision)
 *
 *      Arguments:
 *        double *inb      input array of IEEEE format DP real numbers
 *        unsigned *outb   output array of 32 bit unsigned integers
 *        int count        number of elements in the input array
 *        int precision    integer indicating whether to pack into a
 *                         5 or 6 byte format. If precision != 5 then
 *                         6 is assumed.
 *
 *      Returns the total length in bytes of the elements packed into the
 *      output array.
 *
 *
 *  unpack_doublesc(double *outb, unsigned *inb, int *count, int *precision)
 *
 * 	Arguments:
 * 	  double *outb	   output array of IEEE format DP real numbers
 * 	  unsigned inb     input array of packed DP numbers
 * 	  int count        number of elements in the input array
 * 	  int precision    integer indicating whether to pack into a 
 * 	                   5 or 6 byte format. If precision != 5 then
 * 	                   6 is assumed.
 *
 *
 * Fortran interfaces are also provided where all arguments are pass
 * by reference.
 *
 *  int pack_doublesf(double *inb, unsigned *outb, int *count, int *precision);
 * void unpack_doublesf(double *outb, unsigned *inb, int *count, int *precision);
 * SOURCE
 */


/* External interfaces for C calling convention */
int pack_doublesc (double *inb, unsigned *outb, int count, int precision);
void unpack_doublesc (double *outb, unsigned *inb, int count, int precision);



/* External interfaces for Fortran calling conventions. Underscores are
 * required on some OS but not on others.
 */
int pack_doublesf (double *inb, unsigned *outb, int *count, int *precision);
void unpack_doublesf (double *outb, unsigned *inb, int *count, int *precision);
int pack_doublesf_ (double *inb, unsigned *outb, int *count, int *precision);
void unpack_doublesf_ (double *outb, unsigned *inb, int *count, int *precision);
int pack_doublesf__ (double *inb, unsigned *outb, int *count, int *precision);
void unpack_doublesf__ (double *outb, unsigned *inb, int *count, int *precision);





/* Internal function prototypes */
int pack_doubles_int (double *inb, unsigned *outb, int count, int precision);
void unpack_doubles_int (double *outb, unsigned *inb, int count, int precision);


#define		BIG_ENDIAN	1
#define		LITTLE_ENDIAN	2






int pack_doubles_int (double *inb, unsigned *outb, int count, int precision)
{

    int idx, stop, icount, type = BIG_ENDIAN;
    unsigned u1, *u1ptr, *u2ptr;
    double *tptr1, *tptr2;
    unsigned char *b;

    idx = 1;
    b = (unsigned char *) &idx;
    if (b[0])
        type = LITTLE_ENDIAN;

    u1ptr = (unsigned *) inb;
    u2ptr = (unsigned *) outb;

    switch (type)
    {

    case BIG_ENDIAN:

        if (precision == 5)
        {

            stop = count / 4;
            icount = stop * 20;
            for (idx = 0; idx < stop; idx++)
            {

                u2ptr[0] = u1ptr[0];
                u1 = u1ptr[1] >> 24;
                u2ptr[1] = u1ptr[2];
                u1 += ((u1ptr[3] >> 16) & 0xff00);
                u2ptr[2] = u1ptr[4];
                u1 += ((u1ptr[5] >> 8) & 0xff0000);
                u2ptr[3] = u1ptr[6];
                u1 += (u1ptr[7] & 0xff000000);
                u2ptr[4] = u1;

                u1ptr += 8;
                u2ptr += 5;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 4; idx++)
            {

                tptr2[idx] = tptr1[idx];
                icount += 8;

            }                   /* end for */

        }
        else
        {

            stop = count / 2;
            icount = stop * 12;
            for (idx = 0; idx < stop; idx++)
            {

                u2ptr[0] = u1ptr[0];
                u1 = u1ptr[1] >> 16;
                u2ptr[1] = u1ptr[2];
                u1 += (u1ptr[3] & 0xffff0000);
                u2ptr[2] = u1;

                u1ptr += 4;
                u2ptr += 3;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 2; idx++)
            {

                tptr2[idx] = tptr1[idx];
                icount += 8;

            }                   /* end for */

        }                       /* end if */

        break;

    default:                   /* Let little-endian be the default */

        if (precision == 5)
        {

            stop = count / 4;
            icount = stop * 20;
            for (idx = 0; idx < stop; idx++)
            {

                u2ptr[0] = u1ptr[1];
                u1 = u1ptr[0] & 0xff000000;
                u2ptr[1] = u1ptr[3];
                u1 += ((u1ptr[2] & 0xff000000) >> 8);
                u2ptr[2] = u1ptr[5];
                u1 += ((u1ptr[4] & 0xff000000) >> 16);
                u2ptr[3] = u1ptr[7];
                u1 += ((u1ptr[6] & 0xff000000) >> 24);
                u2ptr[4] = u1;

                u1ptr += 8;
                u2ptr += 5;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 4; idx++)
            {

                tptr2[idx] = tptr1[idx];
                icount += 8;

            }                   /* end for */

        }
        else
        {

            stop = count / 2;
            icount = stop * 12;
            for (idx = 0; idx < stop; idx++)
            {

                u2ptr[0] = u1ptr[1];
                u1 = u1ptr[0] & 0xffff0000;
                u2ptr[1] = u1ptr[3];
                u1 += ((u1ptr[2] & 0xffff0000) >> 16);
                u2ptr[2] = u1;

                u1ptr += 4;
                u2ptr += 3;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 2; idx++)
            {

                tptr2[idx] = tptr1[idx];
                icount += 8;

            }                   /* end for */

        }                       /* end if */

    }                           /* end switch */


    return icount;

}                               /* pack_doubles */


void unpack_doubles_int (double *outb, unsigned *inb, int count, int precision)
{

    int idx, stop, type = BIG_ENDIAN;
    unsigned u1, *u1ptr, *u2ptr;
    double *tptr1, *tptr2;
    unsigned char *b;


    u1 = 1;
    b = (unsigned char *) &u1;
    if (b[0])
        type = LITTLE_ENDIAN;


    u1ptr = (unsigned *) outb;
    u2ptr = (unsigned *) inb;

    switch (type)
    {

    case BIG_ENDIAN:

        if (precision == 5)
        {

            stop = count / 4;
            for (idx = 0; idx < stop; idx++)
            {

                u1 = u2ptr[4];

                u1ptr[0] = u2ptr[0];
                u1ptr[1] = (u1 << 24);
                u1ptr[2] = u2ptr[1];
                u1ptr[3] = ((u1 << 16) & 0xff000000);
                u1ptr[4] = u2ptr[2];
                u1ptr[5] = ((u1 << 8) & 0xff000000);
                u1ptr[6] = u2ptr[3];
                u1ptr[7] = (u1 & 0xff000000);

                u1ptr += 8;
                u2ptr += 5;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 4; idx++)
            {

                tptr1[idx] = tptr2[idx];

            }                   /* end for */

        }
        else
        {

            stop = count / 2;
            for (idx = 0; idx < stop; idx++)
            {

                u1 = u2ptr[2];
                u1ptr[0] = u2ptr[0];
                u1ptr[1] = (u1 << 16);
                u1ptr[2] = u2ptr[1];
                u1ptr[3] = (u1 & 0xffff0000);

                u1ptr += 4;
                u2ptr += 3;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 2; idx++)
            {

                tptr1[idx] = tptr2[idx];

            }                   /* end for */

        }                       /* end if */

        break;

    default:                   /* Let little-endian be the default */

        if (precision == 5)
        {

            stop = count / 4;
            for (idx = 0; idx < stop; idx++)
            {

                u1 = u2ptr[4];

                u1ptr[0] = u1 & 0xff000000;
                u1ptr[1] = u2ptr[0];
                u1ptr[2] = (u1 << 8) & 0xff000000;
                u1ptr[3] = u2ptr[1];
                u1ptr[4] = (u1 << 16) & 0xff000000;
                u1ptr[5] = u2ptr[2];
                u1ptr[6] = (u1 << 24) & 0xff000000;
                u1ptr[7] = u2ptr[3];

                u1ptr += 8;
                u2ptr += 5;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 4; idx++)
            {

                tptr1[idx] = tptr2[idx];

            }                   /* end for */

        }
        else
        {

            stop = count / 2;
            for (idx = 0; idx < stop; idx++)
            {

                u1 = u2ptr[2];
                u1ptr[0] = u1 & 0xffff0000;
                u1ptr[1] = u2ptr[0];
                u1ptr[2] = (u1 << 16) & 0xffff0000;
                u1ptr[3] = u2ptr[1];

                u1ptr += 4;
                u2ptr += 3;

            }                   /* end for */

            /* Cleanup leftovers */
            tptr1 = (double *) u1ptr;
            tptr2 = (double *) u2ptr;
            for (idx = 0; idx < count % 2; idx++)
            {

                tptr1[idx] = tptr2[idx];

            }                   /* end for */

        }                       /* end if */

    }                           /* end switch */


}                               /* end unpack_doubles */



/* Fortran wrappers for different platforms */
int pack_doublesf (double *inb, unsigned *outb, int *count, int *precision)
{

    return pack_doubles_int (inb, outb, *count, *precision);

}                               /* end pack_doublesf_ */

int pack_doublesf_ (double *inb, unsigned *outb, int *count, int *precision)
{

    return pack_doubles_int (inb, outb, *count, *precision);

}                               /* end pack_doublesf */

int pack_doublesf__ (double *inb, unsigned *outb, int *count, int *precision)
{

    return pack_doubles_int (inb, outb, *count, *precision);

}                               /* end pack_doublesf */

/* Fortran wrapper */
void unpack_doublesf (double *outb, unsigned *inb, int *count, int *precision)
{

    unpack_doubles_int (outb, inb, *count, *precision);

}                               /* end pack_doublesf_ */

/* Fortran wrapper */
void unpack_doublesf_ (double *outb, unsigned *inb, int *count, int *precision)
{

    unpack_doubles_int (outb, inb, *count, *precision);

}                               /* end pack_doublesf */

/* Fortran wrapper */
void unpack_doublesf__ (double *outb, unsigned *inb, int *count, int *precision)
{

    unpack_doubles_int (outb, inb, *count, *precision);

}                               /* end pack_doublesf */

/* C wrapper */
int pack_doublesc (double *inb, unsigned *outb, int count, int precision)
{

    return pack_doubles_int (inb, outb, count, precision);

}                               /* end pack_doublesc */


/* C wrapper */
void unpack_doublesc (double *outb, unsigned *inb, int count, int precision)
{

    unpack_doubles_int (outb, inb, count, precision);

}                               /* end pack_doublesc */

/******/
