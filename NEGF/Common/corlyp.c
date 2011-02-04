/************************** SVN Revision Information **************************
 **    $Id: corlyp.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*
 *  * This is a complete rewrite in c of the f90 version of corlyp
 *  * The original code is in comment adjacent to the new code.
 *  * Frisco Rose September 1, 2005
 *  * 
 *
c
c Lee Yang Parr correlation energy functional
c no provisions taken against division by zero
c 
c see e.g.
c C. Lee et al. Phys Rev B 37 785 (1988)
c
c Hartree a.u.
c               
c input
c dp ...... spin up density
c dp1 ..... grad(dp) 
c dp2 ..... laplace(dp)
c dm ...... spin down density
c ...
c ndm ..... 1  input gradients are scalar variables
c           3  input gradients are arrays variables of size 3
c tpot .... T  evaluate correlation energy and potential
c           F  evaluate energy only (a posteriori scheme)
c
c output
c ec ...... correlation energy per electron
c vcp0 .... correlation potential for spin up
c vcm0 .... correlation potential for spin down
c
c Martin Fuchs, FHI 24-01-1996
c
cc    subroutine corlyp_f90(dp,dm,dp1,dm1,dp2,dm2,ec,vcp0,vcm0,ndm,tpot)
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

/*
      subroutine corlyp_f90(dp,dm,dp1,dm1,dp2,dm2,ec,vcp0,vcm0,ndm)
*/
void corlyp (REAL *dp, REAL *dm, REAL *dp1, REAL *dm1, REAL *dp2,
             REAL *dm2, REAL *ec, REAL *vcp0, REAL *vcm0, int *ndm)
{
/*
c
      implicit none
      logical  tpot
      int ndm
      real*8  dp,dm,dp2,dm2,ec,vcp0,vcm0
      real*8  aa,bb,cc,dd,c1,c2,c3,c4,c5,c6,c7,c8,t13,t53,t43,t83,t89,c9
      real*8  zero,d0,dxsq,d2,d0xt13,d0xt53,dpt53,dmt53,z,sc,h,ga,gb
      real*8  gafp,gafm,scf,sc2,chf,hf,hff,h2,zfp,zfm,yz,yz2,z2,ya,yafp
      real*8  yafm,yb,ybfp,ybfm,yb2,d1sq,d1_x_z1
c     real*8, dimension(1:ndm) :: dp1,dm1
c     real*8, dimension(size(dp1)) :: d1,yy1,yz1,z1,yb1
*/
    int i;
    REAL dtmp1, dtmp2, dtmp3;
    int NDM = *ndm;
    REAL DP = *dp, DM = *dm, DP2 = *dp2, DM2 = *dm2, EC = *ec, VCP0 =
        *vcp0, VCM0 = *vcm0;
    REAL *d1, *yy1, *yz1, *z1, *yb1;
    int tpot;
    REAL d0, dxsq, d2, d0xt13, d0xt53, dpt53, dmt53, z, sc, h, ga, gb;
    REAL gafp, gafm, scf, sc2, chf, hf, hff, h2, zfp, zfm, yz, yz2, z2, yafp;
    REAL yafm, yb, ybfp, ybfm, yb2, d1sq, d1_x_z1;

    my_malloc( d1, NDM, REAL );
    my_malloc( yy1, NDM, REAL );
    my_malloc( yz1, NDM, REAL );
    my_malloc( z1, NDM, REAL );
    my_malloc( yb1, NDM, REAL );

/*   parameter(aa=0.04918d0,bb=0.132d0,cc=0.2533d0,dd=0.349d0,
     &         c1=-4*aa,c2=dd,c3=2*bb,c4=cc,c5=4.55779986d0,
     &         c6=1.d0/72.d0,c7=1.d0/18.d0,c8=0.125d0,
     &         t13=1.d0/3.d0,t53=5*t13,t43=4*t13,t83=2*t43,
     &         t89=8.d0/9.d0,c9=t43+t89,zero=0.d0)
c*/
#define	aa  0.049180
#define	bb  0.1320
#define	cc  0.25330
#define	dd  0.3490
#define	c1  (-4*aa)
#define	c2  dd
#define	c3  (2*bb)
#define	c4  cc
#define	c5  4.557799860
#define	c6  (1.0/72.0)
#define	c7  (1.0/18.0)
#define	c8  0.1250
#define	t13  (1.0/3.0)
#define	t53  (5*t13)
#define	t43  (4*t13)
#define	t83  (2*t43)
#define	t89  (8.0/9.0)
#define	c9  (t43+t89)

/*
      tpot = .true.
      d0= dp+ dm
      dxsq= 1.d0/(d0*d0)
      d1= dp1+ dm1
      d1sq=sum(d1*d1)
      d2= dp2+ dm2
      d0xt13= d0**(-t13)
      d0xt53= d0xt13*d0xt13/d0
      dpt53= dp**t53
      dmt53= dm**t53
*/
    tpot = 1;
    d0 = DP + DM;
    dxsq = 1.0 / (d0 * d0);
    d1sq = 0.0;
    for (i = 0; i < NDM; i++)
    {
        d1[i] = dp1[i] + dm1[i];
        d1sq += d1[i] * d1[i];
    }
    d2 = DP2 + DM2;
    d0xt13 = pow (d0, -t13);
    d0xt53 = d0xt13 * d0xt13 / d0;
    dpt53 = pow (DP, t53);
    dmt53 = pow (DM, t53);

/*
c polarization factor
      z= c1*(dp*dm)*dxsq
*/
    z = c1 * (DP * DM) * dxsq;

/*
c scaling function
      sc= 1.d0/(1.d0+ c2*d0xt13)
      if(d0 .lt. 1.d-09) then
         h = 0.d0
      else
          h= c3*d0xt53*exp(-c4*d0xt13)
      endif
*/
    sc = 1.0 / (1.0 + c2 * d0xt13);
    if (d0 < 1.0E-9)
        h = 0.0;
    else
        h = c3 * d0xt53 * exp (-c4 * d0xt13);

/*
c kinetic energy density expansion
      ga= c5*(dp*dpt53+ dm*dmt53)

      gb= c6*(sum(dp1*dp1)-dp*dp2+ sum(dm1*dm1)-dm*dm2) 
     &   +c7*(dp*dp2+ dm*dm2)
     &   +c8*(d0*d2- d1sq)
*/
    ga = c5 * (DP * dpt53 + DM * dmt53);
    dtmp1 = dtmp2 = 0.0;
    for (i = 0; i < NDM; i++)
    {
        dtmp1 += dp1[i] * dp1[i];
        dtmp2 += dm1[i] * dm1[i];
    }
    gb = c6 * (dtmp1 - DP * DP2 + dtmp2 - DM * DM2) + c7 * (DP * DP2 +
                                                            DM * DM2) +
        c8 * (d0 * d2 - d1sq);

/*	
c calculate potential
      if(tpot) then

        gafp= t83*c5*dpt53
        gafm= t83*c5*dmt53

        scf= t13*c2*d0xt13/d0*sc*sc
        sc2= scf*(d2+ 2*(scf/sc- 2*t13/d0)*d1sq)

        chf= t13*(c4*d0xt13 -5)/d0
        hf= chf*h
        hff= h*(chf**2+ t13*(5.d0-4*t13*c4*d0xt13)*dxsq)
        h2= (hf*d2+ hff*d1sq)
  
        zfp= (c1*dm- 2*z*d0)*dxsq
        zfm= (c1*dp- 2*z*d0)*dxsq
        yz= z/c1
        yy1= dp*dm1+dm*dp1
        yz1= (yy1-2*yz*d1*d0)*dxsq
        yz2= (2*yz*d1sq- 2*(sum(yz1*d1) +yz*d2)*d0
     &      -2*sum(d1*yy1)/d0+ (dp*dm2+2*sum(dp1*dm1)+dm*dp2))*dxsq
        z1= c1*yz1
        d1_x_z1=sum(d1*z1)
        z2= c1*yz2
  
        ya= sc*z*d0
        yafp= sc*(d0*zfp+z)+ z*d0*scf
        yafm= sc*(d0*zfm+z)+ z*d0*scf

        yb= sc*z*h
        ybfp= sc*(h*zfp+z*hf)+ z*h*scf
        ybfm= sc*(h*zfm+z*hf)+ z*h*scf
        yb1= sc*(h*z1+z*hf*d1)+ z*h*scf*d1
        yb2= (sc*hf+h*scf)*d1_x_z1+ h*sc*z2
     &     + (sc*d1_x_z1+z*scf*d1sq)*hf+ z*sc*h2
     &     + (z*hf*d1sq+h*d1_x_z1)*scf+ z*h*sc2
*/
    if (tpot)
    {
        gafp = t83 * c5 * dpt53;
        gafm = t83 * c5 * dmt53;

        scf = t13 * c2 * d0xt13 / d0 * sc * sc;
        sc2 = scf * (d2 + 2 * (scf / sc - 2 * t13 / d0) * d1sq);

        chf = t13 * (c4 * d0xt13 - 5) / d0;
        hf = chf * h;
        hff = h * (chf * chf + t13 * (5.0 - 4 * t13 * c4 * d0xt13) * dxsq);
        h2 = (hf * d2 + hff * d1sq);
        zfp = (c1 * DM - 2 * z * d0) * dxsq;
        zfm = (c1 * DP - 2 * z * d0) * dxsq;
        yz = z / c1;

        for (i = 0; i < NDM; i++)
            yy1[i] = DP * dm1[i] + DM * dp1[i];
        for (i = 0; i < NDM; i++)
            yz1[i] = (yy1[i] - 2 * yz * d1[i] * d0) * dxsq;
        dtmp1 = dtmp2 = dtmp3 = 0.0;
        for (i = 0; i < NDM; i++)
        {
            dtmp1 += d1[i] * yz1[i];
            dtmp2 += d1[i] * yy1[i];
            dtmp3 += dp1[i] * dm1[i];
        }
        yz2 =
            (2 * yz * d1sq - 2 * (dtmp1 + yz * d2) * d0 - 2 * dtmp2 / d0 +
             (DP * DM2 + 2 * dtmp3 + DM * DP2)) * dxsq;
        for (i = 0; i < NDM; i++)
            z1[i] = c1 * yz1[i];
        d1_x_z1 = 0.0;
        for (i = 0; i < NDM; i++)
        {
            d1_x_z1 += d1[i] * z1[i];
        }
        z2 = c1 * yz2;

        /*ya = sc*z*d0; */
        yafp = sc * (d0 * zfp + z) + z * d0 * scf;
        yafm = sc * (d0 * zfm + z) + z * d0 * scf;

        yb = sc * z * h;
        ybfp = sc * (h * zfp + z * hf) + z * h * scf;
        ybfm = sc * (h * zfm + z * hf) + z * h * scf;
        for (i = 0; i < NDM; i++)
            yb1[i] = sc * (h * z1[i] + z * hf * d1[i]) + z * h * scf * d1[i];
        yb2 =
            (sc * hf + h * scf) * d1_x_z1 + h * sc * z2 + (sc * d1_x_z1 +
                                                           z * scf * d1sq) *
            hf + z * sc * h2 + (z * hf * d1sq + h * d1_x_z1) * scf +
            z * h * sc2;

/*
c collect contributions
        vcp0= yafp+ ybfp*(ga+gb)
     &      + yb*(gafp+2*c8*(c9*dp2+2*dm2))
     &      + 2*c8*(c9*sum(dp1*yb1)+2*sum(dm1*yb1))
     &      + yb2*c8*(t43*dp+dm)

        vcm0= yafm+ ybfm*(ga+gb)
     &      + yb*(gafm+2*c8*(c9*dm2+2*dp2))
     &      + 2*c8*(c9*sum(dm1*yb1)+2*sum(dp1*yb1))
     &      + yb2*c8*(t43*dm+dp)
*/
        dtmp1 = dtmp2 = 0.0;
        for (i = 0; i < NDM; i++)
        {
            dtmp1 += dp1[i] * yb1[i];
            dtmp2 += dm1[i] * yb1[i];
        }
        VCP0 =
            yafp + ybfp * (ga + gb) + yb * (gafp +
                                            2 * c8 * (c9 * DP2 + 2 * DM2)) +
            2 * c8 * (c9 * dtmp1 + 2 * dtmp2) + yb2 * c8 * (t43 * DP + DM);
        VCM0 =
            yafm + ybfm * (ga + gb) + yb * (gafm +
                                            2 * c8 * (c9 * DM2 + 2 * DP2)) +
            2 * c8 * (c9 * dtmp2 + 2 * dtmp1) + yb2 * c8 * (t43 * DM + DP);
    }
    else
    {
        VCP0 = 0.0;
        VCM0 = 0.0;
    }
/*  
c correlation energy per electron
      ec= z*sc*(d0+ h*(ga+gb))/d0
*/
    EC = z * sc * (d0 + h * (ga + gb)) / d0;

/*
      return
      end
*/
    *dp = DP;
    *dm = DM;
    *dp2 = DP2;
    *dm2 = DM2;
    *ec = EC;
    *vcp0 = VCP0;
    *vcm0 = VCM0;
    my_free (d1);
    my_free (yy1);
    my_free (yz1);
    my_free (z1);
    my_free (yb1);

    return;
}
