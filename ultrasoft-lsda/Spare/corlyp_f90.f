c************************** SVN Revision Information **************************
c **    $Id$    **
c******************************************************************************
 
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
      subroutine corlyp_f90(dp,dm,dp1,dm1,dp2,dm2,ec,vcp0,vcm0,ndm)

c
      implicit none
      logical  tpot
      integer ndm
      real*8  dp,dm,dp2,dm2,ec,vcp0,vcm0
      real*8  aa,bb,cc,dd,c1,c2,c3,c4,c5,c6,c7,c8,t13,t53,t43,t83,t89,c9
      real*8  zero,d0,dxsq,d2,d0xt13,d0xt53,dpt53,dmt53,z,sc,h,ga,gb
      real*8  gafp,gafm,scf,sc2,chf,hf,hff,h2,zfp,zfm,yz,yz2,z2,ya,yafp
      real*8  yafm,yb,ybfp,ybfm,yb2,d1sq,d1_x_z1
      real*8, dimension(1:ndm) :: dp1,dm1
      real*8, dimension(size(dp1)) :: d1,yy1,yz1,z1,yb1

      parameter(aa=0.04918d0,bb=0.132d0,cc=0.2533d0,dd=0.349d0,
     &         c1=-4*aa,c2=dd,c3=2*bb,c4=cc,c5=4.55779986d0,
     &         c6=1.d0/72.d0,c7=1.d0/18.d0,c8=0.125d0,
     &         t13=1.d0/3.d0,t53=5*t13,t43=4*t13,t83=2*t43,
     &         t89=8.d0/9.d0,c9=t43+t89,zero=0.d0)
c
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

c polarization factor
      z= c1*(dp*dm)*dxsq

c scaling function
      sc= 1.d0/(1.d0+ c2*d0xt13)
      if(d0 .lt. 1.d-09) then
         h = 0.d0
      else
          h= c3*d0xt53*exp(-c4*d0xt13)
      endif

c kinetic energy density expansion
      ga= c5*(dp*dpt53+ dm*dmt53)

      gb= c6*(sum(dp1*dp1)-dp*dp2+ sum(dm1*dm1)-dm*dm2) 
     &   +c7*(dp*dp2+ dm*dm2)
     &   +c8*(d0*d2- d1sq)

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
  
c collect contributions
        vcp0= yafp+ ybfp*(ga+gb)
     &      + yb*(gafp+2*c8*(c9*dp2+2*dm2))
     &      + 2*c8*(c9*sum(dp1*yb1)+2*sum(dm1*yb1))
     &      + yb2*c8*(t43*dp+dm)

        vcm0= yafm+ ybfm*(ga+gb)
     &      + yb*(gafm+2*c8*(c9*dm2+2*dp2))
     &      + 2*c8*(c9*sum(dm1*yb1)+2*sum(dp1*yb1))
     &      + yb2*c8*(t43*dm+dp)
  
      else

        vcp0 = zero
        vcm0 = zero

      endif
  
c correlation energy per electron
      ec= z*sc*(d0+ h*(ga+gb))/d0

      return
      end
