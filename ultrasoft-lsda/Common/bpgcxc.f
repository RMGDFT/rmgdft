c************************** SVN Revision Information **************************
c **    $Id$    **
c******************************************************************************
 
c $Header:$
c**********************************************************************c
c
c  Becke exchange for a spin-unpolarized electronic system
c
c  Gradient-corrected exchange energy based on
c     [A.D. Becke, J.Chem.Phys.96, 2155, 1992].
c  The LSDA energy functional, obtained as E{n+,n-}=(E{2n+}+E{2n-})/2,
c     and the functional derivative formula are given by
c     [J.P. Perdew , PRB 33, 8800, 1986].
c     [J.P. Perdew , PRB 34, 7406, 1986].
c  see also [G.Ortiz ...,PRB 43, 6376 (1991)] eq. (A2)
c
c  Hartree a.u.
c
c  input
c  d            density
c  s            abs(grad d)/(2kf*d)
c  u            (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
c           >>  grad(abs(grad d) has mixed derivatives ! <<
c  v            (laplacian d)/(d*(2*kf)**2)
c
c  output
c  ex           exchange energy per electron
c  vx           exchange potential
c
c author -- M. Fuchs, FHI Berlin, 02-1993
c**********************************************************************
c
      subroutine xbecke(d,s,u,v,ex,vx)
c
      implicit real*8 (a-h,o-z)
      data c / .779555417944150792d1/
      data b / .42d-2/
      data bb/-.451357747124625192d-2/
      data ax/-.738558766382022406d0/
      data thrd,thrd4/.333333333333333333d0,.1333333333333333333d1/
c
c exchange enhancement factor f
      x  = c*s
      y1 = 1.d0/sqrt(1.d0+x*x)
      y0 = log(x+1.d0/y1)
      y2 = -x*y1*y1*y1
      ddi= 1.d0/(1.d0 + 6.d0*b*x*y0)
      dd1= 6.d0*b*(y0+x*y1)
      g  = 1.d0 - 0.5d0*x*dd1*ddi
      fs = -2.d0*bb*c*c*ddi
      g1 = -3.d0*b*(y0+x*(3.d0*y1+x*y2-dd1*dd1*ddi/(6.d0*b)))
      fss= fs*c*(g1 - g*dd1)*ddi
      fs = fs*g
      f  = 1.d0 - bb*x*x*ddi
c
c LDA only
      fac= ax*d**thrd
c
c energy
      ex = fac*f
c
c potential
      vx = fac*(thrd4*f-(u-thrd4*s*s*s)*fss-v*fs)
c
      return
      end
ce
c
c***********************************************************************
c
c  gradient-correction to correlation energy from
c
c  [J.P.Perdew, PRB 33, 8822 (1986) and PRB 34, 7406 (1986)]
c
c  to correlation part of the Becke-Perdew gradient-corrected
c  xc-functional
c
c  input
c  d1,d2 : up/down spindensity
c  dp12..: grad(d1)*grad(d2) {* == vector product}, MUST NOT == 0
c  uu    : (grad d)*grad(abs(grad d)) , d = d1 + d2
c  vv    : laplacian d
c
c  output
c  dec  : correction to correlation energy per electron
c  dvcup :         - "" -            potential maj. spin
c  dvcdn :         - "" -            potential min. spin
c
c  Hartree a.u.
c
c  Martin Fuchs, FHI 24.06.1993
c                FHI 28.07.1993
c
c***********************************************************************
c
      subroutine corga86(d1,d2,dp11,dp22,dp12,uu,vv,dec,dvcup,dvcdn)
c
      implicit real*8 (a-h,o-z)
c
      data a1,a2,a3,a4,a5,a6,a7/ 2.568d-3
     &,                           1.443307452d-2
     &,                           2.843543831d-6
     &,                           5.411317331d0
     &,                           1.816419933d-1
     &,                           1.763993811d-2
     &,                           8.12908d-4/
      data t13,t23,t43,t53,t76/	.33333333333333333d0
     &,				.66666666666666667d0
     &,				.13333333333333333d1
     &,				.16666666666666667d1
     &,				.11666666666666667d1/
      data crt2/.1587401052d1/
c
      d    = d1+d2
      dm13 = 1.d0/d**t13
      d43  = d**t43
c
c gradient expansion coefficient
      c1 = a1 + dm13*(a2+dm13*a3)
      c2 = 1.d0 + dm13*(a4+dm13*(a5+dm13*a6))
      c  = 1.667d-3 + c1/c2
c
      dpnorm = sqrt(dp11+dp22+2.d0*dp12)
      dpnorm2= dpnorm*dpnorm
c
      fi = a7*dpnorm/(c*d**t76)
c
c spin interpolation
      zet= (d1-d2)/d
      dd = sqrt(.5d0*((1.d0+zet)**t53 + (1.d0-zet)**t53))
c
c dC(n)/dn
      cp = -t13/(d43*c2*c2)
     &    *(c2*(a2+2.d0*a3*dm13)-c1*(a4+dm13*(2.d0*a5+dm13*3.d0*a6)))
c
c spin-independent terms
      www=( (fi-1.d0)*(cp/c-t43/d)+(fi-2.d0)*fi*(cp/c+t76/d) )*dpnorm2
      uuu=(3.d0-fi)*fi*uu/dpnorm
      vvv=(fi-2.d0)*vv
c
c spin dependent term
      zzz=crt2*.5d0*t53/(dd*d43)**2*(d1**t23-d2**t23)
      zz1= zzz*((1.d0-fi)*d2*dpnorm2-(2.d0-fi)*d*(dp12+dp22))
      zz2=-zzz*((1.d0-fi)*d1*dpnorm2-(2.d0-fi)*d*(dp12+dp11))
c
      dec=c/(dd*exp(fi)*d43)
c
      dvcup=dec*(www+uuu+vvv+zz1)
      dvcdn=dec*(www+uuu+vvv+zz2)
      dec=dec*dpnorm2/d
c
      return
      end
