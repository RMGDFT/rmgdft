c************************** SVN Revision Information **************************
c **    $Id: pbe.f 1004 2008-08-09 18:37:21Z froze $    **
c******************************************************************************
 
c $Header:$
c==========================================================================
c c original version by K Burke
c Perdew-Burke-Ernzerhof GGA 
c see: J.P. Perdew, K. Burke, M. Ernzerhof, Phys Rev Lett 77, 3865 (1996).
c this collection contains the older PW91 and Becke exchange as well, 
c everything not needed for the PBE GGA is commented out w/ "c--"
c
c Martin Fuchs, FHI der MPG, Berlin, 11-1996
c==========================================================================
c
c--c PBE alpha2.1:
c--c Perdew-Burke-Ernzerhof generalized gradient approximation to the
c--c density functional for exchange-correlation energy of a many-electron
c--c system.
c--c  --------------------------------------------------------------------
c--c |WARNING!  PBE is a simplification of PW91, which yields almost      |
c--c |identical numerical results with simpler formulas from a simpler    |
c--c |derivation.  If you should find significant DIFFERENCES between     |
c--c |PBE and PW91 results, please consult kieron@merlin.phy.tulane.edu   |
c--c |or perdew@mailhost.tcs.tulane.edu.  Thank you.                      |
c--c  --------------------------------------------------------------------
c--c Note: Neglects small grad (zeta) contributions to the correlation
c--c energy.
c--c
c--c Programs implement functional in PBE paper, July 1996 version.
c--c 
c--c----------------------------------------------------------------------
c--c Main program testing PBE subroutines for exchange-correlation energy
c--c and potentials, by application to unequal exponential
c--c up- and down-spin densities, showing that the functional derivative
c--c (exchange-correlation potential) correctly generates the energy change
c--c due to a small change in the density.
c--c Kieron Burke, July 2, 1996.
c--C Atomic units are used, so all energies are in hartrees and all 
c--c distances in bohrs.  
c--c 1 hartree=27.2116eV=627.5kcal/mol; 1bohr=0.529E-10m.
c--c The output should be:
c--c Fup Fdn Zup Zdn             Exc           CHNG1          CHNG
c--c 1.0  .0 1.0  .5  -.311916530229   .000000000000   .0000000000
c--c 1.0  .2 1.0  .5  -.336377065446  -.053880102756  -.0538804290
c--c 1.0  .4 1.0  .5  -.369084886012  -.120463328976  -.1204642921
c--c 1.0  .6 1.0  .5  -.406088525151  -.193370595518  -.1933723422
c--c 1.0  .8 1.0  .5  -.446305936853  -.271252139632  -.2712547575
c--c 1.0 1.0 1.0  .5  -.489150144888  -.353405855349  -.3534094042
c--c 1.0 1.0  .5  .5  -.341059977353  -.316599687356  -.3166037653
c--c 1.0 1.0 1.0 1.0  -.653407740519  -.309758886707  -.3097606837
c--c 1.0 1.0 1.5 1.5  -.962039224827  -.307820467953  -.3078216918
c--c 1.0 1.0 2.0 2.0 -1.269410948459  -.307021487395  -.3070225637
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--      IMPLICIT REAL*8 (A-H,O-Z)
c--      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
c--      pi=4.d0*datan(1.d0)
c--      CONF=(3.D0*PI*pi)**THRD
c--      CONRS=(3.D0/(4.D0*PI))**THRD
c--      write(6,*)'Fup Fdn Zup Zdn             Exc'
c--     1,'           CHNG1          CHNG'
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--C BEGIN THE LOOP THAT SELECTS A TRIAL DENSITY
c--c spin-densities are of the form
c--c          rho(r)=f*(Z**3/pi)*dexp(-2*Z*r)
c--c delzdn=small change in zdn to test potentials
c--c jdens=counter for which density
c--      DO JDENS = 1,10
c--        FUP=1.D0
c--        FDN=0.2D0*(JDENS-1)
c--        ZUP=1.D0
c--        ZDN=0.5D0
c--        IF(JDENS.GT.6)then
c--	  FDN=1.D0
c--          ZUP=0.5D0+0.5D0*(JDENS-7)
c--          ZDN=ZUP
c--	endif
c--        DELZDN=1D-5
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--C BEGIN THE LOOP THAT INCREMENTS THE DENSITY DIFFERENTIALLY
c--c kdif=1=>density as above
c--c kdif=2=>Zdn changed by DELZDN
c--        DO KDIF=1,2
c--          IF(KDIF.EQ.2)ZDN=ZDN+DELZDN
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--C BEGIN THE RADIAL LOOP
c--c sumexc=integrated exchange-correlation energy 
c--c chng1=integrated xc energy change, based on vxc
c--c nr=number of points in radial loop
c--c rf=final value of r in integrals
c--c dr=change in r
c--c wt=weight of r in trapezoidal rule
c--c dup=up density
c--c agrup=|grad up|
c--c delgrup=(grad up).(grad |grad up|) 
c--c uplap=grad^2 up=Laplacian of up
c--c dn,agrdn,delgrdn,dnlap=corresponding down quantities
c--c d=up+dn
c--c agrad=|grad rho|
c--c delgrad=(grad rho).(grad |grad rho|) 
c--          sumexc=0.0D0
c--          CHNG1=0.0D0
c--	  nr=10000
c--	  rf=20.d0
c--	  dr=rf/real(nr)
c--          DO I=1,nr
c--            R=I*dr
c--            WT=4.d0*PI*R*R*dr
c--            DUP=FUP*(ZUP**3/PI)*DEXP(-2.D0*ZUP*R)
c--            DDN=FDN*(ZDN**3/PI)*DEXP(-2.D0*ZDN*R)
c--            ZDNNU=ZDN+DELZDN
c--            DELDDN=FDN*(ZDNNU**3/PI)*DEXP(-2.D0*ZDNNU*R)-DDN
c--	    agrup=2.d0*zup*dup
c--	    delgrup=8.d0*(zup**3)*dup*dup
c--	    uplap=4.d0*zup*dup*(zup-1.d0/r)
c--	    agrdn=2.d0*zdn*ddn
c--	    delgrdn=8.d0*(zdn**3)*ddn*ddn
c--	    dnlap=4.d0*zdn*ddn*(zdn-1.d0/r)
c--            D=DUP+DDN
c--            agrad=2.d0*(ZUP*DUP+ZDN*DDN)
c--	    delgrad=4.d0*agrad*(ZUP**2*DUP+ZDN**2*DDN)
c--            call easypbe(dup,agrup,delgrup,uplap,ddn,agrdn,delgrdn,
c--     1           dnlap,agrad,delgrad,1,1,
c--     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
c--     1           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
c--     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
c--	    sumexc=sumexc+d*(expbe+ecpbe)*wt
c--            CHNG1=CHNG1+(vxdnpbe+vcdnpbe)*DELDDN*WT/DELZDN
c--	  enddo
c--          IF(KDIF.EQ.1)then
c--	    sumEXCO=sumEXC
c--	  endif
c--        enddo
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--C  CHNG: DIRECT XC ENERGY INCREMENT
c--C  IF THE FUNCTIONAL DERIVATIVE IS CORRECT, THEN CHNG1=CHNG
c--        CHNG=(sumEXC-sumEXCO)/DELZDN
c--        PRINT 200,FUP,FDN,ZUP,ZDN,sumEXC,CHNG1,chng
c--      enddo
c--      STOP
c--  200 FORMAT(4f4.1,2f16.12,f14.10)
c--      END
c--c----------------------------------------------------------------------
c--c######################################################################
c--c----------------------------------------------------------------------
c--      subroutine easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
c--     1           agr,delgr,lcor,lpot,
c--     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
c--     1           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
c--     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--c EASYPBE is a driver for the PBE subroutines, using simple inputs
c--c K. Burke, May 14, 1996.
c--c inputs: up=up density
c--c	: agrup=|grad up|
c--c	: delgrup=(grad up).(grad |grad up|) 
c--c	: uplap=grad^2 up=Laplacian of up
c--c	: dn,agrdn,delgrdn,dnlap=corresponding down quantities
c--c	: agr=|grad rho|
c--c	: delgr=(grad rho).(grad |grad rho|) 
c--c	: lcor=flag to do correlation(=0=>don't)
c--c	: lpot=flag to do potential(=0=>don't)
c--c outputs: exlsd=LSD exchange energy density, so that
c--c		ExLSD=int d^3r rho(r) exlsd(r)
c--c	 : vxuplsd=up LSD exchange potential
c--c	 : vxdnlsd=down LSD exchange potential
c--c        : exclsd=LSD exchange-correlation energy density
c--c	 : vxcuplsd=up LSD exchange-correlation potential
c--c	 : vxcdnlsd=down LSD exchange-correlation potential
c--c        : expw91,vxuppw91,vxdnpw91,ecpw91,etc.=PW91 quantities
c--c        : expbe,vxuppbe,vxdnpbe,ecpbe,etc.=PBE quantities
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--c needed constants:
c--c pi32=3 pi**2
c--c alpha=(9pi/4)**thrd
c--      implicit real*8(a-h,o-z)
c--      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
c--      parameter(pi32=29.608813203268075856503472999628d0)
c--      parameter(pi=3.1415926535897932384626433832795d0)
c--      parameter(alpha=1.91915829267751300662482032624669d0)
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--c PBE exchange
c--c use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
c--c do up exchange
c--c fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3) 
c--c s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
c--c u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
c--c v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
c--      rho2=2.d0*up
c--      if(rho2.gt.1.d-15)then
c--        fk=(pi32*rho2)**thrd
c--        s=2.d0*agrup/(2.d0*fk*rho2)
c--        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
c--        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
c--        call exchpbe(rho2,s,u,v,0,lpot,exuplsd,vxuplsd)
c--        call exchpw91(rho2,s,u,v,exuppw91,vxuppw91)
c--        call exchpbe(rho2,s,u,v,1,lpot,exuppbe,vxuppbe)
c--      else
c--	exuplsd=0.d0
c--	vxuplsd=0.d0
c--	exuppw91=0.d0
c--	vxuppw91=0.d0
c--	exuppbe=0.d0
c--	vxuppbe=0.d0
c--      endif
c--c repeat for down
c--      rho2=2.d0*dn
c--      if(rho2.gt.1.d-15)then
c--        fk=(pi32*rho2)**thrd
c--        s=2.d0*agrdn/(2.d0*fk*rho2)
c--        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
c--        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
c--        call exchpbe(rho2,s,u,v,0,lpot,exdnlsd,vxdnlsd)
c--        call exchpw91(rho2,s,u,v,exdnpw91,vxdnpw91)
c--        call exchpbe(rho2,s,u,v,1,lpot,exdnpbe,vxdnpbe)
c--      else
c--	exdnlsd=0.d0
c--	vxdnlsd=0.d0
c--	exdnpw91=0.d0
c--	vxdnpw91=0.d0
c--	exdnpbe=0.d0
c--	vxdnpbe=0.d0
c--      endif
c--10    continue 
c--c construct total density and contribution to ex
c--      rho=up+dn
c--      if(rho .gt. 1.d-15) then
c--        exlsd=(exuplsd*up+exdnlsd*dn)/rho
c--        expw91=(exuppw91*up+exdnpw91*dn)/rho
c--        expbe=(exuppbe*up+exdnpbe*dn)/rho
c--      else
c--        exlsd=0.d0
c--        expw91=0.d0
c--        expbe=0.d0
c--      endif
c--      if(lcor.eq.0)return
c--c----------------------------------------------------------------------
c--c----------------------------------------------------------------------
c--c Now do correlation
c--c zet=(up-dn)/rho
c--c g=phi(zeta)
c--c rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
c--c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c--c twoksg=2*Ks*phi
c--c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c--c uu=delgrad/(rho^2*twoksg^3)
c--c rholap=Laplacian
c--c vv=Laplacian/(rho*twoksg^2)
c--c ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
c--c ec=lsd correlation energy
c--c vcup=lsd up correlation potential
c--c vcdn=lsd down correlation potential
c--c h=gradient correction to correlation energy
c--c dvcup=gradient correction to up correlation potential
c--c dvcdn=gradient correction to down correlation potential
c--      if(rho.lt.1.d-18)return
c--      zet=(up-dn)/rho
c--      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
c--      fk=(pi32*rho)**thrd
c--      rs=alpha/fk
c--      sk=sqrt(4.d0*fk/pi)
c--      twoksg=2.d0*sk*g
c--      t=agr/(twoksg*rho)
c--      uu=delgr/(rho*rho*twoksg**3)
c--      rholap=uplap+dnlap
c--      vv=rholap/(rho*twoksg**2)
c--      ww=(agrup**2-agrdn**2-zet*agr**2)/(rho*rho*twoksg**2)
c--      call CORPBE(RS,ZET,T,UU,VV,WW,1,lpot,ec,vcup,vcdn,
c--     1                  H,DVCUP,DVCDN)
c--      eclsd=ec
c--      ecpbe=ec+h
c--      vcuplsd=vcup
c--      vcdnlsd=vcdn
c--      vcuppbe=vcup+dvcup
c--      vcdnpbe=vcdn+dvcdn
c--      call CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
c--      call CORPW91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,DVCUP,DVCDN)
c--      ecpw91=ec+h
c--      vcuppw91=vcup+dvcup
c--      vcdnpw91=vcdn+dvcdn
c--      return
c--      end
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
	ex=exunif
        vx=ex*thrd4
	return
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c calculate potential from [b](24) 
      VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,
     1                  H,DVCUP,DVCDN)
c----------------------------------------------------------------------
c  Official PBE correlation code. K. Burke, May 14, 1996.
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
c       :  UU,VV,WW, only needed for PBE potential
c       : lgga=flag to do gga (0=>LSD only)
c       : lpot=flag to do potential (0=>energy only)
c  output: ec=lsd correlation energy from [a]
c        : vcup=lsd up correlation potential
c        : vcdn=lsd dn correlation potential
c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c        : dvcup=nonlocal correction to vcup
c        : dvcdn=nonlocal correction to vcdn
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
c thrd*=various multiples of 1/3
c numbers for use in LSD energy spin-interpolation formula, [c](9).
c      GAM= 2^(4/3)-2
c      FZZ=f''(0)= 8/(9*GAM)
c numbers for construction of PBE
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at 
c          |zeta|=1.
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find LSD energy contributions, using [c](10) and Table I[c].
c EU=unpolarized LSD correlation energy
c EURS=dEU/drs
c EP=fully polarized LSD correlation energy
c EPRS=dEP/drs
c ALFM=-spin stiffness, [c](3).
c ALFRSM=-dalpha/drs
c F=spin-scaling factor from [c](9).
c construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,rtrs,EU,EURS)
      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,rtRS,EP,EPRS)
      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)
c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c PBE correlation energy
c G=phi(zeta), given after [a](3)
c DELT=bet/gamma
c B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      if(lpot.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm-
     1((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
c slimmed down version of GCOR used in PW91 routines, to interpolate
c LSD correlation energy, as given by (10) of
c J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c K. Burke, May 11, 1996.
      IMPLICIT REAL*8 (A-H,O-Z)
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END
c--c----------------------------------------------------------------------
c--c######################################################################
c--c----------------------------------------------------------------------
c--      SUBROUTINE EXCHPW91(D,S,U,V,EX,VX)
c--C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c--C  INPUT D : DENSITY
c--C  INPUT S:  ABS(GRAD D)/(2*KF*D)
c--C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
c--C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
c--C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c--      IMPLICIT REAL*8 (A-H,O-Z)
c--      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.15084D0,a4=100.d0)
c--      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.004d0)
c--      parameter(thrd=0.333333333333D0,thrd4=1.33333333333D0)
c--c for Becke exchange, set a3=b1=0
c--      FAC = AX*D**THRD
c--      S2 = S*S
c--      S3 = S2*S
c--      S4 = S2*S2
c--      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
c--      P1 = DLOG(A*S+1.D0/P0)
c--      P2 = DEXP(-A4*S2)
c--      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
c--      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
c--      F = P3*P4
c--      EX = FAC*F
c--C  LOCAL EXCHANGE OPTION
c--C     EX = FAC
c--C  ENERGY DONE. NOW THE POTENTIAL:
c--      P5 = B1*S2-(A2-A3*P2)
c--      P6 = A1*S*(P1+A*S*P0)
c--      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
c--      FS = P3*(P3*P5*P6+P7)
c--      P8 = 2.D0*S*(B1-A3*A4*P2)
c--      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
c--      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
c--      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
c--      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
c--      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
c--C  LOCAL EXCHANGE OPTION:
c--C     VX = FAC*THRD4
c--      RETURN
c--      END
c--c----------------------------------------------------------------------
c--c######################################################################
c--c----------------------------------------------------------------------
c--      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
c--C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
c--C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
c--C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
c--C     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
c--C  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
c--      IMPLICIT REAL*8 (A-H,O-Z)
c--      parameter(gam=0.5198421D0,fzz=1.709921D0)
c--      parameter(thrd=0.333333333333D0,thrd4=1.333333333333D0)
c--      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
c--      CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
c--     1    0.49294D0,1.00D0,RS,EU,EURS)
c--      CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
c--     1    0.62517D0,1.00D0,RS,EP,EPRS)
c--      CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
c--     1    0.49671D0,1.00D0,RS,ALFM,ALFRSM)
c--C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
c--      ALFC = -ALFM
c--      Z4 = ZET**4
c--      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
c--C  ENERGY DONE. NOW THE POTENTIAL:
c--      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
c--      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
c--      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
c--     1        -(1.D0-Z4)*ALFM/FZZ)
c--      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
c--      VCUP = COMM + ECZET
c--      VCDN = COMM - ECZET
c--      RETURN
c--      END
c--c----------------------------------------------------------------------
c--c######################################################################
c--c----------------------------------------------------------------------
c--      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
c--C  CALLED BY SUBROUTINE CORLSD
c--      IMPLICIT REAL*8 (A-H,O-Z)
c--      P1 = P + 1.D0
c--      Q0 = -2.D0*A*(1.D0+A1*RS)
c--      RS12 = DSQRT(RS)
c--      RS32 = RS12**3
c--      RSP = RS**P
c--      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
c--      Q2 = DLOG(1.D0+1.D0/Q1)
c--      GG = Q0*Q2
c--      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
c--      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
c--      RETURN
c--      END
c--c----------------------------------------------------------------------
c--c######################################################################
c--c----------------------------------------------------------------------
c--      SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,
c--     1                   DVCUP,DVCDN)
c--C  pw91 CORRELATION, modified by K. Burke to put all arguments 
c--c  as variables in calling statement, rather than in common block
c--c  May, 1996.
c--C  INPUT RS: SEITZ RADIUS
c--C  INPUT ZET: RELATIVE SPIN POLARIZATION
c--C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
c--C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
c--C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
c--C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
c--C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c--C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
c--      IMPLICIT REAL*8 (A-H,O-Z)
c--      parameter(xnu=15.75592D0,cc0=0.004235D0,cx=-0.001667212D0)
c--      parameter(alf=0.09D0)
c--      parameter(c1=0.002568D0,c2=0.023266D0,c3=7.389D-6,c4=8.723D0)
c--      parameter(c5=0.472D0,c6=7.389D-2,a4=100.D0)
c--      parameter(thrdm=-0.333333333333D0,thrd2=0.666666666667D0)
c--      BET = XNU*CC0
c--      DELT = 2.D0*ALF/BET
c--      G3 = G**3
c--      G4 = G3*G
c--      PON = -DELT*EC/(G3*BET)
c--      B = DELT/(DEXP(PON)-1.D0)
c--      B2 = B*B
c--      T2 = T*T
c--      T4 = T2*T2
c--      T6 = T4*T2
c--      RS2 = RS*RS
c--      RS3 = RS2*RS
c--      Q4 = 1.D0+B*T2
c--      Q5 = 1.D0+B*T2+B2*T4
c--      Q6 = C1+C2*RS+C3*RS2
c--      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
c--      CC = -CX + Q6/Q7
c--      R0 = 0.663436444d0*rs
c--      R1 = A4*R0*G4
c--      COEFF = CC-CC0-3.D0*CX/7.D0
c--      R2 = XNU*COEFF*G3
c--      R3 = DEXP(-R1*T2)
c--      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
c--      H1 = R3*R2*T2
c--      H = H0+H1
c--C  LOCAL CORRELATION OPTION:
c--C     H = 0.0D0
c--C  ENERGY DONE. NOW THE POTENTIAL:
c--      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
c--      RSTHRD = RS/3.D0
c--      R4 = RSTHRD*CCRS/COEFF
c--      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
c--      FAC = DELT/B+1.D0
c--      BG = -3.D0*B2*EC*FAC/(BET*G4)
c--      BEC = B2*FAC/(BET*G3)
c--      Q8 = Q5*Q5+DELT*Q4*Q5*T2
c--      Q9 = 1.D0+2.D0*B*T2
c--      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
c--      H0RS = -RSTHRD*H0B*BEC*ECRS
c--      FACT0 = 2.D0*DELT-6.D0*B
c--      FACT1 = Q5*Q9+Q4*Q9*Q9
c--      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
c--      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
c--      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
c--      H0T = 2.*BET*G3*Q9/Q8
c--      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
c--      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
c--      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
c--      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
c--      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
c--      FACT4 = 2.D0-R1*T2
c--      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
c--      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
c--      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
c--      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
c--      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
c--      HRS = H0RS+H1RS
c--      HRST = H0RST+H1RST
c--      HT = H0T+H1T
c--      HTT = H0TT+H1TT
c--      HZ = H0Z+H1Z
c--      HZT = H0ZT+H1ZT
c--      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
c--      PREF = HZ-GZ*T2*HT/G
c--      FACT5 = GZ*(2.D0*HT+T*HTT)/G
c--      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
c--      DVCUP = COMM + PREF
c--      DVCDN = COMM - PREF
c--C  LOCAL CORRELATION OPTION:
c--C     DVCUP = 0.0D0
c--C     DVCDN = 0.0D0
c--      RETURN
c--      END
c--c----------------------------------------------------------------------
