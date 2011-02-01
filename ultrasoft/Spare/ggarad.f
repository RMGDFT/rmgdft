c************************** SVN Revision Information **************************
c **    $Id$    **
c******************************************************************************
 
c $Header:$
c***********************************************************************
c
c exchange potential & energy density
c spherical symmetry
c LSDA - GGA
c
c Hartree a.u.
c
c input
c    mode = 1    LSDA
c    mode = 2    GGA-X Becke
c    mode = 3    GGA-X Perdew
c    r           radius
c    rh()        spin up/down density
c    rhp()       1st radial derivative of rh()
c    rhpp()      2nd radial derivative of rh()
c
c output
c    xpot()      up/down exchange potential
c    xen         exchange energy density, i.e. per electron
c
c convention
c    yy(1) = spin up, yy(2) = spin down
c
c author -- M. Fuchs, FHI Berlin, 07-1992
c***********************************************************************
c
      subroutine ggaxrad(mode,r,rh,rhp,rhpp,xpot,xen)
c
      implicit real*8 (a-h,o-z)
c
      dimension rh(2),rhp(2),rhpp(2),xpot(2)
      data thrd/.333333333333333333d0/
      data pisq3/.296088132032680740d2/
      data eps /1.d-15/
c
c exchange GGA, loop for up & down spin
c
      xen=0.d0
      do 20 i=1,2
        if(rh(i) .le. eps) then
          xpot(i)=0.d0
        else
          d=2.d0*rh(i)
          if(mode .eq. 1) then
            call xlda(d,vx,ex)
          else if(mode .eq. 2 .or. mode .eq. 3) then
            dp  = 2.d0*rhp(i)
            dpp = 2.d0*rhpp(i)
            fk    = (pisq3*d)**thrd
            fk2   = 2.d0*fk
            fk2sq = fk2**2
            s = abs(dp)/(fk2*d)
            u = abs(dp)*dpp/(d*d*fk2*fk2sq)
            v = (dpp+2.d0*dp/r)/(d*fk2sq)
            if(mode .eq. 2) then
              call xbecke(d,s,u,v,ex,vx)
            else if(mode .eq. 3) then
              call exch(d,s,u,v,ex,vx)
            endif
          else
            stop 'ggaxrad : mode improper'
          endif
          xpot(i) = vx
          xen = xen+rh(i)*ex
        endif
 20   continue
      xen = xen/max(rh(1)+rh(2),eps)

      return
      end
c
c**********************************************************************
c
c correlation potential & energy density
c spherical symmetry
c LSDA - GGA
c Hartree a.u.
c
c input
c    mode = 1    LSDA
c    mode = 2    GGA-C Perdew 91
c    mode = 3    GGA-C Perdew 86
c    r           radius
c    rh()        spin up/down density
c    rhp()       1st derivative of rh
c    rhpp()      2nd derivative of rh
c
c output
c    zet         spin polarization
c    cpot()       - " -  correlation potential
c    cen         correlation energy density
c
c convention
c    yy(1) = spin up, yy(2) = spin down
c
c author -- M. Fuchs, FHI Berlin, 07-1992
c**********************************************************************
c
      subroutine ggacrad(mode,r,rh,rhp,rhpp,cpot,cen)
c
      implicit real*8 (a-h,o-z)
      common /gas/fk,sk,g,ec,ecrs,eczet
      dimension rh(2),rhp(2),rhpp(2),cpot(2)
      data thrd,pi /.333333333333333333d0,.314159265358979312d1/
      data pisq3,thrd2 /.296088132032680740d2,.666666666666666666d0/
      data crs,eps /1.91915829267751281d0,1.d-15/
c
c LSDA
      d = rh(1)+rh(2)
      cen=0.d0
      if(d .le. eps) then
        cen=0.d0
        cpot(1)=0.d0
        cpot(2)=0.d0
      else
        zet=(rh(1)-rh(2))/d
        fk=(pisq3*d)**thrd
        rs=crs/fk
        call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
c
c GGA correction to LSDA
        if(mode .eq. 1) then
          dvcup=0.d0
          dvcdn=0.d0
          h=0.d0
        else if(mode .eq. 2) then
          dp=rhp(1)+rhp(2)
          dpp=rhpp(1)+rhpp(2)
          ztp=(rhp(1)-rhp(2)-zet*dp)/d
          sk=2.d0*dsqrt(fk/pi)
          g=((1.d0+zet)**thrd2 + (1.d0-zet)**thrd2) /2.d0
          gks2=2.d0*sk*g
          gks2sq=gks2*gks2
          t=abs(dp)/(d*gks2)
          uu=abs(dp)*dpp/(d*d*gks2sq*gks2)
          vv=(dpp +2.d0*dp/r)/(d*gks2sq)
          ww=dp*ztp/(d*gks2sq)
          call corgga(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn)
        else if(mode .eq. 3) then
          dp  =rhp(1)+rhp(2)
          dpp =rhpp(1)+rhpp(2)
          uu  =abs(dp)*dpp
          vv  =dpp +2.d0*dp/r
          dp11=rhp(1)*rhp(1)
          dp22=rhp(2)*rhp(2)
          dp12=rhp(1)*rhp(2)
          if(rhp(1) .ne. 0.d0 .or. rhp(2) .ne. 0.d0)
     &      call corga86(rh(1),rh(2),dp11,dp22,dp12,uu,vv,h,dvcup,dvcdn)
        else
          stop 'ggacrad : mode improper'
        endif
        cpot(1)=vcup+dvcup
        cpot(2)=vcdn+dvcdn
        cen=ec+h
      endif
c
      return
      end
ce
