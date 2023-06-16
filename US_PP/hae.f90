!
! Copyright (c) 1989-201r by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module haesratom
use iso_c_binding
IMPLICIT NONE
PRIVATE
SAVE
PUBLIC :: aeo, aii, sratom, lschfb, tfapot

CONTAINS

 subroutine sratom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&           rr,vi,zz,mmax,iexc,etot,ierr,uu)
!&           rr,vi,zz,mmax,iexc,etot,ierr,srel) bind(c, name="sratom")

! self-consistent scalar-relativistic all-electron atom
! calculation using log mesh (non-relativistic when srel=.false.)

!na  principal quantum number array, dimension ncv
!la  angular-momenta
!ea  eigenvalues (output)
!fa  occupancies
!rpk  radius of outermost peak of wave function
!nc  number of core states
!ncv  number of core+valence states
!it  number of iterations (output)
!rr  log radial mesh
!vi  all-electron potential (output)
!zz  atomic number
!mmax  size of log grid
!iexc  exchange-correlation function to be used
!etot  all-electron total energy (output)
!ierr  error flag
!uu   atomic wavefunctions
!srel  .true. for scalar-relativistic, .false. for non-relativistic

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 
 integer :: mmax,iexc,nc,ncv,idx
 integer :: na(ncv),la(ncv)
 real(dp) :: zz
 real(dp) :: fa(ncv),rr(mmax)
! logical :: srel

!Output variables
 integer :: it,ierr
 real(dp) :: etot
 real(dp) :: ea(ncv),rpk(ncv)
 real(dp) :: rho(mmax),rhoc(mmax),vi(mmax),uu(mmax,ncv)

!Local function
! real(dp) :: tfapot

!Local variables
 integer :: nin,mch
 real(dp) :: amesh,al
 real(dp) :: dr,eeel,eexc,et,rl,rl1,sd,sf,sn,eeig
 real(dp) :: thl,vn,zion,vofz
 integer :: ii,jj
 logical :: convg
 logical :: srel

 real(dp), allocatable :: u(:),up(:)
 real(dp), allocatable :: vo(:),vi1(:),vo1(:),vxc(:)

! blend parameter for Anderson iterative potential mixing
 real(dp), parameter ::  bl=0.5d0

 allocate(u(mmax),up(mmax))
 allocate(vo(mmax),vi1(mmax),vo1(mmax),vxc(mmax))

! why all this is necessary is unclear, but it seems to be
 u(:)=0.d0; up(:)=0.d0; vo(:)=0.d0; vi1(:)=0.d0; vo1(:)=0.d0; vxc(:)=0.d0
 dr=0.d0; eeel=0.d0; eexc=0.d0; et=0.d0; rl=0.d0; rl1=0.d0
 sd=0.d0; sf=0.d0; sn=0.d0; eeig=0.d0; thl=0.d0; vn=0.d0; zion=0.d0 
 nin=0; mch=0
 srel = .true.

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

 do ii=1,mmax
  vi(ii)=tfapot(rr(ii),zz)
 end do

! starting approximation for energies
 sf=0.0d0
 do ii=1,ncv
   sf=sf+fa(ii)
   zion=zz+1.0d0-sf
   ea(ii)=-0.5d0*(zion/na(ii))**2
   if(ea(ii)>vi(mmax)) ea(ii)=2.0d0*vi(mmax)
 end do

! big self  self-consietency loop

 do it=1,500
   convg=.true.

   rhoc(:) = 0.0d0
   rho(:)=0.0d0

! solve for bound states in turn
   eeig=0.0d0
   do ii=1,ncv
     et=ea(ii)
     ierr = 0
     call lschfb(na(ii),la(ii),ierr,et, &
&                rr,vi,u,up,zz,mmax,mch,srel)
! EMIL for C++
     do jj=1,mmax
       uu(jj,ii) = u(jj)
     end do
     if(ierr>0) then
       write(6,'(/a,3i4)') 'sratom: lschfb convergence ERROR n,l,iter=', &
&       na(ii),la(ii),it
       stop
     end if

! overall convergence criterion based on eps within lschfb
     if(ea(ii)/= et) convg=.false.
     ea(ii)=et

! accumulate charge and eigenvalues
     eeig = eeig + fa(ii) * ea(ii)
     rho(:)=rho(:) + fa(ii)*(u(:)/rr(:))**2
     if(ii<=nc) then
       rhoc(:)=rhoc(:) + fa(ii)*(u(:)/rr(:))**2
     end if


! find outermost peak of wavefunction
     do jj=mch-1,1,-1
       if(up(jj)*up(jj+1)<0.0d0) then
         rpk(ii)=rr(jj)
         exit
       end if
     end do

   end do

   if(ierr/=0) then
    exit
   end if


! output potential
   call vout(0,rho,rhoc,vo,vxc,sf-zz,eeel,eexc, &
&            rr,mmax,iexc)

! generate next iteration using d. g. anderson''s
! method
   thl=0.0d0
   if(it>1) then
     sn=0.0d0
     sd=0.0d0
     do ii=1,mmax
       rl=vo(ii)-vi(ii)
       rl1=vo1(ii)-vi1(ii)
       dr=rl-rl1
       sn=sn + rl*dr*rr(ii)**2
       sd=sd + dr*dr*rr(ii)**2
     end do
     thl=sn/sd
   end if

   do ii=1,mmax
     vn=(1.0d0-bl)*((1.0d0-thl)*vi(ii) + thl*vi1(ii)) &
  &   + bl*((1.0d0-thl)*vo(ii) + thl*vo1(ii))
     vi1(ii)=vi(ii)
     vo1(ii)=vo(ii)
!     vi(ii)=vn
     vi(ii)=0.3*(vofz(zz, rr(ii)) + vo(ii)) + 0.7*vi(ii)
   end do

   if(convg) exit

   if(it==500 .and. .not. convg) then
     write(6,'(/a)') 'sratom: WARNING failed to converge'
   end if

 end do !it

! EMIL all electon stuff
! do idx=1,mmax
!     write(6,*)'AEPLOT  ',rr(idx),'  ', vo(idx),'  ',vi(idx)
! end do

 if(.not. convg .and. ierr==0) then
   ierr=100
 end if

! total energy output

! output potential for e-e interactions

 call vout(0,rho,rhoc,vo,vxc,sf,eeel,eexc, &
&          rr,mmax,iexc)

 etot =  eeig + eexc - 0.5d0*eeel

 deallocate(u,up)
 deallocate(vo,vi1,vo1,vxc)
 return

 end subroutine sratom


!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! tfapot
 function tfapot(rr,zz) bind(c, name="tfapot")

! generalized Thomas-Fermi atomic potential

!...to an article of N. H. March ( "The Thomas-Fermi Approximation in 
! Quantum Mechanics", Adv. Phys. 6, 1 (1957)). He has the formula,
! but it is not the result of his work. The original publication is: 
!     R. Latter, Phys. Rev. 99, 510 (1955).
! He says that it''s an analytic fit to an improved calculation of the 
! potential distribution of a Thomas-Fermi atom without exchange first 
! performed by Miranda (C. Miranda, Mem. Acc. Italia 5, 285 (1934)).
!                                 Alexander Seidl, TU Munich

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr,zz
 
!Output function
 real(dp) :: tfapot

!Local variables
 real(dp) :: bb,tt,xx,xs

 bb=(0.69395656d0/zz)**.33333333d0
 xx=rr/bb
 xs=sqrt(xx)

 tt=zz/(1.0d0+xs*(0.02747d0 - xx*(0.1486d0 - 0.007298d0*xx)) &
&   + xx*(1.243d0 + xx*(0.2302d0 + 0.006944d0*xx)))

 if(tt<1.0d0) tt=1.0d0
 tfapot=-tt/rr

 return
 end function tfapot


!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! calculates Hartree and exchange-correlation potentials and total-energy
! term to add to eigenvalue sum

 subroutine vout(mode,rho,rhoc,vo,vxc,zion,eeel,eexc, &
&                rr,mmax,iexc)

!mode  1=> add rhoc to rho, 0 => don't
!rho  total charge density or valence/pseudovalence charge density
!rhoc  core charge density, full or model
!vo  output total potential
!vc  output exchange-correlation potential
!zion  charge of screened ion
!eeel  electron-electron interaction energy
!eexc  exchange-correlation correction to eigenvalue sum for total energy
!rr  log radial mesh
!output electrostatic and exchange-correlation potential

 implicit none

 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0_dp*pi

!Input vaiables
 integer :: iexc,mmax,mode
 real(dp) :: rho(mmax),rhoc(mmax),rr(mmax)
 real(dp) :: zion

!Output variables
 real(dp) :: vo(mmax),vxc(mmax)
 real(dp) :: eeel,eexc

!Local variables
 integer ii
 real(dp) :: al,tv

!Function
! real(dp) :: aii

 real(dp), allocatable :: rvp(:),rv(:),excca(:),difxc(:)
 real(dp), allocatable :: exca(:),vxcd(:),rhot(:)
 real(dp), allocatable :: foutv(:),foutx(:),foutfc(:)

 allocate(rvp(mmax),rv(mmax),excca(mmax),difxc(mmax))
 allocate(exca(mmax),vxcd(mmax),rhot(mmax))
 allocate(foutv(mmax),foutx(mmax),foutfc(mmax))


 al = 0.01d0 * dlog(rr(101) / rr(1))

 foutfc(:)=0.0d0

! integration for electrostatic potential
 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**3
 end do

 rv(mmax)=zion
 rv(mmax-1)=zion
 rv(mmax-2)=zion

 do ii=mmax-2,2,-1
   rv(ii-1)=rv(ii)+aii(rvp,ii)
 end do

 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**2
 end do

 tv=0.0d0
 do ii=mmax-2,2,-1
   tv=tv+aii(rvp,ii)
   rv(ii-1)=rv(ii-1)-rr(ii-1)*tv
 end do

 do ii=1,mmax
   vo(ii)=rv(ii)/rr(ii)
 end do

! electron-electron interaction for total energy
 do ii=1,mmax
  foutv(ii)=rho(ii)*vo(ii)*rr(ii)**3
 end do

 eeel=(9.0d0*foutv(1) + 28.0d0*foutv(2) &
&   + 23.0d0*foutv(3))/24.0d0
 do ii=4,mmax
   eeel=eeel + foutv(ii)
 end do
 eeel=al*eeel + foutv(1)/3.0d0

 if(mode .eq. 0) then
   do ii = 1, mmax
     rhot(ii) = rho(ii)
   end do
 else if(mode .eq. 1) then
   do ii = 1, mmax
     rhot(ii) = rho(ii) + rhoc(ii)
   end do
 else
   write(6,'(/a,i4)') 'vout: ERROR bad input mode =',mode
   stop
 end if

! exchange-correlation potential added

 if(iexc .eq. 1) then
   call excwig(rhot,vxc,exca,mmax)
 else if(iexc .eq. 2) then
   call exchdl(rhot,vxc,exca,mmax)
 else if(iexc .eq. 3) then
   call excpzca(rhot,vxc,exca,mmax)
 else if(iexc .eq. 4) then
   call excggc(rhot,vxc,exca,rr,mmax)
 else if (iexc < 0) then
   call exc_libxc(iexc,al,rhot,vxc,exca,rr,mmax)
 else
   write(6,'(/a,i4)') 'vout: ERROR bad input iexc =',iexc
   stop
 end if

 do ii = 1, mmax
   difxc(ii)=exca(ii) - vxc(ii)
   vo(ii)=vo(ii) + vxc(ii)
 end do

! exchange-correlation correction for total energy

 do ii=1,mmax
  foutx(ii)=rho(ii)*difxc(ii)*rr(ii)**3
 end do
 eexc=(9.0d0*foutx(1) + 28.0d0*foutx(2) &
&   + 23.0d0*foutx(3))/24.0d0
 do ii=4,mmax
   eexc=eexc + foutx(ii)
 end do

 if(mode .eq. 1 .and. rhoc(1) .ne. 0.0d0) then
   if(iexc .eq. 1) then
     call excwig(rhoc,vxcd,excca,mmax)
   else if(iexc .eq. 2) then
     call exchdl(rhoc,vxcd,excca,mmax)
   else if(iexc .eq. 3) then
     call excpzca(rhoc,vxcd,excca,mmax)
   else if(iexc .eq. 4) then
     call excggc(rhoc,vxcd,excca,rr,mmax)
   else if (iexc < 0) then
     call exc_libxc(iexc,al,rhoc,vxcd,excca,rr,mmax)
   else
     write(6,'(/a,i4)') 'vout: ERROR bad input iexc =',iexc
     stop
   end if

   do ii=1,mmax
    foutfc(ii)=rhoc(ii)*(exca(ii) - excca(ii))*rr(ii)**3
   end do
   eexc=eexc + (9.0d0*foutfc(1) + 28.0d0*foutfc(2) &
&     + 23.0d0*foutfc(3))/24.0d0
   do ii=4,mmax
     eexc=eexc + foutfc(ii)
   end do
 end if

 eexc=al*eexc + foutx(1)/3.0d0


 if(mode .eq. 1) then
   eexc=eexc + foutfc(1)/3.0d0
 end if

 deallocate(rvp,rv,excca,difxc)
 deallocate(exca,vxcd,rhot)
 deallocate(foutv,foutx,foutfc)
 return

 end subroutine vout


!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine exchdl(rho,vxc,exc,mmax)

!calculates Hedin-Lundquist exchange-correlation potential and energy
!density

 implicit none

 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4
 real(dp), parameter :: thrd=1.0d0/3.0d0

!Inpupt variables
 integer :: mmax
 real(dp) :: rho(mmax)

!Output variables
 real(dp) :: vxc(mmax),exc(mmax)

!Local variables
 integer :: ii
 real(dp) :: conrs,convx,conex
 real(dp) :: rs,ecp,aln,xx,rh

 conrs = (3.d0/(4.d0*pi))**thrd
 convx = (1.5d0/pi)**(2.0d0/3.0d0)
 conex = 0.75d0*convx

! Hedin-Lundqvist correlation

 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
     rh=rho(ii)*pi4i
     rs=conrs/rh**thrd
     xx=rs/21.0d0
     aln=dlog(1.0d0 + 1.0d0/xx)
     ecp = aln+(xx**3*aln-xx*xx)+0.5d0*xx-thrd
!    ecp = aln+(xx**3*aln-xx*xx)+xx/2-1.0d0/3.0d0
!    exc(ii)=-0.458175d0/rs - 0.0225d0*ecp
!    vxc(ii)=-0.6109d0/rs - 0.0225d0*aln
     exc(ii)=-conex/rs - 0.0225d0*ecp
     vxc(ii)=-convx/rs - 0.0225d0*aln
   else
     vxc(ii)=0.0d0 ; exc(ii)=0.0d0
   end if
 end do

 return
 end subroutine exchdl
!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University, and M. Verstratte, University of Liege
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 
 
subroutine derivs(mmax, rho, al, rr, dpr, dppr, dlap)
! NB: presumes incoming rho is properly scaled without 4pi or r**2 factors.
 
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 
 integer, intent(in) :: mmax
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: rho(mmax)
 real(dp), intent(out) :: dpr(mmax), dppr(mmax), dlap(mmax)
 
! local vars
 integer :: i
 real(dp) :: dpn(mmax), dppn(mmax)
 real(dp) :: c11,c12,c13,c14,c15
 real(dp) :: c21,c22,c23,c24,c25
 
 c11 =   2.0d0 / 24.0d0
 c12 = -16.0d0 / 24.0d0
 c13 =   0.0d0 / 24.0d0
 c14 =  16.0d0 / 24.0d0
 c15 =  -2.0d0 / 24.0d0
!
 c21 =   -1.0d0 / 12.0d0
 c22 =   16.0d0 / 12.0d0
 c23 =  -30.0d0 / 12.0d0
 c24 =   16.0d0 / 12.0d0
 c25 =   -1.0d0 / 12.0d0
!
! n derivatives of d
!     
 i=1
 dpn(i) = -25.d0/12.d0*rho(i) +4.d0*rho(i+1) -3.d0*rho(i+2) &
&         +4.d0/3.d0*rho(i+3) -1.d0/4.d0*rho(i+4)
 dppn(i) = 15.d0/4.d0*rho(i) -77.d0/6.d0*rho(i+1) +107.d0/6.d0*rho(i+2) &
&         -13.d0*rho(i+3) +61.d0/12.d0*rho(i+4) -5.d0/6.d0*rho(i+5)
 i=2
 dpn(i) = -25.d0/12.d0*rho(i) +4.d0*rho(i+1) -3.d0*rho(i+2)  &
&         +4.d0/3.d0*rho(i+3) -1.d0/4.d0*rho(i+4)
 dppn(i) = 15.d0/4.d0*rho(i) -77.d0/6.d0*rho(i+1) +107.d0/6.d0*rho(i+2) &
&         -13.d0*rho(i+3) +61.d0/12.d0*rho(i+4) -5.d0/6.d0*rho(i+5)
 
 do i = 3, mmax - 2
   dpn(i) =  c11*rho(i-2) + c12*rho(i-1) + c14*rho(i+1) + c15*rho(i+2)
   dppn(i) = c21*rho(i-2) + c22*rho(i-1) + c23*rho(i)   + c24*rho(i+1) &
&           +c25*rho(i+2)
 end do
 
 i=mmax-1
 dpn(i) = +25.d0/12.d0*rho(i) -4.d0*rho(i-1) +3.d0*rho(i-2) &
&         -4.d0/3.d0*rho(i-3) +1.d0/4.d0*rho(i-4)
 dppn(i) = -15.d0/4.d0*rho(i) +77.d0/6.d0*rho(i-1) -107.d0/6.d0*rho(i-2) &
&          +13.d0*rho(i-3) -61.d0/12.d0*rho(i-4) +5.d0/6.d0*rho(i-5)
 i=mmax
 dpn(i) = +25.d0/12.d0*rho(i) -4.d0*rho(i-1) +3.d0*rho(i-2) &
&         -4.d0/3.d0*rho(i-3) +1.d0/4.d0*rho(i-4)
 dppn(i) = -15.d0/4.d0*rho(i) +77.d0/6.d0*rho(i-1) -107.d0/6.d0*rho(i-2) &
&          +13.d0*rho(i-3) -61.d0/12.d0*rho(i-4) +5.d0/6.d0*rho(i-5)
 
!
! r derivatives of d
!
 do i = 1, mmax
   dpr(i) = dpn(i) / (al * rr(i))
   dppr(i) = (dppn(i) - al * dpn(i)) / (al * rr(i))**2
   dlap(i) = (dppn(i) + al * dpn(i)) / (al * rr(i))**2
 end do
 
end subroutine derivs
!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University, and M. Verstratte, University of Liege
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 
subroutine exc_libxc(iexc,al,rho,vxc,exc,rr,mmax)
!
!  dummy substitute for interface to libxc library to allow successful
!  build without this library
!
 implicit none
 
 integer, parameter :: dp=kind(1.0d0)
 
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),rr(mmax)
 real(dp), intent(out) :: vxc(mmax),exc(mmax)
 integer, intent(in) :: mmax, iexc
 
 write(6,'(a,i8,a)') 'exc_libxc_stub: ERROR iexc = ',iexc,' requires libxc'
 write(6,'(a)') 'The present oncvpsp executable was built without libxc.'
 write(6,'(a)') 'Program will stop.'
 
 stop
 return
end subroutine exc_libxc
!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine excpzca(rho,vxc,exc,mmax)

!calculates Perdew-Zunger-Ceperly_Alder exchange-correlation potential 
! and energy density J. P. Perdew and A. Zunger, Phys. Rev. B23, 5048 (1981)

!rho  charge density
!vxc  exchange-correlation potential
!exc  exchange-correlation energy density
!mmax  dimension of log radial grid

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4
 real(dp), parameter :: thrd=1.0d0/3.0d0

!Input variables
 integer :: mmax
 real(dp) :: rho(mmax)

!Output variables
 real(dp) :: vxc(mmax),exc(mmax)

!Local variables
 integer :: ii
 real(dp) :: rs,rh,sqrs,rsl,den

 real(dp), parameter :: conrs = (3.d0/(4.d0*pi))**thrd
 real(dp), parameter :: convx = (1.5d0/pi)**(2.0d0/3.0d0)
 real(dp), parameter ::conex = 0.75d0*convx


!rs>1, Ceperley-Alder fit to QMC
 real(dp), parameter :: gam = -0.1423d0
 real(dp), parameter :: bt1 =  1.0529d0
 real(dp), parameter :: bt2 =  0.3334d0

!rs<1, Perdew-Zunger leading terms of RPA high-density expansion
!AA and BB known from Gell-Mann & Bruckner
!PZ adjust CC and DD for contiunity of exc and vxc at rs=1
 real(dp), parameter :: AA =   0.0311d0
 real(dp), parameter :: BB =  -0.0480d0
 real(dp), parameter :: CC =   0.0020d0
 real(dp), parameter :: DD =  -0.0116d0

!PZ modified coeffients from libxc
!real(dp), parameter :: CC =   0.0020191519406228d0
!real(dp), parameter :: DD =  -0.0116320663789130d0

 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
      rh=rho(ii)/pi4
      rs=conrs*rh**(-thrd)
      if(rs > 1.0d0) then
        sqrs=dsqrt(rs)
        den=1.0d0 + bt1*sqrs +bt2*rs
        exc(ii)=-conex/rs + gam/den
        vxc(ii)=-convx/rs + gam*(1.0d0 + (3.5d0*bt1*sqrs &
   &            + 4.0d0*bt2*rs)*thrd)/den**2
      else
        rsl=dlog(rs)
        exc(ii)=-conex/rs + AA*rsl + BB + CC*rs*rsl + DD*rs
        vxc(ii)=-convx/rs + AA*rsl + (BB-thrd*AA) + 2.0d0*thrd*CC*rs*rsl &
   &            + thrd*(2.0d0*DD-CC)*rs
      end if
   else
      vxc(ii)=0.0d0 ; exc(ii)=0.0d0
   end if
 end do

 return
 end subroutine excpzca
!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine excwig(rho,vxc,exc,mmax)

! Wigner ingerpolation formula for exchange-correlation potential and
! energy density

!rho  charge density
!vxc  exchange-correlation potential
!exc  exchange-correlation energy density
!mmax  dimension of log radial grid

 implicit none

 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4
 real(dp), parameter :: thrd=1.0d0/3.0d0

!Input variables
 real(dp) ::  rho(mmax)
 integer :: mmax

!Output variables
 real(dp) :: vxc(mmax),exc(mmax)

!Local vafriables
 real(dp) ::  rh,rs
 real(dp) :: conrs,convx,conex
 integer :: ii

 conrs = (3.d0/(4.d0*pi))**thrd
 convx = (1.5d0/pi)**(2.0d0/3.0d0)
 conex = 0.75d0*convx

! Wigner interpolation formula, E. Wigner, Phys. Rev. 46, 1002 (1934).

 do ii=1,mmax
  if(rho(ii)>=1.0d-15) then
   rh=rho(ii)*pi4i
   rs=conrs*rh**(-thrd)
   vxc(ii)=-convx/rs - 0.44d0*(4.0d0*thrd*rs+7.8d0)/(rs+7.8d0)**2
   exc(ii)=-conex/rs - 0.44d0/(rs+7.8d0)
  end if
 end do
 return
 end subroutine excwig


!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine lschfb(nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,srel) bind(c, name="lschfb")

!Finds bound states of an all-electron atomic potential using
!Pauli-type  scalar-relativistic Schroedinger equation

!nn  principal quantum number
!ll  angular-momentum quantum number
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!rr  log radial mesh
!vv  local atomic potential
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations
!srel .true. for scalar-relativistic, .false. for non-relativistic

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr(mmax),vv(mmax)
 real(dp) :: zz
 integer :: nn,ll
 integer :: mmax
 logical :: srel

!Output variables
 real(dp) :: uu(mmax),up(mmax)
 real(dp) :: ee
 integer :: ierr,mch

!Local variables

! real(dp) :: aei,aeo,aii,aio,als !functions in aeo.f90
 real(dp) als
 real(dp) :: de,emax,emin
 real(dp) :: eps,fss,tfss,gamma,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap
 real(dp) :: amesh,al,xx
 integer :: ii,it,nint,node,nin

 real(dp), allocatable :: upp(:),cf(:),dv(:),fr(:),frp(:)
 allocate(upp(mmax),cf(mmax),dv(mmax),fr(mmax),frp(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
 eps=1.0d-10
 ierr = 60

! relativistic - non-relativistic switch
 if(srel) then
  fss=(1.0d0/137.036d0)**2
 else
  fss=1.0d-12
 end if

 if(ll==0) gamma=dsqrt(1.0d0-fss*zz**2)
 if(ll>0) gamma=(ll*dsqrt(ll**2-fss*zz**2) + &
& (ll+1)*dsqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

 sls=ll*(ll+1)

 emax=vv(mmax)+0.5d0*sls/rr(mmax)**2
 emin=0.0d0
 do ii=1,mmax
   emin=dmin1(emin,vv(ii)+0.5d0*sls/rr(ii)**2)
 end do
 if(ee>emax) ee=1.25d0*emax
 if(ee<emin) ee=0.75d0*emin
 if(ee>emax) ee=0.5d0*(emax+emin)

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0

 als=al**2

! return point for bound state convergence
 do nint=1,60
  
! coefficient array for u in differential eq.
   do ii=1,mmax
     cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
   end do
  
! calculate dv/dr for darwin correction
   dv(:)=0.0d0

   dv(1)=(-50.d0*vv(1)+96.d0*vv(2)-72.d0*vv(3)+32.d0*vv(4) &
&         -6.d0*vv(5))/(24.d0*al*rr(1))
   dv(2)=(-6.d0*vv(1)-20.d0*vv(2)+36.d0*vv(3)-12.d0*vv(4) &
&         +2.d0*vv(5))/(24.d0*al*rr(2))
  
   do ii=3,mmax-2
     dv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
&          -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
   end do
  
!  relativistic coefficient arrays for u (fr) and up (frp).
   do ii=1,mmax
     tfss=fss
     fr(ii)=als*(rr(ii)**2)*(-tfss*(vv(ii)-ee)**2 + 0.5d0*tfss*dv(ii)/ &
&     (rr(ii)*(1.0d0+0.5d0*tfss*(ee-vv(ii)))))
     frp(ii)=-al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0+0.5d0*tfss*(ee-vv(ii)))
   end do
  
! find classical turning point for matching
   mch=0
   do ii=mmax,2,-1
     if(cf(ii-1)<=0.d0 .and. cf(ii)>0.d0) then
       mch=ii
       exit
     end if
   end do

   if(mch==0) then
    ierr=-1
    exit
   end if
  
! start wavefunction with series
  
   do ii=1,4
     uu(ii)=rr(ii)**gamma
     up(ii)=al*gamma*rr(ii)**gamma
     upp(ii)=(al+frp(ii))*up(ii)+(cf(ii)+fr(ii))*uu(ii)
   end do
  
! outward integration using predictor once, corrector
! twice
   node=0
  
   do ii=4,mch-1
     uu(ii+1)=uu(ii)+aeo(up,ii)
     up(ii+1)=up(ii)+aeo(upp,ii)
     do it=1,2
       upp(ii+1)=(al+frp(ii+1))*up(ii+1)+(cf(ii+1)+fr(ii+1))*uu(ii+1)
       up(ii+1)=up(ii)+aio(upp,ii)
       uu(ii+1)=uu(ii)+aio(up,ii)
     end do
     if(uu(ii+1)*uu(ii) .le. 0.0d0) node=node+1
   end do
  
   uout=uu(mch)
   upout=up(mch)
  
  
   if(node-nn+ll+1==0) then
  
! start inward integration at 10*classical turning
! point with simple exponential
  
     nin=mch+2.3d0/al
     if(nin+4>mmax) nin=mmax-4
     xkap=dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee))
  
     do ii=nin,nin+4
       uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       up(ii)=-rr(ii)*al*xkap*uu(ii)
       upp(ii)=(al+frp(ii))*up(ii)+(cf(ii)+fr(ii))*uu(ii)
     end do
  
! integrate inward
  
     do ii=nin,mch+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2
         upp(ii-1)=(al+frp(ii-1))*up(ii-1)+(cf(ii-1)+fr(ii-1))*uu(ii-1)
         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do
  
! scale outside wf for continuity
  
     sc=uout/uu(mch)
  
     do ii=mch,nin
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do
  
     upin=up(mch)
  
! perform normalization sum
  
     ro=rr(1)/dsqrt(amesh)
     sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)
  
     do ii=1,nin-3
       sn=sn+al*rr(ii)*uu(ii)**2
     end do
  
     sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2)**2 &
&              + 28.0d0*rr(nin-1)*uu(nin-1)**2 &
&              +  9.0d0*rr(nin  )*uu(nin  )**2)/24.0d0
  
! normalize u
  
     cn=1.0d0/dsqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin
  
     do ii=1,nin
       up(ii)=cn*up(ii)
       uu(ii)=cn*uu(ii)
     end do
     do ii=nin+1,mmax
       uu(ii)=0.0d0
     end do
  
! perturbation theory for energy shift
  
     de=0.5d0*uout*(upout-upin)/(al*rr(mch))
  
! convergence test and possible exit
  
     if(dabs(de)<dmax1(dabs(ee),0.2d0)*eps) then
       ierr = 0
       exit
     end if
  
     if(de>0.0d0) then 
       emin=ee
     else
       emax=ee
     end if
     ee=ee+de
     if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)
  
   else if(node-nn+ll+1<0) then
! too few nodes
     emin=ee
     ee=0.5d0*(emin+emax)
  
   else
! too many nodes
     emax=ee
     ee=0.5d0*(emin+emax)
   end if
  
 end do

!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
 end if

 deallocate(upp,cf,dv,fr,frp)
 return

 end subroutine lschfb

 function aeo(yy,jj)

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aeo,yy(*)
 integer :: jj

 aeo=(4.16666666667d-2)*(55.0d0*yy(jj)-59.0d0*yy(jj-1) &
& +37.0d0*yy(jj-2)-9.0d0*yy(jj-3))
 return
 end function aeo

 function aio(yy,jj)

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aio,yy(*)
 integer :: jj

 aio=(4.16666666667d-2)*(9.0d0*yy(jj+1)+19.0d0*yy(jj) &
& -5.0d0*yy(jj-1)+yy(jj-2))
 return
 end function aio

 function aei(yy,jj)

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aei,yy(*)
 integer :: jj

 aei=-(4.16666666667d-2)*(55.0d0*yy(jj)-59.0d0*yy(jj+1) &
& +37.0d0*yy(jj+2)-9.0d0*yy(jj+3))
 return
 end function aei

 function aii(yy,jj)

 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aii,yy(*)
 integer :: jj

 aii=-(4.16666666667d-2)*(9.0d0*yy(jj-1)+19.0d0*yy(jj) &
& -5.0d0*yy(jj+1)+yy(jj+2))
 return
 end function aii

end module haesratom
