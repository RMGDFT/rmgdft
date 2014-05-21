subroutine commutp (P0,P1,Om,N,thrs,maxiter,niter,errmax,iprint)
!use utilmod
      implicit none 
      integer :: N, Nsq,Nsq2
      real*8  :: P0(N,N,2),P1(N,N,2), Om(N,N)
      real*8,allocatable, dimension(:,:,:) :: dP, C 
      real*8  :: alpha, thrs,err ,errmax
      integer :: iter, maxiter, ierr ,i,j,niter,iprint
      logical ::  tConv
      
!
! P = P0 +(dt/ih) [F,P]+ (dt/ih)^2/2! [F[F,P] + (dt/ih)^3/3! [F,F[F,P]] +... 
!
! P(k)  <-- P(k-1) + 1/(ihk) *[Om,dP(k-1)] = 
!         = P(k-1) +1/(hk)  *[Om,A ] -i/(hk) *[Om,S ]  
!   
!   where dP = S+i*A 
!  
!  Om = (Fav*dt)/I
!

      allocate (dP(N,N,2), C(N,N,2))
      
!      write(*,*)'.......Entering commut:'
!      write(*,*)'N,thrs,maxiter =',N,thrs,maxiter
!      write(*,*)'****** P0 = '; call dumpm(P0)
!      write(*,*)'****** Om = '; call dumpm(Om)
     
      P1=P0
      dP=P0
       
      alpha=1.0d0
      iter=1
      Nsq = N*N
      Nsq2 =2*Nsq
      errmax =0.d0

      lp_bch: do while (iter  <= maxiter) 
        alpha = 1.0d0 /iter  
        call  commatr(Om,dP(1,1,2),C(1,1,1),alpha,N)   ! real contrib correction 
        call  commstr(Om,dP(1,1,1),C(1,1,2),-alpha,N)  ! imag contrib correction
        dP =C 
        P1 =P1 +dP
!--v1:
!         ierr= idamax(Nsq2,dP,1)
!         err = dP(ierr) 
!         if (abs(err) <= thrs)  exit lp_bch
!-v2:
        call tstconv(dP,2*Nsq,N,thrs,ierr,err,tconv)
         errmax=max(errmax, abs(err))
        if (tConv) exit lp_bch 
        iter= iter +1 
      end do lp_bch
      write(*,*)'Converged: Niter, errmax =',  iter,errmax

      niter= iter
      if(iprint >0) write(*,*)'bch:  errmax= ', errmax,  ' Niter =',iter


end 
!c-------TSTSCONV
subroutine  tstconv(C,M,N,thrs,ierr,err,tconv)
      implicit  none
      real*8,intent(in)   :: C(M)
      integer  :: N,M, ix,iy,im,inc,ierr
      real*8   :: err ,thrs
      logical  ::  tconv
      integer :: idamax
      inc=1
      ierr = idamax(M,C,inc)
      err= C(ierr)
      im = ierr/(N*N)+1
      !Iy = ierr/N    +1    -(Im-1)*N
      !Ix = modulo(ierr,N)
      Iy = (ierr-1)/N    +1    -(Im-1)*N
      Ix = modulo(ierr-1,N)
      write(*,1000)ix,iy,im,err,thrs
1000  format('Err: dP(i,j,k)=',I6,I6,I6, E16.6,'   Thrs =',E12.4)
      tconv=.false.
      if (abs(err) <= thrs) tconv= .true.
end subroutine tstconv

!----------
!
!
!
subroutine commstr(A,B,C,alpha,Nb) 
! [S1,S2] =S1*S2-*S2*S1= S1*S2 -transpose(S1*S2)
        implicit none 
        real*8  :: A(Nb,Nb), B(Nb,Nb),C(Nb,Nb)
        real*8, allocatable,dimension(:,:) ::  W 
        real*8 ::  alpha , beta 
        integer :: Nb
        !alpha=1.0d0 
        beta=0.0d0
        allocate( W(Nb,Nb))
        call dgemm('N','N',Nb,Nb,Nb,alpha,A,Nb,B,Nb,beta,C,Nb)   
        !call transpose(C,W,Nb)
        W =transpose(C)
        C =C -W 
        deallocate (W)

end subroutine  commstr
!
!
!
subroutine commatr(A,B,C,alpha,Nb) 
!  [S,A] = S*A-A*S =  S*A +transpose(S*A)
        implicit none 
        real*8  :: A(Nb,Nb), B(Nb,Nb),C(Nb,Nb)
        real*8, allocatable,dimension(:,:) ::  W 
        real*8 ::  alpha , beta 
        integer :: Nb
        !alpha=1.0d0 
        beta=0.0d0
        allocate( W(Nb,Nb))
        call dgemm('N','N',Nb,Nb,Nb,alpha,A,Nb,B,Nb,beta,C,Nb)   
        !call transpose(C,W,Nb)
        W =transpose(C)
        C =C +W 
        deallocate (W)

end subroutine commatr
