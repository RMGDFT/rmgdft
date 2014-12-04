subroutine diagev (P0,P1,Fdt,N,iprint)
   use utilmod
      implicit none 
      integer :: N, Nsq
      real*8  :: P0(N,N,2),P1(N,N,2), Fdt(N,N)
      real*8,allocatable, dimension(:,:) ::  C 
      real*8,allocatable, dimension(:)   ::  eig,work 
      integer :: Lwork, info,iprint
      integer :: i,j 
      real*8, allocatable, dimension(:,:) :: Ure, Uim,Wre,Wim
      real*8 :: expre, expim,one, zero
      integer :: idamax

      one=1.0d0 
      zero= 0.0d0
   allocate( C(N,N), eig(N), work(5*N*N) )
   C = Fdt

   lwork =5*N*N
   call   DSYEV( 'Vectors', 'Upper', N, C, N, eig, Work, LWORK, INFO ) 
    !write(*,*)  'Eig-max =', eig(idamax(N,eig,1)) 
    if (iprint>0) write(*,1000)( eig(i),i=1,N)
1000 format(10(' ',E12.6))


  deallocate(work)

!   write(*,*)"**** Fdt" ; call  dumpm(Fdt)
!   write(*,*)"**** eig" ; call  dumpm(eig)
!   write(*,*)"**** C  " ; call  dumpm(C)

!
!  time evolution operator
! U  = CC * exp(   F*dt/i*hbar) * CC' 
!    = CC * exp( -i*Fdt/hbar  ) * CC'  = C*expRe*C'  +  i* C*expIm*C'
!   
  allocate(Ure(N,N),Uim(N,N))
  allocate(Wre(N,N),Wim(N,N))
  Ure=0.0d0 
  Uim=0.0d0 
  Wre=0.0d0
  Wim=0.0d0

!
!  W=C*exp  
!
  do i=1,N
    expRe =  cos(-eig(i))
    expIm =  sin(-eig(i))
    call daxpy(N,expRe,C(1,i),1,Wre(1,i),1)   !  Ur(1,i),1)
    call daxpy(N,expIm,C(1,i),1,Wim(1,i),1)   !  Ui(1,i),1)
  enddo

!
!  U =W*C^T
!
  call  dgemm('N','T',n,n,n, one, Wre,n,C,n,zero,Ure,n)
  call  dgemm('N','T',n,n,n, one, Wim,n,C,n,zero,Uim,n)
  
!
! now propagate:  P1 = U* P0 *U' 
!       where U=(Ur+i*Ui)  and  U'=(Ur   -i*Ui^T)
!
!   P1 = (Ur+i*Ui)*(Pr +i*Pi) *( Ur -i*Ui^T) 
!
! dgem: C := alpha*op( A )*op( B ) + beta*C, 
!  DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 
!
  call  dgemm('N','N',n,n,n, one,Ure,n,P0(1,1,1),n,zero,Wre,      n)      !  Wre  =   Ur*P0r 
  call  dgemm('N','N',n,n,n, one,Wre,n,Ure      ,n,zero,P1(1,1,1),n)      !  P1r  =   UrP0r*Ur 
  call  dgemm('N','T',n,n,n,-one,Wre,n,Uim      ,n,zero,P1(1,1,2),n)      !  P1i  = -i*UrP0r*Ur 

  call  dgemm('N','N',n,n,n, one,Uim,n,P0(1,1,1),n,zero,Wim,      n)      !  Wim  =  i*Ui*Pr 
  call  dgemm('N','N',n,n,n, one,Wim,n,Ure      ,n,one, P1(1,1,2),n)      !  P1i +=  i*UiPr*Ur
  call  dgemm('N','T',n,n,n,+one,Wim,n,Uim      ,n,one, P1(1,1,1),n)      !  P1r +=  i*UiPr*(-i*Ui^t) = +UiPr*UiT

  call  dgemm('N','N',n,n,n, one,Ure,n,P0(1,1,2),n,zero,Wim,      n)      !  Wim  =  i*Ur*Pi   
  call  dgemm('N','N',n,n,n, one,Wim,n,Ure      ,n,one ,P1(1,1,2),n)      !  P1i +=  i*UrPi*Ur  
  call  dgemm('N','T',n,n,n,+one,Wim,n,Uim      ,n,one ,P1(1,1,1),n)      !  P1r +=  i*UrPi*(-i*Ui^t) = +UrPi*UiT
  
  call  dgemm('N','N',n,n,n,-one,Uim,n,P0(1,1,2),n,zero,Wre,      n)      !  Wre  = (i*i)*Ui*Pi       =   -UiPi
  call  dgemm('N','N',n,n,n, one,Wre,n,Ure      ,n,one ,P1(1,1,1),n)      !  P1r += (-UiPi)*Ur 
  call  dgemm('N','T',n,n,n,-one,Wre,n,Uim      ,n,one ,P1(1,1,2),n)      !  P1i += (-UiPi)*(-i*Ui^t)


! write(*,*)'***FINAL: P1'; call dumpm(P1)     


end
