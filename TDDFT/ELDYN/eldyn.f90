

!-----------------------------------------------

subroutine eldyn(N,S,F,Pn0,Pn1,Ieldyn,iprint)
use  timing
use  utilmod
use  initmd

!  S,Pn0,Pn1, F
implicit none 
integer  ::  N, iprint, Nsq, Maxiter,niter, ieldyn
real*8 :: S(N,N),F(N,N), Pn0(N,N,2), Pn1(N,N,2)
!real*8, allocatable, dimension(:,:)   :: S, Si, Fn, Fo, U, Ui, F
!real*8, allocatable, dimension(:,:,:) :: Pn0, Pn1, Po0,Po1, W, P
real*8, allocatable, dimension(:,:)   :: Si, Fn, Fo, U, Ui 
real*8, allocatable, dimension(:,:,:) :: Po0,Po1, W, P
real*8 ::   thrs, errmax
integer :: iPn2Po, iPo2Pn, iFo2Fn, iFn2Fo
real :: t0,t1
logical  :: tCommutp,  tDiagp

write(*,*) "enter eldyn"

!
! select a propagator
! 
  if (Ieldyn.eq.0) then 
     tCommutp =.true.
     tDiagp   =.true.
  else if( Ieldyn.eq.1) then 
     tCommutp =.true.
     tDiagp   =.false.
  else if( Ieldyn.eq.1) then 
     tCommutp =.false.
     tDiagp   =.true.
  else  
    write(*,*)' Incorrect ieldyn. Should  be 0, 1 or 2. Exiting.'
    stop
  endif
         
!
! set threshold and  max number of  iteration  for commutator expansion
!

!write(*,*)'-- eldyn -1'
  Nsq= N*N
  thrs =  1.0d-8
  Maxiter = 100

!
!  generate  matrices
  call start_time(1)
  call start_time(2)
  
  allocate(Si(n,n),U(n,n), Ui(n,n))
  allocate(Fo(n,n),Po0(n,n,2),Po1(n,n,2))
!  call init(S,P,F,n)
  call stop_time(2)

!
! set o/n transformation flags:
    iPo2Pn = 1  ; iPn2Po = 2  ; iFo2Fn = 3  ; iFn2Fo = 4

!
!  cholesky decomposition and  O/N transformation matrices U,Ui
!
 
!write(*,*)'-- eldyn -2'
   call dochol(S,Si,U,Ui,N)   ! includes timing inside dochol (split on Sinverse and U,Ui

!
!  O/N transformation Pn-Po
!
!write(*,*)'-- eldyn -3'
   call start_time(5)
   call  ontrans(U,Ui,F         , Fo        ,n,iFn2Fo)
   call  ontrans(U,Ui,Pn0(1,1,1), Po0(1,1,1),n,iPn2Po)
   call  ontrans(U,Ui,Pn0(1,1,2), Po0(1,1,2),n,iPn2Po)
   call stop_time(5)


!write(*,*)'-- eldyn -4'

   if(iprint.gt.0) then ; 
     write(*,*)'*******  S :' ; call dumpm(S,iprint)  !iprint)
     write(*,*)'******* Pn0:' ; call dumpm(Pn0,iprint)
     write(*,*)'******* F  :' ; call dumpm(F,iprint)
     write(*,*)'***chol: U :' ; call dumpm(U,iprint)
     write(*,*)'***chol: Ui:' ; call dumpm(Ui,iprint)
     write(*,*)'***chol: Si:' ; call dumpm(Si,iprint)
     write(*,*)'***chol: Fo:' ; call dumpm(Fo,iprint)
     write(*,*)'***chol: Po0:'; call dumpm(Po0,iprint)

     endif


!
!  now commut propagation
!
!write(*,*)'-- eldyn -5'
   if (tCommutp) then 
     call  start_time(6)
     call commutp (Po0,Po1,Fo,N,thrs,maxiter,niter,errmax,iprint)
     call  stop_time(6)
     if (iprint.gt.0)  then 
        write(*,*)"=============================Commut expansion:========"
        write(*,*)"**** completed,  Po1:" ;  call dumpm(Po1,iprint)
        endif
    endif 
     
!write(*,*)'-- eldyn -6'
!
!  now diag propagation
!
   if (tDiagp) then 
      call start_time(7)
      call diagev (Po0,Po1,Fo,N,iprint)
      call stop_time(7)
      if (iprint.gt.0)  then 
         write(*,*)"=============================diagev ================="
         write(*,*)"**** completed,  Po1:" ;  call dumpm(Po1,iprint)
         endif
   endif 

!write(*,*)'-- eldyn -7'
!
!--  now back n/o transformation: Po->Pn
!
   call  start_time(5)
!   call  ontrans(U,Ui,Fo       , Fn        ,n,iFo2Fn)
   call  ontrans(U,Ui,Po1(1,1,1), Pn1(1,1,1),n,iPo2Pn)
   call  ontrans(U,Ui,Po1(1,1,2), Pn1(1,1,2),n,iPo2Pn)
   call  stop_time(5)

   if (iprint.gt.0)  then 
      write(*,*)"**** o/n transformed,  Pn1:" ;  call dumpm(Pn1,iprint)
      endif

!write(*,*)'-- eldyn -8'

!----- stop timer
   call stop_time(1)
end subroutine

