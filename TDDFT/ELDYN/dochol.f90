subroutine dochol(S,Si,U,Ui,N)
use utilmod
use timing
implicit none
  real*8 :: S(n,n),Si(n,n),U(n,n), Ui(n,n)
  real*8, allocatable,dimension(:,:) ::  W
  integer  :: n, nsq,info,lt
  real :: t0,t1
  allocate (W(n,n))

  Nsq= N*N
  lt =1 ; 

!
!-- get Si  directly from  U (wihtout Ui)
!
!  write(*,*)' dochol -0 '
!  write(*,*)'Nsq,N = ',Nsq, N
!  write(*,*)'**:S'
!  write(*,100) s
100 format(9f8.4)
  !call  dumpm(S)
  !write(*,*)'**:Si',Si
  
  call   dcopy(Nsq,S,1,Si,1)
! write(*,*)' dochol -1 '
  call   start_time(3)
  call   dpotrf('U',N,Si,N,info)   !  cholesky 
! write(*,*)' dochol -2 '
  if (info.ne.0)  then ; write(*,*) 'Chol-1-err';stop; endif
  call   dpotri('U',N,Si,N,info)   ! inverse  of S using  cholesky 
  if (info.ne.0)  then ; write(*,*) 'Chol-2-err';stop; endif
  call  filllt(Si,N,lt)
! write(*,*)' dochol -3 '
  call  stop_time(3)
!  write(*,*)' **** chol: Si='  ;  call dumpm(Si)

!
!-- transf.  matrices
!
  
  call  start_time(4)
  call  dcopy(Nsq,S,1,U,1)
! write(*,*)' dochol -4 '
  call  dpotrf('U',N,U,N,info)   !  cholesky 
  if (info.ne.0)  then ; write(*,*) 'Chol-3-err';stop; endif
  call  dcopy(Nsq,U,1,Ui,1) 
  call  dtrtri('U','N',n,Ui,n,info)
  if (info.ne.0)  then ; write(*,*) 'Chol-4-err';stop; endif
  call  cleanlt(U,N,lt)
  call  cleanlt(Ui,N,lt)
  call  stop_time(4)
!  write(*,*)' **** chol: U =' ;  call dumpm(U)
!  write(*,*)' **** chol: Ui=' ;  call dumpm(Ui)


!  dpptrf
!  dtptri 

! write(*,*)' dochol -5 '

end subroutine

!---------------------------------       
subroutine cleanLT(A,N,lt)
implicit none
 real *8 :: A(N,N), zero
 data zero /0.0d0/
 integer :: i,j,N,lt
 if (lt.eq.1) then
   do i=1,N
     do j=1,i-1
       A(i,j) = zero
      enddo
   enddo
 elseif (lt.eq.2) then
   do i=1,N
     do j=1,i-1
         A(j,i) = zero
      enddo
    enddo
 endif
end  subroutine cleanlt
!---------------------------------       
subroutine fillLT(A,N,lt)
!  lt=1 : fill L part by copying from U
!  lt=2 : fill U part by copying from L
implicit none
 real *8 :: A(N,N), zero
 data zero /0.0d0/
 integer :: i,j,N,lt
 if (lt.eq.1) then
   do i=1,N
     do j=1,i-1
       A(i,j) = A(j,i)
      enddo
   enddo
 elseif (lt.eq.2) then
   do i=1,N
     do j=1,i-1
         A(j,i) = A(i,j)
      enddo
    enddo
 endif
end  subroutine filllt

subroutine  ontrans (U,Ui, Ain, Aout,N,ifon)
 !-------------------------------------------------
 !    Sn = L*U        ;  Sni = Ui*Li
 !    So = Li*Sn*Ui   ; 
 !    Sn = L* I* U    ; 
 !    Fn = L*Fo* U    ;  Fo = Li*Fn*Ui  
 !    Pn = Ui*Po*Li   ;  Po = U *Pn*L
 !------------------------------------------------
 ! ifon=1 :   Po->  Pn       :     Ui * Po * (Ui^t)
 !      2 :   Pn->  Po       :     U  * Pn * (U^t)
 !      3 :   Fo->  Fn       :  (U^t) * Fo *  U
 !      4 :   Fn->  Fo       : (Ui^t) * Fn *  Ui
 !--------------------------------------------------
  implicit none
  integer, intent(in)                 :: ifon,N
  real*8, dimension(n,n), intent(in)  :: U, Ui
  real*8  :: Aout(*),Ain(*)
  real*8, allocatable, dimension(:,:) :: W
  integer ::  Nsq ,i 
  real*8 :: zero,one 


  zero = 0.0d0 
  one  = 1.0d0

  Nsq=N*N

  allocate ( W(n,n) )

  if (ifon.eq.1)  then              !  Po->Pn  :  Pn= Ui * Po * (Ui^t)
    call  dgemm('N','N',n,n,n,one,Ui,n,Ain,n,zero,W   ,n)    
    call  dgemm('N','T',n,n,n,one,W ,n,Ui ,n,zero,Aout,n)    

  elseif (ifon.eq.2) then           !  Pn->Po  :  Po= U  * Pn * (U^t)
    call  dgemm('N','N',n,n,n,one,U ,n,Ain,n,zero,W   ,n)    
    call  dgemm('N','T',n,n,n,one,W ,n,U  ,n,zero,Aout,n)    
   
  elseif (ifon.eq.3) then           !  Fo->Fn  :  Fn= (U^t) * Fo *  U 
    call  dgemm('T','N',n,n,n,one,U ,n,Ain,n,zero,W   ,n)    
    call  dgemm('N','N',n,n,n,one,W ,n,U  ,n,zero,Aout,n)    
   
  elseif (ifon.eq.4) then           !  Fn->Fo  :  Fo= (Ui^t) * Fn *  Ui
    call  dgemm('T','N',n,n,n,one,Ui,n,Ain,n,zero,W   ,n)    
    call  dgemm('N','N',n,n,n,one,W ,n,Ui ,n,zero,Aout,n)    
100   format(6(' ', E16.8))
  else
    write(*,*) 'Wrong  o/n transformation'
    stop
  endif  
  deallocate(W) 
end
