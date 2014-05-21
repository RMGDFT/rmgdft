!--------------------------
module utilmod 
    implicit none
interface  dumpm
   module procedure dumpm_d      ! print matrix
   module procedure dumpm_s
   module procedure dumpm2_d
   module procedure dumpm2_s
   module procedure dumpv_d     ! print vector
   module procedure dumpv_s
   module procedure dumpv2_d
   module procedure dumpv2_s
   module procedure dump2m_d    !print 2 or more matrices (whole)
   module procedure dump2m2_d   !print 2 or more matrices (block only)
   module procedure dump2m_s   
   module procedure dump2m2_s  
   
end  interface dumpm

 
contains 
  
  !@
  !@  dump  vector  A, double precision
  !@
  subroutine   dumpv_d(A) 
    integer :: i,i0,i1,N 
    double precision, dimension(:) :: A
    N=size(A)
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = ubound(A,1) 
    write(*,100) (A(i),i=i0,i1) 
   !write(*,*)  'size', size(A)
  end subroutine   dumpv_d 

  
  !@
  !@  dump  vector  A, single  precision
  !@
  subroutine   dumpv_s(A) 
    integer :: i,i0,i1,N 
    real, dimension(:) :: A
    N=size(A)
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = ubound(A,1) 
    write(*,100) (A(i),i=i0,i1) 
   !write(*,*)  'size', size(A)
  end subroutine   dumpv_s 

  
  !@
  !@  dump  vector  A, double precision
  !@
  subroutine   dumpv2_d(A,m) 
    integer :: i,i0,i1,N,m 
    double precision, dimension(:) :: A
    N=size(A)
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = min(ubound(A,1),i0+m-1) 
    write(*,100) (A(i),i=i0,i1) 
   !write(*,*)  'size', size(A)
  end subroutine   dumpv2_d 

  
  !@
  !@  dump  vector  A, single  precision
  !@
  subroutine   dumpv2_s(A,m) 
    integer :: i,i0,i1,N,m
    real, dimension(:) :: A
    N=size(A)
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = min(ubound(A,1),i0+m-1) 
    write(*,100) (A(i),i=i0,i1) 
   !write(*,*)  'size', size(A)
  end subroutine   dumpv2_s 





  !@
  !@  dump  whole matrix A, double precision
  !@
  subroutine   dumpm_d(A) 
    integer :: i,i0,i1, j,j0,j1 
    double precision, dimension(:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = ubound(A,1) 
    j0 = lbound(A,2)  ; j1 = ubound(A,2) 
    do i=i0,i1 
      write(*,100) (A(i,j),j=j0,j1) 
     enddo 
   !write(*,*)  'size', size(A)
  end subroutine   dumpm_d 

  !@
  !@  dump  whole matrix A, single precision
  !@
  subroutine   dumpm_s(A) 
    integer :: i,i0,i1, j,j0,j1 
    real, dimension(:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = ubound(A,1) 
    j0 = lbound(A,2)  ; j1 = ubound(A,2) 
    do i=i0,i1 
      write(*,100) (A(i,j),j=j0,j1) 
     enddo 
   !write(*,*)  'size', size(A)
  end subroutine   dumpm_s  

  !@
  !@ dump a small block size (m*m)  of large A  , double precision
  !@   
  subroutine   dumpm2_d(A,m) 
    integer :: i,i0,i1, j,j0,j1,m 
    double precision, dimension(:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = min(ubound(A,1), i0+m-1) 
    j0 = lbound(A,2)  ; j1 = min(ubound(A,2), j0+m-1) 
    do i=i0,i1 
      write(*,100) (A(i,j),j=j0,j1) 
     enddo 
   !write(*,*)  'size', size(A)
  end subroutine   dumpm2_d 

  !@
  !@ dump a small block size m*m  of large A  , single precision
  !@   
  subroutine   dumpm2_s(A,m) 
    integer :: i,i0,i1, j,j0,j1,m 
    real, dimension(:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = min(ubound(A,1), i0+m-1) 
    j0 = lbound(A,2)  ; j1 = min(ubound(A,2), j0+m-1) 
    do i=i0,i1 
      write(*,100) (A(i,j),j=j0,j1) 
     enddo 
   !write(*,*)  'size', size(A)
  end subroutine   dumpm2_s  


  !@
  !@  dump  whole matrix A, double precision
  !@  print several matrices (supposedly 2)
  !@ 
  subroutine   dump2m_d(A) 
    integer :: i,i0,i1, j,j0,j1, k,k0,k1 
    double precision, dimension(:,:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = ubound(A,1) 
    j0 = lbound(A,2)  ; j1 = ubound(A,2) 
    k0 = lbound(A,3)  ; k1 = ubound(A,3) 
    do k=k0,k1
    do i=i0,i1 
      write(*,100) (A(i,j,k),j=j0,j1) 
     enddo 
     write(*,*)
   enddo
   !write(*,*)  'size', size(A)
  end subroutine   dump2m_d 


  !@
  !@ dump a small block size (m*m)  of large A  , double precision
  !@  print several matrices (supposedly 2)
  !@   
  subroutine   dump2m2_d(A,m) 
    integer :: i,i0,i1, j,j0,j1,k,k0,k1,m 
    double precision, dimension(:,:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = min(ubound(A,1), i0+m-1) 
    j0 = lbound(A,2)  ; j1 = min(ubound(A,2), j0+m-1) 
    k0 = lbound(A,3)  ; k1 = ubound(A,3) 
    do k=k0,k1
    do i=i0,i1 
      write(*,100) (A(i,j,k),j=j0,j1) 
     enddo 
     write(*,*)
    enddo
   !write(*,*)  'size', size(A)
  end subroutine   dump2m2_d 

!--

  !@
  !@  dump  whole matrix A, double precision
  !@  print several matrices (supposedly 2)
  !@ 
  subroutine   dump2m_s(A) 
    integer :: i,i0,i1, j,j0,j1, k,k0,k1 
    real,  dimension(:,:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = ubound(A,1) 
    j0 = lbound(A,2)  ; j1 = ubound(A,2) 
    k0 = lbound(A,3)  ; k1 = ubound(A,3) 
    do k=k0,k1
    do i=i0,i1 
      write(*,100) (A(i,j,k),j=j0,j1) 
     enddo 
     write(*,*)
   enddo
   !write(*,*)  'size', size(A)
  end subroutine   dump2m_s 


  !@
  !@ dump a small block size (m*m)  of large A  , double precision
  !@  print several matrices (supposedly 2)
  !@   
  subroutine   dump2m2_s(A,m) 
    integer :: i,i0,i1, j,j0,j1,k,k0,k1,m 
    real, dimension(:,:,:) :: A
    100 format (100(E16.8,' '))
    i0 = lbound(A,1)  ; i1 = min(ubound(A,1), i0+m-1) 
    j0 = lbound(A,2)  ; j1 = min(ubound(A,2), j0+m-1) 
    k0 = lbound(A,3)  ; k1 = ubound(A,3) 
    do k=k0,k1
    do i=i0,i1 
      write(*,100) (A(i,j,k),j=j0,j1) 
     enddo 
     write(*,*)
    enddo
   !write(*,*)  'size', size(A)
  end subroutine   dump2m2_s 
!----


  end module  utilmod 

