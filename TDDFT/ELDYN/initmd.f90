module initmd 
contains
subroutine init(S,P,F,N)
 use utilmod
 implicit none
 real*8  :: S(n,n), F(n,n),  P(n,n,2)  
 real*8, allocatable, dimension(:,:)   :: AA,B, SS
 integer :: N,Nsq,i
 real*8  :: tmp

 
  Nsq =n*n

  allocate(AA(n,n),B(n,n),SS(n,n))

!
!   overlap
!
   call random_number(S)
   B =  S -0.5d0
   tmp=1.0d0/N
   write(*,*) tmp
   S=matmul(B,transpose(B)) *tmp
   do i=1,N
     S(i,i) = S(i,i) + 1.0d0
   enddo

!
! symmetric   part P  
   call random_number(SS)
   B =  SS -0.5d0
   tmp=1.0d0/N
   write(*,*) tmp
   SS=matmul(B,transpose(B)) *tmp
   do i=1,N
     SS(i,i) = SS(i,i) + 1.0d0
   enddo


!
!  antisymmetric part of P
   call random_number(AA)
   B =  AA -0.5d0
   AA =B-transpose(B)

  call dcopy(Nsq,SS,1,P(1,1,1),1)
  call dcopy(Nsq,AA,1,P(1,1,2),1)

!
! random  F
!
   call random_number(F)
   B =  F -0.5d0
   F= B+transpose(B)
   
  F = F /sqrt(25.0d0*real(N))


 
 end subroutine init
end module  initmd
