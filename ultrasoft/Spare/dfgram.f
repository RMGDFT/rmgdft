c************************** SVN Revision Information **************************
c **    $Id$    **
c******************************************************************************
 
 
      subroutine dfgram (c, vel, numst, maxst, numpt, maxpt, dr)
      implicit none
 
c ## gram-schmidt orthogonalization of complex wavefunctions
c ## only for use on the T3E
 
      integer numst, maxst
      integer numpt, maxpt
      double complex c(maxpt, numst), ct1, zr1
      real*8 dznrm2
      double complex zdotc
      real*8 vel, rt1
 
      integer st, st1
      integer info
      integer length, idx, idj
      real*8 tmp
      double complex dr(maxst)
      double complex onec, nonec, zeroc
 
      onec = dcmplx(1.0 , 0.0)
      nonec = cmplx(-1.0 , 0.0)
      zeroc = dcmplx(0.0 , 0.0)
 
c ##########
 
c ## compute the lower-triangular part of the overlap matrix
 
 
      do st = 1, numst
 
c ##     normalize this wavefunction
         rt1 = dznrm2(numpt, c(1, st), 1)
         rt1 = rt1 * rt1 * vel
 
c ##     next we have to sum over all processors           
         call global_sums(rt1, 1) 
         rt1 = 1.0 / dsqrt(rt1)
 
c ##     scale the eigenfunction
         call zdscal(numpt, rt1, c(1, st), 1)
 
         do st1 = st + 1, numst
c ##       compute the projection along the remaining vectors
 
           dr(st1) = zdotc(numpt, c(1, st), 1, c(1, st1), 1)
 
         enddo

         call global_sums(dr, 2*maxst)

         zr1 = cmplx(-vel, 0.0)
         call zdscal(maxst, zr1, dr, 1)

                                                                                 
         do st1 = st + 1, numst
 
            call zaxpy(numpt, dr(st1), c(1, st), 1, c(1, st1), 1)
 
         enddo
 
      enddo                     ! st
 
 
      end                                                                       
