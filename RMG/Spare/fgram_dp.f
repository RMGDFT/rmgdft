c************************** SVN Revision Information **************************
c **    $Id$    **
c******************************************************************************
 
 
 
      subroutine fgram (c, vel, numst, maxst, numpt, maxpt, dr)
      implicit none
 
c ## gram-schmidt orthogonalization of complex wavefunctions
c ## only for use on the T3E
 
      integer numst, maxst
      integer numpt, maxpt
      complex c(maxpt, numst), ct1
      complex scnrm2, cdotc
      real vel, rt1
 
      integer st, st1
      integer info
      integer length, idx, idj
      real tmp
      complex dr(maxst)
      complex onec, nonec, zeroc
 
      onec = cmplx(1.0 , 0.0)
      nonec = cmplx(-1.0 , 0.0)
      zeroc = cmplx(0.0 , 0.0)
 
c ##########
 
c ## compute the lower-triangular part of the overlap matrix
 
 
      do st = 1, numst
 
c ##     normalize this wavefunction
         rt1 = scnrm2(numpt, c(1, st), 1)
         rt1 = rt1 * rt1 * vel
 
c ##     next we have to sum over all processors           
         dr(1) = rt1
         call global_sums(dr(1), 1) 
         rt1 = 1.0 / sqrt(dr(1))
 
c ##     scale the eigenfunction
         call csscal(numpt, rt1, c(1, st), 1)
 
         do st1 = st + 1, numst
c ##       compute the projection along the remaining vectors
 
           dr(st1) = cdotc(numpt, c(1, st), 1, c(1, st1), 1)
 
           dr(st1) = -vel * dr(st1)
 
         enddo
 
         call global_sums(dr, 2*maxst)
                                                                                 
         do st1 = st + 1, numst
 
            call caxpy(numpt, dr(st1), c(1, st), 1, c(1, st1), 1)
 
         enddo
 
      enddo                     ! st
 
 
      end                                                                       
