#include "vdW_params.h"
MODULE vdW_kernel
use iso_c_binding
USE kinds,             ONLY : dp
USE constants

! No implicit variables

IMPLICIT NONE
include 'mpif.h'


! ----------------------------------------------------------------------
! By default everything is private

PRIVATE


! ----------------------------------------------------------------------
! Save all objects in this module

SAVE

PUBLIC  :: generate_vdw_kernel

INTEGER, PARAMETER       :: Nr_points = VDW_NRPOINTS
! The number of radial points (also the number of k points) used in the
! formation of the kernel functions for each pair of q values.
! Increasing this value will help in case you get a run-time error
! saying that you are trying to use a k value that is larger than the
! largest tabulated k point since the largest k point will be 2*pi/r_max
! * Nr_points. Memory usage of the vdW_DF piece of PWSCF will increase
! roughly linearly with this variable.

REAL(DP), PARAMETER      :: r_max     = VDW_RMAX
! The value of the maximum radius to use for the real-space kernel
! functions for each pair of q values. The larger this value is the
! smaller the smallest k value will be since the smallest k point value
! is 2*pi/r_max. Be careful though, since this will also decrease the
! maximum k point value and the vdW_DF code will crash if it encounters
! a g-vector with a magnitude greater than 2*pi/r_max *Nr_points.

REAL(DP), PARAMETER      :: dr = r_max/Nr_points, dk = 2.0D0*pi/VDW_RMAX
! Real space and k-space spacing of grid points.

REAL(DP), PARAMETER      :: q_min = 1.0D-5, q_cut = 5.0D0
! The maximum and minimum values of q. During a vdW run, values of q0
! found larger than q_cut will be saturated (SOLER equation 5) to q_cut.

INTEGER, PARAMETER       :: Nqs    = VDW_NQS
! The above two parameters define the q mesh to be used in the vdW_DF
! code. These are perhaps the most important to have set correctly.
! Increasing the number of q points will DRAMATICALLY increase the
! memory usage of the vdW_DF code because the memory consumption depends
! quadratically on the number of q points in the mesh. Increasing the
! number of q points may increase accuracy of the vdW_DF code, although,
! in testing it was found to have little effect. The largest value of
! the q mesh is q_cut. All values of q0 (DION equation 11) larger than
! this value during a run will be saturated to this value using equation
! 5 of SOLER. In testing, increasing the value of q_cut was found to
! have little impact on the results, although it is possible that in
! some systems it may be more important. Always make sure that the
! variable Nqs is consistent with the number of q points that are
! actually in the variable q_mesh. Also, do not set any q value to 0.
! This will cause an infinity in the Fourier transform.

INTEGER,  PARAMETER      :: Nintegration_points = VDW_INTEGRATION_POINTS

REAL(DP), PARAMETER      :: a_min = 0.0D0, a_max = VDW_AMAX
! Min/max values for the a and b integration in DION equation 14.

REAL(DP) :: kernel( 0:Nr_points, Nqs, Nqs ), d2phi_dk2( 0:Nr_points, Nqs, Nqs )
! Matrices holding the Fourier transformed kernel function  and its
! second derivative for each pair of q values. The ordering is
! kernel(k_point, q1_value, q2_value).

INTEGER                  :: idx, idx1, ierr
! Indexing variable and MPI error return

! ----------------------------------------------------------------------
! General variables

INTEGER                  :: inlc            = 1
! The non-local correlation

INTEGER                  :: vdW_DF_analysis = 0
! vdW-DF analysis tool as described in PRB 97, 085115 (2018)

REAL(DP) :: W_ab( Nintegration_points, Nintegration_points )
! Defined in DION equation 16.

REAL(DP) :: a_points( Nintegration_points ), a_points2( Nintegration_points )
! The values of the "a" points (DION equation 14) and their squares.

CONTAINS

  FUNCTION h_function(y)

     IMPLICIT NONE
     REAL(DP)             :: y, y2, y4, h_function
     REAL(DP), PARAMETER  :: g1 = fpi/9.0D0                                     ! vdW-DF1/2
     REAL(DP), PARAMETER  :: a3 = 0.94950D0, g3 = 1.12D0, g32 = g3*g3           ! vdW-DF3-opt1
     REAL(DP), PARAMETER  :: a4 = 0.28248D0, g4 = 1.29D0, g42 = g4*g4           ! vdW-DF3-opt2
     REAL(DP), PARAMETER  :: a5 = 2.01059D0, b5 = 8.17471D0, g5 = 1.84981D0, &  ! vdW-DF-C6
                             AA = ( b5 + a5*(a5/2.0D0-g5) ) / ( 1.0D0+g5-a5 )   !


     y2 = y*y

     IF ( inlc == 1 .OR. inlc == 2 ) THEN

        h_function = 1.0D0 - EXP( -g1*y2 )

     ELSE IF ( inlc == 3 ) THEN

        y4 = y2*y2
        h_function = 1.0D0 - 1.0D0 / ( 1.0D0 + g3*y2 + g32*y4 + a3*y4*y4 )

     ELSE IF ( inlc == 4 ) THEN

        y4 = y2*y2
        h_function = 1.0D0 - 1.0D0 / ( 1.0D0 + g4*y2 + g42*y4 + a4*y4*y4 )

     ELSE IF ( inlc == 5 ) THEN

        y4 = y2*y2
        h_function = 1.0D0 - ( 1.0D0 + ( (a5-g5)*y2 + AA*y4 ) / ( 1.0D0+AA*y2 ) ) * EXP( -a5*y2 )

     END IF

  END FUNCTION


! ----------------------------------------------------------------------
! Public functions



  ! ####################################################################
  !                           |                 |
  !                           | GENERATE_KERNEL |
  !                           |_________________|
  !
  ! The original definition of the kernel function is given in DION
  ! equations 14-16. The Soler method makes the kernel function a
  ! function of only 1 variable (r) by first putting it in the form
  ! phi(q1*r, q2*r). Then, the q-dependence is removed by expanding the
  ! function in a special way (see SOLER equation 3). This yields a
  ! separate function for each pair of q points that is a function of r
  ! alone. There are (Nqs^2+Nqs)/2 unique functions, where Nqs is the
  ! number of q points used. In the Soler method, the kernel is first
  ! made in the form phi(d1, d2) but this is not done here. It was found
  ! that, with q's chosen judiciously ahead of time, the kernel and the
  ! second derivatives required for interpolation could be tabulated
  ! ahead of time for faster use of the vdW-DF functional. Through
  ! testing we found no need to soften the kernel and correct for this
  ! later (see SOLER eqations 6-7).
  !
  ! The algorithm employed here is "embarrassingly parallel," meaning
  ! that it parallelizes very well up to (Nqs^2+Nqs)/2 processors,
  ! where, again, Nqs is the number of q points chosen. However,
  ! parallelization on this scale is unnecessary. In testing the code
  ! runs in under a minute on 16 Intel Xeon processors.
  !
  ! IMPORTANT NOTICE: Results are very sensitive to compilation details.
  ! In particular, the usage of FMA (Fused Multiply-and-Add)
  ! instructions used by modern CPUs such as AMD Interlagos (Bulldozer)
  ! and Intel Ivy Bridge may affect quite heavily some components of the
  ! kernel (communication by Ake Sandberg, Umea University). In practice
  ! this should not be a problem, since most affected elements are the
  ! less relevant ones.
  !
  ! Some of the algorithms here are somewhat modified versions of those
  ! found in:
  !
  !    Numerical Recipes in C; William H. Press, Brian P. Flannery, Saul
  !    A. Teukolsky, and William T. Vetterling. Cambridge University
  !    Press (1988).
  !
  ! hereafter referred to as NUMERICAL_RECIPES. The routines were
  ! translated to Fortran, of course and variable names are generally
  ! different.
  !
  ! For the calculation of the kernel we have benefited from access to
  ! earlier vdW-DF implementation into PWscf and ABINIT, written by Timo
  ! Thonhauser, Valentino Cooper, and David Langreth. These codes, in
  ! turn, benefited from earlier codes written by Maxime Dion and Henrik
  ! Rydberg.

  SUBROUTINE generate_vdw_kernel(kernel_in, d2phi_dk2_in, q_mesh, intra_image_comm, mpime, nproc)
  use iso_c_binding
  IMPLICIT NONE

  INTEGER  :: a_i, b_i, q1_i, q2_i, r_i
  ! Indexing variables.

  REAL(c_double), INTENT(OUT), dimension((Nr_points+1) * Nqs * Nqs ) :: kernel_in
  REAL(c_double), INTENT(OUT), dimension((Nr_points+1) * Nqs * Nqs ) :: d2phi_dk2_in
  REAL(c_double), INTENT(IN), dimension(Nqs) :: q_mesh
  !REAL(DP), INTENT(INOUT) :: d2phi_dk2_in((Nr_points+1) * Nqs * Nqs )
  !REAL(DP), INTENT(IN)    :: q_mesh(Nqs)
  INTEGER  :: intra_image_comm
  INTEGER  :: mpime
  INTEGER  :: nproc

  REAL(DP) :: weights( Nintegration_points )
  ! Array to hold dx values for the Gaussian-Legendre integration of the kernel.

  REAL(DP) :: sin_a( Nintegration_points ), cos_a( Nintegration_points )
  ! Sine and cosine values of the aforementioned points a.

  REAL(DP) :: d1, d2, d, integral
  ! Intermediate values.

  ! --------------------------------------------------------------------
  ! The following variables control the parallel environment.

  INTEGER :: my_start_q, my_end_q, Ntotal
  ! Starting and ending q value for each  processor, also the total
  ! number of calculations to do, i.e. (Nqs^2 + Nqs)/2.

  REAL(DP), ALLOCATABLE :: phi(:,:), phi_deriv(:,:)
  ! Arrays to store the kernel functions and their second derivatives.
  ! They are stored as phi(radial_point, idx).

  INTEGER, ALLOCATABLE  :: indices(:,:), proc_indices(:,:)
  ! Indices holds the values of q1 and q2 as partitioned out to the
  ! processors. It is an Ntotal x 2 array stored as indices(index of
  ! point number, q1:q2). Proc_indices holds the section of the indices
  ! array that is assigned to each processor. This is a Nproc x 2
  ! array, stored as proc_indices(processor_number,
  ! starting_index:ending_index)

  INTEGER :: Nper, Nextra, start_q, end_q
  ! Baseline number of jobs per processor, number of processors that
  ! get an extra job in case the number of jobs doesn't split evenly
  ! over the number of processors, starting index into the indices
  ! array, ending index into the indices array.

  ! Make these are passed from the C code now
  !INTEGER :: nproc, mpime
  ! Number or procs, rank of current processor.

  INTEGER :: proc_i, my_Nqs

  ! --------------------------------------------------------------------
  ! The total number of phi_alpha_beta functions that have to be
  ! calculated.

  Ntotal = (Nqs**2 + Nqs)/2
  ALLOCATE ( indices(Ntotal, 2) )


  ! --------------------------------------------------------------------
  ! This part fills in the indices array. It just loops through the q1
  ! and q2 values and stores them. Sections of this array will be
  ! assigned to each of the processors later.

  idx = 1

  DO q1_i = 1, Nqs
     DO q2_i = 1, q1_i
        indices(idx, 1) = q1_i
        indices(idx, 2) = q2_i
        idx = idx + 1
     END DO
  END DO


  ! --------------------------------------------------------------------
  ! Figure out the baseline number of functions to be calculated by each
  ! processor and how many processors get one extra job.

!  nproc  = mp_size( intra_image_comm )
!  mpime  = mp_rank( intra_image_comm )
  Nper   = Ntotal/nproc
  Nextra = MOD(Ntotal, nproc)

  ALLOCATE( proc_indices(nproc, 2) )

  start_q = 0
  end_q   = 0


  ! --------------------------------------------------------------------
  ! Loop over all the processors and figure out which section of the
  ! indices array each processor should do. All processors figure this
  ! out for every processor so there is no need to communicate results.

  DO proc_i = 1, nproc

     start_q = end_q + 1
     end_q   = start_q + (Nper - 1)
     IF (proc_i <= Nextra) end_q = end_q + 1

     ! This is to prevent trouble if number of processors exceeds Ntotal.
     IF ( proc_i > Ntotal ) THEN
        start_q    = Ntotal
        end_q      = Ntotal
     END IF

     IF ( proc_i == (mpime+1) ) THEN
        my_start_q = start_q
        my_end_q   = end_q
     END IF

     proc_indices(proc_i, 1) = start_q
     proc_indices(proc_i, 2) = end_q

  END DO


  ! --------------------------------------------------------------------
  ! Store how many jobs are assigned to me.

  my_Nqs    = my_end_q - my_start_q + 1
  ALLOCATE( phi( 0:Nr_points, my_Nqs ), phi_deriv( 0:Nr_points, my_Nqs ) )

  phi       = 0.0D0
  phi_deriv = 0.0D0
  kernel    = 0.0D0
  d2phi_dk2 = 0.0D0


  ! --------------------------------------------------------------------
  ! Find the integration points we are going to use in the
  ! Gaussian-Legendre integration.

  CALL prep_gaussian_quadrature( weights )


  ! --------------------------------------------------------------------
  ! Get a, a^2, sin(a), cos(a) and the weights for the Gaussian-Legendre
  ! integration.

  DO a_i=1, Nintegration_points
     a_points (a_i) = TAN( a_points(a_i) )
     a_points2(a_i) = a_points(a_i)**2
     weights(a_i)   = weights(a_i)*( 1 + a_points2(a_i) )
     cos_a(a_i)     = COS( a_points(a_i) )
     sin_a(a_i)     = SIN( a_points(a_i) )
  END DO


  ! --------------------------------------------------------------------
  ! Calculate the value of the W function defined in DION equation 16
  ! for each value of a and b.

  DO a_i = 1, Nintegration_points
  DO b_i = 1, Nintegration_points
     W_ab(a_i, b_i) = 2.0D0 * weights(a_i)*weights(b_i) * (           &
        (3.0D0-a_points2(a_i))*a_points(b_i) *sin_a(a_i)*cos_a(b_i) + &
        (3.0D0-a_points2(b_i))*a_points(a_i) *cos_a(a_i)*sin_a(b_i) + &
        (a_points2(a_i)+a_points2(b_i)-3.0D0)*sin_a(a_i)*sin_a(b_i) - &
        3.0D0*a_points(a_i)*a_points(b_i)*cos_a(a_i)*cos_a(b_i) )   / &
        (a_points(a_i)*a_points(b_i) )
  END DO
  END DO


  ! --------------------------------------------------------------------
  ! vdW-DF analysis tool as described in PRB 97, 085115 (2018).

  IF      ( vdW_DF_analysis == 1 ) THEN

     DO a_i = 1, Nintegration_points
     DO b_i = 1, Nintegration_points
        W_ab(a_i, b_i) = weights(a_i)*weights(b_i) *                  &
           a_points(a_i)*a_points(b_i)*sin_a(a_i)*sin_a(b_i)
     END DO
     END DO

  ELSE IF ( vdW_DF_analysis == 2 ) THEN

     DO a_i = 1, Nintegration_points
     DO b_i = 1, Nintegration_points
        W_ab(a_i, b_i) = W_ab(a_i, b_i) - weights(a_i)*weights(b_i) *  &
           a_points(a_i)*a_points(b_i)*sin_a(a_i)*sin_a(b_i)
     END DO
     END DO

  END IF


  ! --------------------------------------------------------------------
  ! Now, we loop over all the pairs (q1,q2) that are assigned to us and
  ! perform our calculations.

  DO idx = 1, my_Nqs

     ! -----------------------------------------------------------------
     ! First, get the value of phi(q1*r, q2*r) for each r and the
     ! particular values of q1 and q2 we are using.

     DO r_i = 1, Nr_points
        d1  = q_mesh( indices(idx+my_start_q-1, 1) ) * dr * r_i
        d2  = q_mesh( indices(idx+my_start_q-1, 2) ) * dr * r_i
        phi(r_i, idx) = phi_value(d1, d2)
     END DO


     ! -----------------------------------------------------------------
     ! Now, perform a radial FFT to turn our phi_alpha_beta(r) into
     ! phi_alpha_beta(k) needed for SOLER equation 8.

     CALL radial_fft( phi(:,idx) )


     ! -----------------------------------------------------------------
     ! Determine the spline interpolation coefficients for the Fourier
     ! transformed kernel function.

     CALL set_up_splines( phi(:, idx), phi_deriv(:, idx) )

  END DO


  ! --------------------------------------------------------------------
  ! Finally, we collect the results after letting everybody catch up.
  kernel = 0.0;
  d2phi_dk2 = 0.0;

  CALL MPI_Barrier( intra_image_comm, ierr )
  DO proc_i = 0, nproc-1

!     IF ( proc_i >= Ntotal ) EXIT

!     CALL mp_get ( phi      , phi      , mpime, 0, proc_i, 0, intra_image_comm )
!     CALL mp_get ( phi_deriv, phi_deriv, mpime, 0, proc_i, 0, intra_image_comm )

!     IF ( mpime == 0 ) THEN
     IF ( proc_i == mpime ) THEN

        DO idx = proc_indices(proc_i+1,1), proc_indices(proc_i+1,2)
        !DO idx = proc_indices(mpime+1,1), proc_indices(mpime+1,2)
           q1_i = indices(idx, 1)
           q2_i = indices(idx, 2)
           kernel    (:, q1_i, q2_i) = phi       (:, idx - proc_indices(proc_i+1,1) + 1)
           d2phi_dk2 (:, q1_i, q2_i) = phi_deriv (:, idx - proc_indices(proc_i+1,1) + 1)
           kernel    (:, q2_i, q1_i) = kernel    (:, q1_i, q2_i)
           d2phi_dk2 (:, q2_i, q1_i) = d2phi_dk2 (:, q1_i, q2_i)
        END DO

     END IF

  END DO

  !CALL mp_bcast ( kernel(1,1,1) , 0, intra_image_comm )
  !CALL mp_bcast ( d2phi_dk2, 0, intra_image_comm )
 CALL MPI_ALLREDUCE(MPI_IN_PLACE, kernel, Nqs*Nqs*(Nr_points+1), &
      MPI_DOUBLE, MPI_SUM, intra_image_comm, ierr);
 CALL MPI_ALLREDUCE(MPI_IN_PLACE, d2phi_dk2, Nqs*Nqs*(Nr_points+1), &
      MPI_DOUBLE, MPI_SUM, intra_image_comm, ierr);

  idx1 = 1
  do idx = 0, Nr_points
      do q1_i = 1, Nqs
          do q2_i = 1, Nqs
              kernel_in(idx1) = kernel(idx, q2_i, q1_i)
              d2phi_dk2_in(idx1) = d2phi_dk2(idx, q2_i, q1_i)
              idx1 = idx1 + 1
          end do
      end do
  end do
  ! --------------------------------------------------------------------
  ! Keep the lines below for testing and combatibility with the old
  ! kernel file reading/writing method.
  !
  ! Writing the calculated kernel.
  !
  ! IF ( ionode ) THEN
  !    WRITE(stdout,'(/ / A)') "     vdW-DF kernel table calculated and written to file."
  !    OPEN(UNIT=21, FILE='kernel_table', STATUS='replace', FORM='formatted', ACTION='write')
  !    WRITE(21, '(2i5,f13.8)') Nqs, Nr_points
  !    WRITE(21, '(1p4e23.14)') r_max
  !    WRITE(21, '(1p4e23.14)') q_mesh
  !    DO q1_i = 1, Nqs
  !       DO q2_i = 1, q1_i
  !          WRITE(21, '(1p4e23.14)') kernel(:, q1_i, q2_i)
  !       END DO
  !    END DO
  !    DO q1_i = 1, Nqs
  !       DO q2_i = 1, q1_i
  !          WRITE(21, '(1p4e23.14)') d2phi_dk2(:, q1_i, q2_i)
  !       END DO
  !    END DO
  !    CLOSE (21)
  ! END IF
  !
  !
  ! Reading the kernel from an old kernel file.
  !
  ! IF (ionode) WRITE(stdout,'(/ / A)') "     vdW-DF kernel read from file."
  ! OPEN(UNIT=21, FILE='vdW_kernel_table', STATUS='old', FORM='formatted', ACTION='read')
  ! read(21, '(/ / / / / /)')
  ! DO q1_i = 1, Nqs
  !    DO q2_i = 1, q1_i
  !       READ(21, '(1p4e23.14)') kernel(:, q1_i, q2_i)
  !       kernel(:, q2_i, q1_i) = kernel(:, q1_i, q2_i)
  !    END DO
  ! END DO
  ! DO q1_i = 1, Nqs
  !    DO q2_i = 1, q1_i
  !       READ(21, '(1p4e23.14)')    d2phi_dk2(:, q1_i, q2_i)
  !       d2phi_dk2(:, q2_i, q1_i) = d2phi_dk2(:, q1_i, q2_i)
  !    END DO
  ! END DO
  ! CLOSE (21)


  DEALLOCATE( indices, proc_indices, phi, phi_deriv )

  END SUBROUTINE generate_vdw_kernel








  ! ####################################################################
  !                    |                            |
  !                    |  PREP_GAUSSIAN_QUADRATURE  |
  !                    |____________________________|
  !

  SUBROUTINE prep_gaussian_quadrature( weights )

  !! Routine to calculate the points and weights for the
  !! Gaussian-Legendre integration. This routine is modeled after the
  !! routine GAULEG from NUMERICAL RECIPES.
  
  REAL(DP), INTENT(INOUT) :: weights(:)
  ! The points and weights for the Gaussian-Legendre integration.

  INTEGER  :: Npoints
  ! The number of points we actually have to calculate. The rest will
  ! be obtained from symmetry.

  REAL(DP) :: poly_1, poly_2, poly_3
  ! Temporary storage for Legendre polynomials.

  INTEGER  :: i_point, i_poly
  ! Indexing variables.

  REAL(DP) :: root, dp_dx, last_root
  ! The value of the root of a given Legendre polynomial, the derivative
  ! of the polynomial at that root and the value of the root in the last
  ! iteration (to check for convergence of Newton's method).

  real(dp) :: midpoint, length
  ! The middle of the x-range and the length to that point.




  Npoints  = (Nintegration_points + 1)/2
  midpoint = 0.5D0 * ( ATAN(a_min) + ATAN(a_max) )
  length   = 0.5D0 * ( ATAN(a_max) - ATAN(a_min) )

  DO i_point = 1, Npoints
     ! -----------------------------------------------------------------
     ! Make an initial guess for the root.

     root = COS(DBLE(pi*(i_point - 0.25D0)/(Nintegration_points + 0.5D0)))

     DO
        ! --------------------------------------------------------------
        ! Use the recurrence relations to find the desired polynomial,
        ! evaluated at the approximate root. See NUMERICAL_RECIPES.

        poly_1 = 1.0D0
        poly_2 = 0.0D0

        DO i_poly = 1, Nintegration_points

           poly_3 = poly_2
           poly_2 = poly_1
           poly_1 = ((2.0D0 * i_poly - 1.0D0)*root*poly_2 - (i_poly-1.0D0)*poly_3)/i_poly

        END DO


        ! --------------------------------------------------------------
        ! Use the recurrence relations to find the desired polynomial.
        ! Find the derivative of the polynomial and use it in Newton's
        ! method to refine our guess for the root.

        dp_dx = Nintegration_points * (root*poly_1 - poly_2)/(root**2 - 1.0D0)

        last_root = root
        root      = last_root - poly_1/dp_dx


        ! --------------------------------------------------------------
        ! Check for convergence.

        IF (abs(root - last_root) <= 1.0D-14) EXIT

     END DO


     ! -----------------------------------------------------------------
     ! Fill in the array of evaluation points.

     a_points(i_point) = midpoint - length*root
     a_points(Nintegration_points + 1 - i_point) = midpoint + length*root


     ! -----------------------------------------------------------------
     ! Fill in the array of weights.

     weights(i_point) = 2.0D0 * length/((1.0D0 - root**2)*dp_dx**2)
     weights(Nintegration_points + 1 - i_point) = weights(i_point)

  END DO

  END SUBROUTINE prep_gaussian_quadrature








  ! ####################################################################
  !                            |             |
  !                            |  PHI_VALUE  |
  !                            |_____________|
  !
  
  REAL(DP) FUNCTION phi_value(d1, d2)

  !! This function returns the value of the kernel calculated via DION
  !! equation 14.
  
  REAL(DP), INTENT(IN) :: d1, d2
  ! The point at which to evaluate the kernel. d1 = q1*r and d2 = q2*r.

  REAL(DP) :: w, x, y, z, T
  ! Intermediate values.

  REAL(DP) :: nu(Nintegration_points), nu1(Nintegration_points)
  ! Defined in the discussio below equation 16 of DION.

  INTEGER  :: a_i, b_i
  ! Indexing variables.




  ! --------------------------------------------------------------------
  ! Loop over all integration points and calculate the value of the nu
  ! functions defined in the discussion below equation 16 in DION.

  DO a_i = 1, Nintegration_points
     nu(a_i)  = a_points2(a_i)/( 2.0D0 * h_function( a_points(a_i)/d1 ))
     nu1(a_i) = a_points2(a_i)/( 2.0D0 * h_function( a_points(a_i)/d2 ))
  END DO


  ! --------------------------------------------------------------------
  ! Carry out the integration of DION equation 13.

  phi_value = 0.0D0

  DO a_i = 1, Nintegration_points
     w = nu(a_i)
     y = nu1(a_i)
     DO b_i = 1, Nintegration_points
        x = nu(b_i)
        z = nu1(b_i)
        T = (1.0D0/(w+x) + 1.0D0/(y+z))*(1.0D0/((w+y)*(x+z)) + 1.0D0/((w+z)*(y+x)))
        phi_value = phi_value + T * W_ab(a_i, b_i)
     END DO
  END DO

  phi_value = 1.0D0/pi**2*phi_value

  END FUNCTION phi_value








  ! ####################################################################
  !                            |              |
  !                            |  RADIAL_FFT  |
  !                            |______________|
  !
  
  SUBROUTINE radial_fft(phi)
  
  !! This subroutine performs a radial Fourier transform on the
  !! real-space kernel functions.
  ! Basically, this is just:
  !            int(4*pi*r^2*phi*sin(k*r)/(k*r))dr
  ! integrated from 0 to r_max.
  ! That is, it is the kernel function phi integrated with the 0^th spherical
  ! Bessel function radially, with a 4*pi assumed from angular
  ! integration since we have spherical symmetry. The spherical symmetry
  ! comes in because the kernel function depends only on the magnitude
  ! of the vector between two points. The integration is done using the
  ! trapezoid rule.
  
  REAL(DP), INTENT(INOUT) :: phi(0:Nr_points)
  ! On input holds the real-space function phi_q1_q2(r).
  ! On output hold the reciprocal-space function phi_q1_q2(k).

  REAL(DP) :: phi_k(0:Nr_points)
  ! Temporary storage for phi_q1_q2(k).

  INTEGER  :: k_i, r_i
  ! Indexing variables.

  REAL(DP) :: r, k
  ! The real and reciprocal space points.




  phi_k = 0.0D0

  ! --------------------------------------------------------------------
  ! Handle the k=0 point separately.

  DO r_i = 1, Nr_points
     r        = r_i * dr
     phi_k(0) = phi_k(0) + phi(r_i)*r**2
  END DO


  ! --------------------------------------------------------------------
  ! Subtract half of the last value off because of the trapezoid rule.

  phi_k(0) = phi_k(0) - 0.5D0 * (Nr_points*dr)**2 * phi(Nr_points)


  ! --------------------------------------------------------------------
  ! Integration for the rest of the k-points.

  DO k_i = 1, Nr_points
     k = k_i * dk
     DO r_i = 1, Nr_points
        r          = r_i * dr
        phi_k(k_i) = phi_k(k_i) + phi(r_i) * r * SIN(k*r) / k
     END DO
     phi_k(k_i) = phi_k(k_i) - 0.5D0 * phi(Nr_points) * r * SIN(k*r) / k
  END DO


  ! --------------------------------------------------------------------
  ! Add in the 4*pi and the dr factor for the integration.

  phi = 4.0D0 * pi * phi_k * dr

  END SUBROUTINE radial_fft








  ! ####################################################################
  !                          |                  |
  !                          |  SET UP SPLINES  |
  !                          |__________________|
  !
  
  SUBROUTINE set_up_splines(phi, D2)
  
  !! This subroutine accepts a function (phi) and finds at each point the
  !! second derivative (D2) for use with spline interpolation. This
  !! function assumes we are using the expansion described in SOLER
  !! equation 3.
  ! That is, the derivatives are those needed to interpolate
  ! Kronecker delta functions at each of the q values. Other than some
  ! special modification to speed up the algorithm in our particular
  ! case, this algorithm is taken directly from NUMERICAL_RECIPES.

  REAL(DP), INTENT(IN)    :: phi(0:Nr_points)
  ! The k-space kernel function for a particular q1 and q2.

  REAL(DP), INTENT(INOUT) :: D2(0:Nr_points)
  ! The second derivatives to be used in the interpolation expansion
  ! (SOLER equation 3).

  REAL(DP), ALLOCATABLE   :: temp_array(:)         ! Temporary storage.
  REAL(DP)                :: temp_1, temp_2

  INTEGER  :: r_i
  ! Indexing variable.




  ALLOCATE( temp_array(0:Nr_points) )

  D2         = 0
  temp_array = 0

  DO r_i = 1, Nr_points - 1
     temp_1  = DBLE(r_i - (r_i - 1))/DBLE( (r_i + 1) - (r_i - 1) )
     temp_2  = temp_1 * D2(r_i-1) + 2.0D0
     D2(r_i) = (temp_1 - 1.0D0)/temp_2
     temp_array(r_i) = ( phi(r_i+1) - phi(r_i))/DBLE( dk*((r_i+1) - r_i) ) - &
          ( phi(r_i) - phi(r_i-1))/DBLE( dk*(r_i - (r_i-1)) )
     temp_array(r_i) = (6.0D0*temp_array(r_i)/DBLE( dk*((r_i+1) - (r_i-1)) )-&
          temp_1*temp_array(r_i-1))/temp_2
  END DO

  D2(Nr_points) = 0.0D0
  DO  r_i = Nr_points-1, 0, -1
     D2(r_i) = D2(r_i)*D2(r_i+1) + temp_array(r_i)
  END DO

  DEALLOCATE( temp_array )

  END SUBROUTINE set_up_splines

END MODULE vdw_kernel
