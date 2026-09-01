! Copyright (C) 2001-2009 Quantum ESPRESSO group
! Copyright (C) 2015 Brian Kolb, Timo Thonhauser
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! ----------------------------------------------------------------------


module vdw_splines
use iso_c_binding
implicit none
public  :: spline_interpolation, initialize_spline_interpolation

CONTAINS

  ! ####################################################################
  !                       |                        |
  !                       |  SPLINE_INTERPOLATION  |
  !                       |________________________|
  !
  ! This routine is modeled after an algorithm from "Numerical Recipes
  ! in C" by Cambridge University press, page 97.  It was adapted for
  ! Fortran, of course and for the problem at hand, in that it finds the
  ! bin a particular x value is in and then loops over all the P_i
  ! functions so we only have to find the bin once.

  SUBROUTINE spline_interpolation (x, Nx, evaluation_points, Ngrid_points, values, d2y_dx2)

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
  INTEGER, PARAMETER :: i4b = selected_int_kind(9)

  real(dp), intent(in) :: x(1:Nx), evaluation_points(1:Ngrid_points)     ! Input variables. The x values used to
                                                         ! form the interpolation (q_mesh in this
                                                         ! case) and the values of q0 for which we
                                                         ! are interpolating the function.

  complex(dp), intent(inout) :: values(1:Ngrid_points*Nx)              ! An output array (allocated outside this
                                                         ! routine) that stores the interpolated
                                                         ! values of the P_i (SOLER equation 3)
                                                         ! polynomials. The format is
                                                         ! values(grid_point, P_i).

  integer :: Ngrid_points, Nx                            ! Total number of grid points to evaluate
                                                         ! and input x points.

  real(dp), intent(in) :: d2y_dx2(1:Nx*Nx)                   ! The second derivatives required to do
                                                         ! the interpolation.

  integer :: i_grid, lower_bound, upper_bound, idx, P_i  ! Some indexing variables.

  real(dp), allocatable :: y(:)                          ! Temporary variables needed for the
  real(dp) :: a, b, c, d, dx                             ! interpolation.



  ! --------------------------------------------------------------------
  ! Allocate the temporary array.
!WRITE(6,'(5X,A,I3,A,I6,A,F8.3)' ) "Nx    = ",Nx,"  Npoints = ", Ngrid_points," X[1] = ",x(1)
  allocate( y(Nx) )


  ! --------------------------------------------------------------------
  ! If this is the first time this routine has been called we need to
  ! get the second derivatives (d2y_dx2) required to perform the
  ! interpolations. So we allocate the array and call
  ! initialize_spline_interpolation to get d2y_dx2.

  do i_grid=1, Ngrid_points

     lower_bound = 1
     upper_bound = Nx

     do while ( (upper_bound - lower_bound) > 1 )

        idx = (upper_bound+lower_bound) / 2

        if ( evaluation_points(i_grid) > x(idx) ) then
           lower_bound = idx
        else
           upper_bound = idx
        end if

     end do

     dx = x(upper_bound)-x(lower_bound)

     a = (x(upper_bound) - evaluation_points(i_grid))/dx
     b = (evaluation_points(i_grid) - x(lower_bound))/dx
     c = ((a**3-a)*dx**2)/6.0D0
     d = ((b**3-b)*dx**2)/6.0D0

     do P_i = 1, Nx

        y = 0
        y(P_i) = 1

        !values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
        values((P_i-1)*Ngrid_points + i_grid) = a*y(lower_bound) + b*y(upper_bound) &
             + (c*d2y_dx2(P_i + (lower_bound-1)*Nx) + d*d2y_dx2(P_i +  (upper_bound-1)*Nx))
     end do

  end do

  deallocate( y )

  END SUBROUTINE spline_interpolation







  ! ####################################################################
  !                  |                                   |
  !                  |  INITIALIZE_SPLINE_INTERPOLATION  |
  !                  |___________________________________|
  !
  ! This routine is modeled after an algorithm from "Numerical Recipes
  ! in C" by Cambridge University Press, pages 96-97. It was adapted
  ! for Fortran and for the problem at hand.

  SUBROUTINE initialize_spline_interpolation (x, Nx, d2y_dx2)

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
  INTEGER, PARAMETER :: i4b = selected_int_kind(9)


  real(dp), intent(in)  :: x(:)                 ! The input abscissa values.
  real(dp), intent(inout) :: d2y_dx2(:,:)       ! The output array (allocated outside this routine)
                                                ! that holds the second derivatives required for
                                                ! interpolating the function.

  integer :: Nx, P_i, idx                       ! The total number of x points and some indexing
                                                ! variables.

  real(dp), allocatable :: temp_array(:), y(:)  ! Some temporary arrays required. y is the array
                                                ! that holds the funcion values (all either 0 or
                                                ! 1 here).

  real(dp) :: temp1, temp2                      ! Some temporary variables required.




  allocate( temp_array(Nx), y(Nx) )

  do P_i=1, Nx

     ! -----------------------------------------------------------------
     ! In the Soler method, the polynomicals that are interpolated are Kroneker
     ! delta funcions at a particular q point. So, we set all y values to 0
     ! except the one corresponding to the particular function P_i.

     y = 0.0D0
     y(P_i) = 1.0D0

     d2y_dx2(P_i,1) = 0.0D0
     temp_array(1) = 0.0D0

     do idx = 2, Nx-1

        temp1 = (x(idx)-x(idx-1))/(x(idx+1)-x(idx-1))
        temp2 = temp1 * d2y_dx2(P_i,idx-1) + 2.0D0
        d2y_dx2(P_i,idx) = (temp1-1.0D0)/temp2

        temp_array(idx) = (y(idx+1)-y(idx))/(x(idx+1)-x(idx)) &
             - (y(idx)-y(idx-1))/(x(idx)-x(idx-1))
        temp_array(idx) = (6.0D0*temp_array(idx)/(x(idx+1)-x(idx-1)) &
             - temp1*temp_array(idx-1))/temp2

     end do

     d2y_dx2(P_i,Nx) = 0.0D0

     do idx=Nx-1, 1, -1

        d2y_dx2(P_i,idx) = d2y_dx2(P_i,idx) * d2y_dx2(P_i,idx+1) + temp_array(idx)

     end do

  end do
  deallocate( temp_array, y)

  END SUBROUTINE initialize_spline_interpolation

end module vdw_splines
