
SUBROUTINE dftd3_driver(natoms, nSpecies, coords, atnum, latt, functype, d3_version, edisp, grads, stress, func_len) 
    use dftd3_api
    implicit none

    ! Working precision
    integer, parameter :: wp = kind(1.0d0)


    integer, intent(in) :: natoms

    integer, intent(in) :: nSpecies
    integer, intent(in) :: func_len

    ! Coordinates in Angstrom as found in dna.xyz/dna.poscar
    ! They must be converted to Bohr before passed to dftd3
    real(wp), intent(in) :: coords(3,nAtoms)
    integer, intent(in) :: atnum(nAtoms)
    real(wp), intent(in) :: latt(3,3)
    character(func_len), intent(in) :: functype
    integer, intent(in) :: d3_version

    real(wp), intent(out) :: edisp
    real(wp), intent(out) :: grads(3,nAtoms)
    real(wp), intent(out) :: stress(3,3)

    type(dftd3_input) :: input
    type(dftd3_calc) :: dftd3
    !ele = ['C ', 'H ', 'H ', 'H ']

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize input
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! You can set input variables if you like, or just leave them on their
    ! defaults, which are the same as the dftd3 program uses.

!! Threebody interactions (default: .false.)
    !input%threebody = .true.
    !
!! Numerical gradients (default: .false.)
    !input%numgrad = .false.
    !
    !! Cutoffs (below you find the defaults)
    !input%cutoff = sqrt(9000.0_wp)
!input%cutoff_cn = sqrt(1600.0_wp)

    ! Initialize dftd3
    call dftd3_init(dftd3, input)

    ! Choose functional. Alternatively you could set the parameters manually
    ! by the dftd3_set_params() function.
!call dftd3_set_functional(dftd3, func=functype, version=d3_version, tz=.false.)
    call dftd3_set_functional(dftd3, func=functype, version=d3_version, tz=.false.)

    ! Calculate dispersion and gradients for periodic case
    call dftd3_pbc_dispersion(dftd3, coords, atnum, latt, edisp, grads, stress)
    !write(*, "(A)") "*** Dispersion for periodic case"
    !write(*, "(A,ES20.12)") "Energy [au]:", edisp
    !write(*, "(A)") "Gradients [au]:"
    !write(*, "(3ES20.12)") grads
    !write(*, "(A)") "Stress [au]:"
    !write(*, "(3ES20.12)") stress
    end SUBROUTINE
