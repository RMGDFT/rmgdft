!
! Copyright (C) 2001-2011 Quantum ESPRESSO  group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE irreducible_BZ( nrot, s, nsym, minus_q, magnetic_sym, at, bg, &
                           npk, nks, xk, wk, t_rev )
  !-----------------------------------------------------------------------
  !! This routine finds the special points in the irreducible wedge
  !! of the true point group (or small group of q) of the crystal, 
  !! starting from the points in the irreducible BZ wedge of the 
  !! point group of the Bravais lattice.
  !
  USE kinds,   ONLY: DP
  USE iso_c_binding
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: nrot
  !! order of the parent point group
  INTEGER,  INTENT(IN) :: nsym
  !! order of the subgroup
  INTEGER,  INTENT(IN) :: npk
  !! maximum number of special points
  INTEGER,  INTENT(IN) :: s(3,3,48)
  !! symmetry matrices, in crystal axis
  INTEGER,  INTENT(IN) :: t_rev(48)
  !! time reversal operation
  REAL(DP), INTENT(IN) :: at(3,3)
  !! basis vectors of the Bravais lattice
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! basis vectors of the reciprocal lattice
  INTEGER,  INTENT(IN) :: minus_q
  !! it is .TRUE. if symmetries q=-q+G are acceptable
  INTEGER,  INTENT(IN) :: magnetic_sym
  !! magnetic_sym = noncolin .AND. domag
  INTEGER,  INTENT(INOUT) :: nks
  !! number of special points
  REAL(DP), INTENT(INOUT) :: xk(3,npk)
  !! special points
  REAL(DP), INTENT(INOUT) :: wk(npk)
  !! weights for special points
  !
  ! ... local variables
  !
  INTEGER :: table(48,48), invs(3,3,48), irg(48)
  ! table: multiplication table of the group
  ! invs:  contains the inverse of each rotation
  ! irg:   gives the correspondence of symmetry operations forming
  !        a n-th coset
  INTEGER :: isym, jsym
  LOGICAL :: sym(48)
  !
  !
  ! ... We compute the multiplication table of the group.
  !
  CALL multable( nrot, s, table )
  !
  ! ... And we set the matrices of the inverse.
  !
  DO isym = 1, nrot
     DO jsym = 1, nrot
        IF (table(isym,jsym) == 1) invs(:,:,isym) = s(:,:,jsym)
     ENDDO
  ENDDO
  !
  ! ... Find the coset in the point group of the Bravais lattice.
  !
  IF ( magnetic_sym .eq. 1 ) THEN
     CALL irrek_nc( at, bg, nrot, invs, nsym, irg, npk, nks, xk, &
                    wk, t_rev )
  ELSE
     sym(1:nsym) = .TRUE.
     sym(nsym+1:) = .FALSE.
     CALL coset( nrot, table, sym, nsym, irg )
     !
     ! ... here we set the k-points in the irreducible wedge of the point grou
     !     of the crystal.
     !
     CALL irrek( at, bg, nrot, invs, nsym, irg, minus_q, npk, nks, xk, &
                 wk, t_rev )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE irreducible_BZ 
!
!
!-----------------------------------------------------------------------
SUBROUTINE irrek( at, bg, nrot, invs, nsym, irg, minus_q, npk, &
                  nks, xk, wk, t_rev )
  !-----------------------------------------------------------------------
  !! Given a set of special points in the Irreducible Wedge of some
  !! group, finds the equivalent special points in the IW of one of
  !! its subgroups.
  !
  USE kinds,   ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npk
  !! maximum number of special points
  INTEGER, INTENT(IN) :: nrot
  !! order of the parent point group
  INTEGER, INTENT(IN) :: nsym
  !! order of the subgroup
  INTEGER, INTENT(IN) :: invs(3,3,48)
  !! inverse of the elements of the symmetry group
  INTEGER, INTENT(IN) :: irg(nrot)
  !! partition of the elements of the symmetry group into left cosets
  !! as given by subroutine COSET.
  INTEGER, INTENT(INOUT) :: nks
  !! number of special points
  INTEGER, INTENT(IN) :: t_rev(48)
  !! time reversal operation
  REAL(DP), INTENT(IN) :: at(3,3)
  !! basis vectors of the Bravais lattice
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! basis vectors of the reciprocal lattice
  REAL(DP), INTENT(INOUT) :: xk(3,npk)
  !! special points
  REAL(DP), INTENT(INOUT) :: wk(npk)
  !! weights for special points
  INTEGER, INTENT(IN) :: minus_q
  ! .TRUE. if symmetries q=-q+G are acceptable
  !
  ! ... local variables
  !
  INTEGER :: nks0, jk, kpol, irot, jrot, ncos, jc, ic, isym
  ! nks0: used to save the initial number of k-points
  ! ncos: total number of cosets
  REAL(DP) :: xkg(3), xks(3,48), w(48), sw, one
  ! coordinates of the k point in crystal axis
  ! coordinates of the rotated k point
  ! weight of each coset
  ! buffer which contains the weight of k points
  ! total weight of k-points
  LOGICAL :: latm, satm
  ! true if a k-point is equivalent to a previous one
  ! true if equivalent point found
  !
  !
  nks0 = nks
  !
  DO jk = 1, nks0
     !
     ! ... The k point is first computed in crystal axis
     !
     DO kpol = 1, 3
        ! xkg are the components ofx k in the crystal RL base
        xkg(kpol) = at(1,kpol) * xk(1,jk) + &
                    at(2,kpol) * xk(2,jk) + &
                    at(3,kpol) * xk(3,jk)
     ENDDO
     !
     ! ... Then it is rotated with each symmetry of the global group. Note that
     !     the irg vector is used to divide all the rotated vector in cosets
     !
     DO irot = 1, nrot
        jrot = irg(irot)
        DO kpol = 1, 3
           ! the rotated of xkg with respect to the group operations
           xks(kpol,irot) = invs(kpol,1,jrot) * xkg(1) + &
                            invs(kpol,2,jrot) * xkg(2) + &
                            invs(kpol,3,jrot) * xkg(3)
        ENDDO
        IF (t_rev(jrot) == 1)  xks(:,irot) = -xks(:,irot)
     ENDDO
     !
     !    For each coset one point is tested with all the preceding
     !
     ncos = nrot / nsym
     DO ic = 1, ncos
        irot = (ic - 1) * nsym + 1
        latm = .FALSE.
        !
        ! ... latm = .TRUE. if the present k-vector is equivalent to some previous
        !
        DO jc = 1, ic-1
           DO isym = 1, nsym
              !
              ! ... satm = .TRUE. if the present symmetry operation makes 
              !     the ir and ik k-vectors equivalent ...
              !
              jrot = (jc - 1) * nsym + isym
              satm = ABS(  xks(1,irot) - xks(1,jrot) - &
                     NINT( xks(1,irot) - xks(1,jrot) ) ) < 1.0d-5 .AND. &
                     ABS(  xks(2,irot) - xks(2,jrot) - &
                     NINT( xks(2,irot) - xks(2,jrot) ) ) < 1.0d-5 .AND. &
                     ABS(  xks(3,irot) - xks(3,jrot) - &
                     NINT( xks(3,irot) - xks(3,jrot) ) ) < 1.0d-5
              !
              !  .... or equivalent to minus each other when minus_q=.t.
              !
              IF (minus_q .eq. 1) satm = satm .OR. &
                   ABS(  xks(1,irot) + xks(1,jrot) - &
                   NINT( xks(1,irot) + xks(1,jrot) ) ) < 1.0d-5 .AND. &
                   ABS(  xks(2,irot) + xks(2,jrot) - &
                   NINT( xks(2,irot) + xks(2,jrot) ) ) < 1.0d-5 .AND. &
                   ABS(  xks(3,irot) + xks(3,jrot) - &
                   NINT( xks(3,irot) + xks(3,jrot) ) ) < 1.0d-5
              latm = latm .OR. satm
              IF (satm .AND. w(jc) /= 0.d0) THEN
                 w(jc) = w(jc) + 1.d0
                 GOTO 100
              ENDIF
           ENDDO
           !
        ENDDO
        !
100     CONTINUE
        !
        IF (latm) THEN
           w(ic) = 0.d0
        ELSE
           w(ic) = 1.d0
        ENDIF
        !
     ENDDO
     !
     ! ... here the k-point list is updated
     !
     sw = wk(jk) / SUM( w(1:ncos) )
     wk(jk) = sw * w(1)
     DO ic = 2, ncos
        irot = (ic - 1) * nsym + 1
        IF (w(ic) /= 0.d0) THEN
           nks = nks + 1
           IF (nks > npk) CALL errore( 'irrek', 'too many k-points', nks )
           wk (nks) = sw * w(ic)
           DO kpol = 1, 3
              xk(kpol, nks) = bg(kpol,1) * xks(1,irot) + &
                              bg(kpol,2) * xks(2,irot) + &
                              bg(kpol,3) * xks(3,irot)
           ENDDO
        ENDIF
     ENDDO
     !
  ENDDO
  !
  ! normalize weights to one
  !
  one = SUM( wk(1:nks) )
  IF ( one > 0.d0 ) wk(1:nks) = wk(1:nks) / one
  !
  RETURN
  !
END SUBROUTINE irrek
!
!
!-----------------------------------------------------------------------
SUBROUTINE irrek_nc( at, bg, nrot, invs, nsym, irg, npk, &
                     nks, xk, wk, t_rev )
  !-----------------------------------------------------------------------
  !! Given a set of special points in the Irreducible Wedge of some
  !! group, finds the equivalent special points in the IW of one of
  !! its subgroups.
  !
  USE kinds,    ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npk
  !! maximum number of special points
  INTEGER, INTENT(IN) :: nrot
  !! order of the parent point group
  INTEGER, INTENT(IN) :: nsym
  !! order of the subgroup
  INTEGER, INTENT(IN) :: invs(3,3,48)
  !! inverse of the elements of the symmetry group
  INTEGER, INTENT(IN) :: irg(nrot)
  !! partition of the elements of the symmetry group into left cosets
  !! as given by SUBROUTINE COSET
  INTEGER, INTENT(INOUT) :: nks
  !! number of special points
  INTEGER, INTENT(IN) :: t_rev(48)
  !! time reversal operation
  REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
  ! basis vectors of the Bravais and reciprocal lattice
  REAL(DP), INTENT(INOUT) :: xk(3,npk), wk(npk)
  ! special points and weights
  !
  ! ... local variables
  !
  INTEGER :: nks0, jk, kpol, irot, jrot, isym, ik, iks, start_k
  ! nks0: used to save the initial number of k-points
  ! ncos: total number of cosets
  REAL(DP) :: xkg(3), xks(3), xkn(3), one, xk_new(3,npk), wk_new(npk), &
              xk_cart(3)
  ! coordinates of the k point in crystal axis
  ! coordinates of the rotated k point
  ! weight of each coset
  ! buffer which contains the weight of k points
  ! total weight of k-points
  LOGICAL :: satm
  ! true if equivalent point found
  !
  !
  nks0 = nks
  nks = 0
  start_k = 0
  !
  DO jk = 1, nks0
     !
     ! ... The k point is first computed in crystal axis
     !
     ! xkg are the components of xk in the crystal base
     xkg(:) = at(1,:) * xk(1,jk) + &
              at(2,:) * xk(2,jk) + &
              at(3,:) * xk(3,jk)
     !
     ! ... Then it is rotated with each symmetry of the global group. 
     !
     DO irot = 1, nrot
        xks(:) = invs(:,1,irot) * xkg(1) + &
                 invs(:,2,irot) * xkg(2) + &
                 invs(:,3,irot) * xkg(3)
        ! 
        !  ... Now check if there is an operation of the subgroup that
        !      makes xks equivalent to some other already found k point
        !
        DO jrot = 1, nsym
           xkn(:) = invs(:,1,jrot) * xks(1) + &
                    invs(:,2,jrot) * xks(2) + &
                    invs(:,3,jrot) * xks(3)
           IF (t_rev(jrot) == 1) xkn = -xkn
           DO ik = start_k+1, nks
              satm = ABS(  xk_new(1,ik) - xkn(1) - &
                     NINT( xk_new(1,ik) - xkn(1) ) ) < 1.0d-5 .AND. &
                     ABS(  xk_new(2,ik) - xkn(2) - &
                     NINT( xk_new(2,ik) - xkn(2) ) ) < 1.0d-5 .AND. &
                     ABS(  xk_new(3,ik) - xkn(3) - &
                     NINT( xk_new(3,ik) - xkn(3) ) ) < 1.0d-5
              IF ( satm ) THEN
                 wk_new(ik) = wk_new(ik) + wk(jk)
                 GOTO 100
              ENDIF   
           ENDDO
        ENDDO
        nks = nks+1
        IF (nks > npk) CALL errore( 'irrek_nc','too many k points', 1 )
        xk_new(:,nks) = xks
        wk_new(nks) = wk(jk)
100     CONTINUE
     ENDDO
     start_k=nks
  ENDDO
  !
  !  The order of the original k points is preserved
  !
  iks = nks0
  DO ik = 1, nks
     !
     !    for each new k point found, check if it was in the original list
     !
     DO jk = 1, nks0
        !
        xkg(:) = at(1,:) * xk(1,jk) + &
                 at(2,:) * xk(2,jk) + &
                 at(3,:) * xk(3,jk)
        satm = ABS(  xk_new(1,ik) - xkg(1) - &
               NINT( xk_new(1,ik) - xkg(1) ) ) < 1.0d-5 .AND. &
               ABS(  xk_new(2,ik) - xkg(2) - &
               NINT( xk_new(2,ik) - xkg(2) ) ) < 1.0d-5 .AND. &
               ABS(  xk_new(3,ik) - xkg(3) - &
               NINT( xk_new(3,ik) - xkg(3) ) ) < 1.0d-5
               !
        IF (satm) THEN
           ! If it was, just update the weight
           wk(jk) = wk_new(ik)
           GOTO 200
        ENDIF
        !
     ENDDO
     !
     ! If it was not, bring xk_new in cartesian coodinates and copy it in the
     ! first free place available
     !
     iks = iks + 1
     xk_cart(:) = bg(:,1) * xk_new(1,ik) + &
                  bg(:,2) * xk_new(2,ik) + &
                  bg(:,3) * xk_new(3,ik)
     xk(:,iks) = xk_cart(:)
     wk(iks) = wk_new(ik)
     !
200  CONTINUE
     !
  ENDDO
  !
  IF (iks /= nks ) CALL errore( 'irrek_nc', 'Internal problem with k points', 1 )
  !
  ! normalize weights to one
  !
  one = SUM( wk(1:nks) )
  IF ( one > 0.d0 ) wk(1:nks) = wk(1:nks) / one
  !
  RETURN
  !
END SUBROUTINE irrek_nc


!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE multable( nsym, s, table )
  !-----------------------------------------------------------------------
  !! Checks that {S} is a group and calculates multiplication table
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nsym
  !! the number of symmetry operations
  INTEGER, INTENT(IN) :: s(3,3,nsym)
  !! rotation matrix (in crystal axis, represented by integers)
  INTEGER, INTENT(OUT) :: table(48,48)
  ! multiplication table:  S(n)*S(m) = S( table(n,m) )
  !
  ! ... local variables
  !
  INTEGER :: isym, jsym, ksym, ss(3,3)
  LOGICAL :: found, smn
  !
  !
  DO isym = 1, nsym
     DO jsym = 1, nsym
        ! 
        ss = MATMUL( s(:,:,jsym), s(:,:,isym) )
        !
        ! ... here we check that the input matrices really form a group
        !     and we set the multiplication table
        !
        found = .FALSE.
        DO ksym = 1, nsym
           smn =  ALL( s(:,:,ksym) == ss(:,:) )
           IF (smn) THEN
              IF (found) CALL errore( 'multable', 'Not a group', 1 )
              found = .TRUE.
              table (jsym,isym) = ksym
           ENDIF
        ENDDO
        IF ( .NOT. found) CALL errore( 'multable', ' Not a group', 2 )
        !
     ENDDO
  ENDDO
  !
  !
  RETURN
  !
END SUBROUTINE multable


!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE coset( nrot, table, sym, nsym, irg )
  !-----------------------------------------------------------------------
  !!  Divides the elements of a given group into left cosets of one
  !!  of its subgroups.
  !
  !!  The input is the array sym which is true only for the
  !!  operations of the subgroup, the output is nsym, and the array irg,
  !!  which contains as its first elements the indices of the subgroup,
  !!  and then its right cosets.
  !
  !!  Revised layout 1 may 1995 by A. Dal Corso
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER :: nrot
  !! input: order of the group
  INTEGER :: table(48, 48)
  !! input: multiplication table of the group
  INTEGER :: nsym
  !! output: order of the subgroup
  INTEGER :: irg(48)
  !! output: gives the correspondence of symme
  !! operations forming a n-th coset
  LOGICAL :: sym(48)
  !! input: flag indicating if an operation
  !! belongs to the subgroup
  !
  ! ... local variables
  !
  LOGICAL :: done(48)
  ! if true the operation has been already ch
  !
  INTEGER :: irot, ncos, isym, nc, nelm
  ! counter on rotations
  ! number of cosets (=nrot/nsym)
  ! counter on symmetries
  ! counter on cosets
  ! counter on the number of elements
  !
  !    here we count the elements of the subgroup and set the first part o
  !    irg which contain the subgroup
  !
  nsym = 0
  DO irot = 1, nrot
     done (irot) = sym (irot)
     IF ( sym (irot) ) THEN
        nsym = nsym + 1
        irg (nsym) = irot
     ENDIF
  ENDDO
  !
  ! ... we check that the order of the subgroup is a divisor of the order
  ! total group. ncos is the number of cosets
  !
  IF ( nsym == 0 ) CALL errore( 'coset', 'nsym == 0', 1 ) 
  !
  ncos = nrot / nsym
  IF ( ncos * nsym /= nrot ) CALL errore( 'coset', &
   'The order'//' of the group is not a multiple of that of the subgroup', 1 )
  !
  ! ... here we set the other elements of irg, by using the multiplication
  !
  nelm = nsym
  DO nc = 2, ncos
     DO irot = 1, nrot
        IF ( .NOT.done(irot) ) THEN
           DO isym = 1, nsym
              nelm = nelm + 1
              irg (nelm) = table (irot, irg (isym) )
              done (irg (nelm) ) = .TRUE.
           ENDDO
        ENDIF
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE coset
