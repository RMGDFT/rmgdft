SUBROUTINE infomsg( routine, message )
  !----------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an info message
  ! ... from a given routine to output.
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*) :: routine, message
  ! the name of the calling routine
  ! the output message
  !
!  IF ( ionode ) THEN   !if not ionode it is redirected to /dev/null anyway
     !
     WRITE( * , '(5X,"Message from routine ",A,":")' ) routine
     WRITE( * , '(5X,A)' ) message
     !
!  END IF
  !
  RETURN
  !
END SUBROUTINE infomsg

