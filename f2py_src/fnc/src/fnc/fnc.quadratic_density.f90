

SUBROUTINE quadratic_density(r, p, y)
!############################################################
! LENNARD JONES FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
IF(r .GT. p(1))THEN
  y = 0.0D0
ELSE
  y = (r - d)**2
END IF
END SUBROUTINE quadratic_density

SUBROUTINE quadratic_density_v(r, p, y)
!############################################################
! LENNARD JONES VECTOR FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL quadratic_density(r(n), p,  y(n))
END DO
END SUBROUTINE quadratic_density_v

