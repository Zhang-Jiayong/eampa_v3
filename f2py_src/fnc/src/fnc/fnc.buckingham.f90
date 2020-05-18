

SUBROUTINE buckingham(r, p, y)
!############################################################
! BUCKINGHAM POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
y = p(1) * exp(-1 * p(2) * r) - p(3) / r**6
END SUBROUTINE buckingham




SUBROUTINE buckingham_v(r, p, y)
!############################################################
! BUCKINGHAM POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
y(:) = p(1) * exp(-1 * p(2) * r(:)) - p(3) / r(:)**6
END SUBROUTINE buckingham_v





