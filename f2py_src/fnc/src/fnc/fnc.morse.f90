

SUBROUTINE morse(r, p, y)
!############################################################
! MORSE FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
y = p(1) * (exp(-2.0D0 * p(2) * (r - p(3))) - 2.0D0 * exp(-p(2)*(r - p(3))))
END SUBROUTINE morse




SUBROUTINE morse_v(r, p, y)
!############################################################
! MORSE FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
y(:) = p(1) * (exp(-2.0D0 * p(2) * (r(:) - p(3))) - 2.0D0 * exp(-p(2)*(r(:) - p(3))))
END SUBROUTINE morse_v





