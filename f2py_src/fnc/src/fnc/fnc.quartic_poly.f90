! V-Ti-Cr Fu, Li, Johansson, Zhao

SUBROUTINE quartic_poly(r, p, y)
!############################################################
! Quartic Poly (r - rc)^2(c0 + c1 r + c2 r^2)
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
y = (r - p(1))**2 * (p(2) + p(3) * r + p(4) * r**2)
END SUBROUTINE quartic_poly

SUBROUTINE quartic_poly_v(r, p, y)
!############################################################
! Quartic Poly (r - rc)^2(c0 + c1 r + c2 r^2)
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL quartic_poly(r(n), p,  y(n))
END DO
END SUBROUTINE quartic_poly_v
