! Two Band Modelling Fe-Cr Olsson, Wallenius

! sum (ai (r - ri)^3H(ri - r)

SUBROUTINE summed_spline(r, p, y)
!############################################################
! Quartic Poly (r - rc)^2(c0 + c1 r + c2 r^2)
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: H
REAL(kind=DoubleReal) :: rc
REAL(kind=DoubleReal) :: a
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0
DO n = 1, SIZE(p) / 2
  rc = p(2 * n - 1)
  a = p(2 * n)
  IF((rc - r) < 0.0D0)THEN
    H = 0.0D0
  ELSE
    H = 1.0D0
  END IF
  y = y + a * (r - rc)**3 * H
END DO
END SUBROUTINE summed_spline


SUBROUTINE summed_spline_v(r, p, y)
!############################################################
! SUMMED SPLINE VECTOR FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL summed_spline(r(n), p,  y(n))
END DO
END SUBROUTINE summed_spline_v



