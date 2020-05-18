! Two Band Modelling Fe-Cr Olsson, Wallenius

! sum (ai (r - ri)^3H(ri - r)

SUBROUTINE simple_spline(r, p, y)
!############################################################
! Simple spline ri = 0.5, 1.0, 1.5, 2.0, 3.0, 4.0
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:6)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: H
REAL(kind=DoubleReal) :: rc(1:6)
INTEGER(kind=StandardInteger) :: n
!############################################################
rc(1) = 0.5D0
rc(2) = 1.0D0
rc(3) = 1.5D0
rc(4) = 2.0D0
rc(5) = 3.0D0
rc(6) = 4.0D0
!#
y = 0.0D0
DO n = 1, 6
  IF((rc(n) - r) < 0.0D0)THEN
    H = 0.0D0
  ELSE
    H = 1.0D0
  END IF
  y = y + p(n) * (r - rc(n))**3 * H
END DO
END SUBROUTINE simple_spline

SUBROUTINE simple_spline_v(r, p, y)
!############################################################
! SIMPLE SPLINE VECTOR FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:6)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL simple_spline(r(n), p,  y(n))
END DO
END SUBROUTINE simple_spline_v
