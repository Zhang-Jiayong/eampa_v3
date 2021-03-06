

SUBROUTINE mendelev_embedding(r, p, y)
!############################################################
! f(x) = -sqrt(rho) + A*rho**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
y = -1.0D0 * sqrt(r) + p(1) * r**2
END SUBROUTINE mendelev_embedding


SUBROUTINE mendelev_embedding_v(r, p, y)
!############################################################
! f(x) = -sqrt(rho) + A*rho**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL mendelev_embedding(r(n), p,  y(n))
END DO
END SUBROUTINE mendelev_embedding_v



