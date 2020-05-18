! FINNIS SINCLAIR EMBEDDING

SUBROUTINE fs_embedding(r, p, y)
!############################################################
! f(x) = -A sqrt(rho)
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
y = -1.0D0 * p(1) * sqrt(r)
END SUBROUTINE fs_embedding


SUBROUTINE fs_embedding_v(r, p, y)
!############################################################
! f(x) = -A sqrt(rho)
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL fs_embedding(r(n), p,  y(n))
END DO
END SUBROUTINE fs_embedding_v



