MODULE dummy_interface
   INTERFACE
      SUBROUTINE curvp1(n, x, y, p, yp, temp, sigma, ierr) BIND(C, NAME='c_curvp1')
         IMPLICIT NONE
         INTEGER(INT), INTENT(IN) :: n
         REAL(DOUBLE), INTENT(IN), DIMENSION(n) :: x, y
         REAL(DOUBLE), INTENT(IN) :: p, sigma
         REAL(DOUBLE), INTENT(IN), DIMENSION(n*2) :: temp
         REAL(DOUBLE), INTENT(OUT), DIMENSION(n) :: yp
         INTEGER(INT), INTENT(OUT) :: ierr
      END SUBROUTINE curvp1
   END INTERFACE
END MODULE dummy_interface

MODULE DUMMY
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: curvp1

CONTAINS
   SUBROUTINE curvp1(n, x, y, p, yp, temp, sigma, ierr) BIND(C, NAME='c_curvp1')
      USE, INTRINSIC :: ISO_C_BINDING
      USE fitpack_interf, ONLY: f_curvp1 => curvp1
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: n
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: x, y
      REAL(C_DOUBLE), INTENT(IN) :: p, sigma
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(n*2) :: temp
      REAL(C_DOUBLE), INTENT(OUT), DIMENSION(n) :: yp
      INTEGER(C_INT), INTENT(OUT) :: ierr

      CALL f_curvp1(n, x, y, p, yp, temp, sigma, ierr)
   END SUBROUTINE curvp1
END MODULE DUMMY
