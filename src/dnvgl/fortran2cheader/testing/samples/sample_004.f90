MODULE fitpack_interf
   INTERFACE
      SUBROUTINE curv1(n, x, y, slp1, slpn, islpsw, yp, temp, sigma, ierr)
         IMPLICIT NONE
         INTEGER(C_INT), INTENT(IN) :: n
         REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: x
         REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: y
         REAL(C_DOUBLE), INTENT(IN) :: slp1
         REAL(C_DOUBLE), INTENT(IN) :: slpn
         INTEGER(C_INT), INTENT(IN) :: islpsw
         REAL(C_DOUBLE), INTENT(OUT), DIMENSION(n) :: yp
         REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: temp(n)
         REAL(C_DOUBLE), INTENT(IN) :: sigma
         INTEGER(C_INT), INTENT(OUT):: ierr
      END SUBROUTINE curv1
      FUNCTION curv2(t, n, x, y, yp, sigma) result(res) 
         IMPLICIT NONE
         REAL(C_DOUBLE) res
         INTEGER(C_INT) n
         REAL(C_DOUBLE) t, x(n), y(n), yp(n), sigma
         REAL(C_DOUBLE) result

         res = f_curv2(t, n, x, y, yp, sigma)
      END FUNCTION curv2
   END INTERFACE
END MODULE fitpack_interf

MODULE fitpack_c
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: curv1, curv2

CONTAINS

   SUBROUTINE curv1(n, x, y, slp1, slpn, islpsw, yp, temp, sigma, ierr) &
        & BIND(C, NAME='c_curv1')
      USE, INTRINSIC :: ISO_C_BINDING
      USE fitpack_interf, ONLY: f_curv1 => curv1

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: n
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: x
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: y
      REAL(C_DOUBLE), INTENT(IN) :: slp1
      REAL(C_DOUBLE), INTENT(IN) :: slpn
      INTEGER(C_INT), INTENT(IN) :: islpsw
      REAL(C_DOUBLE), INTENT(OUT), DIMENSION(n) :: yp
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: temp(n)
      REAL(C_DOUBLE), INTENT(IN) :: sigma
      INTEGER(C_INT), INTENT(OUT):: ierr

      CALL f_curv1(n, x, y, slp1, slpn, islpsw, yp, temp, sigma, ierr)
   END SUBROUTINE curv1
   FUNCTION curv2(t, n, x, y, yp, sigma) result(res) BIND(C, NAME='c_curv2')
      USE, INTRINSIC :: ISO_C_BINDING
      USE fitpack_interf, ONLY: f_curv2 => curv2

      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: res
      INTEGER(C_INT), INTENT(IN) :: n
      REAL(C_DOUBLE), INTENT(IN) :: t, sigma
      REAL(C_DOUBLE), INTENT(IN), DIMENTSION(n) :: x, y, yp
      REAL(C_DOUBLE) result

      res = f_curv2(t, n, x, y, yp, sigma)
   END FUNCTION curv2
END MODULE fitpack_c
