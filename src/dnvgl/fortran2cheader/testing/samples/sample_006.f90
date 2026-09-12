MODULE dummy_interface
   INTERFACE
      SUBROUTINE surf1(m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, &
           & zxymn, islpsw, zp, temp, sigma, ierr)
         IMPLICIT NONE
         INTEGER m, n, iz, islpsw, ierr
         DOUBLE PRECISION x(m), y(n), z(iz, n), zx1(n), zxm(n), zy1(m), zyn(m), &
              & zxy11, zxym1, zxy1n, zxymn, zp(m, n, 3), temp(1), sigma
      END SUBROUTINE surf1
   END INTERFACE
END MODULE dummy_interface

MODULE DUMMY
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: surf1

CONTAINS
   SUBROUTINE surf1(m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
        & islpsw, zp, temp, sigma, ierr) BIND(C, NAME='c_surf1')
      USE, INTRINSIC :: ISO_C_BINDING
      USE dummy_interf, ONLY: f_surf1 => surf1

      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: m, n, iz, islpsw
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(m) :: x, zy1, zyn
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(iz, n) :: z
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(n) :: y, zx1, zxm
      REAL(C_DOUBLE), INTENT(IN) :: zxy11, zxym1, zxy1n, zxymn, sigma
      REAL(C_DOUBLE), INTENT(IN), DIMENSION(n + n + m) :: temp
      REAL(C_DOUBLE), INTENT(OUT), DIMENSION(m, n, 3) :: zp
      INTEGER(C_INT), INTENT(OUT):: ierr

      CALL f_surf1(m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
        & islpsw, zp, temp, sigma, ierr)
   END SUBROUTINE surf1
END MODULE DUMMY
