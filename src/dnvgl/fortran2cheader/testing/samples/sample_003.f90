FUNCTION pstr(s) RESULT(x) BIND(c, name='pstr')
   use iso_c_binding ! C bindings
   COMPLEX(C_DOUBLE_COMPLEX) :: s
   COMPLEX(C_DOUBLE_COMPLEX) :: x
END SUBROUTINE pstr
