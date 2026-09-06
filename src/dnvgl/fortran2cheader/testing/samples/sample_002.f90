subroutine pstr(s) bind(c,name='pstr')
  use iso_c_binding ! C bindings
  COMPLEX(C_DOUBLE_COMPLEX) :: s
end subroutine pstr
