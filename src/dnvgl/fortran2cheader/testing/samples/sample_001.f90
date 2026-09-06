subroutine pstr(s) bind(c,name='pstr')
  use iso_c_binding ! C bindings
  character(kind=c_char,len=1), intent(in) :: &
     & s(*)
end subroutine pstr
