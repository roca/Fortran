program get_mantissa_exponent
  ! real(kind=16) :: x = 127.23402
  real :: x = 127.23402
  integer :: exp_val
  real :: frac_val

  exp_val = EXPONENT(x)
  frac_val = FRACTION(x)

  print *, "Exponent:", exp_val
  print *, "Fraction:", frac_val
  print *, "Mantissa:", frac_val * 2.0**exp_val
end program get_mantissa_exponent   
