program Test
  implicit none
  integer(kind=4):: i

  print *, "Hello, World!"

  do while (i <= 5)
    print *, "This is iteration number ", i
    i = i + 1
  end do
end program Test
