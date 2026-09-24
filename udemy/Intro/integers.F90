program Integers
  ! max values if kind=1 is 127, kind=2 is 32767, kind=4 is 2147483647
  integer(kind=1):: small_num = 127
  integer(kind=2):: medium_num = 32767
  integer(kind=4):: large_num = 2147483647
  integer(kind=8):: very_large_num = 9223372036854775807_8
  PRINT *, small_num, medium_num, large_num, very_large_num

end program Integers
