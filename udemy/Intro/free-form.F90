  program FreeForm
  implicit none
      real(kind=4), parameter :: pi = 3.14159
      INTEGER :: x=10, y, z; y=10; z = x+y

      print *, "x: ", x, " y: ", y, " z: ", z
      PRINT '(i5, i5, i5, F10.5)', x, y, z, pi
  end program FreeForm
