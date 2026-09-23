C     This is an example of a Fortran 77 program that demonstrates the use of old-format input and output statements. The program reads two integers from the user, adds them together, and prints the result.

      PROGRAM OLD_FORMAT
      INTEGER A, B, SUM

      PRINT *, 'Enter two integers:'
      READ *, A, B

      SUM = A + B
C     The following line uses the old-format output statement to print the
C     result of the addition. The '1' in the second line indicates that the output should continue on the next line.
      PRINT *, 'The sum of', A,
     1 'and', B, 'is', SUM
      END PROGRAM OLD_FORMAT
