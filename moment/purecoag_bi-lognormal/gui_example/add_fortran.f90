module main
  real :: x, y
  real :: z
contains
  subroutine addTwoNumbers
    print *, "x = ", x
    print *, "y = ", y
    z = x + y
  end subroutine addTwoNumbers
end module main
