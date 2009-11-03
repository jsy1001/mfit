!$Id: search.f90,v 1.1 2009/11/03 17:26:27 jsy1001 Exp $

module Search
  
  implicit none

  public
  
contains

  !============================================================================

  function locate(xx, x)

    !return index i such that xx(i) < x < xx(i+1)
    !xx must be monotonic

    !function arguments
    real, dimension(:) :: xx
    real :: x
    integer :: locate

    !local variables
    integer n, lower, mid, upper

    n = size(xx)
    lower = 0
    upper = n + 1
    do while (upper - lower > 1)
       mid = (upper + lower)/2
       if ((xx(n) > xx(1)).eqv.(x > xx(mid))) then
          lower = mid
       else
          upper = mid
       end if
    end do
    locate = lower

  end function locate

  !============================================================================

end module Search
