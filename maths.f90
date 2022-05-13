! Copyright (C) 2003-2018 John Young, Matthew Worsley
!
! This file is part of mfit.
!
! Mfit is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/ .

module Maths

  implicit none

  public

  !public subroutines contained:
  !
  !bessel(nu,x,n,output) - returns array output of values for bessel function
  !                         with real argument x and
  !                         orders nu, nu+1 ... nu+(n-1)
  !                         NB 0 <= nu < 1
  !inv_mat(A) - inverts the matrix A
  !invdet_mat(A, det) - inverts the matrix A, also calculates determinant
  !laguerre(n,x, output) - returns real array of laguerre polynomial values for
  !                        argument x. Returns polynomials of order 0 up to n.
  !                        Returned array is double precision with elements 0:n
  !
  !public functions contained:
  !
  !argument(z) -   returns the argument of double complex z
  !bess0(x) -      returns 0th order bessel function of double precision x
  !bess1(x) -      returns 1st order bessel function of double precision x
  !comb(n,p) -     returns combinatorial nCp for integers n, p
  !fact(x) -       returns factorial of integer x
  !gamma(x) -      returns gamma function of double precision x
  !machine_max() - returns maximum machine number
  !machine_min() - returns minimum machine number
  !machine_precision() - returns the machine precsion
  !modulus(z) -    returns the modulus of complex z
  !modx(a,b) -     difference mod 360 of a and b

  !public parameters and conversions
  double precision, parameter :: pi = 3.1415926535897932D0
  double precision, parameter :: mas2rad = pi/6.48D8
  double precision, parameter :: rad2mas = 6.48D8/pi
  double precision, parameter :: deg2rad = pi/180D0
  double precision, parameter :: rad2deg = 180D0/pi

contains

  !============================================================================

  subroutine bessel(nu, x, num, output)

    !function arguments
    double precision :: nu, x
    integer :: num
    double precision, dimension(:), allocatable :: output

    !local variables
    integer :: ncalc

    allocate(output(num))

    call RJBESL(x, nu, num, output, ncalc)
    if (ncalc < 0) then
       stop 'maths error: bessel: an argument to RJBESL is out of range'
    else if (ncalc < num) then
       stop 'maths error: bessel: not all requested function values could be calculated accurately'
    end if

  end subroutine bessel

  !============================================================================

  subroutine inv_mat(A)

    !inverts matrix A in situ

    !subroutine arguments
    double precision, dimension(:,:) :: A

    !local variables
    integer :: ifail
    integer, dimension(:), allocatable :: ipvt
    double precision, dimension(:), allocatable :: work
    double precision, dimension(2) :: dummy !not referenced

    ifail = 0
    allocate(ipvt(size(A, 2)))
    allocate(work(size(A, 2)))

    call PDA_DGEFA(A, size(A, 1), size(A, 2), ipvt, ifail)
    if ( ifail /= 0 ) then
       stop 'maths error: inv_mat: zero-valued element in factorisation; cannot invert matrix'
    else
       call PDA_DGEDI(A, size(A, 1), size(A, 2), ipvt, dummy, work, 1)
    end if

    !clean up
    deallocate(ipvt)
    deallocate(work)

  end subroutine inv_mat

  !==============================================================================

  subroutine invdet_mat(A, det)

    !inverts matrix A in situ, also calculates determinant

    !subroutine arguments
    double precision, dimension(:,:) :: A
    double precision, intent(out) :: det

    !local variables
    integer :: ifail
    integer, dimension(:), allocatable :: ipvt
    double precision, dimension(:), allocatable :: work
    double precision, dimension(2) :: det0

    ifail = 0
    allocate(ipvt(size(A, 2)))
    allocate(work(size(A, 2)))

    call PDA_DGEFA(A, size(A, 1), size(A, 2), ipvt, ifail)
    if ( ifail /= 0 ) then
       stop 'maths error: inv_mat: zero-valued element in factorisation; cannot invert matrix'
    else
       call PDA_DGEDI(A, size(A, 1), size(A, 2), ipvt, det0, work, 11)
       det = det0(1) * 10.0**det0(2)
    end if

    !clean up
    deallocate(ipvt)
    deallocate(work)

  end subroutine invdet_mat

  !============================================================================

  subroutine laguerre(n, x, output)

    !subroutine arguments
    double precision :: x
    double precision, dimension(:), allocatable :: output
    integer :: n

    !local variables
    integer :: i
    double precision :: a, b, c

    !fatal erros
    if (n < 0) stop 'maths error: laguerre: n < 0'

    allocate(output(0:n))

    !generate laguerre sequence using recurrence relation
    !that tL[t](x) = (2t-1-x)L[t-1](x) - (t-1)L[t-2](x)
    !noting that L[0](x) = 1 and L[1](x) = 1-x

    a = 1D0
    b = 1D0-x
    output(0) = a
    if (n > 0) output(1) = b

    do i = 2, n
       c = ((((2D0*dble(i))-1D0-x)*b)-((dble(i)-1D0)*a))/(dble(i))
       a = b
       b = c
       output(i) = c
    end do

  end subroutine laguerre

  !============================================================================

  function argument(z)

    !returns argument of complex number z in the range -pi<arg<=pi

    !function arguments
    double complex :: z
    double precision :: argument

    !find argument
    argument = atan2(aimag(z), dble(z))

  end function argument

  !============================================================================

  function bess0(x)

    !function arguments
    double precision :: bess0, x

    !local variables
    integer :: ncalc
    double precision, dimension(1) :: output

    call RJBESL(x, 0D0, 1, output, ncalc)
    if (ncalc < 0) then
       stop 'maths error: bessel: an argument to RJBESL is out of range'
    else if (ncalc < 1) then
       stop 'maths error: bessel: not all requested function values could be calculated accurately'
    end if
    bess0 = output(1)

  end function bess0

  !============================================================================

  function bess1(x)

    !function arguments
    double precision :: bess1, x

    !local variables
    integer :: ncalc
    double precision, dimension(2) :: output
    integer :: ifail

    call RJBESL(x, 0D0, 2, output, ncalc)
    if (ncalc < 0) then
       stop 'maths error: bess1: an argument to RJBESL is out of range'
    else if (ncalc < 2) then
       stop 'maths error: bess1: not all requested function values could be calculated accurately'
    end if
    bess1 = output(2)

  end function bess1

  !============================================================================

  function comb(n, p)

    !function arguments
    double precision :: comb
    integer :: n, p

    !fatal errors
    if (n < 0) stop 'maths error: comb: n < 0'
    if (n > 170) stop 'maths error: comb: n > 170'
    if (p < 0) stop 'maths error: comb: p < 0'
    if (p > n) stop 'maths error: comb: p > n'

    comb = fact(n)/(fact(p)*fact(n-p))

  end function comb

  !============================================================================

  function fact(x)

    !function arguments
    double precision :: fact
    integer :: x

    !local variables
    integer :: i

    !fatal errors
    if (x < 0) stop 'maths error: fact: x < 0'
    if (x > 170) stop 'maths error: fact: x > 170'

    fact = 1D0
    do i = 1, x
       fact = fact * dble(i)
    end do

  end function fact

  !============================================================================

  function gamma(x)

    !function arguments
    double precision :: gamma, DGAMMA, x

    gamma = DGAMMA(x)

  end function gamma

  !============================================================================

  function machine_max()

    !returns the machine max number (largest magnitude)

    !function arguments
    double precision :: machine_max

    !local variables
    double precision :: PDA_D1MACH

    machine_max = PDA_D1MACH(2)

  end function machine_max

  !============================================================================

  function machine_min()

    !returns the machine min number (smallest magnitude)

    !function arguments
    double precision :: machine_min

    !local variables
    double precision :: PDA_D1MACH

    machine_min = PDA_D1MACH(1)

  end function machine_min

  !============================================================================

  function machine_precision()

    !returns the machine precision

    !function arguments
    double precision :: machine_precision

    !local variables
    double precision :: PDA_D1MACH

    machine_precision = PDA_D1MACH(3)

  end function machine_precision

  !============================================================================

  function modulus(z)

    !returns modulus of complex z, uses method from
    !  Wilkinson J H and Reinsch C (1971) "Handbook for Automatic Computation II,
    !  Linear Algebra", Springer-Verlag
    !for highest precision and to
    !avoid underflow errors if im(z)<<re(z) or im(z)>>re(z)

    !function arguments
    double complex :: z
    double precision :: modulus

    !local variables
    double precision :: re, imag, a, b

    re = abs(dble(z))
    imag = abs(aimag(z))

    if (re == 0D0 .and. imag == 0D0) then
       modulus = 0D0
       return
    end if

    if (re > imag) then
       a = re
       b = imag
    else
       a = imag
       b = re
    end if

    modulus = a*sqrt(1D0 + (b/a)**2D0)

  end function modulus

  !============================================================================

  function modx(a,b)

    !returns mod 360 difference of a and b

    !function arguments
    double precision :: a, b, modx

    !local variables

    modx = modulo(a-b,360D0)
    if (modx>180D0) modx = 360D0 - modx

  end function modx

  !============================================================================

end module Maths
