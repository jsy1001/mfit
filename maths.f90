module Maths

!subroutines contained
!
!bessel(nu,x,y,n,output) - returns complex array output of bessel function
!                          values for bessel function with complex argument
!                          x + iy and orders nu, nu+1 ... nu+(n-1)
!inv_rsis(A) - returns the inverse of real symmetric indefinite square matrix A
!laguerre(n,x) - returns real array of laguerre polynomial values for argument
!                x. returns polynomials of order 0 up to n. array is dble
!                precision with elemens 0:n
!rnd_initialise() - initialise random number generator 
!
!functions contained
!
!argument(z) - returns the argument of double complex z
!bess1(x) - returns 1st order bessel function of double precision x
!comb(n,p) - returns combinatorial nCp for integers n, p
!fact(x) - returns factorial of integer x
!gamma(x) - returns gamma function of double precision x
!machine_max() - returns maximum machine number
!machine_min() - returns minimum machine number
!machine_precision() - returns the machine precsion
!modulus(z) - returns the modulus of complex z
!modx(a,b) - difference mod 360 of a and b
!rnd_gauss(mu,sigma) - returns random value selected from N(mu,sigma)
!rnd_unit - returns random value selected in range 0<rnd_unit<1

implicit none

!parameters and conversions
double precision, parameter :: pi = 3.1415926535897932D0
double precision, parameter :: mas2rad = pi/6.48D8
double precision, parameter :: rad2mas = 6.48D8/pi
double precision, parameter :: deg2rad = pi/180D0
double precision, parameter :: rad2deg = 180D0/pi

contains

!==============================================================================
subroutine bessel(nu, x, y, num, output)

  !function arguments
  double precision :: nu, x, y
  integer :: num
  double complex :: z
  double complex, dimension(:), allocatable :: output
  
  !local variables
  integer :: ifail, nz, i
  character(len=1) :: scale

  !fatal errors - leave to NAG
  
  !allocate final output
  allocate(output(num))

  !find complex input
  z = dcmplx(x, y)

  !call NAG bessel function subroutine
  ifail = 0
  scale = 'u'
  call S17DEF(nu, z, num, scale, output, nz, ifail)
  if (nz /= 0) stop 'maths error: bessel: underflow in output elements' 
  
end subroutine bessel

!==============================================================================

subroutine inv_rsis(A)

  !subroutine arguments
  double precision, dimension(:,:) :: A

  !local variables
  integer :: n, info, nmax, ifail, i, j
  integer, dimension(:), allocatable :: ipiv, work

  !maximum matrix size nmax x nmax
  nmax = 10

  !check size of A and fix dimensionality stuff
  if (size(A,1)/=size(A,2)) stop 'maths error: inv_rsis: non-square matrix'
  n = size(A,1)
  if (n>nmax) stop 'maths error: inv_rsis: matrix too large'
  allocate(ipiv(n))
  allocate(work(n**3))

  !factorise A
  call F07MDF('L', n, A, n, ipiv, work, n**3, info)

  !compute inverse
  if (info==0) then
     call F07MJF('L', n, A, n, ipiv, work, info)
  else
     stop 'maths error: inv_rsis: singular D, matrix inversion impossible'
  end if

  !write upper right triangle back
  do i = 2, n
     do j = 1, i-1
        A(j,i) = A(i,j)
     end do
  end do

  !clean up
  deallocate(ipiv)
  deallocate(work)

end subroutine inv_rsis

!==============================================================================

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

!==============================================================================

subroutine rnd_initialise()

  !initialises random number generation
  !this must be done prior to any looping over random value generation

  call G05CCF

end subroutine rnd_initialise

!==============================================================================

function argument(z)

  !returns argument of complex number z in the range -pi<arg<=pi

  !function arguments
  double complex :: z
  double precision :: argument
  
  !find argument
  argument = atan2(dimag(z), dble(z))

end function argument

!==============================================================================

function bess1(x)

  !function arguments
  double precision :: bess1, S17AFF, x

  !local variables
  integer :: ifail

  !fatal errors - leave to NAG

  !call NAG gamma function routine
  ifail = 0
  bess1 = S17AFF(x, ifail)

end function bess1

!==============================================================================

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

!==============================================================================

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

!==============================================================================

function gamma(x)

  !function arguments
  double precision :: gamma, S14AAF, x

  !local variables
  integer :: ifail

  !fatal errors - leave to NAG

  !call NAG gamma function routine
  ifail = 0
  gamma = S14AAF(x, ifail)

end function gamma

!==============================================================================

function machine_max()

  !returns the machine max number using NAG

  !function arguments
  double precision :: machine_max

  !local variables 
  double precision :: X02ALF

  machine_max = X02ALF()

end function machine_max

!==============================================================================

function machine_min()

  !returns the machine min number using NAG

  !function arguments
  double precision :: machine_min

  !local variables 
  double precision :: X02AKF

  machine_min = X02AKF()

end function machine_min

!==============================================================================

function machine_precision()

  !returns the machine precision using NAG machine constants routine

  !function arguments
  double precision :: machine_precision

  !local variables
  double precision :: X02AJF

  machine_precision = X02AJF()

end function machine_precision

!==============================================================================

function modulus(z)

  !returns modulus of complex z, uses NAG for highest precision and to
  !avoid underflow errors if im(z)<<re(z) or im(z)>>re(z)

  !function arguments
  double complex :: z
  double precision :: modulus

  !local variables
  integer :: ifail
  double precision :: A02ABF

  !fatal errors - leave to NAG

  !call NAG modulus routine
  ifail = 0
  modulus = A02ABF(dble(z), dimag(z))

end function modulus

!==============================================================================

function modx(a,b)

  !returns mod 360 difference of a and b

  !subroutine arguments
  double precision :: a, b, modx
  
  !local variables

  modx = modulo(a-b,360D0)
  if (modx>180D0) modx = 360D0 - modx

end function modx

!==============================================================================

function rnd_gauss(mu,sigma)

  !returns random number selected from gaussian distribution with
  !mean mu and standard deviation sigma

  !function arguments
  double precision :: rnd_gauss, mu, sigma

  !local variables
  double precision :: G05DDF

  !call NAG normal dist random generator routine
  rnd_gauss = G05DDF(mu, sigma)

end function rnd_gauss

!==============================================================================

function rnd_unit()

  !returns random value selected from range 0 - 1
  !dummy is dummy argument to NAG

  !function arguments
  double precision :: rnd_unit

  !local variables
  double precision :: G05CAF, dummy

  !call NAG unit dist random number generator routine
  dummy = 1D0
  rnd_unit = G05CAF(dummy)

end function rnd_unit

!==============================================================================

end module Maths



































