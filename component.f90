module Component

  use Maths

  implicit none

  !public functions contained:
  !
  !return transform function (on origin i.e. real) as
  !function of rho and d_lambda
  !all arguments to double precision where possible
  !
  !uniform           happy - thesis ok
  !thin-ring         JAG added
  !taylor            happy - thesis ok
  !square_root       not happy - db thesis seems incorrect
  !gaussian          happy - thesis ok
  !hestroffer        happy - agrees with paper
  !gauss_hermite     not happy - db thesis seems incorrect
  !two_layer         happy-ish can take a long time to run
  !clvvis            tested

contains

  !============================================================================

  function uniform(a, rho)

    !function arguments
    double precision :: uniform, a, rho

    !evaluate visibility
    uniform = (2D0*bess1(pi*a*rho))/(pi*a*rho)

  end function uniform

  !============================================================================

  function thin_ring(a, rho)

    !function arguments
    double precision :: thin_ring, a, rho

    !evaluate visibility
    thin_ring = bess0(pi*a*rho)

  end function thin_ring

  !============================================================================

  function taylor(a, rho, nmax, alpha)

    !function arguments, note alpha is an array of coefficients
    double precision :: taylor, a, rho
    integer :: nmax
    double precision, dimension(0:10) :: alpha

    !local variables
    double complex :: x1, x2, y1, y2, z1, z2, z3, z4, z5, z6
    integer :: i, p, n
    double precision, dimension(:), allocatable :: even_output
    double precision, dimension(:), allocatable :: odd_output
    double precision, dimension(:), allocatable :: bessel_data

    !set up bessel array space
    allocate(bessel_data(0:nmax))

    !bessel() can only give separate odd/even bessel data arrays
    !bessel call is done in bulk for highest efficiency
    !no need for odd output call if nmax=0 case
    call bessel(0D0, pi*a*rho, 2+floor(dble(nmax)/2D0), &
         even_output)
    if (nmax > 0) call bessel(0.5D0, pi*a*rho, 1+ceiling(dble(nmax)/2D0), &
         odd_output)

    !interlace odd/even arrays to get array of J[(p/2)+1](pi*a*rho)
    !values that are referenced by p value
    do i = 2, size(even_output)
       bessel_data((2*i)-4) = even_output(i)
    end do
    do i = 2, size(odd_output)
       bessel_data((2*i)-3) = odd_output(i)
    end do

    !sigma over n
    x1 = 0
    x2 = 0
    do n = 0, nmax

       !sigma over p
       y1 = 0
       y2 = 0
       do p = 0, n
          z1 = dble((-1D0)**(p))
          z2 = (2D0)**(dble(p)/2D0)
          z3 = comb(n,p)
          z4 = gamma((dble(p)/2D0)+1D0)
          z5 = bessel_data(p) !J[p/2+1](pi*a*rho)
          z6 = (pi*a*rho)**((dble(p)/2D0)+1D0)
          y1 = y1 + (z1*z2*z3*z4*(z5/z6))
          y2 = y2 + ((z1*z3)/(dble(p)+2D0))
       end do

       x1 = x1 + ((alpha(n))*(y1))
       x2 = x2 + ((alpha(n))*(y2))

    end do

    !deallocate
    deallocate(bessel_data)
    deallocate(even_output)
    if (allocated(odd_output)) deallocate(odd_output)

    taylor = (x1/x2)

  end function taylor

  !============================================================================

  function square_root(a, rho, alpha, beta)

    !subroutine arguments, alpha and beta are single parameters
    !alpha is (1-mu) coeff and beta is (1-root(mu)) coeff
    double precision :: square_root, a, rho, alpha, beta

    !local variables
    double precision :: x1, x2, x3, x4
    double precision, dimension(:), allocatable :: bs1, bs2, bs3

    !get bessel function values
    call bessel(0D0, pi*a*rho, 2, bs1)
    call bessel(0.5D0, pi*a*rho, 2, bs2)
    call bessel(0.25D0, pi*a*rho, 2, bs3)

    !calculate visibility
    x1 = ((1D0-alpha-beta)*(bs1(2)))/(pi*a*rho)
    x2 = (alpha*(sqrt(pi/2D0))*bs2(2))/((pi*a*rho)**(1.5D0))
    x3 = (beta*((2D0)**(0.25D0))*(gamma(1.25D0))*bs3(2))/((pi*a*rho)**(1.25D0))
    x4 = (0.5D0)-((1D0/6D0)*alpha)-((1D0/10D0)*beta)

    square_root = (x1+x2+x3)/x4

    !deallocate bessel data arrays
    deallocate(bs1)
    deallocate(bs2)
    deallocate(bs3)

  end function square_root

  !============================================================================

  function gaussian(a, rho)

    !subroutine arguments
    double precision :: gaussian, a, rho

    gaussian = exp(-((pi*a*rho)**2D0)/(4D0*log(2D0)))

  end function gaussian

  !============================================================================

  function hestroffer(a, rho, alpha)

    !subroutine arguments, note alpha is a single parameter
    double precision :: hestroffer, a, rho, alpha

    !local variables
    integer :: num
    double precision :: order, frac, x1, x2, x3
    double precision, dimension(:), allocatable :: bessel_data

    order = (alpha/2D0)+1D0
    frac = order - floor(order)
    num = floor(order) + 1
    call bessel(frac, pi*a*rho, num, bessel_data)

    x1 = gamma((alpha/2D0)+2D0)
    x2 = bessel_data(num)
    x3 = ((pi*a*rho)/2D0)**((alpha/2D0)+1D0)

    !deallocate bessel array
    deallocate(bessel_data)

    hestroffer = x1*(x2/x3)

  end function hestroffer

  !============================================================================

  function gauss_hermite(a, rho, nmax, alpha)

    !subroutine arguments, note alpha is an array of coefficients
    double precision :: gauss_hermite, a, rho
    integer :: nmax
    double precision, dimension(0:10) :: alpha

    !local variables
    double precision, dimension(:), allocatable :: laguerre_data
    double precision :: x1, x2, y1, y2, z1, z2, z3, z4, z5, z6, z7
    integer :: n, t

    !get laguerre coefficients
    call laguerre(nmax, ((pi*rho*a)**2D0)/(4D0*log(2D0)), laguerre_data)

    !sigma over n
    x1 = 0
    x2 = 0
    do n = 0, nmax  

       !sigma over t
       y1 = 0
       y2 = 0
       do t = 0, n        
          z1 = dble((-1D0)**t)
          z2 = gamma(dble(t)+1D0)
          z3 = (2D0)**(3D0*dble(t))
          z4 = fact(2*t)
          z5 = fact(n-t)
          z6 = exp(-((pi*rho*a)**2D0)/(4D0*log(2D0)))
          z7 = laguerre_data(t)
          y1 = y1 + ((z1*z2*z3*z6*z7)/(z4*z5))
          y2 = y2 + ((z1*z2*z3)/(z4*z5))
       end do

       x1 = x1 + (alpha(n)*y1*(dble((-1D0)**n))*(fact(2*n)))
       x2 = x2 + (alpha(n)*y2*(dble((-1D0)**n))*(fact(2*n)))

    end do

    deallocate(laguerre_data)

    gauss_hermite = (x1/x2)

  end function gauss_hermite

  !============================================================================

  function two_layer (a, rho, lambda, alpha)

    !subroutine arguments, note alpha is an array of parameters
    !alpha = <steps> <temp1> <temp2> <r2/r1> <tau>
    double precision :: two_layer, a, rho, lambda
    double precision :: alpha(0:10) !though we only care about 1:4

    !constants
    double precision :: h = 6.63D-34  !plancks constant
    double precision :: c = 3.0D+08   !speed of light
    double precision :: k = 1.38D-23  !boltzmann constant

    !function variables
    integer :: i, steps
    double precision :: step_size, cosine, r, sample, norm
    double precision :: p1, p2
    steps = int (alpha (1))

    !calculate planck functions
    p1 = 2 * h * c * c / lambda**5 / (exp (h * c / lambda / k / alpha (2)) -1)
    p2 = 2 * h * c * c / lambda**5 / (exp (h * c / lambda / k / alpha (3)) -1)

    !calculate intensities and multiply by J0(2PIqr)r as we integrate
    step_size = a * alpha (4) / steps
    two_layer = 0
    norm = 0
    do i = 1, steps-1, 2
       r = i * step_size
       sample = calc_sample (i) * r
       norm = norm + sample
       two_layer = two_layer + sample * bess0 (2 * pi * rho * r)
    end do
    two_layer = two_layer * 2
    norm = norm * 2
    do i = 2, steps-1, 2
       r = i * step_size
       sample = calc_sample (i) * r
       norm = norm + sample
       two_layer = two_layer + sample * bess0 (2 * pi * rho * r)
    end do
    two_layer = two_layer * 2
    norm = norm * 2

    r = steps * step_size
    sample = calc_sample (steps) * r
    norm = norm + sample
    two_layer = two_layer + sample * bess0 (2 * pi * rho * r)

    two_layer = two_layer / norm

  contains

    function calc_sample (step)
      integer :: step
      double precision :: calc_sample, cosine

      cosine = sqrt (1 - (dble (step) / (steps+1))**2)
      if (r < a) then
         calc_sample = exp (-alpha (5) / cosine)
         calc_sample = p1 * calc_sample + p2 * (1 - calc_sample)
      else
         calc_sample = p2 * (1 - exp (-2 * alpha (5) / cosine))
      end if
    end function calc_sample

  end function two_layer

  !============================================================================

end module Component
