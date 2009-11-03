!$Id: visibility.f90,v 1.16 2009/11/03 17:26:27 jsy1001 Exp $

module Visibility

  use Maths
  use Search
  use Model

  implicit none

  private

  !public function contained:
  !cmplx_vis - returns complex visibility for the model
  !            by summing over model components present it
  !            calls the functions below in the summation.

  public :: cmplx_vis

  !private functions contained:
  !
  !return transform function (on origin i.e. real) as
  !function of rho and d_lambda
  !all arguments to double precision where possible
  !
  !uniform           happy - thesis ok
  !taylor            happy - thesis ok
  !square_root       not happy - db thesis seems incorrect
  !gaussian          happy - thesis ok
  !hestroffer        happy - agrees with paper
  !gauss_hermite     not happy - db thesis seems incorrect
  !two_layer         happy-ish can take a long time to run
  !clvvis            tested

contains

  !===========================================================================
  
  !! Return complex visibility for multiple component model
  !!
  !! Model component details are supplied in 
  !! the spec and param arrays. u, v must be supplied in metres, lambda in nm.
  !! mjd argument gives epoch of observation as Modified Julian Day.
  function cmplx_vis(spec, param, lambda, delta_lambda, u, v, mjd)

    !function arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    double precision, intent(in) :: lambda, delta_lambda, u, v, mjd
    double complex :: cmplx_vis

    !parameters
    double precision, parameter :: sig = 0.1D0  !waveband must match to sig nm

    !local variables
    integer :: i, iwave, navg, iavg, num_comps, order, iwb, ipar, relto
    double precision :: r, theta_t0, t0, dtheta_dt, theta
    double precision :: B, a, phi, B_total, lambda1
    double precision :: epsilon, rho, F, x1, x2, x3
    double precision :: alpha(0:10)
    double complex :: sumvis

    !set number of components
    num_comps = size(spec,1)

    !trap zero baseline case
    if ((u == 0D0) .and. (v == 0D0)) then 
       cmplx_vis = cmplx(1D0,0D0)
       B_total = 1D0

    else

       !find index into model_wb
       if (nwave > 1) then
          if (.not. allocated(model_wb)) &
               stop 'Wavebands not supplied to read_model'
          iwb = -1
          do iwave = 1, size(model_wb, 1)
             if (lambda  >= model_wb(iwave, 1)-sig &
                  .and. lambda <= model_wb(iwave, 1)+sig &
                  .and. delta_lambda >= model_wb(iwave, 2)-sig &
                  .and. delta_lambda <= model_wb(iwave, 2)+sig) then
                iwb = iwave
                exit
             end if
          end do
          if (iwb == -1) stop 'Datum has no matching waveband'
       else
          iwb = 1
       end if

       !clear real and imag vis
       cmplx_vis = 0D0
       B_total = 0D0

       !choose number of wavelengths (within band) to average over
       navg = 1D2*delta_lambda/lambda
       if (navg < 2) navg = 2

       !sum over components
       do i = 1, num_comps

          !get position r and theta, brightness I
          !major axis a, orientation phi, ellipticity epsilon
          !ld order order, ld parameters alpha
          relto = model_pos_relto(i)
          if (relto >= 1 .and. relto <= num_comps) then
             r = mas2rad*param(i, 2+4*model_wldep(1)*(iwb-1)) &
                  *param(relto, 2+4*model_wldep(1)*(iwb-1))
             theta_t0 = deg2rad*param(i, 3+4*model_wldep(1)*(iwb-1)) &
                  *param(relto, 3+4*model_wldep(1)*(iwb-1))
          else
             r = mas2rad*param(i, 2+4*model_wldep(1)*(iwb-1))
             theta_t0 = deg2rad*param(i, 3+4*model_wldep(1)*(iwb-1))
          end if
          t0 = param(i, 4+4*model_wldep(1)*(iwb-1))
          dtheta_dt = deg2rad*param(i, 5+4*model_wldep(1)*(iwb-1))
          if (dtheta_dt /= 0D0 .and. mjd < 0D0) &
               stop 'Need observation time for time-dependent model'
          theta = theta_t0 + (mjd - t0)*dtheta_dt
          B = param(i, 6+4*model_wldep(1)*(nwave-1)+model_wldep(2)*(iwb-1))

          ipar = 7 + (4*model_wldep(1)+model_wldep(2))*(nwave-1) &
               + 3*model_wldep(3)*(iwb-1)
          a = mas2rad*param(i, ipar)
          phi = deg2rad*param(i, ipar+1)
          epsilon = param(i, ipar+2)
          order = int(param(i, 1)) !ld order *for this cpt*

          !sum up brightnesses (need later to normalise whole model)
          B_total = B_total + B

          !configure array of ld parameters
          !nb alpha(0) set differently for different models
          if (trim(spec(i,3)) == 'two-layer') then
             ipar = 10 + (4*model_wldep(1) + model_wldep(2) &
                  + 3*model_wldep(3))*(nwave-1)
             alpha(1:4) = param(i, ipar:ipar+3)
             alpha(5) = param(i, ipar+4+model_wldep(4)*(iwb-1))
          else
             ipar = 10 + (4*model_wldep(1) + model_wldep(2) &
                  + 3*model_wldep(3))*(nwave-1) &
                  + max_order*model_wldep(4)*(iwb-1)
             alpha(1:order) = param(i, ipar:ipar+order-1)
          end if

          !loop over wavelengths within bandpass (will normalise later)
          sumvis = 0D0
          do iavg = 1, navg
             lambda1 = lambda - delta_lambda/2D0 + &
                  (iavg-1)*delta_lambda/(navg-1)

             !calculate parameter rho
             x1 = (epsilon**2)*(((u*cos(phi))-(v*sin(phi)))**2)
             x2 = (((v*cos(phi))+(u*sin(phi)))**2)
             rho = (sqrt(x1 + x2))/(lambda1*(1D-9))

             !deal with point shape, everything else is elliptical
             !(disc case already set to have ellipticity 1.0)
             if ((trim(spec(i,2)) == 'point') .or. (a==0)) then
                F = 1D0
             else
                select case (trim(spec(i,3)))
                case ('uniform')
                   F = uniform(a, rho)
                case ('taylor')
                   !alpha is array of taylor coeffs 0-20 (0th defined as -1)
                   alpha(0) = -1D0
                   F = taylor(a, rho, order, alpha)
                case ('square-root')
                   !alpha(1) and alpha(2) are sqroot coeffs alpha and beta
                   F = square_root(a, rho, alpha(1), alpha(2))
                case ('gaussian')
                   F = gaussian(a, rho)
                case ('hestroffer')
                   !alpha(1) is hestroffer power law parameter
                   F = hestroffer(a, rho, alpha(1))
                case ('gauss-hermite')
                   !alpha is array of gauss-hermite coeffs -20 (0th def'd as 1)
                   alpha(0) = 1D0
                   F = gauss_hermite(a, rho, order, alpha)
                case ('two-layer')
                   !alpha is array of two-layer atmosphere model parameters
                   F = two_layer(a,rho,lambda1*1D-9,alpha)
                case default
                   !must be numerical CLV
                   F = clvvis(a, rho, iwb)
                end select
             end if

             !calculate real and imaginary parts (incorporate off-centre
             !position of component r, theta and brightness B)
             x1 = (2D0*pi*r*((u*sin(theta))+(v*cos(theta))))/(lambda1*(1D-9))
             x2 = B*F*cos(x1)
             x3 = B*F*sin(x1)

             !add this component visibility to sum
             sumvis = sumvis + cmplx(x2,x3)

          end do

          !add normalised sum over wavelengths to sum over components
          cmplx_vis = cmplx_vis + sumvis/navg

       end do

    end if

    !normalise whole model visiblity
    if (B_total /= 0D0) cmplx_vis = cmplx_vis/B_total

    return

  end function cmplx_vis

  !============================================================================

  function uniform(a, rho)

    !function arguments
    double precision :: uniform, a, rho

    !evaluate visibility
    uniform = (2D0*bess1(pi*a*rho))/(pi*a*rho)

  end function uniform

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

  function clvvis(a, rho, iwb)

    !subroutine arguments
    double precision :: clvvis, a, rho
    integer :: iwb

    !local variables
    real bas_scaled, frac
    integer ifind

    bas_scaled = a*rho*1D-6/(mas2rad*clv_mdiam) !scaled baseline in Mega-lambda

    ifind = locate(clv_mbase, bas_scaled)
    if (ifind == 0) then
       clvvis = clv_mvis(1, iwb)
    else if (ifind == nxsiz+1) then
       clvvis = clv_mvis(nxsiz+1, iwb)
    else
       frac = (bas_scaled - clv_mbase(ifind)) / &
            (clv_mbase(ifind+1) - clv_mbase(ifind))
       clvvis = clv_mvis(ifind, iwb) &
            + frac*(clv_mvis(ifind+1, iwb) - clv_mvis(ifind, iwb))
    end if

  end function clvvis

  !============================================================================

end module Visibility






