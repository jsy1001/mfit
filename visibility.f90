!$Id: visibility.f90,v 1.9 2005/01/06 18:45:12 jsy1001 Exp $

module Visibility

!callable function contained
!cmplx_vis - returns complex visibility for the model
!            by summing over model components present it
!            calls the functions below in the summation.  
!
!local functions contained
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
!clvvis            tested

use Maths !picks up pi parameter from here
use Model

implicit none

contains
!==============================================================================

function cmplx_vis(spec, param, lambda, delta_lambda, u, v)

  !returns complex visibility real and imaginary parts for
  !multiple component model. model component details stored in 
  !the model_spec and model_param arrays.
  !u, v must be supplied in metres, lambda in nm
  
  !subroutine arguments
  character(len=128), dimension(:,:), intent(in) :: spec
  double precision, dimension(:,:), intent(in) :: param
  double precision, intent(in) :: lambda, delta_lambda, u, v
  double complex :: cmplx_vis
  
  !parameters
  double precision, parameter :: sig = 0.1D0  !waveband must match to sig nm

  !local variables
  integer :: i, iwave, nwave, num_comps, nmax, iwb
  double precision :: r, theta, B, a, phi, B_total, lambda1
  double precision :: epsilon, rho, F, x1, x2, x3
  double precision, dimension(0:10) :: alpha
  double complex :: sumvis

  !set number of components
  num_comps = size(spec,1)
  
  !trap zero baseline case
  if ((u == 0D0) .and. (v == 0D0)) then 
     cmplx_vis = cmplx(1D0,0D0)
     B_total = 1D0

  else
  
     !clear real and imag vis
     cmplx_vis = 0D0
     B_total = 0D0
  
     !choose number of wavelengths to average over
     nwave = 1D2*delta_lambda/lambda
     if (nwave < 2) nwave = 2
     
     !sum over components
     do i = 1, num_comps
        
        !get position r and theta, brightness I
        !major axis a, orientation phi, ellipticity epsilon
        !ld order nmax, ld parameters alpha
        r = mas2rad*param(i,2)
        theta = deg2rad*param(i,3)
        B = param(i,4)

        a = mas2rad*param(i,5)
        phi = deg2rad*param(i,6)
        epsilon = param(i,7)
        nmax = int(param(i,1))
        
        !sum up brightnesses (need later to normalise whole model)
        B_total = B_total + B

        !configure array of ld parameters
        !nb alpha(0) set differently for different models
        alpha(1:10) = param(i,8:nmax+7)

        !find index into clv_wb
        if (spec(i,3)(1:1) == '<') then
           if (size(clv_inten, 2) > 1) then
              !wavelength-dependent CLV
              iwb = -1
              do iwave = 1, size(clv_wb, 1)
                 if (lambda  >= clv_wb(iwave, 1)-sig &
                      .and. lambda <= clv_wb(iwave, 1)+sig &
                      .and. delta_lambda >= clv_wb(iwave, 2)-sig &
                      .and. delta_lambda <= clv_wb(iwave, 2)+sig) then
                    iwb = iwave
                    exit
                 end if
              end do
              if (iwb == -1) stop 'Datum has no matching waveband'
           else
              iwb = 1
           end if
        end if
        
        !loop over wavelengths within bandpass (will normalise later)
        sumvis = 0D0
        do iwave = 1, nwave
           lambda1 = lambda - delta_lambda/2D0 + (iwave-1)*delta_lambda/(nwave-1)
           
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
                 F = taylor(a, rho, nmax, alpha)
              case ('square-root')
                 !alpha(1) and alpha(2) are sqroot coeffs alpha and beta
                 F = square_root(a, rho, alpha(1), alpha(2))
              case ('gaussian')
                 F = gaussian(a, rho)
              case ('hestroffer')
                 !alpha(1) is hestroffer power law parameter
                 F = hestroffer(a, rho, alpha(1))
              case ('gauss-hermite')
                 !alpha is array of gauss-hermite coeffs -20 (0th defined as 1)
                 alpha(0) = 1D0
                 F = gauss_hermite(a, rho, nmax, alpha)
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
        cmplx_vis = cmplx_vis + sumvis/nwave

     end do
   
  end if 

  !normalise whole model visiblity
  if (B_total /= 0D0) cmplx_vis = cmplx_vis/B_total

  return

end function cmplx_vis

!==============================================================================

function uniform(a, rho)

  !function arguments
  double precision :: uniform, a, rho

  !evaluate visibility
  uniform = (2D0*bess1(pi*a*rho))/(pi*a*rho)
  
end function uniform

!==============================================================================

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

!==============================================================================

function square_root(a, rho, alpha, beta)

  !subroutine arguments, alpha and beta are single parameters
  !alpha is (1-mu) coeff and beta is (1-root(mu)) coeff
  double precision :: square_root, a, rho, alpha, beta

  !local variables
  double precision :: x1, x2, x3, x4
  double precision, dimension(:), allocatable :: bs1, bs2, bs3

  !get bessel function values
  call bessel(0D0, pi*a*rho, 1, bs1)
  call bessel(0.5D0, pi*a*rho, 1, bs2)
  call bessel(0.25D0, pi*a*rho, 1, bs3)

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

!==============================================================================

function gaussian(a, rho)

  !subroutine arguments
  double precision :: gaussian, a, rho

  gaussian = exp(-((pi*a*rho)**2D0)/(4D0*log(2D0)))

end function gaussian

!==============================================================================

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

!==============================================================================

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

!==============================================================================

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

!==============================================================================

end module Visibility






