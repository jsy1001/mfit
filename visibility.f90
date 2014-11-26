!$Id: visibility.f90,v 1.18 2013/02/18 13:18:16 jsy1001 Exp $

module Visibility

  use Maths, only: mas2rad, deg2rad
  use Search, only: locate
  use Model
  use Component

  implicit none

  !public function contained:
  !cmplx_vis - returns complex visibility for the model
  !            by summing over model components present it
  !            calls the functions below in the summation.

  !private function contained:
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
		case ('thin-ring')
                   F = thin_ring(a, rho)
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

  function clvvis(a, rho, iwb)

    !subroutine arguments
    double precision :: clvvis, a, rho
    integer :: iwb, j

    !local variables
    real bas_scaled, frac
    integer ifind

    bas_scaled = a*rho*1D-6/(mas2rad*clv_mdiam) !scaled baseline in Mega-lambda
    if(size(clv_mvis, 2) == 1) then
       j = 1  !wavelength-independent CLV
    else
       j = iwb
    end if

    ifind = locate(clv_mbase, bas_scaled)
    if (ifind == 0) then
       clvvis = clv_mvis(1, j)
    else if (ifind == nxsiz+1) then
       clvvis = clv_mvis(nxsiz+1, j)
    else
       frac = (bas_scaled - clv_mbase(ifind)) / &
            (clv_mbase(ifind+1) - clv_mbase(ifind))
       clvvis = clv_mvis(ifind, j) &
            + frac*(clv_mvis(ifind+1, j) - clv_mvis(ifind, j))
    end if

  end function clvvis

  !============================================================================

end module Visibility






