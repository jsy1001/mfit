!$Id: fit.f90,v 1.21 2006/08/17 13:45:40 jsy1001 Exp $

module Fit

use Maths
use Visibility
use Model

implicit none

private

!callable subroutines contained
!
!minimiser - wrapper to the minimising routine. this subroutine preps
!            the arrays into true parameter space and fixed model space
!            before calling the actual minimising alogrithm
!free_fit - frees storage allocated by minimiser
!likelihood - returns negative log likelihood of data given data, model
!prior - returns negative log prior given model 
!
!local subroutines contained
!
!posterior - returns negative log posterior of model and data

public :: minimiser, free_fit, gof, likelihood, prior

!public module variables contained:

public :: vis_data, triple_data, sel_wavebands
public :: fit_param, x_pos, x_bound

!data arrays
double precision, dimension(:,:), allocatable :: vis_data, triple_data

!wavebands included in data (for wavelength-dependent models)
double precision, dimension(:, :), allocatable :: sel_wavebands

!storage for model parameters (both variable and non-variable) during fit
double precision, dimension(:,:), allocatable :: fit_param

!positions of variable parameters in model_param and fit_param arrays
integer, dimension(:,:), allocatable :: x_pos

!1st and 2nd columns hold lower and upper bounds on the
!value of the parameter in question - used to prevent negative/other illegal
!values causing problems (esp in visibility calculations)
double precision, dimension(:,:), allocatable :: x_bound


!private module variables contained:

!typical length scales for variables
double precision, dimension(:), allocatable :: x_scale

contains

!==============================================================================
  
  subroutine minimiser(info, force_symm, sol, flag, desc, &
       hes, cov, cor, chisqrd, nlposterior, nlevidence)

    !On exit, flag holds success state of the minimisation:
    ! flag =  0:  Optimal solution found
    ! flag =  1:  Terminated with gradient small,
    !             parameter values probably optimal
    ! flag =  2:  Terminated with step size small,
    !             parameter values probably optimal
    ! flag =  3:  Lower point cannot be found,
    !             parameter values probably optimal
    ! flag =  4:  Iteration limit (150) exceeded
    ! flag =  5:  Too many large steps,
    !             function may be unbounded
    ! flag =  6:  position found by minimiser appears not to be a minimum
    !
    !Corresponding message text will be in info,
    !parameter values will be in fit_param
    !
    !force_symm is a logical that forces centrosymmetric only models

    !subroutine arguments
    character(len=128), intent(out) :: info
    character(len=55), dimension(:), allocatable, intent(out) :: desc
    double precision, dimension(:,:), allocatable, intent(out) :: sol, hes, &
         cov, cor
    double precision, intent(out) :: chisqrd, nlposterior, nlevidence
    integer, intent(out) :: flag
    logical, intent(in) :: force_symm

    !local variables
    integer :: i, j, k, iwave, ipar, n, lwork
    character(len=128) :: name, comp
    double precision :: P, P_l, P_u, P_i, P_j, deltai, deltaj, diff, det
    double precision, dimension(:), allocatable ::  x, temp_x, grad
    double precision, dimension(:), allocatable :: work
    logical :: illegal
    !for E04XAF
    logical :: hes_fail
    integer :: iwarn, ifail
    integer, dimension(:), allocatable :: pinf
    integer, dimension(1) :: iuser
    double precision :: objf
    double precision, dimension(:), allocatable :: x_us, hforw, hcntrl
    double precision, dimension(1) :: user

    interface
       subroutine PDA_UNCMND(n, x0, fcn, x, f, pinf, w, lw)
         integer, intent(in) :: n, lw
         double precision, dimension(n) :: x0, x
         external fcn
         double precision, intent(out) :: f
         integer, intent(out) :: pinf
         double precision, dimension(lw), intent(out) :: w
       end subroutine PDA_UNCMND

       subroutine E04XAF(msglvl, n, epsrf, x, mode, objfun, lhes, hforw, &
            objf, objgrd, hcntrl, hesian, &
            iwarn, work, iuser, user, pinf, ifail)
         integer, intent(in) :: msglvl, n, mode, lhes
         integer, intent(inout) :: ifail
         integer, intent(out) :: iwarn
         integer, dimension(n), intent(out) :: pinf
         integer, dimension(*), intent(in) :: iuser
         double precision, intent(in) :: epsrf
         double precision, dimension(n), intent(in) :: x
         double precision, dimension(n), intent(inout) :: hforw
         double precision, intent(out) :: objf
         double precision, dimension(n), intent(out) :: objgrd, hcntrl
         double precision, dimension(lhes, *), intent(out) :: hesian
         double precision, dimension(*), intent(out) :: work
         double precision, dimension(*), intent(in) :: user
         external objfun
       end subroutine E04XAF
   end interface

    !go through model arrays and count number of variable parameters
    n = 0
    do i = 1, size(model_param,1)
       do j = 1, size(model_param,2)
          if (model_prior(i,j) /= 0D0) n=n+1
       end do
    end do

    !allocate arrays associated with fit variables
    !x is a vector holding the SCALED values of the model parameters
    !that are allowed to vary,
    !x_pos holds the positions of the values in the model_param array
    !that defines the model.
    allocate(x(n))
    allocate(temp_x(n))
    allocate(grad(n))
    allocate(x_pos(n,2))
    allocate(desc(n))
    k = 0
    do i = 1, size(model_param,1) !loop over cpts
       do j = 1, size(model_param,2)
          if (model_prior(i,j) /= 0D0) then
             k = k+1
             x_pos(k,1) = i
             x_pos(k,2) = j
             desc(k) = model_desc(i,j)
          end if
       end do
    end do

    !check for illegal types of minimisation:
    illegal = .false.
    !must have at least 1 freedom to minimise with
    if (n==0) illegal = .true.
    !for single component model cannot vary r/theta (change in position
    !only constant phase offset) or B (flux is normalised anyway) 
    if (size(model_param,1)==1) then
       do iwave = 1, nwave
          if (model_prior(1, 2+4*model_wldep(1)*(iwave-1)) /= 0D0) &
               illegal = .true.
          if (model_prior(1, 3+4*model_wldep(1)*(iwave-1)) /= 0D0) &
               illegal = .true.
          if (model_prior(1, &
               6+4*model_wldep(1)*(nwave-1)+model_wldep(2)*(iwave-1)) /= 0D0) &
               illegal = .true.
       end do
    end if
    !for any component cannot vary theta if r is zero and not free to vary
    !(position angle has no effect if position radius fixed at zero)
    do i = 1, size(model_param,1)
       do iwave = 1, nwave
          ipar = 2+4*model_wldep(1)*(iwave-1)
          if ((model_param(i,ipar)==0D0) .and. (model_prior(i,ipar)==0D0) &
               .and. (model_prior(i,ipar+1)/=0D0)) illegal = .true.
       end do
    end do
    !for any component cannot vary phi if epsilon is unity and not free to vary
    !(orientation is meaningless if ellipse reduced to circular disc)
    do i = 1, size(model_param,1)
       do iwave = 1, nwave
          ipar = 7 + (4*model_wldep(1)+model_wldep(2))*(nwave-1) &
               + model_wldep(3)*(iwave-1)
          if ((model_param(i,ipar+2)==1D0) .and. (model_prior(i,ipar+2)==0D0) &
               .and. (model_prior(i,ipar+1)/=0D0)) illegal = .true.
       end do
    end do
    !doesn't make sense to have wavelength-dependent numerical CLV with
    !wavelength-dependent diameter
    if (size(clv_mvis, 2) > 1 .and. model_wldep(3) == 1) illegal = .true.
    if (illegal) goto 90

    !if centrosymmetric model is forced then must have 
    !eccentricity epsilon fixed to be unity and position radius fixed at zero
    if (force_symm .and. .not. symm) goto 92

    !allocate fit_param array
    !fit_param is identical to model_param, it holds the def'n of the model
    !with the difference that fit_param is designed to vary as the model is
    !adjusted to fit. x_pos is used to update the variable parameters in 
    !the array with their values from x (after undoing the scaling)
    if (allocated(fit_param)) deallocate(fit_param)
    allocate(fit_param(size(model_param,1),size(model_param,2)))
    fit_param = model_param

    !compile x_scale array
    !contains typical length scales for variables - based on prior width
    allocate(x_scale(n))
    do i = 1, n
       x_scale(i) = 1D-4*model_prior(x_pos(i,1),x_pos(i,2))
    end do

    !assign initial values of fit variables (=SCALED variable model parameters)
    do i = 1, n
       x(i) = model_param(x_pos(i,1),x_pos(i,2))/x_scale(i)
    end do

    !compile x_bound array
    !1st and 2nd columns hold lower and upper bounds on the
    !value of the parameter in question - used to prevent negative/other
    !illegal values causing problems (esp in visibility calculations)
    allocate(x_bound(n,2))
    !preset limits to values from the limits array used in model input
    do i = 1, n
       x_bound(i,1) = model_limits(x_pos(i,1),x_pos(i,2),1)
       x_bound(i,2) = model_limits(x_pos(i,1),x_pos(i,2),2)
       !alpha(1) in hestroffer model must be > 0
       if (trim(model_spec(x_pos(i,1),3)) == 'hestroffer') then
          if (model_wldep(4) == 1) then
             do iwave = 1, nwave
                if (x_pos(i,2) == (10+(4*model_wldep(1)+model_wldep(2) &
                     +3*model_wldep(3))*(nwave-1) &
                     + max_order*model_wldep(4)*(iwave-1))) x_bound(i,1) = 0D0
             end do
          else
             if (x_pos(i,2) == (10+(4*model_wldep(1)+model_wldep(2) &
                  +3*model_wldep(3))*(nwave-1))) x_bound(i,1) = 0D0
          end if
       end if
    end do

    !allocate sol, hes, cov, cor arrays
    allocate(sol(n,2))
    allocate(hes(n,n))
    allocate(cov(n,n))
    allocate(cor(n,n))
    sol = 0D0
    hes = 0D0
    cov = 0D0
    cor = 0D0

    print *, 'calculating initial goodness-of-fit...'
    !calculate goodness-of-fit statistic chi squared
    call gof(model_spec, fit_param, chisqrd)
    print *,' '
    print *,'initial chi squared =',real(chisqrd) 

    !call minimising algorithm
    info = ''
    lwork = n*(n+10)
    allocate(work(lwork))
    temp_x = x !starting point
    call PDA_UNCMND(n, temp_x, posterior, x, nlposterior, flag, work, lwork)
    deallocate(work)
    ! scaled solution is in x
    do i = 1, n
       sol(i,1) = x(i)*x_scale(i)
    end do

    !set return message
    select case(flag)
    case(-1)
       stop 'Insufficient workspace for PDA_UNCMND'
    case(0)
       info = 'Optimal solution found'
    case(1)
       info = 'Terminated with gradient small, parameter values probably optimal'
    case(2)
       info = 'Terminated with step size small, parameter values probably optimal'
    case(3)
       info = 'Lower point cannot be found, parameter values probably optimal'
    case(4)
       info = 'Iteration limit exceeded'
    case(5)
       info = 'Too many large steps, likelihood function may be unbounded'
    case default
       stop 'Unexpected status value from PDA_UNCMND'
    end select

    !XXX don't trust flag equal to 2 or 3 either?
    !XXX but almost never get 0 or 1, even when it works!
    if (flag < 0 .or. flag > 3) return ! not a minimum, so skip hessian calc.

    !find error in x position more rigorously by considering hessian
    !hessian matrix Aij has components (d2/(dyi*dyj))P({y}) where yi and yj
    !are parameters i and j, P({y}) is the neg log posterior for model with
    !parameters {y} = y1, y2, ... yi, ... yj, ... yN
    !estimate each element Aij via taylor 2nd order symmetric numeric 
    !differentiation method, for off-diagonals this is
    ![P(i+di,j+dj)+P(i-di,j-dj)-P(i+di,j-dj)-P(i-di,j+dj)]/[4di*dj], for diagonals
    ![P(i+di)+P(i-di)-2P(i)]/[di^2]
    !increments deltai, deltaj are chosen to minimise
    !uncertainty in second derivative in presence of roundoff error

    !Find hessian elements (nb Aij=Aji symmetric) by finite differences. Steps
    !are uniform in SCALED parameters.
    !We rescale the matrix elements to get the true Hessian that
    !relates to the unscaled parameters
    deltai = 1D0
    deltaj = 1D0
    do i = 1, n
       do j = 1, i
          temp_x = x
          if (j == i) then
             !P(i)
             call posterior(n, temp_x, P)
             !P(i-di)
             temp_x(i) = x(i)-deltai
             call posterior(n, temp_x, P_l)
             !P(i+di)
             temp_x(i) = x(i)+deltai
             call posterior(n, temp_x, P_u)
             !check we are at a minimum w.r.t. this parameter
             grad(i) = (P_u - P_l)/(2D0*deltai)
             if (abs(grad(i)) > 1D0) then
                flag = 6
                info = 'Final position not a local minimum'
                return
             end if
             !hessian by numeric method
             diff = P_u + P_l - (2D0*P)
             if (diff/=0D0) hes(i,i) = diff / ((deltai**2)*(x_scale(i)**2))
          else
             !P(i+di,j+dj)
             temp_x(i) = x(i)+deltai
             temp_x(j) = x(j)+deltaj
             call posterior(n, temp_x, P_u)
             !P(i+di,j-dj)
             temp_x(j) = x(j)-deltaj
             call posterior(n, temp_x, P_i)   
             !P(i-di,j-dj)
             temp_x(i) = x(i)-deltai
             call posterior(n, temp_x, P_l)
             !P(i-di,j+dj)
             temp_x(j) = x(j)+deltaj
             call posterior(n, temp_x, P_j)  
             !hessian by numeric method
             diff = P_u + P_l - P_i - P_j
             if (diff/=0D0) hes(i,j) = &
                  diff / ((4D0*deltai*deltaj)*(x_scale(i)*x_scale(j)))
             hes(j,i) = hes(i,j) !since hessian is symmetric
          end if
       end do
    end do

    !calculate Hessian using NAg routine
    allocate(pinf(n))
    allocate(x_us(n))
    allocate(hforw(n))
    allocate(hcntrl(n))
    allocate(work(n*(n+1)))
    !work with unscaled parameters, so we calculate true Hessian
    do i = 1, n
       x_us(i) = sol(i, 1)
       hforw(i) = x_scale(i)
    end do
    ifail = 1 !silent
    call E04XAF(1, n, 1D-9, x_us, 2, objfun, size(hes,1), hforw, &
         objf, grad, hcntrl, hes, &
         iwarn, work, iuser, user, pinf, ifail)
    if (ifail == 1) then
       stop 'E04XAF returned with ifail=1, indicates programming error'
    else if (ifail == 2) then
       print *, 'E04XAF had a problem with parameter(s):'
       hes_fail = .false.
       do i = 1, n
          if (pinf(i) /= 0) then
             print *, 'Parameter', i, ': diagnostic code = ', pinf(i)
          if (pinf(i) < 4) hes_fail = .true.
          end if
       end do
       if (hes_fail) then
          info = 'Problem(s) calculating Hessian'
          return
       end if
    end if
    do i = 1, n
       print *, hcntrl(i), hcntrl(i)/x_scale(i)
    end do
    deallocate(pinf)
    deallocate(hforw)
    deallocate(hcntrl)
    deallocate(work)

    !calculate covariance matrix (is inverse of the hessian)
    cov = hes
    call invdet_mat(cov, det) !inversion - also calculate determinant

    !calculate correlation matrix
    !correlation defined as: corr(x,y) = cov(x,y) / sqroot(var(x)*var(y))
    !noting that var(x) = cov(x,x)
    do i = 1, n
       do j = 1, i
          if (cov(i,i)*cov(j,j)>0D0) then
             cor(i,j) = cov(i,j)/sqrt(cov(i,i)*cov(j,j))
          else
             flag = 1
          end if
          cor(j,i) = cor(i,j)
       end do
    end do

    !calculate estimated error on solution
    !ideal case: error on a single parameter x is taken as the sqroot(cov(x,x))
    !negative cov(x,x) case: estimated error reported as zero (as occurs for
    !minimisation failure too)
    !also need to undo parameter scaling
    do i = 1, n
       if (cov(i,i)>0) then
          sol(i,2) = sqrt(cov(i,i))
       end if
    end do

    !calculate goodness-of-fit statistic chi squared
    do i = 1, n
       fit_param(x_pos(i,1),x_pos(i,2)) = sol(i,1)
    end do
    call gof(model_spec, fit_param, chisqrd)

    !estimate -log(evidence)
    !if det -ve returns zero
    if (det > 0) then
       nlevidence = nlposterior - 0.5D0*n*log(2D0*pi) + 0.5D0*log(det)
    else
       nlevidence = 0D0
    end if

    !clean-up and return
200 continue

    !deallocate local storage
    if (allocated(x)) deallocate(x)
    if (allocated(temp_x)) deallocate(temp_x)
    if (allocated(grad)) deallocate(grad)

    return

    !error trapping 
90  info = 'illegal freedom(s) in model'
    flag = -2
    goto 200
92  info = 'for vis/nvis data must have guaranteed symmetric model'
    flag = -3
    goto 200

  end subroutine minimiser

!==============================================================================

  subroutine gof(spec, param, chisqrd)

    !Goodness of fit information calculation
    !
    !Computes the chisqrd value of a model fit and data set
    !chisqrd = sum[i] { (data[i]-theory[i])^2/data_error[i]^2 }
    !chisqrd is thus almost identical to likelihood
    !could save code here and significantly combine this function
    !with the likelihood one but maybe less risky to keep them separated...
    !Could in future add code here to look up critical chi squared values
    !for num of data points and num of free parameters

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    double precision, intent(out) :: chisqrd

    !local variables
    integer :: i
    double precision :: lambda, delta_lambda, u, v, data_vis, data_vis_err
    double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
    double precision :: data_phase_err, model_vis, model_amp, model_phase, mjd
    double complex :: vis, vis1, vis2, vis3

    chisqrd = 0D0

    !sum over the visibility data points
    do i = 1, size(vis_data,1)
       
       !extract points
       lambda = vis_data(i,1)
       delta_lambda = vis_data(i,2)
       u = vis_data(i,3)
       v = vis_data(i,4)
       data_vis = vis_data(i,5)
       data_vis_err = vis_data(i,6)
       mjd = vis_data(i,7)

       !ignore if -ve or zero error:
       if (data_vis_err>0D0) then

          !compute model-predicted visibility amplitude squared
          vis = cmplx_vis(spec, param, lambda, delta_lambda, u, v, mjd)
          model_vis = (modulus(vis))**2D0

          !compute contribution to chisqrd
          chisqrd = chisqrd + &
               (((data_vis-model_vis)/data_vis_err)**2D0)

       end if

    end do

    !sum over triple product amplitude and phase data points
    do i = 1, size(triple_data,1)

       !extract points
       lambda = triple_data(i,1)
       delta_lambda = triple_data(i,2)
       u1 = triple_data(i,3)
       v1 = triple_data(i,4)
       u2 = triple_data(i,5)
       v2 = triple_data(i,6)
       data_amp = triple_data(i,7)
       data_amp_err = triple_data(i,8)
       data_phase = triple_data(i,9)
       data_phase_err = triple_data(i,10)
       mjd = triple_data(i,11)

       if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
          !compute model-predicted triple amplitude
          vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
          vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
          vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
               -(u1+u2), -(v1+v2), mjd)
          vis = vis1*vis2*vis3
          model_amp = modulus(vis)
          model_phase = rad2deg*argument(vis)
       end if

       !compute rmsdev and chisqrd contributions
       !nb phase calculations modulo 360 degrees
       if (data_amp_err>0D0) then
          chisqrd = chisqrd + &
               (((data_amp-model_amp)/(data_amp_err))**2D0)
       end if
       if (data_phase_err>0D0) then
          chisqrd = chisqrd + &
               (modx(data_phase,model_phase)/(data_phase_err))**2D0

       end if

    end do

  end subroutine gof

!==============================================================================

subroutine objfun(mode, n, x_us, objf, objgrd, nstate, iuser, user)

! Calculate -ln(posterior), given arguments supplied by E04XAF

!subroutine arguments - as specified by E04XAF
integer, intent(in) :: mode, n, nstate
double precision, dimension(n), intent(in) :: x_us
double precision, intent(out) :: objf
double precision, dimension(n), intent(out) :: objgrd !not used if mode /= 1
!we don't use iuser or user
integer, dimension(*), intent(in) :: iuser
double precision, dimension(*), intent(in) :: user

!local variables
double precision :: lhd, pri
integer :: i
logical :: violation

!check for range violations in current position
!if position is out of range then return very large energy 
violation = .false.
do i = 1, n !over variable parameters
   if (x_us(i) < x_bound(i,1)) violation = .true.
   if (x_us(i) > x_bound(i,2)) violation = .true.
end do
if (violation) then
   objf = 1.0d10
   return
end if 

!update fit_param with current position x_us
do i = 1, n
   fit_param(x_pos(i,1),x_pos(i,2)) = x_us(i)
end do

!form E = neg log posterior = neg log likelihood + neg log prior
lhd = likelihood(vis_data, triple_data, model_spec, fit_param)
pri = prior(x_pos, fit_param, model_param, model_prior)
objf = lhd + pri

end subroutine objfun

!==============================================================================

subroutine posterior(n, x, post)

! Black box minimisers require this function have specific arguments
! Everything else needed comes from module variables

!subroutine arguments
integer, intent(in) :: n
double precision, dimension(n), intent(in) :: x
double precision, intent(out) :: post

!local variables
double precision :: lhd, pri
integer :: i
logical :: violation

!check for range violations in current position
!if position is out of range then return very large energy 
violation = .false.
do i = 1, n !over variable parameters
   if (x(i)*x_scale(i) < x_bound(i,1)) violation = .true.
   if (x(i)*x_scale(i) > x_bound(i,2)) violation = .true.
end do
if (violation) then
   post = 1.0d10
   return
end if 

!update fit_param with current position x, undoing scaling
do i = 1, n
   fit_param(x_pos(i,1),x_pos(i,2)) = x(i)*x_scale(i)
end do

!form E = neg log posterior = neg log likelihood + neg log prior
lhd = likelihood(vis_data, triple_data, model_spec, fit_param)
pri = prior(x_pos, fit_param, model_param, model_prior)
post = lhd + pri

end subroutine posterior

!==============================================================================

function likelihood(vis_data, triple_data, model_spec, param)

!Computes the negative log likelihood of data given model

!subroutine arguments
double precision :: likelihood
double precision, dimension(:,:), intent(in) :: vis_data, triple_data
character(len=128), dimension(:,:), intent(in) :: model_spec
double precision, dimension(:,:) :: param

!local variables
integer :: i
double precision :: lambda, delta_lambda, u, v, data_vis, data_vis_err
double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
double precision :: data_phase_err, model_vis, model_amp, model_phase, mjd
double complex :: vis, vis1, vis2, vis3

likelihood = 0D0

!sum over the visibility data points
do i = 1, size(vis_data,1)
   
   !extract points
   lambda = vis_data(i,1)
   delta_lambda = vis_data(i,2)
   u = vis_data(i,3)
   v = vis_data(i,4)
   data_vis = vis_data(i,5)
   data_vis_err = vis_data(i,6)
   mjd = vis_data(i,7)

   if (data_vis_err>0D0) then
      
      !compute model-predicted visibility amplitude squared
      vis = cmplx_vis(model_spec, param, lambda, delta_lambda, u, v, mjd)
      model_vis = (modulus(vis)**2D0)
      
      !compute contribution to likelihood
      likelihood = likelihood + &
           (((data_vis-model_vis)**2D0)/(2D0*(data_vis_err**2D0)))

   end if

end do

!sum over triple product amplitude and phase data points
do i = 1, size(triple_data,1)

   !extract points
   lambda = triple_data(i,1)
   delta_lambda = triple_data(i,2)
   u1 = triple_data(i,3)
   v1 = triple_data(i,4)
   u2 = triple_data(i,5)
   v2 = triple_data(i,6)
   data_amp = triple_data(i,7)
   data_amp_err = triple_data(i,8)
   data_phase = triple_data(i,9)
   data_phase_err = triple_data(i,10)
   mjd = triple_data(i,11)

   if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
      !compute model-predicted triple product
      vis1 = cmplx_vis(model_spec, param, lambda, delta_lambda, u1, v1, mjd)
      vis2 = cmplx_vis(model_spec, param, lambda, delta_lambda, u2, v2, mjd)
      vis3 = cmplx_vis(model_spec, param, lambda, delta_lambda, &
           -(u1+u2), -(v1+v2), mjd)
      vis = vis1*vis2*vis3
      model_amp = modulus(vis)
      model_phase = rad2deg*argument(vis)
   end if

   !compute likelihood contributions
   !nb phase calculations modulo 360 degrees
   if (data_amp_err>0D0) then
      likelihood = likelihood + &
           (((data_amp-model_amp)**2D0)/(2D0*(data_amp_err**2D0)))
   end if
   if (data_phase_err>0D0) then
      likelihood = likelihood + &
           (((modx(data_phase,model_phase))**2D0)/(2D0*(data_phase_err**2D0)))
   end if
end do

end function likelihood

!==============================================================================

function prior(x_pos, fit_param, model_param, model_prior)

!function arguments
double precision :: prior
integer, dimension(:,:), intent(in) :: x_pos
double precision, dimension(:,:), intent(in) :: fit_param, model_param
double precision, dimension(:,:), intent(in) :: model_prior

!local variables
integer :: i
double precision :: value, value0, err

prior = 0D0

!loop over the variable parameters
do i = 1, size(x_pos,1)

   value0 = model_param(x_pos(i,1),x_pos(i,2)) !parameter original value
   err = model_prior(x_pos(i,1),x_pos(i,2)) !parameter prior width
   value = fit_param(x_pos(i,1),x_pos(i,2)) !parameter current value as it is varied

   !calculate normalised contribution to the prior
   prior = prior + (((value-value0)**2D0)/(2D0*(err**2D0))) + log(err*sqrt(2D0*pi))

end do

end function prior

!==============================================================================

subroutine free_fit()

  !Deallocate fitting storage
  if (allocated(fit_param)) deallocate(fit_param)
  if (allocated(x_pos)) deallocate(x_pos)
  if (allocated(x_bound)) deallocate(x_bound)
  if (allocated(x_scale)) deallocate(x_scale)

end subroutine free_fit

!==============================================================================

end module Fit
