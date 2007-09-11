!$Id: fit.f90,v 1.17.2.1 2007/09/11 11:09:28 jsy1001 Exp $

module Fit

!callable subroutines contained
!
!minimiser - wrapper to the minimising routine. this subroutine preps
!            the arrays into true parameter space and fixed model space
!            before calling the actual minimising alogrithm
!likelihood - returns negative log likelihood of data given data, model
!prior - returns negative log prior given model 
!
!local subroutines contained
!
!posterior - returns negative log posterior of model and data

use Maths
use Visibility
use Model

implicit none

!module variables contained:

!data arrays
double precision, dimension(:,:), allocatable :: vis_data, triple_data

!wavebands included in data (for wavelength-dependent models)
double precision, dimension(:, :), allocatable :: sel_wavebands

!storage for model parameters (both variable and non-variable) during fit
double precision, dimension(:,:), allocatable :: fit_param

!positions of variable parameters in model_param and fit_param arrays
integer, dimension(:,:), allocatable :: x_pos

!1st column is of typical length scales for variables and set equal to the
!prior width. 2nd and 3rd columns hold lower and upper bounds on the
!value of the parameter in question - used to prevent negative/other illegal
!values causing problems (esp in visibility calculations)
double precision, dimension(:,:), allocatable :: x_info

contains

!==============================================================================
  
  subroutine minimiser(info, force_symm, sol, ifail, desc, &
       hes, cov, cor, chisqrd, nlposterior, nlevidence)

    !On exit, ifail holds success state:
    ! ifail = -3 : vis/nvis format data & model not guaranteed symmetric
    ! ifail = -2 : illegal freedom(s) in model
    ! ifail = -1 : problem with Hessian calculation
    ! ifail =  0 : success
    !Other positive values indicate minimiser failure -
    !refer to E04JYF NAg documentation
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
    integer, intent(out) :: ifail
    logical, intent(in) :: force_symm

    !local variables
    integer :: i, j, k, iwave, ipar, n, liwork, lwork
    character(len=128) :: name, comp
    double precision :: P, P_l, P_u, P_i, P_j, deltai, deltaj, diff, det
    double precision, dimension(:), allocatable ::  x, temp_x, work
    integer, dimension(:), allocatable :: iwork
    logical :: illegal
    integer, dimension(1) :: iuser !dummy
    double precision, dimension(1) :: user !dummy

    interface
       subroutine E04JYF(n, ibound, funct1, bl, bu, x, f, iw, liw, w, lw, &
            iuser, user, ifail)
         integer, intent(in) :: n, ibound, liw, lw
         integer, dimension(liw), intent(out) :: iw
         integer, dimension(1), intent(in) :: iuser
         integer, intent(inout) :: ifail
         external funct1
         double precision, dimension(n), intent(inout) :: bl, bu, x
         double precision, intent(out) :: f
         double precision, dimension(lw), intent(out) :: w
         double precision, dimension(1), intent(in) :: user
       end subroutine E04JYF
    end interface

    !go through model arrays and count number of variable parameters
    n = 0
    do i = 1, size(model_param,1)
       do j = 1, size(model_param,2)
          if (model_prior(i,j) /= 0D0) n=n+1
       end do
    end do

    !allocate variable array and read in the variables
    !x is a vector holding the values of the model parameters that are allowed
    !to vary, x_pos holds the positions of the values in the model_param array
    !that defines the model.
    allocate(x(n))
    allocate(temp_x(n))
    allocate(x_pos(n,2))
    allocate(desc(n))
    k = 0
    do i = 1, size(model_param,1) !loop over cpts
       do j = 1, size(model_param,2)
          if (model_prior(i,j) /= 0D0) then
             k = k+1
             x_pos(k,1) = i
             x_pos(k,2) = j
             x(k) = model_param(i,j)
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
    !the array with their values from x
    if (allocated(fit_param)) deallocate(fit_param)
    allocate(fit_param(size(model_param,1),size(model_param,2)))
    fit_param = model_param

    !compile x_info array
    !1st column is of typical length scales for variables and set equal to the
    !prior width. 2nd and 3rd columns hold lower and upper bounds on the
    !value of the parameter in question - used to prevent negative/other
    !illegal values causing problems (esp in visibility calculations)
    allocate(x_info(n,3))
    do i = 1, n
       x_info(i,1) = model_prior(x_pos(i,1),x_pos(i,2))
    end do
    !preset limits to values from the limits array used in model input
    do i = 1, n
       x_info(i,2) = model_limits(x_pos(i,2),1)
       x_info(i,3) = model_limits(x_pos(i,2),2)
       !alpha(1) in hestroffer model must be > 0
       if (trim(model_spec(x_pos(i,1),3)) == 'hestroffer') then
          if (model_wldep(4) == 1) then
             do iwave = 1, nwave
                if (x_pos(i,2) == (10+(4*model_wldep(1)+model_wldep(2) &
                     +3*model_wldep(3))*(nwave-1) &
                     + max_order*model_wldep(4)*(iwave-1))) x_info(i,2) = 0D0
             end do
          else
             if (x_pos(i,2) == (10+(4*model_wldep(1)+model_wldep(2) &
                  +3*model_wldep(3))*(nwave-1))) x_info(i,2) = 0D0
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
    lwork = n*(n-1)/2+12*n
    if (lwork < 13) lwork = 13
    allocate(work(lwork))
    liwork = n+2
    allocate(iwork(liwork))
    ifail = -1
    call E04JYF(n, 0, posterior, x_info(:,2), x_info(:,3), x, nlposterior, &
         iwork, liwork, work, lwork, iuser, user, ifail)
    deallocate(work)
    deallocate(iwork)
    !solution is in x
    sol(:,1) = x

    !set return message
    select case(ifail)
    case(1)
       stop 'E04JYF returned with ifail=1, indicates programming error'
    case(0)
       info = 'Optimal solution found'
    case(2)
       info = 'Failed - iteration limit exceeded'
    case(3)
       info = 'Lower point cannot be found, but some condition(s) for minimum not met'
    case(4)
       info = 'Failed - overflow occurred'
    case(9)
       info = 'Failed - modulus of a model parameter has become large'
    case(10)
       info = 'Failed with ifail=10 - try again with different initial model'
    case default
       if (ifail >= 5 .and. ifail <= 8) &
            write (info, *) 'Possible minimum found - ifail=', ifail
    end select

    if (ifail /= 0 .and. ifail /= 3) return ! not a minimum, so skip hessian calc.
    
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

    !find hessian elements (nb Aij=Aji symmetric)
    deltai = 5D-3
    deltaj = 5D-3
    do i = 1, n
       do j = 1, i
          temp_x = sol(:,1)
          if (j == i) then
             !P(i)
             call posterior(n, temp_x, P, iuser, user)
             !P(i-di)
             temp_x(i) = sol(i,1)-deltai
             call posterior(n, temp_x, P_l, iuser, user)
             !P(i+di)
             temp_x(i) = sol(i,1)+deltai
             call posterior(n, temp_x, P_u, iuser, user) 
             !hessian by numeric method
             diff = P_u + P_l - (2D0*P)
             if (diff/=0D0) hes(i,i) = diff/(deltai**2)
          else
             !P(i+di,j+dj)
             temp_x(i) = sol(i,1)+deltai
             temp_x(j) = sol(j,1)+deltaj
             call posterior(n, temp_x, P_u, iuser, user)
             !P(i+di,j-dj)
             temp_x(j) = sol(j,1)-deltaj
             call posterior(n, temp_x, P_i, iuser, user)   
             !P(i-di,j-dj)
             temp_x(i) = sol(i,1)-deltai
             call posterior(n, temp_x, P_l, iuser, user)
             !P(i-di,j+dj)
             temp_x(j) = sol(j,1)+deltaj
             call posterior(n, temp_x, P_j, iuser, user)  
             !hessian by numeric method
             diff = P_u + P_l - P_i - P_j
             if (diff/=0D0) hes(i,j) = diff/(4D0*deltai*deltaj)
             hes(j,i) = hes(i,j) !since hessian is symmetric
          end if
       end do
    end do

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
             ifail = -1
          end if
          cor(j,i) = cor(i,j)
       end do
    end do

    !calculate estimated error on solution
    !ideal case: error on a single parameter x is taken as the sqroot(cov(x,x))
    !negative cov(x,x) case: estimated error reported as zero (as occurs for
    !minimisation failure too)
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

    return

    !error trapping 
90  info = 'illegal freedom(s) in model'
    ifail = -2
    goto 200
92  info = 'for vis/nvis data must have guaranteed symmetric model'
    ifail = -3
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
    double precision :: data_phase_err, model_vis, model_amp, model_phase
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

       !ignore if -ve or zero error:
       if (data_vis_err>0D0) then

          !compute model-predicted visibility amplitude squared
          vis = cmplx_vis(spec, param, lambda, delta_lambda, u, v)
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

       if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
          !compute model-predicted triple amplitude
          vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1)
          vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2)
          vis3 = cmplx_vis(spec, param, lambda, delta_lambda, -(u1+u2), -(v1+v2))
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

subroutine posterior(n, x, post, iuser, user)

! Black box minimisers require this function have specific arguments
! Everything else needed comes from module variables

!subroutine arguments
integer, intent(in) :: n
double precision, dimension(n), intent(in) :: x
double precision, intent(out) :: post
integer, dimension(*) :: iuser !required by E04JYF, not used
double precision, dimension(*) :: user !required by E04JYF, not used

!local variables
double precision :: lhd, pri
integer :: i
logical :: violation

!check for range violations in current position
!if position is out of range then return very large energy
!** not needed as NAg minimiser handles constraints
!violation = .false.
!do i = 1, n !over variable parameters
!   if (x(i) < x_info(i,2)) violation = .true.
!   if (x(i) > x_info(i,3)) violation = .true.
!end do
!if (violation) then
!   post = 1.0d10
!   return
!end if 

!update fit_param with current position x
do i = 1, n
   fit_param(x_pos(i,1),x_pos(i,2)) = x(i)
end do

!form E = neg log posterior = neg log likelihood + neg log prior
lhd = likelihood(vis_data, triple_data, model_spec, fit_param)
pri = prior(x, x_pos, model_param, model_prior)
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
double precision :: data_phase_err, model_vis, model_amp, model_phase
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

   if (data_vis_err>0D0) then
      
      !compute model-predicted visibility amplitude squared
      vis = cmplx_vis(model_spec, param, lambda, delta_lambda, u, v)
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

   if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
      !compute model-predicted triple product
      vis1 = cmplx_vis(model_spec, param, lambda, delta_lambda, u1, v1)
      vis2 = cmplx_vis(model_spec, param, lambda, delta_lambda, u2, v2)
      vis3 = cmplx_vis(model_spec, param, lambda, delta_lambda, -(u1+u2), -(v1+v2))
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

function prior(x, x_pos, model_param, model_prior)

!function arguments
double precision :: prior
double precision, dimension(:), intent(in) :: x
integer, dimension(:,:), intent(in) :: x_pos
double precision, dimension(:,:), intent(in) :: model_param
double precision, dimension(:,:), intent(in) :: model_prior

!local variables
integer :: i
double precision :: value, value0, err

prior = 0D0

!loop over the variable parameters
do i = 1, size(x_pos,1)

   value0 = model_param(x_pos(i,1),x_pos(i,2)) !parameter original value
   err = model_prior(x_pos(i,1),x_pos(i,2)) !parameter prior width
   value = x(i) !parameter current value as it is varied

   !calculate normalised contribution to the prior
   prior = prior + (((value-value0)**2D0)/(2D0*(err**2D0))) + log(err*sqrt(2D0*pi))

end do

end function prior

!==============================================================================

subroutine free_fit()

  !Deallocate fitting storage
  if (allocated(fit_param)) deallocate(fit_param)
  if (allocated(x_pos)) deallocate(x_pos)
  if (allocated(x_info)) deallocate(x_info)

end subroutine free_fit

!==============================================================================

end module Fit
