!$Id: fit.f90,v 1.8 2003/05/29 12:40:52 jsy1001 Exp $

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
    !
    !Corresponding message text will be in info,
    !parameter values will be in fit_param
    !
    !force_symm is a logical that forces centrosymmetric only models

    !subroutine arguments
    character(len=128), intent(out) :: info
    character(len=35), dimension(:), allocatable, intent(out) :: desc
    double precision, dimension(:,:), allocatable, intent(out) :: sol, hes, cov, cor
    double precision, intent(out) :: chisqrd, nlposterior, nlevidence
    integer, intent(out) :: flag
    logical, intent(in) :: force_symm

    !local variables
    integer :: i, j, k, n, lwork
    character(len=128) :: name, comp
    character(len=2), dimension(10) :: numbers
    double precision :: P, Pl, Pu, Pi, Pj, deltai, deltaj, eta, diff, det
    double precision, dimension(:), allocatable ::  x, temp_x, work
    logical :: illegal

    interface
       subroutine PDA_UNCMND(n, x0, fcn, x, f, info, w, lw)
         integer, intent(in) :: n, lw
         double precision, dimension(n) :: x0, x
         external fcn
         double precision, intent(out) :: f
         integer, intent(out) :: info
         double precision, dimension(lw), intent(out) :: w
       end subroutine PDA_UNCMND
    end interface

    !go through model arrays and count number of variable parameters
    n = 0
    do i = 1, size(model_param,1)
       do j = 1, 17
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
    k = 0
    do i = 1, size(model_param,1)
       do j = 1, 17
          if (model_prior(i,j) /= 0D0) then
             k = k+1
             x_pos(k,1) = i
             x_pos(k,2) = j
             x(k) = model_param(i,j)
          end if
       end do
    end do

    !assign names to the variables listed in the sol solution array
    allocate(desc(n))
    name = ''
    numbers = (/' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10'/)
    do i = 1, n
       comp = adjustl(numbers(x_pos(i,1)))
       select case(x_pos(i,2))
       case (2)
          name = 'position radius (mas)'
       case (3)
          name = 'position angle (deg)'
       case (4)
          name = 'flux (arb units)'
       case (5)
          name = 'major axis (mas)'
       case (6)
          name = 'orientation (deg)'
       case (7)
          name = 'ellipticity'
       case default
          name = 'LD parameter ' // adjustl(numbers(x_pos(i,2)-7))
       end select
       desc(i) = 'component ' // trim(comp) // ', ' // trim(name)
    end do

    !check for illegal types of minimisation:
    illegal = .false.
    !must have at least 1 freedom to minimise with
    if (n==0) illegal = .true.
    !for single component model cannot vary r/theta (change in position
    !only contstant phase offset) or B (flux is normalised anyway) 
    if (size(model_param,1)==1) then
       do i = 1, size(x_pos,1)
          if (x_pos(i,2)==2) illegal = .true.
          if (x_pos(i,2)==3) illegal = .true.
          if (x_pos(i,2)==4) illegal = .true.
       end do
    end if
    !for any component cannot vary theta if r is zero and not free to vary
    !(position angle has no effect if position radius fixed at zero)
    do i = 1, size(model_param,1)
       if ((model_param(i,2)==0D0) .and. (model_prior(i,2)==0D0) &
            .and. (model_prior(i,3)/=0D0)) illegal = .true.
    end do
    !for any component cannot vary phi if epsilon is unity and not free to vary
    !(orientation is meaningless is ellipse reduced to circular disc)
    do i = 1, size(model_param,1)
       if ((model_param(i,7)==1D0) .and. (model_prior(i,7)==0D0) &
            .and. (model_prior(i,6)/=0D0)) illegal = .true.
    end do
    if (illegal) goto 90

    !if centrosymmetric model is forced then must have 
    !eccentricity epsilon fixed to be unity and position radius fixed at zero
    if (force_symm .and. .not. symm) goto 92

    !allocate fit_param array
    !fit_param is identical to model_param, it holds the definition of the model
    !with the difference that fit_param is designed to vary as the model is
    !adjusted to fit. x_pos is used to update the variable parameters in 
    !the array with their values from x
    if (allocated(fit_param)) deallocate(fit_param)
    allocate(fit_param(size(model_param,1),size(model_param,2)))
    fit_param = model_param

    !compile x_info array
    !1st column is of typical length scales for variables and set equal to the
    !prior width. 2nd and 3rd columns hold lower and upper bounds on the
    !value of the parameter in question - used to prevent negative/other illegal
    !values causing problems (esp in visibility calculations)
    allocate(x_info(n,3))
    do i = 1, n
       x_info(i,1) = model_prior(x_pos(i,1),x_pos(i,2))
    end do
    !preset limits to values from the limits array used in model input
    do i = 1, n
       x_info(i,2) = model_limits(x_pos(i,2),1)
       x_info(i,3) = model_limits(x_pos(i,2),2)
       !alpha(1) in hestroffer model must be > 0
       if (x_pos(i,2)==8) then
          if (trim(model_spec(x_pos(i,1),3))=='hestroffer') x_info(i,2) = 0D0
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
    call gof(chisqrd)
    print *,' '
    print *,'initial chi squared =',real(chisqrd) 

    !call minimising algorithm
    info = ''
    lwork = n*(n+10)
    allocate(work(lwork))
    temp_x = x !starting point
    call PDA_UNCMND(n, temp_x, posterior, x, nlposterior, flag, work, lwork)
    deallocate(work)
    ! solution is in x

    !put best-fit values in fit_param & sol
    do i = 1, n
       fit_param(x_pos(i,1),x_pos(i,2)) = x(i)
       sol(i,1) = x(i)
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

    if (flag < 0 .or. flag > 3) return ! not a minimum, so skip hessian calc.
    
    !find error in x position more rigourously by considering hessian
    !hessian matrix Aij has components (d2/(dyi*dyj))P({y}) where yi and yj
    !are parameters i and j, P({y}) is the neg log posterior for model with
    !parameters {y} = y1, y2, ... yi, ... yj, ... yN
    !estimate each element Aij via taylor 2nd order symmetric numeric 
    !differentiation method, for off-diagonals this is
    ![P(i+di,j+dj)+P(i-di,j-dj)-P(i+di,j-dj)-P(i-di,j+dj)]/[4di*dj], for diagonals
    ![P(i+di)+P(i-di)-2P(i)]/[di^2]
    !eta is scaling factor in selecting increment in numeric differentiation,
    !increment is eta * upper limit on parameter
    eta = 1D-4
    !find hessian elements (nb Aij=Aji symmetric)
    do i = 1, n
       deltai = eta * x_info(i,3)
       do j = 1, i
          deltaj = eta * x_info(j,3)
          temp_x = sol(:,1)
          if (j == i) then
             !P(i)
             call posterior(n, temp_x, P)
             !P(i-di)
             temp_x(i) = sol(i,1)-deltai
             call posterior(n, temp_x, Pl)
             !P(i+di)
             temp_x(i) = sol(i,1)+deltai
             call posterior(n, temp_x, Pu) 
             !hessian by numeric method
             diff = Pu + Pl - (2D0*P)
             if (diff/=0D0) hes(i,i) = diff/(deltai**2)
          else
             !P(i+di,j+dj)
             temp_x(i) = sol(i,1)+deltai
             temp_x(j) = sol(j,1)+deltaj
             call posterior(n, temp_x, Pu)
             !P(i+di,j-dj)
             temp_x(j) = sol(j,1)-deltaj
             call posterior(n, temp_x, Pi)   
             !P(i-di,j-dj)
             temp_x(i) = sol(i,1)-deltai
             call posterior(n, temp_x, Pl)
             !P(i-di,j+dj)
             temp_x(j) = sol(j,1)+deltaj
             call posterior(n, temp_x, Pj)  
             !hessian by numeric method
             diff = Pu + Pl - Pi - Pj
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
             flag = 1
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
    call gof(chisqrd)

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
    goto 200
91  info = 'minimisation routine error: ' // info
    goto 200
92  info = 'for vis/nvis data must have guaranteed symmetric model'
    goto 200

  end subroutine minimiser

!==============================================================================

subroutine gof(chisqrd)

!Goodness of fit information calculation
!
!Computes the chisqrd value of a model fit and data set
!chisqrd = sum[i] { (data[i]-theory[i])^2/data_error[i]^2 }
!chisqrd is thus almost identical to likelihood
!could save code here and significantly combine this function
!with the likeihood one but maybe less risky to keep them separated...
!Could in future add code here to look up critical chi sqaured values
!for num of data points and num of free parameters

!subroutine arguments
double precision, intent(out) :: chisqrd

!local variables
integer :: i
double precision :: lambda, u, v, data_vis, data_vis_err, model_vis
double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
double precision :: data_phase_err, model_amp, model_phase
double complex :: vis, vis1, vis2, vis3

chisqrd = 0D0

!sum over the visibility data points
do i = 1, size(vis_data,1)
   
   !extract points
   lambda = vis_data(i,1)
   u = vis_data(i,3)
   v = vis_data(i,4)
   data_vis = vis_data(i,5)
   data_vis_err = vis_data(i,6)

   !ignore if -ve or zero error:
   if (data_vis_err>0D0) then
      
      !compute model-predicted visibility amplitude squared
      vis = cmplx_vis(model_spec, fit_param, lambda, u, v)
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
      vis1 = cmplx_vis(model_spec, fit_param, lambda, u1, v1)
      vis2 = cmplx_vis(model_spec, fit_param, lambda, u2, v2)
      vis3 = cmplx_vis(model_spec, fit_param, lambda, -(u1+u2), -(v1+v2))
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

subroutine posterior(n, x, post)

! Black box minimisers require this function have specific arguments
! Everything else needed comes from module variables

!subroutine arguments
double precision, dimension(n), intent(in) :: x
integer, intent(in) :: n
double precision, intent(out) :: post

!local variables
double precision :: lhd, pri
integer :: i
logical :: violation

!check for range violations in current position
!if position is out of range then return very large energy 
violation = .false.
do i = 1, n !over variable parameters
   if (x(i) < x_info(i,2)) violation = .true.
   if (x(i) > x_info(i,3)) violation = .true.
end do
if (violation) then
   post = 1.0d10
   return
end if 

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
double precision :: lambda, u, v, data_vis, data_vis_err, model_vis
double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
double precision :: data_phase_err, model_amp, model_phase
double complex :: vis, vis1, vis2, vis3

likelihood = 0D0

!sum over the visibility data points
do i = 1, size(vis_data,1)
   
   !extract points
   lambda = vis_data(i,1)
   u = vis_data(i,3)
   v = vis_data(i,4)
   data_vis = vis_data(i,5)
   data_vis_err = vis_data(i,6)

   if (data_vis_err>0D0) then
      
      !compute model-predicted visibility amplitude squared
      vis = cmplx_vis(model_spec, param, lambda, u, v)
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
      vis1 = cmplx_vis(model_spec, param, lambda, u1, v1)
      vis2 = cmplx_vis(model_spec, param, lambda, u2, v2)
      vis3 = cmplx_vis(model_spec, param, lambda, -(u1+u2), -(v1+v2))
      vis = vis1*vis2*vis3
      model_amp = modulus(vis)
      model_phase = rad2deg*argument(vis)
   end if

   !compute normalised likelihood contributions
   !nb phase calculations modulo 360 degrees
   if (data_amp_err>0D0) then
      likelihood = likelihood + &
           (((data_amp-model_amp)**2D0)/(2D0*(data_amp_err**2D0))) + &
           log(data_amp_err*sqrt(2D0*pi))
   end if
   if (data_phase_err>0D0) then
      likelihood = likelihood + &
           (((modx(data_phase,model_phase))**2D0)/(2D0*(data_phase_err**2D0))) + &
           log(data_phase_err*sqrt(2D0*pi))
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
