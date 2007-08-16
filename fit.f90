!$Id: fit.f90,v 1.23 2007/08/16 16:40:32 jsy1001 Exp $

module Fit

  use Maths
  use Bayes
  use Wrap
  use Model

  implicit none

  private

  !public subroutines contained
  !
  !minimiser - wrapper to the minimising routine (can call multiple times to
  !            implement a minimisation strategy)
  !err_est - check for minimum, estimate parameter error bars at solution point
  !free_fit - frees storage allocated by minimiser

  public :: minimiser, err_est, free_fit


  !private module variables contained:
  
  type(allparam), save :: fitpar !! Used by minimiser, err_est, posterior

  double precision, allocatable :: work(:) !! Workspace
  integer :: lwork = 0 !! Size of work array

  double precision :: post2_nlnorm !! Scaling in post2

contains

  !============================================================================
  
  !! (Re-)allocate workspace
  subroutine alloc_work(lwork_req)

    !subroutine arguments
    integer :: lwork_req !! Minimum size required for work array

    if (lwork_req > lwork) then
       if (allocated(work)) deallocate(work)
       lwork = lwork_req
       allocate(work(lwork))
    end if

  end subroutine alloc_work

  !============================================================================
  
  !! Run minimiser and calculate goodness-of-fit at stopping point
  subroutine minimiser(inpar, sol, chisqrd, nlposterior, success, info)

    !subroutine arguments
    type(allparam), intent(inout) :: inpar !! Model parameters
    !! Variable parameter values at minimum
    double precision, intent(out) :: sol(:)
    double precision, intent(out) :: chisqrd    !! chi-squared at sol
    double precision, intent(out) :: nlposterior!! -ln(posterior) at sol
    !! True if minimiser completed (but may not have stopped at a minimum)
    logical, intent(out) :: success
    character(len=*), intent(out) :: info   !! Error/warning message

    !local variables
    integer :: n, i, flag
    double precision :: minval
    double precision :: x_scale(size(sol,1)), x0(size(sol,1)), x(size(sol,1))

    interface
       subroutine PDA_UNCMND(n, x0, fcn, x, f, flag, w, lw)
         integer, intent(in) :: n, lw
         double precision, intent(in) :: x0(n)
         external fcn
         double precision, intent(out) :: x(n)
         double precision, intent(out) :: f
         integer, intent(out) :: flag
         double precision, intent(out) :: w(lw)
       end subroutine PDA_UNCMND
    end interface

    !copy model parameters to private module variable
    call allparam_copy(inpar, fitpar)

    !set variable scaling, and
    !assign initial values of fit variables (=SCALED variable model parameters)
    n = inpar%nvar !number of fit variables
    do i = 1, n
       x_scale(i) = 1D-3*model_prior(inpar%var_pos(i,1), inpar%var_pos(i,2))
       x0(i) = &
            model_param(inpar%var_pos(i,1), inpar%var_pos(i,2)) / x_scale(i)
    end do
    call allparam_setscale(fitpar, x_scale)

    !set default values for intent(out) arguments
    info = ''
    success = .false.
    chisqrd = 0D0
    nlposterior = 0D0

    !call minimising algorithm
    call alloc_work(n*(n+10))
    call PDA_UNCMND(n, x0, posterior, x, minval, flag, work, lwork)
    nlposterior = minval
    sol = x*x_scale !x contains scaled solution

    !restart from where it finished, this time minimising
    !1/(posterior probability) rather than -log(prob)
    if (flag /= -1 .and. flag <= 3) then
       x0 = x
       post2_nlnorm = nlposterior
       call PDA_UNCMND(n, x0, post2, x, minval, flag, work, lwork)
       print *, x-(sol/x_scale) !is shift significant?
       print *, nlposterior
       nlposterior = log(minval) + post2_nlnorm
       sol = x*x_scale !x contains scaled solution
    end if

    !set error/warning message
    select case(flag)
    case(-1)
       stop 'Insufficient workspace for PDA_UNCMND'
    case(0)
       info = 'Optimal solution found'
       success = .true.
    case(1)
       info = 'Terminated with gradient small, parameter values probably optimal'
       success = .true.
    case(2)
       info = 'Terminated with step size small, parameter values probably optimal'
       success = .true. !flag==2 usually means it worked
    case(3)
       info = 'Lower point cannot be found, parameter values probably optimal'
       success = .true. !flag==3 usually means it worked
    case(4)
       info = 'Iteration limit exceeded'
    case(5)
       info = 'Too many large steps, likelihood function may be unbounded'
    case default
       stop 'Unexpected status value from PDA_UNCMND'
    end select

    !calculate final goodness-of-fit statistic chi squared
    call allparam_setvar(fitpar, x) !fitpar set up to use scaled variables
    call gof(model_spec, fitpar%param, chisqrd)

  end subroutine minimiser

  !===========================================================================
  
  !! Check for minimum, estimate parameter error bars at solution point
  !! Need to call minimiser() before invoking this routine
  subroutine err_est(n, sol, &
       found_min, hes_valid, hes, cov, cor, err, nlevidence, info)

    !subroutine arguments
    integer, intent(in) :: n !! Number of variable parameters
    !! Variable parameter values at minimum
    double precision, intent(in) :: sol(n)
    logical, intent(out) :: found_min !! If true, found a local minimum
    !! If true, Hessian estimate in hes (and hence subsequent args) valid
    logical, intent(out) :: hes_valid
    double precision, intent(out) :: hes(n,n) !! Estimate of Hessian matrix
    double precision, intent(out) :: cov(n,n) !! Estimate of covariance matrix
    double precision, intent(out) :: cor(n,n) !! Estimate of correlation matrix
    !! Parameter errors from covariance  matrix
    double precision, intent(out) :: err(n)
    double precision, intent(out) :: nlevidence !! -ln(evidence) at sol
    character(len=*), intent(out) :: info   !! Error/warning message

    !parameters
    double precision, parameter :: maxgrad = 1D-1 !! Max. acceptable gradient

    !local variables
    integer :: i, j
    double precision :: nlposterior
    double precision :: P, P_l, P_u, P_i, P_j, deltai, deltaj, diff, det
    double precision :: x_scale(n), x(n), temp_x(n), grad(n)

    !set default values for intent(out) arguments
    info = ''
    found_min = .false.
    hes_valid = .false.
    err = 0D0
    hes = 0D0
    cov = 0D0
    cor = 0D0

    if (.not. fitpar%done_init) then
       info = 'minimiser() not called prior to calling err_est()'
       return
    end if

    !posterior() works with SCALED parameters
    x_scale = fitpar%var_scale
    do i = 1, n
       x(i) = sol(i) / x_scale(i)
    end do
    call posterior(n, x, nlposterior)

    !find error in sol position by considering hessian
    !hessian matrix Aij has components (d2/(dyi*dyj))P({y}) where yi and yj
    !are parameters i and j, P({y}) is the neg log posterior for model with
    !parameters {y} = y1, y2, ... yi, ... yj, ... yN
    !estimate each element Aij via taylor 2nd order symmetric numeric 
    !differentiation method, for off-diagonals this is
    ![P(i+di,j+dj)+P(i-di,j-dj)-P(i+di,j-dj)-P(i-di,j+dj)]/[4di*dj],
    !for diagonals
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
             !hessian by numeric method
             diff = P_u + P_l - (2D0*P)
             hes(i,i) = diff / ((deltai**2)*(x_scale(i)**2))
             !check we are at a minimum w.r.t. this parameter
             grad(i) = (P_u - P_l)/(2D0*deltai)
             print *, 'gradient', i, grad(i)
             if (abs(grad(i)) > maxgrad .or. hes(i,i) < 0D0) then
                info = 'Final position not a local minimum'
                return
             end if
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
             hes(i,j) = diff / ((4D0*deltai*deltaj)*(x_scale(i)*x_scale(j)))
             hes(j,i) = hes(i,j) !since hessian is symmetric
          end if
       end do
    end do
    found_min = .true.

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
          end if
          cor(j,i) = cor(i,j)
       end do
    end do

    !calculate estimated error on solution
    !ideal case: error on a single parameter x is taken as the sqroot(cov(x,x))
    !negative cov(x,x) case: estimated error reported as zero (as occurs for
    !minimisation failure too)
    hes_valid = .true.
    do i = 1, n
       if (cov(i,i) > 0) then
          err(i) = sqrt(cov(i,i))
       else
          print *, 'Diagonal element', i, 'of covariance matrix negative'
          hes_valid = .false.
          info = 'Illegal value(s) in estimated covariance matrix'
       end if
    end do

    !estimate -log(evidence)
    !if det -ve return zero
    if (det > 0) then
       nlevidence = nlposterior - 0.5D0*n*log(2D0*pi) + 0.5D0*log(det)
    else
       nlevidence = 0D0
    end if

  end subroutine err_est

  !============================================================================

  !! Calculate negative log posterior probability.
  !! Black box minimisers require this function have specific arguments;
  !! everything else needed comes from module variables:
  !! - fitpar from this module
  !! - vis_data, triple_data from Bayes
  !! - model_spec, model_param, model_prior from Model
  subroutine posterior(n, x, funcval)

    !subroutine arguments
    integer, intent(in) :: n !! Number of variable parameters
    double precision, intent(in) :: x(n) !! Values of variable params
    double precision, intent(out) :: funcval !! Value of function to minimise
    
    !local variables
    double precision :: lhd, pri
    
    !check for range violations in current position
    !if position is out of range then return very large energy 
    if (.not. allparam_inlimit(fitpar, x)) then
       funcval = 1.0d10
       return
    end if

    !form E = neg log posterior = neg log likelihood + neg log prior
    call allparam_setvar(fitpar, x)
    lhd = likelihood(vis_data, triple_data, model_spec, fitpar%param)
    pri = prior(fitpar%var_pos, fitpar%param, model_param, model_prior)
    funcval = lhd + pri

  end subroutine posterior

  !============================================================================

  !! Calculate 1/(posterior probability).
  !! Black box minimisers require this function have specific arguments;
  !! everything else needed comes from module variables:
  !! - fitpar from this module
  !! - vis_data, triple_data from Bayes
  !! - model_spec, model_param, model_prior from Model
  subroutine post2(n, x, funcval)

    !subroutine arguments
    integer, intent(in) :: n !! Number of variable parameters
    double precision, intent(in) :: x(n) !! Values of variable params
    double precision, intent(out) :: funcval !! Value of function to minimise

    !local variables
    double precision :: nlprob
    logical, save :: set_max_min = .false.
    double precision, save :: log_max
    double precision, save :: log_min

    !set log_max and log_min if unset
    if (.not. set_max_min) then
       set_max_min = .true.
       log_max = floor(log(machine_max()))
       log_min = -floor(-log(machine_min()))
    end if

    call posterior(n, x, nlprob)
    nlprob = nlprob - post2_nlnorm

    if (nlprob > log_max) then
       funcval = exp(log_max)
       !nothing to worry about
       print *, 'Trapped overflow in post2()'
    else if (nlprob < log_min) then
       funcval = exp(log_min)
       print *, 'Trapped underflow in post2()'
    else
       funcval = exp(nlprob)
    end if

  end subroutine post2

  !============================================================================

  !! Deallocate module storage
  subroutine free_fit()
    
    call allparam_free(fitpar)
    if (allocated(work)) deallocate(work)

  end subroutine free_fit

  !============================================================================

end module Fit
