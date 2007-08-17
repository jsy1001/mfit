!$Id: marginalise.F90,v 1.2 2007/08/17 17:50:27 jsy1001 Exp $

module Marginalise

  !public subroutines contained:
  !
  ! marg_post - return marginalised -ln(posterior)
  ! marg_err  - determine parameter error bars from marginalised -ln(posterior)
  !
  !local subroutines contained:
  !
  ! prob, prob1d - return posterior of model and data, called by integrators
  !                used in marg_post
  ! ferr - function marg_err finds roots of

  use Maths
  use Bayes
  use Wrap
  use Model
  use Fit, only: minimiser

  implicit none

  public
  private :: prob, prob1d, ferr

  !private module variables contained:
  private :: mg_par, mg_nlnorm, &
       mgerr_par, mgerr_ivar, mgerr_marg_var, mgerr_minpost

  !! Model parameters for prob()
  type(allparam), save :: mg_par

  !! Used to normalise values calculated by prob, to avoid overflow
  double precision :: mg_nlnorm


  !used by marg_err/ferr:

  type(allparam), save :: mgerr_par
  !! Variable to estimate error bar for
  integer :: mgerr_ivar
  !! Variables to marginalise over (all but mgerr_ivar)
  logical, allocatable :: mgerr_marg_var(:)
  double precision :: mgerr_minpost

contains

  !============================================================================

  !! Evaluate negative log of (posterior probability marginalised over
  !! specified variables)
  function marg_post(info, inpar, marg_var)

    double precision :: marg_post

    !function arguments
    character(len=*), intent(out) :: info !! Error/warning message
    type(allparam) :: inpar !! Point at which to evaluate -ln(marg. post.)
    !! If marg_var(i) true, marginalise over variable i in inpar
    logical, intent(in) :: marg_var(:)

    !local variables
    integer :: iv, iwave, idim, ndim, minpts, maxpts, lenwrk, leniwrk
    integer :: ifail, alpha, indx1, indx2
    double precision :: epsabs, epsrel, acc, abserr, finval, sol0, sig
    double precision :: lhd, pri, delt, lbound, ubound, chisqrd, unmg_post
    double precision, dimension(:), allocatable :: sol(:), alim(:), blim(:)
    double precision, dimension(:), allocatable :: wrk(:)
    integer, allocatable :: iwrk(:)
    type(allparam) :: fitpar
    logical :: fit_ok

    integer :: num_points, i
    double precision :: val

    interface
       subroutine D01FCF(ndim, alim, blim, minpts, maxpts, functn, epsrel, &
            acc, lenwrk, wrk, finval, ifail)
         integer, intent(in) :: ndim, lenwrk
         double precision, dimension(ndim), intent(in) :: alim, blim
         integer, intent(inout) :: minpts
         integer, intent(in) :: maxpts
         double precision, external :: functn
         double precision, intent(in) :: epsrel
         double precision, intent(out) :: acc, finval
         double precision, dimension(lenwrk), intent(out) :: wrk
         integer, intent(inout) :: ifail
       end subroutine D01FCF

       subroutine D01AJF(f, aa, bb, epsabs, epsrel, finval, abserr, &
            wrk, lenwrk, iwrk, leniwrk, ifail)
         integer, intent(in) :: lenwrk, leniwrk
         double precision, external :: f
         double precision, intent(in) :: aa, bb
         double precision, intent(in) :: epsabs, epsrel
         double precision, intent(out) :: abserr, finval
         double precision, dimension(lenwrk), intent(out) :: wrk
         integer, dimension(leniwrk), intent(out) :: iwrk
         integer, intent(inout) :: ifail
       end subroutine D01AJF
    end interface

    !evaluate unmarginalised -ln(postprob) (used for scaling)
    !at minimum found by varying to-marginalise parameters only
    call allparam_copy(inpar, fitpar, (.not. marg_var))
    ndim = fitpar%nvar !no. of dimensions to marginalise over
    if (ndim > 0) then
       allocate(sol(ndim))
       !minimiser sets parameter scaling using model_prior
       call minimiser(fitpar, sol, chisqrd, unmg_post, fit_ok, info)
       if (.not. fit_ok) then
          info = 'Fit failed in marg_post'
          return
       end if
    else
       lhd = likelihood(vis_data, triple_data, model_spec, fitpar%param)
       pri = prior(fitpar%var_pos, fitpar%param, model_param, model_prior)
       unmg_post = lhd + pri
    end if
    
    !special case - no marginalisation to do
    if (ndim == 0) then
       marg_post = unmg_post
       info = ''
       return
    end if

    !choose integration limits
    allocate(alim(ndim))
    allocate(blim(ndim))
    do idim = 1, ndim
       indx1 = fitpar%var_pos(idim,1)
       indx2 = fitpar%var_pos(idim,2)
       sol0 = sol(idim) !centre integration range on peak
       sig = model_prior(indx1, indx2)
       alim(idim) = sol0 - sig
       blim(idim) = sol0 + sig
       !ensure limits in bounds
       lbound = model_limits(indx1, indx2, 1)
       ubound = model_limits(indx1, indx2, 2)
       if (alim(idim) < lbound) alim(idim) = lbound
       if (blim(idim) > ubound) blim(idim) = ubound
       !beware periodicities in theta
       !!CAUSES PROBS for relto theta
       !do iwave = 1, nwave
       !   if (indx2 == 3+4*model_wldep(1)*(iwave-1)) then
       !      alim(idim) = modulo(alim(idim), 360D0)
       !      blim(idim) = modulo(blim(idim), 360D0)
       !      if (alim(idim) >= blim(idim)) alim(idim) = alim(idim) - 360D0
       !   end if
       !end do
    end do

    !set up mg_par module variable, used by prob()
    call allparam_copy(inpar, mg_par, (.not. marg_var))
    call allparam_setnoscale(mg_par)

    !call appropriate NAg routine to perform integral
#ifdef HAVE_NAG
    if (ndim >= 2) then
       minpts = 1
       alpha = 2**ndim + 2*ndim**2 + 2*ndim + 1
       maxpts = 50000*alpha
       lenwrk = (ndim+2)*(1+maxpts/alpha)
       if (lenwrk < (2*ndim + 4)) lenwrk = 2*ndim + 4
       allocate(wrk(lenwrk))
       !need to overspecify epsrel to get sensible results
       !epsrel = 5D-2
       epsrel = 1D-3
       mg_nlnorm = unmg_post + 0.5D0*log(machine_max())
       ifail = -1
       call D01FCF(ndim, alim, blim, minpts, maxpts, prob, epsrel, acc, &
            lenwrk, wrk, finval, ifail)
       !set return message
       select case(ifail)
       case(1)
          stop 'D01FCF returned with ifail=1, indicates programming error'
       case(2)
          info = 'MAXPTS too small to obtain required accuracy for integral'
       case(3)
          info = 'LENWRK too small to obtain required accuracy for integral'
       case (0)
          info = ''
       end select
       !can get zero result if NAg routine only explores zero-valued region
       if (finval == 0D0) then
          info = 'Integration result zero - try narrowing prior'
          !return unmarginalised -ln(postprob)
          marg_post = unmg_post
       else
          marg_post = -log(finval)+mg_nlnorm !useful even if ifail = 2 or 3
          print '(a,f15.3,a,f8.4)', 'marg_post =', -log(finval)+mg_nlnorm, &
               ' +/-', acc
          print *, '  in', minpts, 'function evalutions'
       end if
    else
       !special case of 1d integral
       lenwrk = 2000
       leniwrk = lenwrk/4
       allocate(wrk(lenwrk))
       allocate(iwrk(leniwrk))
       epsabs = -1D0
       epsrel = 5D-2
       mg_nlnorm = unmg_post + 0.5D0*log(machine_max())
       ifail = -1
       call D01AJF(prob1d, alim(1), blim(1), epsabs, epsrel, finval, abserr, &
            wrk, lenwrk, iwrk, leniwrk, ifail)
       !set return message
       select case(ifail)
       case(6)
          stop 'D01AJF returned with ifail=6, indicates programming error'
       case(1)
          info = 'Max. subdivisions reached'
       case(2)
          info = 'Roundoff error prevents requested tolerance being reached'
       case (3)
          info = 'Bad local integrand behaviour'
       case(4)
          info = 'requested tolerance cannot be achieved'
       case(5)
          info = 'Integral slowly convergent or divergent'
       case(0)
          info = ''
       end select
       !can get zero result if D01AJF only explores zero-valued region
       if (finval == 0D0) then
          info = 'Integration result zero - try narrowing prior'
          !return unmarginalised -ln(postprob)
          marg_post = unmg_post
       else
          marg_post = -log(finval)+mg_nlnorm
          print '(a,f15.4,a,f8.4)', 'marg_post =', marg_post, &
               ' +/-', abserr/finval
       end if
       if (.false.) then
          !save integrand values to text file
          write (33, *) '# ', marg_post
          write (33, *) '# '//trim(info)
          write (33, *) '# ', mg_par%param
          num_points = 20
          do i = 1, num_points
             val = alim(1) + ((i-1.)/(num_points-1.))*(blim(1)-alim(1))
             write (33,*) val, prob1d(val)
          end do
       end if
    end if
#endif

    !free storage
    call allparam_free(fitpar)
    call allparam_free(mg_par)
    if (allocated(sol)) deallocate(sol)
    if (allocated(alim)) deallocate(alim)
    if (allocated(blim)) deallocate(blim)
    if (allocated(wrk)) deallocate(wrk)
    if (allocated(iwrk)) deallocate(iwrk)

  end function marg_post

  !============================================================================

  !! Return posterior of model and data, called by the integrators
  !! used in marg_post
  !!
  !! Uses mg_par, mg_nlnorm from this module
  !! Uses vis_data, triple_data from Fit module
  !! Uses model_spec, model_param, model_prior from Model module
  function prob(ndim, z)

    double precision :: prob

    !function arguments
    integer, intent(in) :: ndim !! No. of dimensions being integrated over
    double precision, dimension(ndim), intent(in) :: z !!Param values for these

    !local variables
    logical, save :: set_max_min = .false.
    double precision, save :: log_max
    double precision, save :: log_min
    double precision :: lhd, pri, lprob

    !set log_max and log_min if unset
    if (.not. set_max_min) then
       set_max_min = .true.
       log_max = floor(log(machine_max()))
       log_min = -floor(-log(machine_min()))
    end if

    !copy variable parameters from z
    call allparam_setvar(mg_par, z)

    !neg log posterior = neg log likelihood + neg log prior - neg log evidence
    lhd = likelihood(vis_data, triple_data, model_spec, mg_par%param)
    pri = prior(mg_par%var_pos, mg_par%param, model_param, model_prior)
    lprob = -(lhd + pri - mg_nlnorm)
    !As above, if model very improbable, a small *fractional* reduction in
    !chi-squared from varying marginalised parameters can cause overflow
    if (lprob > log_max) then
       prob = exp(log_max)
       print *, 'Trapped overflow in prob()'
    else if (lprob < log_min) then
       prob = exp(log_min)
       !underflow is nothing to worry about
       !print *, 'Trapped underflow in prob()'
    else
       prob = exp(lprob)
    end if

  end function prob

  !============================================================================

  !! Wrapper for prob(), for use by 1d integrator
  function prob1d(val)

    double precision :: prob1d

    !function arguments
    double precision, intent(in) :: val

    !local variables
    double precision :: val1(1)

    val1(1) = val
    prob1d = prob(1, val1)

  end function prob1d

  !============================================================================

  subroutine marg_err(inpar, errguess, index, err)

    !subroutine arguments
    type(allparam) :: inpar !! Model parameters (variable and non-variable)
    !! Initial guess errors for variable parameters
    double precision, intent(in) :: errguess(:)
    !! Variable parameter to estimate error for
    integer, intent(in) :: index
    !! Final +/- 1-sigma errors:
    double precision, intent(out) :: err(2)

    !local variables
    double precision x0, x1, sigma, rtb
    double precision dx, xmid, fmid
    double precision :: err_quad(2) !! Estimate of errors by quadratic interp.
    integer nb, i, j, ifail
    integer, parameter :: max_iter = 20
    real, parameter :: tol = 0.1
    integer, dimension(2), parameter :: sign = (/1, -1/)
    logical, parameter :: usenag = .true.

    interface
       subroutine C05ADF(a, b, eps, eta, f, x, ifail)
         double precision, intent(in) :: a, b, eps, eta
         double precision, external :: f
         double precision, intent(out) :: x
         integer, intent(inout) :: ifail
       end subroutine C05ADF
    end interface

    !assign to module variables needed by ferr
    call allparam_copy(inpar, mgerr_par)
    mgerr_ivar = index
    allocate(mgerr_marg_var(inpar%nvar))
    mgerr_marg_var = .true.
    !marginalise over all variables except one we want error bar for:
    mgerr_marg_var(mgerr_ivar) = .false.

    !find 2 roots of f = nlmpost - min(nlmpost) - 0.5 = 0
    !(where nlmpost is negative log of marginalised posterior probability)

    !find value at minimum
    x0 = inpar%param(inpar%var_pos(index,1), inpar%var_pos(index,2))
    sigma = errguess(index)
    mgerr_minpost = 0D0 !so we can use ferr() for next line
    mgerr_minpost = ferr(x0) + 0.5D0
    do j = 1, 2
       !find bracketing point other side of root
       do nb = 1, 10
          x1 = x0 + sign(j)*(2**nb)*sigma
          fmid = ferr(x1)
          if (fmid > 0D0) exit
       end do
       if (fmid <= 0D0) stop 'marg_err: cannot find bracketing point'
       !find root by quadratic interpolation
       err_quad(j) = sign(j)*(x1-x0)/sqrt(2D0*fmid+1D0)
       !find root by bisection
#ifdef HAVE_NAG
       ifail = 0
       call C05ADF(x0, x1, tol*sigma, 0D0, ferr, rtb, ifail)
       err(j) = sign(j)*(rtb-x0)
#else
       dx = x1 - x0
       rtb = x0
       bisect: do i = 1, max_iter
          dx = dx/2D0
          xmid = rtb + dx
          fmid = ferr(xmid)
          if (fmid <= 0D0) rtb = xmid
          if (abs(dx) < tol*sigma .or. fmid == 0D0) then
             err(j) = sign(j)*(rtb-x0)
             exit bisect
          end if
       end do bisect
#endif
       !compare estimates
       print *, err(j), err_quad(j), (err(j)-err_quad(j))/err(j)
    end do

    !clean up
    if (allocated(mgerr_marg_var)) deallocate(mgerr_marg_var)

  end subroutine marg_err

  !============================================================================

  function ferr(val)

    double precision :: ferr
    double precision, intent(in) :: val

    !local variables
    character(len=128) :: info

    !assign value for unmarginalised variable:
    call allparam_setone(mgerr_par, mgerr_ivar, val)

    ferr = marg_post(info, mgerr_par, mgerr_marg_var) - mgerr_minpost - 0.5D0
    if (info /= '') print *, trim(info)
    
  end function ferr

  !============================================================================

end module Marginalise
