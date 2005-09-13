!$Id: marginalise.f90,v 1.2 2005/09/13 09:52:51 jsy1001 Exp $

module Marginalise

  !callable subroutines contained:
  !
  ! alloc_mg, free_mg -- allocate/free module storage
  ! marg_post - return marginalised -ln(posterior)
  ! marg_err  - determine parameter error bars from marginalised -ln(posterior)
  !
  !local subroutines contained:
  !
  ! prob, prob1d - return posterior of model and data, called by integrators
  !                used in marg_post
  ! ferr - function marg_err finds roots of

  use Fit
  use Model

  implicit none

  !module variables contained:

  !values of both variable and non-variable model parameters
  double precision, dimension(:, :), allocatable :: mg_param

  !values of all variable model parameters
  !mg_var(i) corresponds to mg_param(Fit::x_pos(i))
  double precision, dimension(:), allocatable :: mg_var

  !if mg_marg(i) true, marginalise over parameter Fit::x_pos(i)
  logical, dimension(:), allocatable :: mg_marg

  double precision :: mg_nlnorm, post

  !used by marg_err/ferr
  integer :: ivar
  double precision :: minpost

contains

!==============================================================================

  subroutine alloc_mg

    allocate(mg_param(size(fit_param, 1), 17))
    allocate(mg_var(size(x_pos, 1)))
    allocate(mg_marg(size(x_pos, 1)))

  end subroutine alloc_mg

!==============================================================================

  subroutine free_mg

    if (allocated(mg_param)) deallocate(mg_param)
    if (allocated(mg_var)) deallocate(mg_var)
    if (allocated(mg_marg)) deallocate(mg_marg)

  end subroutine free_mg

!==============================================================================

  function marg_post(info)

    !function arguments
    double precision :: marg_post
    character(len=128), intent(out) :: info

    !local variables
    integer :: iv, iwave, idim, ndim, minpts, maxpts, lenwrk, leniwrk
    integer :: ifail, alpha
    double precision :: epsabs, epsrel, acc, abserr, finval, pr0, sig
    double precision :: lhd, pri, delt, peak
    double precision, dimension(:), allocatable :: alim, blim, wrk
    integer, dimension(:), allocatable :: iwrk

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

    mg_param = fit_param !ensure non-variables have correct values
    
    !count dimensions to marginalise over
    ndim = 0
    do iv = 1, size(x_pos, 1)
       if (mg_marg(iv)) then
          ndim = ndim + 1
          mg_var(iv) = fit_param(x_pos(iv, 1), x_pos(iv, 2)) !needed for prior
       else
          mg_param(x_pos(iv, 1), x_pos(iv, 2)) = mg_var(iv)
       end if
    end do

    !evaluate unmarginalised -ln(postprob): used for scaling
    !neg log post = neg log likelihood + neg log prior - neg log evidence
    lhd = likelihood(vis_data, triple_data, model_spec, mg_param)
    pri = prior(mg_var, x_pos, model_param, model_prior)
    post = lhd + pri - mg_nlnorm
    !BUT, if model very improbable, a small *fractional* reduction in
    !chi-squared from varying marginalised parameters can still cause overflow

    !special case - no marginalisation to do
    if (ndim == 0) then
       marg_post = post
       info = ''
       return
    end if

    !choose integration limits
    allocate(alim(ndim))
    allocate(blim(ndim))
    idim = 1
    do iv = 1, size(mg_var, 1)
       if (mg_marg(iv)) then
          pr0 = model_param(x_pos(iv,1), x_pos(iv,2))
          sig = model_prior(x_pos(iv,1), x_pos(iv,2))
          alim(idim) = pr0 - sig !XXX
          blim(idim) = pr0 + sig
          !ensure limits in bounds
          if (alim(idim) < x_info(iv, 2)) alim(idim) = x_info(iv, 2)
          if (blim(idim) > x_info(iv, 3)) blim(idim) = x_info(iv, 3)
          !beware periodicities in theta
          !!CAUSES PROBS for relto theta
          !do iwave = 1, nwave
          !   if (x_pos(iv, 2) == 3+4*model_wldep(1)*(iwave-1)) then
          !      alim(idim) = modulo(alim(idim), 360D0)
          !      blim(idim) = modulo(blim(idim), 360D0)
          !      if (alim(idim) >= blim(idim)) alim(idim) = alim(idim) - 360D0
          !   end if
          !end do
          idim = idim + 1
       end if
    end do

    !call appropriate NAg routine to perform integral
    if (ndim >= 2) then
       minpts = 1
       alpha = 2**ndim + 2*ndim**2 + 2*ndim + 1
       maxpts = 50000*alpha
       lenwrk = (ndim+2)*(1+maxpts/alpha)
       if (lenwrk < (2*ndim + 4)) lenwrk = 2*ndim + 4
       allocate(wrk(lenwrk))
       epsrel = 5D-2
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
          marg_post = post
       else
          marg_post = -log(finval)+post !useful even if ifail = 2 or 3
          print *, marg_post, post
       end if
    else
       !special case of 1d integral
       epsabs = -1D0
       epsrel = 5D-2
       lenwrk = 2000
       leniwrk = lenwrk/4
       allocate(wrk(lenwrk))
       allocate(iwrk(leniwrk))
       ifail = -1
       !finval = 0D0
       !delt = (blim(1) - alim(1))/100D0
       !do iv = 1, 100
       !   finval = finval + prob1d(alim(1)+iv*delt)*delt
       !end do
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
          marg_post = post
       else
          marg_post = -log(finval)+post
          print *, marg_post, post
       end if
    end if

    !free storage
    if (allocated(alim)) deallocate(alim)
    if (allocated(blim)) deallocate(blim)
    if (allocated(wrk)) deallocate(wrk)
    if (allocated(iwrk)) deallocate(iwrk)

  end function marg_post

!==============================================================================

  function prob(ndim, z)

    !uses vis_data, triple_data, x_pos from Fit module
    !uses model_spec, model_param, model_prior from Model module

    !function arguments
    double precision :: prob
    integer, intent(in) :: ndim
    double precision, dimension(ndim), intent(in) :: z

    !local variables
    integer :: iv, iz
    double precision :: lhd, pri, lprob

    !copy subset of variable parameters in z to mg_param and mg_var,
    !otherwise use fit_param values (copied in marg_post)
    iz = 1
    do iv = 1, size(mg_var, 1)
       if (mg_marg(iv)) then
          mg_var(iv) = z(iz)
          mg_param(x_pos(iv, 1), x_pos(iv, 2)) = z(iz)
          iz = iz + 1
       end if
    end do

    !neg log posterior = neg log likelihood + neg log prior - neg log evidence
    lhd = likelihood(vis_data, triple_data, model_spec, mg_param)
    pri = prior(mg_var, x_pos, model_param, model_prior) !initial guesses in model_param
    lprob = -(lhd + pri - mg_nlnorm - post)
    !As above, if model very improbable, a small *fractional* reduction in
    !chi-squared from varying marginalised parameters can cause overflow
    if (lprob < 500D0) then
       prob = exp(lprob)
    else
       prob = exp(500D0)
    end if
    !!print '(5f9.1,g11.3)', lhd, pri, mg_nlnorm, post, lprob, prob

  end function prob

!==============================================================================

  function prob1d(val)

    double precision :: prob1d
    double precision, intent(in) :: val

    double precision, dimension(1) :: val1

    val1(1) = val
    prob1d = prob(1, val1)

  end function prob1d

!==============================================================================

  subroutine marg_err(sol, errguess, index, nlevidence, err)

    !subroutine arguments
    double precision, dimension(:), intent(in) :: sol !best fit parameters
    double precision, dimension(:), intent(in) :: errguess !initial guess errs
    integer, intent(in) :: index !index into sol
    double precision, intent(in) :: nlevidence
    !final +/- 1-sigma errors:
    double precision, dimension(2), intent(out) :: err

    !local variables
    double precision x0, x1, sigma, rtb
    double precision dx, xmid, fmid
    integer sign, nb, i, j, ifail
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
    ivar = index
    mg_marg = .true.
    mg_marg(index) = .false.
    do i = 1, size(x_pos, 1)
       mg_var(i) = fit_param(x_pos(i, 1), x_pos(i, 2))
    end do
    mg_param = fit_param
    mg_nlnorm = nlevidence

    !find 2 roots of f = nlmpost - min(nlmpost) - 0.5 = 0
    !(where nlmpost is negative log of marginalised posterior probability)

    !find value at minimum
    x0 = sol(index)
    sigma = errguess(index)
    minpost = 0D0 !so we can use ferr() for next line
    minpost = ferr(x0) + 0.5D0
    do j = 1, 2
       !find bracketing point other side of root
       do nb = 1, 10
          x1 = x0 + sign(j)*(2**nb)*sigma
          fmid = ferr(x1)
          if (fmid > 0D0) exit
       end do
       if (fmid <= 0D0) stop 'marg_err: cannot find bracketing point'
       !find root by quadratic interpolation
       err(j) = sign(j)*(x1-x0)/sqrt(2D0*fmid+1D0)
       !find root by bisection
       if (usenag) then
          ifail = 0
          call C05ADF(x0, x1, tol*sigma, 0D0, ferr, rtb, ifail)
          print *, err(j), sign(j)*(rtb-x0)
          err(j) = sign(j)*(rtb-x0)
       else
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
       end if
    end do

  end subroutine marg_err

!==============================================================================

  function ferr(val)

    double precision :: ferr
    double precision, intent(in) :: val

    !local variables
    character(len=128) :: info

    mg_var(ivar) = val !value for unmarginalised variable
    mg_param(x_pos(ivar, 1), x_pos(ivar, 2)) = val
    ferr = marg_post(info) - minpost - 0.5D0
    if (info /= '') print *, trim(info)
    
  end function ferr

!==============================================================================

end module Marginalise
