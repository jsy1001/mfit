!$Id: postplot.f90,v 1.1 2006/08/31 08:52:52 jsy1001 Exp $

module Postplot
  
  use Plot
  use Bayes
  use Wrap
  use Model
  use Marginalise

  implicit none

  public

  !public subroutines contained:
  !
  !plot_post1d - plot 1d cut through -ln(posterior) or
  !-ln(marginalised posterior)
  !
  !plot_post2d - plot 2d slice through -ln(posterior) or
  !-ln(marginalised posterior)

contains

  !============================================================================

  !! Plot 1d cut through -ln(posterior) or -ln(marginalised posterior)
  subroutine plot_post1d(plotmargd, nlevidence, inpar, indx, &
       x_title, y_title, top_title, uxmin, uxmax, device)

    !subroutine arguments
    logical, intent(in) :: plotmargd
    double precision, intent(in) :: nlevidence
    !! Model parameters (variable and non-variable)
    type(allparam), intent(in) :: inpar
    !! Index into variable parameters, specifying x-axis
    integer, intent(in) :: indx
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in) :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    type(allparam) :: plotpar
    logical :: mg_marg(inpar%nvar)
    real, allocatable :: post_points(:, :), mpost_points(:, :)
    double precision :: val, lhd, pri, lxbound, uxbound
    integer :: nvar, num_points, i, istat
    real :: xmin, xmax, ymin, ymax, pmin
    character(len=128) :: info

    !functions
    integer :: pgopen

    !return if nothing to plot
    nvar = inpar%nvar
    if (indx < 1 .or. indx > nvar) return
    lxbound = model_limits(inpar%var_pos(indx, 1), inpar%var_pos(indx, 2), 1)
    uxbound = model_limits(inpar%var_pos(indx, 1), inpar%var_pos(indx, 2), 2)
    if (uxmax < lxbound) return
    if (uxmin > uxbound) return

    !adjust x range if necessary - need to keep within bound limits
    if (uxmin < lxbound) then
       xmin = lxbound
    else
       xmin = uxmin
    end if
    if (uxmax > uxbound) then
       xmax = uxbound
    else
       xmax = uxmax
    end if

    !copy model parameters
    call allparam_copy(inpar, plotpar)

    !calculate grid of points
    if (plotmargd) then
       num_points = 20
       allocate(mpost_points(num_points, 2))
       mg_marg = .true.
       mg_marg(indx) = .false.
    else
       num_points = 50
    end if
    allocate(post_points(num_points, 2))
    do i = 1, num_points
       val = xmin + ((i-1.)/(num_points-1.))*(xmax-xmin)
       post_points(i, 1) = val
       call allparam_setone(plotpar, indx, val)
       if (plotmargd) then
          mpost_points(i, 1) = val
          mpost_points(i, 2) = marg_post(info, plotpar, mg_marg)
          if (info /= '') print *, trim(info)
       end if
       lhd = likelihood(vis_data, triple_data, model_spec, plotpar%param)
       pri = prior(plotpar%var_pos, plotpar%param, model_param, model_prior)
       post_points(i, 2) = lhd + pri - nlevidence !fix normalisation later
    end do

    !calculate y range
    if (plotmargd) then
       ymin = minval(mpost_points(:, 2))
       ymax = maxval(mpost_points(:, 2))
       pmin = minval(post_points(:, 2))
       post_points(:, 2) = post_points(:, 2) + (ymin - pmin)
    else
       ymin = minval(post_points(:, 2))
       ymax = maxval(post_points(:, 2))
    end if

    !plot posterior
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    if (.not. plotmargd) then
       call pgsci(1)
    else
       call pgsci(7)
       call pgsls(2)
    end if
    call plot_pts(x_title, '-ln(postprob)', num_points, post_points, -1, &
         'post1d.dat')
    call pgsci(1)
    call pgsls(1)
    if (plotmargd) &
         call plot_pts(x_title, y_title, num_points, mpost_points, -1, &
         'mpost1d.dat')

    call allparam_free(plotpar)
    if (allocated(post_points)) deallocate(post_points)
    if (allocated(mpost_points)) deallocate(mpost_points)

  end subroutine plot_post1d

  !============================================================================

  !! Plot 2d slice through -ln(posterior) or -ln(marginalised posterior)
  subroutine plot_post2d(plotmargd, nlposterior, inpar, indx, &
       x_title, y_title, top_title, uxmin, uxmax, uymin, uymax, device)

    !subroutine arguments
    logical, intent(in) :: plotmargd
    !! Minimum of negative log posterior (in case this is not
    !! sufficiently close to a grid point), needed to plot correct n-sigma
    !! contours in unmarginalised case
    double precision, intent(in) :: nlposterior
    !! Model parameters (variable and non-variable)
    type(allparam), intent(in) :: inpar
    !! Indices into variable parameters, specifying axes for plot
    integer, intent(in) :: indx(2)
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in) :: uxmin, uxmax
    double precision, intent(in) :: uymin, uymax
    character(len=*), intent(in), optional :: device

    !local variables
    type(allparam) :: plotpar
    logical :: mg_marg(inpar%nvar)
    integer, parameter :: iunit = 12
    character(len=20) :: save_filename
    real, allocatable :: post_points(:, :)
    double precision :: xval, yval, lhd, pri
    double precision :: lxbound, uxbound, lybound, uybound
    double precision :: sol(2)
    integer :: nvar, num_points, i, j, istat
    real :: xmin, xmax, ymin, ymax
    character(len=128) :: info
    logical :: invgray
    real :: tr(6)
    real :: fg, bg
    integer, parameter :: nlevs_max = 3
    !delta (chi^2)/2 values for 2-parameter 1-, 2-, 3-sigma confidence regions
    real :: levs(nlevs_max) = (/1.15,3.08,5.90/)

    !functions
    integer :: pgopen

    !return if nothing to plot
    nvar = inpar%nvar
    if (indx(1) < 1 .or. indx(1) > nvar) return
    if (indx(2) < 1 .or. indx(2) > nvar) return
    lxbound = &
         model_limits(inpar%var_pos(indx(1), 1), inpar%var_pos(indx(1), 2), 1)
    uxbound = &
         model_limits(inpar%var_pos(indx(1), 1), inpar%var_pos(indx(1), 2), 2)
    if (uxmax < lxbound) return
    if (uxmin > uxbound) return
    lybound = &
         model_limits(inpar%var_pos(indx(2), 1), inpar%var_pos(indx(2), 2), 1)
    uybound = &
         model_limits(inpar%var_pos(indx(2), 1), inpar%var_pos(indx(2), 2), 2)
    if (uymax < lybound) return
    if (uymin > uybound) return

    !adjust ranges if necessary - need to keep within bound limits
    if (uxmin < lxbound) then
       xmin = lxbound
    else
       xmin = uxmin
    end if
    if (uxmax > uxbound) then
       xmax = uxbound
    else
       xmax = uxmax
    end if
    if (uymin < lybound) then
       ymin = lybound
    else
       ymin = uymin
    end if
    if (uymax > uybound) then
       ymax = uybound
    else
       ymax = uymax
    end if

    !copy model parameters
    call allparam_copy(inpar, plotpar)
    if (plotmargd) then
       mg_marg = .true.
       mg_marg(indx(1)) = .false.
       mg_marg(indx(2)) = .false.
    end if
    sol(1) = inpar%param(inpar%var_pos(indx(1),1), inpar%var_pos(indx(1),2))
    sol(2) = inpar%param(inpar%var_pos(indx(2),1), inpar%var_pos(indx(2),2))

    !calculate grid of points
    if (plotmargd) then
       num_points = 20
       save_filename = 'mpost2d.dat'
    else
       num_points = 40
       save_filename = 'post2d.dat'
    end if
    allocate(post_points(num_points, num_points))
    do i = 1, num_points
       xval = xmin + ((i-1.)/(num_points-1.))*(xmax-xmin)
       call allparam_setone(plotpar, indx(1), xval)
       do j = 1, num_points
          yval = ymin + ((j-1.)/(num_points-1.))*(ymax-ymin)
          call allparam_setone(plotpar, indx(2), yval)
          if (plotmargd) then
             post_points(i, j) = marg_post(info, plotpar, mg_marg)
             if (info /= '') print *, trim(info)
          else
             lhd = likelihood(vis_data, triple_data, model_spec, plotpar%param)
             pri = prior(plotpar%var_pos, plotpar%param, &
                  model_param, model_prior)
             post_points(i, j) = lhd + pri - nlposterior
          end if
       end do
    end do

    !plot posterior
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, -2)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(1)

    !mapping from array elements to sky coordinates
    tr(2) = (xmax-xmin)/real(num_points-1)
    tr(6) = (ymax-ymin)/real(num_points-1)
    tr(1) = xmin - tr(2)
    tr(4) = ymin - tr(6)
    tr(3) = 0.
    tr(5) = 0.

    print *, minval(post_points)
    if (plotmargd) then
       !(marginalisation can shift minimum)
       post_points = post_points - minval(post_points)
    end if

    invgray = .true. !make minimum bright
    if (invgray) then
       fg = minval(post_points)
       bg = 0.6*maxval(post_points)
    else
       fg = maxval(post_points)
       bg = minval(post_points)
    endif
    call pggray(post_points, num_points, num_points, &
         1, num_points, 1, num_points, fg, bg, tr)
      
    write (51, *) (levs(i), i=1, size(levs,1))
    call pgsci(2)
    call pgcont(post_points, num_points, num_points, &
         1, num_points, 1, num_points, levs, size(levs,1), tr)
    call pgsci(1)
    call pgbox('ABCNST', 0., 0, 'ABCNST', 0., 0)
    call pgsci(3)
    call pgpt1(real(sol(1)), real(sol(2)), 12) !mark pos'n found by minimiser

    !write plotted points to text file
    open (unit=iunit, file=save_filename, status='replace', action='write', &
         err=91)
    write (iunit, '(a)') &
         '# -ln(norm. posterior probability) points from last such mfit plot'
    write (iunit, '(a, 3a25)') '# ', trim(x_title), trim(y_title), &
         trim(top_title)
    do j = 1, num_points
       yval = ymin + ((j-1.)/(num_points-1.))*(ymax-ymin)
       do i = 1, num_points
          xval = xmin + ((i-1.)/(num_points-1.))*(xmax-xmin)
          write (iunit, '(2x, 3f25.6)') xval, yval, post_points(i,j)
       end do
       write (iunit, *)
    end do
    close (iunit)

    call allparam_free(plotpar)
    if (allocated(post_points)) deallocate(post_points)
    return

91  print *, 'Cannot open file '//trim(save_filename)

  end subroutine plot_post2d

  !============================================================================

end module Postplot
