!$Id: plot.f90,v 1.19 2006/08/16 09:20:22 jsy1001 Exp $

module Plot
  
  use Model
  use Fit
  use Marginalise
  use Visibility

  !subroutines contained:
  !
  !plot_pts - call PGPLOT to plot points; save points to file
  !
  !plot_model_pts - wrapper for plot_pts, used for model points
  !
  !plot_data_pts - call PGPLOT to plot data points; optionally save to file
  !
  !plot_triple_phase_bas - plot closure phase against longest proj. baseline
  !in triangle
  !
  !plot_triple_amp_bas - plot triple product amplitude against longest proj.
  !baseline in triangle
  !
  !plot_vis_bas - plot squared visibility against projected baseline
  !
  !plot_uv - plot uv coverage
  !
  !plot_post - plot 1d cut through -ln(posterior) or
  !-ln(marginalised posterior)
  !
  !plot_post2d - plot 2d slice through -ln(posterior) or
  !-ln(marginalised posterior)
  !
  !plot_vis - plot squared visibility against specified data column
  !
  !plot_triple_phase - plot closure phase against specified data column
  !
  !plot_triple_amp - plot triple product amplitude specified data column

  implicit none

contains

  !============================================================================

  subroutine plot_pts(x_title, y_title, n, pts, symbol, save_filename)

    !subroutine arguments
    character(len=*), intent(in) :: x_title, y_title, save_filename
    integer, intent(in) :: n, symbol
    real, dimension(:, :), intent(in) :: pts

    !local variables
    integer, parameter :: iunit = 12
    integer :: i

    !plot points/line
    if (symbol > 0) then
       call pgpt(n, pts(:, 1), pts(:, 2), symbol)
    else
       call pgline(n, pts(:, 1), pts(:, 2))
    end if

    !save points to text file, suitable for use with e.g. gnuplot
    open (unit=iunit, file=save_filename, status='replace', action='write', &
         err=91)
    write (iunit, '(a)') '# Model points from last mfit plot'
    write (iunit, '(a, 2a25)') '# ', trim(x_title), trim(y_title)
    do i = 1, n
       write (iunit, '(2x, 2f25.6)') pts(i, 1), pts(i, 2)
    end do
    close (iunit)

    return

91  print *, 'Cannot open file '//trim(save_filename)

  end subroutine plot_pts

  !============================================================================

  subroutine plot_model_pts(x_title, y_title, n, pts, symbol)

    !subroutine arguments
    character(len=*), intent(in) :: x_title, y_title
    integer, intent(in) :: n, symbol
    real, dimension(:, :), intent(in) :: pts

    !local variables
    character(len=*), parameter :: save_filename = 'model_pts.dat'

    call plot_pts(x_title, y_title, n, pts, symbol, save_filename)

  end subroutine plot_model_pts

  !============================================================================

  subroutine plot_data_pts(x_title, y_title, n, pts, symbol, savepts)

    !subroutine arguments
    character(len=*), intent(in) :: x_title, y_title
    integer, intent(in) :: n, symbol
    real, dimension(:, :), intent(in) :: pts
    logical, intent(in) :: savepts

    !local variables
    integer, parameter :: iunit = 12
    character(len=*), parameter :: save_filename = 'data_pts.dat'
    integer :: i

    !plot points/line
    if (symbol > 0) then
       call pgpt(n, pts(:, 1), pts(:, 2), symbol)
       call pgerry(n, pts(:, 1), pts(:, 3), pts(:, 4), 1.0)
    else
       call pgline(n, pts(:, 1), pts(:, 2))
       call pgerry(n, pts(:, 1), pts(:, 3), pts(:, 4), 1.0)
    end if

    if (savepts) then
       !save points to text file, suitable for use with e.g. gnuplot
       open (unit=iunit, file=save_filename, status='replace', &
            action='write', err=91)
       write (iunit, '(a)') '# Data points from last mfit plot'
       write (iunit, '(a, 3a25)') '# ', trim(x_title), trim(y_title), 'Error'
       do i = 1, n
          write (iunit, '(2x, 3f25.6)') pts(i, 1), pts(i, 2), &
               0.5*abs(pts(i, 4) - pts(i, 3))
       end do
       close (iunit)
    end if

    return

91  print *, 'Cannot open file '//trim(save_filename)

  end subroutine plot_data_pts

  !============================================================================

  subroutine plot_triple_phase_bas(spec, param, x_title, y_title, top_title, &
       uxmin, uxmax, device)

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, model_points
    real :: xmin, xmax, ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax
    double precision :: lambda, delta_lambda, u1, v1, u2, v2, mjd
    double complex :: vis1, vis2, vis3
    real :: phase, phase_err, model_phase
    real, dimension(3) :: bas
    integer :: num_data, num_flagged, num_model, i, istat

    !functions
    integer :: pgopen

    ! make up real arrays for unflagged & flagged data points plus model points
    ! columns are x, y, (y+delta, y-delta)
    allocate(data_points(size(triple_data, 1), 4))
    allocate(flagged_points(size(triple_data, 1), 4))
    allocate(model_points(size(triple_data, 1), 2))
    num_model = 0
    num_data = 0
    num_flagged = 0
    do i = 1, size(triple_data, 1)
       lambda = triple_data(i, 1)
       delta_lambda = triple_data(i, 2)
       u1 = triple_data(i, 3)
       v1 = triple_data(i, 4)
       u2 = triple_data(i, 5)
       v2 = triple_data(i, 6)
       mjd = triple_data(i, 11)
       bas(1) = 1000.*sqrt(u1**2. + v1**2.)/lambda
       bas(2) = 1000.*sqrt(u2**2. + v2**2.)/lambda
       bas(3) = 1000.*sqrt((u1+u2)**2. + (v1+v2)**2.)/lambda
       if (present(uxmin) .and. maxval(bas) .lt. uxmin) cycle!out of plot range
       if (present(uxmax) .and. maxval(bas) .gt. uxmax) cycle!out of plot range
       phase = modulo(triple_data(i, 9), 360D0)
       if (phase > 180.) phase = phase - 360.
       phase_err = triple_data(i, 10)
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
            -(u1+u2), -(v1+v2), mjd)
       num_model = num_model + 1
       model_points(num_model, 1) = maxval(bas)
       model_phase = modulo(rad2deg*argument(vis1*vis2*vis3), 360D0)
       if (model_phase > 180.) model_phase = model_phase - 360.
       model_points(num_model, 2) = model_phase
       if (phase_err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = maxval(bas)
          flagged_points(num_flagged, 2) = phase
          flagged_points(num_flagged, 3) = phase - phase_err
          flagged_points(num_flagged, 4) = phase + phase_err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = maxval(bas)
          data_points(num_data, 2) = phase
          data_points(num_data, 3) = phase + phase_err
          data_points(num_data, 4) = phase - phase_err
       end if
    end do

    !calculate y range required
    dymin = minval(data_points(:num_data, 4))
    dymax = maxval(data_points(:num_data, 3))
    fymin = minval(flagged_points(:num_flagged, 4))
    fymax = maxval(flagged_points(:num_flagged, 3))
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if (num_data > 0) then
       !exclude flagged extrema by default
       ymin = min(dymin, mymin)
       ymax = max(dymax, mymax)
    else if (num_flagged > 0) then
       ymin = min(fymin, mymin)
       ymax = max(fymax, mymax)
    else
       !nothing to plot
       return
    end if
    !change y range
    if (ymax > 180.) then
       ymax = 180.
    else if (ymax < -180.) then
       ymax = -180.
    end if

    !calculate x range required
    if (present(uxmin)) then
       xmin = uxmin
    else
       xmin = 0.9*minval(model_points(:num_model, 1))
    end if
    if (present(uxmax)) then
       xmax = uxmax
    else
       xmax = 1.1*maxval(model_points(:num_model, 1))
    end if

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 13)

    if (allocated(data_points)) deallocate(data_points)
    if (allocated(flagged_points)) deallocate(flagged_points)
    if (allocated(model_points)) deallocate(model_points)

  end subroutine plot_triple_phase_bas
  !============================================================================

  subroutine plot_triple_amp_bas(spec, param, x_title, y_title, top_title, &
       uxmin, uxmax, device)

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, model_points
    real :: xmin, xmax, ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax
    double precision :: lambda, delta_lambda, u1, v1, u2, v2, mjd
    double complex :: vis1, vis2, vis3
    real :: amp, amp_err
    real, dimension(3) :: bas
    integer :: num_data, num_flagged, num_model, i, istat

    !functions
    integer :: pgopen
    ! make up real arrays for unflagged & flagged data points plus model points
    ! columns are x, y, (y+delta, y-delta)
    allocate(data_points(size(triple_data, 1), 4))
    allocate(flagged_points(size(triple_data, 1), 4))
    allocate(model_points(size(triple_data, 1), 2))
    num_model = 0
    num_data = 0
    num_flagged = 0
    do i = 1, size(triple_data, 1)
       lambda = triple_data(i, 1)
       delta_lambda = triple_data(i, 2)
       u1 = triple_data(i, 3)
       v1 = triple_data(i, 4)
       u2 = triple_data(i, 5)
       v2 = triple_data(i, 6)
       mjd = triple_data(i, 11)
       bas(1) = 1000.*sqrt(u1**2. + v1**2.)/lambda
       bas(2) = 1000.*sqrt(u2**2. + v2**2.)/lambda
       bas(3) = 1000.*sqrt((u1+u2)**2. + (v1+v2)**2.)/lambda
       if (present(uxmin) .and. maxval(bas) .lt. uxmin) cycle!out of plot range
       if (present(uxmax) .and. maxval(bas) .gt. uxmax) cycle!out of plot range
       amp = triple_data(i, 7)
       amp_err = triple_data(i, 8)
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
            -(u1+u2), -(v1+v2), mjd)
       num_model = num_model + 1
       model_points(num_model, 1) = maxval(bas)
       model_points(num_model, 2) = modulus(vis1*vis2*vis3)
       if (amp_err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = maxval(bas)
          flagged_points(num_flagged, 2) = amp
          flagged_points(num_flagged, 3) = amp - amp_err
          flagged_points(num_flagged, 4) = amp + amp_err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = maxval(bas)
          data_points(num_data, 2) = amp
          data_points(num_data, 3) = amp + amp_err
          data_points(num_data, 4) = amp - amp_err
       end if
    end do

    !calculate y range required
    dymin = minval(data_points(:num_data, 4))
    dymax = maxval(data_points(:num_data, 3))
    fymin = minval(flagged_points(:num_flagged, 4))
    fymax = maxval(flagged_points(:num_flagged, 3))
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if (num_data > 0) then
       !exclude flagged extrema by default
       ymin = min(dymin, mymin)
       ymax = max(dymax, mymax)
    else if (num_flagged > 0) then
       ymin = min(fymin, mymin)
       ymax = max(fymax, mymax)
    else
       !nothing to plot
       return
    end if
    !change y range
    if (ymin > 0) then
       ymin = ymin*0.9
    else if (ymin < 0) then !yes, may be -ve
       ymin = ymin*1.2
    end if
    ymax = ymax*1.1

    !calculate x range required
    if (present(uxmin)) then
       xmin = uxmin
    else
       xmin = 0.9*minval(data_points(:num_data, 1))
    end if
    if (present(uxmax)) then
       xmax = uxmax
    else
       xmax = 1.1*maxval(data_points(:num_data, 1))
    end if

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 7)

    if (allocated(data_points)) deallocate(data_points)
    if (allocated(flagged_points)) deallocate(flagged_points)
    if (allocated(model_points)) deallocate(model_points)

  end subroutine plot_triple_amp_bas

  !============================================================================

  subroutine plot_vis_bas(spec, param, mod_line, x_title, y_title, top_title, &
       uxmin, uxmax, device)

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    logical, intent(in) :: mod_line !plot continuous line for model?
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, &
         model_points
    double precision :: u, v, lambda, delta_lambda, mjd
    integer :: num_data, num_flagged, num_model, i, istat
    real :: bas, vsq, err, xmin, xmax
    real :: ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax
    integer, parameter :: grid_size = 200

    !functions
    integer :: pgopen

    ! make up real arrays for unflagged & flagged data points
    ! columns are x, y, y+delta, y-delta
    allocate(data_points(size(vis_data, 1), 4))
    allocate(flagged_points(size(vis_data, 1), 4))
    num_data = 0
    num_flagged = 0
    do i = 1, size(vis_data, 1)
       lambda = vis_data(i, 1)
       delta_lambda = vis_data(i, 2)
       u = vis_data(i, 3)
       v = vis_data(i, 4)
       bas = 1000.*sqrt(u**2. + v**2.)/lambda
       if (present(uxmin) .and. bas .lt. uxmin) cycle !out of plot range
       if (present(uxmax) .and. bas .gt. uxmax) cycle !out of plot range
       vsq = vis_data(i, 5)
       err = vis_data(i, 6)
       if (err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = bas
          flagged_points(num_flagged, 2) = vsq
          flagged_points(num_flagged, 3) = vsq - err
          flagged_points(num_flagged, 4) = vsq + err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = bas
          data_points(num_data, 2) = vsq
          data_points(num_data, 3) = vsq + err
          data_points(num_data, 4) = vsq - err
       end if
    end do

    !calculate x range required
    if (present(uxmin)) then
       xmin = uxmin
    else
       xmin = 0.9*minval(data_points(:num_data, 1))
    end if
    if (present(uxmax)) then
       xmax = uxmax
    else
       xmax = 1.1*maxval(data_points(:num_data, 1))
    end if

    if (mod_line) then
       !XXX can't have time-dependent model?
       !calculate grid of model points
       num_model = grid_size
       allocate(model_points(num_model, 2))
       lambda = vis_data(1, 1)
       delta_lambda = vis_data(1, 2)
       do i = 1, num_model
          u = xmin*lambda/1000. + ((i-1.)/(num_model-1.))*(xmax-xmin)*lambda/1000.
          model_points(i, 1) = 1000.*u/lambda
          model_points(i, 2) = modulus(cmplx_vis(spec, param, lambda, delta_lambda, u, 0D0, mjd))**2.
       end do
    else
       !calculate model points corresponding to plotted data points
       allocate(model_points(size(vis_data, 1), 2))
       num_model = 0
       do i = 1, size(vis_data, 1)
          lambda = vis_data(i, 1)
          delta_lambda = vis_data(i, 2)
          u = vis_data(i, 3)
          v = vis_data(i, 4)
          mjd = vis_data(i, 7)
          bas = 1000.*sqrt(u**2. + v**2.)/lambda
          if (bas .ge. xmin .and. bas .le. xmax) then
             num_model = num_model + 1
             model_points(num_model, 1) = bas
             model_points(num_model, 2) = modulus(cmplx_vis(spec, param, lambda, delta_lambda, u, v, mjd))**2.
          end if
       end do
    end if

    !calculate y range
    dymin = minval(data_points(:num_data, 4))
    dymax = maxval(data_points(:num_data, 3))
    fymin = minval(flagged_points(:num_flagged, 4))
    fymax = maxval(flagged_points(:num_flagged, 3))
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if (num_data > 0) then
       !exclude flagged extrema by default
       ymin = min(dymin, mymin)
       ymax = max(dymax, mymax)
    else if (num_flagged > 0) then
       ymin = min(fymin, mymin)
       ymax = max(fymax, mymax)
    else if (num_model > 0) then
       !user x range excludes all data, ensure model points in range
       ymin = mymin
       ymax = mymax
    else
       !nothing to plot
       return
    end if
    !change y range
    if (ymin > 0) then
       ymin = ymin*0.9
    else if (ymin < 0) then !yes, may be -ve
       ymin = ymin*1.2
    end if
    ymax = ymax*1.1

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    if (mod_line) then
       call plot_model_pts(x_title, y_title, num_model, model_points, -1)
    else
       call plot_model_pts(x_title, y_title, num_model, model_points, 13)
    end if

    if (allocated(data_points)) deallocate(data_points)
    if (allocated(flagged_points)) deallocate(flagged_points)
    if (allocated(model_points)) deallocate(model_points)

  end subroutine plot_vis_bas

  !============================================================================

  subroutine plot_uv(x_title, y_title, top_title, device)

    !subroutine arguments
    character(len=*), intent(in) :: x_title, y_title, top_title
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points
    real :: u, v, lambda, delta_lambda, bas, max
    integer :: num_data, num_flagged, i, istat

    !functions
    integer :: pgopen

    ! make up real arrays for unflagged & flagged data points
    ! columns are x, y
    allocate(data_points(size(vis_data, 1), 2))
    allocate(flagged_points(size(vis_data, 1), 2))
    num_data = 0
    num_flagged = 0
    max = 0.
    do i = 1, size(vis_data, 1)
       lambda = vis_data(i, 1)
       delta_lambda = vis_data(i, 2)
       u = 1000.*vis_data(i, 3)/lambda
       v = 1000.*vis_data(i, 4)/lambda
       bas = sqrt(u**2. + v**2.)
       if (bas > max) max = bas
       if (vis_data(i, 6) <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = u
          flagged_points(num_flagged, 2) = v
       else
          num_data = num_data + 1
          data_points(num_data, 1) = u
          data_points(num_data, 2) = v
       end if
    end do

    !plot points
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(-max, max, -max, max, 1, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 13)
    call pgpt(num_flagged, -flagged_points(:, 1), -flagged_points(:, 2), 7)
    call pgsci(1)
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 17)
    call pgpt(num_data, -data_points(:, 1), -data_points(:, 2), 22)

    if (allocated(data_points)) deallocate(data_points)
    if (allocated(flagged_points)) deallocate(flagged_points)

  end subroutine plot_uv

  !============================================================================

  subroutine plot_post(plotmargd, nlevidence, indx, &
       x_title, y_title, top_title, uxmin, uxmax, device)

    !indx is into x_pos array of variable parameters
    !x_pos, x_bound must be initialised on entry to this routine i.e. must have
    !called minimiser() from Fit module

    !subroutine arguments
    logical, intent(in) :: plotmargd
    double precision, intent(in) :: nlevidence
    integer, intent(in) :: indx
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in) :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: post_points, mpost_points
    double precision :: val, lhd, pri, sav
    integer :: nvar, num_points, i, istat
    real :: xmin, xmax, ymin, ymax, pmin
    character(len=128) :: info

    !functions
    integer :: pgopen

    !return if nothing to plot
    nvar = size(x_pos, 1)
    if (indx < 1 .or. indx > nvar) return
    if (uxmax < x_bound(indx, 1)) return
    if (uxmin > x_bound(indx, 2)) return

    !adjust x range if necessary - need to keep within x_bound limits
    if (uxmin < x_bound(indx, 1)) then
       xmin = x_bound(indx, 1)
    else
       xmin = uxmin
    end if
    if (uxmax > x_bound(indx, 2)) then
       xmax = x_bound(indx, 2)
    else
       xmax = uxmax
    end if

    !copy model parameters
    sav = fit_param(x_pos(indx, 1), x_pos(indx, 2))

    !calculate grid of points
    if (plotmargd) then
       num_points = 20
       allocate(mpost_points(num_points, 2))
       mg_marg = .true.
       mg_marg(indx) = .false.
       mg_param = fit_param
    else
       num_points = 50
    end if
    allocate(post_points(num_points, 2))
    do i = 1, num_points
       val = xmin + ((i-1.)/(num_points-1.))*(xmax-xmin)
       post_points(i, 1) = val
       if (plotmargd) then
          mpost_points(i, 1) = val
          mg_param(x_pos(indx, 1), x_pos(indx, 2)) = val
          mpost_points(i, 2) = marg_post(info)
          if (info /= '') print *, trim(info)
       end if
       fit_param(x_pos(indx, 1), x_pos(indx, 2)) = val
       lhd = likelihood(vis_data, triple_data, model_spec, fit_param)
       pri = prior(x_pos, fit_param, model_param, model_prior)
       post_points(i, 2) = lhd + pri - nlevidence !fix normalisation later
    end do
    fit_param(x_pos(indx, 1), x_pos(indx, 2)) = sav

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

    if (allocated(post_points)) deallocate(post_points)
    if (allocated(mpost_points)) deallocate(mpost_points)

  end subroutine plot_post

  !============================================================================

  subroutine plot_post2d(plotmargd, nlposterior, indx, &
       x_title, y_title, top_title, uxmin, uxmax, uymin, uymax, device)

    !nlposterior is minimum of negative log posterior (in case this is not
    ! sufficiently close to a grid point), needed to plot correct n-sigma
    ! contours in unmarginalised case
    !indx(i) is into x_pos array of variable parameters
    !x_pos, x_bound must be initialised on entry to this routine i.e. must have
    ! called minimiser() from Fit module

    !subroutine arguments
    logical, intent(in) :: plotmargd
    double precision, intent(in) :: nlposterior
    integer, intent(in), dimension(2) :: indx
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in) :: uxmin, uxmax
    double precision, intent(in) :: uymin, uymax
    character(len=*), intent(in), optional :: device

    !local variables
    integer, parameter :: iunit = 12
    character(len=*), parameter :: save_filename = 'post2d.dat'
    real, dimension(:, :), allocatable :: post_points
    double precision :: xval, yval, lhd, pri
    double precision, dimension(2) :: sav
    integer :: nvar, num_points, i, j, istat
    real :: xmin, xmax, ymin, ymax
    character(len=128) :: info
    logical :: invgray
    real, dimension(6) :: tr
    real fg, bg
    integer, parameter :: nlevs_max=3
    !delta (chi^2)/2 values for 2-parameter 1-, 2-, 3-sigma confidence regions
    real, dimension(nlevs_max) :: levs=(/1.15,3.08,5.90/)

    !functions
    integer :: pgopen

    !return if nothing to plot
    nvar = size(x_pos, 1)
    if (indx(1) < 1 .or. indx(1) > nvar) return
    if (indx(2) < 1 .or. indx(2) > nvar) return
    if (uxmax < x_bound(indx(1), 1)) return
    if (uxmin > x_bound(indx(1), 2)) return
    if (uymax < x_bound(indx(2), 1)) return
    if (uymin > x_bound(indx(2), 2)) return

    !adjust ranges if necessary - need to keep within x_bound limits
    if (uxmin < x_bound(indx(1), 1)) then
       xmin = x_bound(indx(1), 1)
    else
       xmin = uxmin
    end if
    if (uxmax > x_bound(indx(1), 2)) then
       xmax = x_bound(indx(1), 2)
    else
       xmax = uxmax
    end if
    if (uymin < x_bound(indx(2), 1)) then
       ymin = x_bound(indx(2), 1)
    else
       ymin = uymin
    end if
    if (uymax > x_bound(indx(2), 2)) then
       ymax = x_bound(indx(2), 2)
    else
       ymax = uymax
    end if

    !copy model parameters
    sav(1) = fit_param(x_pos(indx(1), 1), x_pos(indx(1), 2))
    sav(2) = fit_param(x_pos(indx(2), 1), x_pos(indx(2), 2))
    if (plotmargd) then
       mg_marg = .true.
       mg_marg(indx(1)) = .false.
       mg_marg(indx(2)) = .false.
       mg_param = fit_param
    end if

    !calculate grid of points
    if (plotmargd) then
       num_points = 20
    else
       num_points = 40
    end if
    allocate(post_points(num_points, num_points))
    do i = 1, num_points
       xval = xmin + ((i-1.)/(num_points-1.))*(xmax-xmin)
       do j = 1, num_points
          yval = ymin + ((j-1.)/(num_points-1.))*(ymax-ymin)
          if (plotmargd) then
             mg_param(x_pos(indx(1), 1), x_pos(indx(1), 2)) = xval
             mg_param(x_pos(indx(2), 1), x_pos(indx(2), 2)) = yval
             post_points(i, j) = marg_post(info)
             if (info /= '') print *, trim(info)
          else
             fit_param(x_pos(indx(1), 1), x_pos(indx(1), 2)) = xval
             fit_param(x_pos(indx(2), 1), x_pos(indx(2), 2)) = yval
             lhd = likelihood(vis_data, triple_data, model_spec, fit_param)
             pri = prior(x_pos, fit_param, model_param, model_prior)
             post_points(i, j) = lhd + pri - nlposterior
          end if
       end do
    end do
    fit_param(x_pos(indx(1), 1), x_pos(indx(1), 2)) = sav(1)
    fit_param(x_pos(indx(2), 1), x_pos(indx(2), 2)) = sav(2)

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
    call pgpt1(real(sav(1)), real(sav(2)), 12) !mark pos'n found by minimiser

    !write plotted points to text file
    open (unit=iunit, file=save_filename, status='replace', action='write', &
         err=91)
    write (iunit, '(a)') &
         '# -ln(norm. posterior probability) points from last such mfit plot'
    write (iunit, '(a, 3a25)') '# ', trim(x_title), trim(y_title), &
         '-ln(postprob/peak)'
    do j = 1, num_points
       yval = ymin + ((j-1.)/(num_points-1.))*(ymax-ymin)
       do i = 1, num_points
          xval = xmin + ((i-1.)/(num_points-1.))*(xmax-xmin)
          write (iunit, '(2x, 3f25.6)') xval, yval, post_points(i,j)
       end do
       write (iunit, *)
    end do
    close (iunit)

    if (allocated(post_points)) deallocate(post_points)
    return

91  print *, 'Cannot open file '//trim(save_filename)

  end subroutine plot_post2d

  !============================================================================

  subroutine plot_vis(xindex, spec, param, x_title, y_title, top_title, &
       xzero, uxmin, uxmax, device)

    !plot squared visibility against specified data column
    !xindex gives index into 2nd axis of vis_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    real, intent(in) :: xzero
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, &
         model_points
    double precision :: u, v, lambda, delta_lambda, mjd
    integer :: num_data, num_flagged, num_model, i, istat
    real :: vsq, err, xmin, xmax
    real :: ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax

    !functions
    integer :: pgopen

    ! make up real arrays for unflagged & flagged data points
    ! columns are x, y, y+delta, y-delta
    allocate(data_points(size(vis_data, 1), 4))
    allocate(flagged_points(size(vis_data, 1), 4))
    num_data = 0
    num_flagged = 0
    do i = 1, size(vis_data, 1)
       !skip if x value outside plot range
       if (present(uxmin) .and. vis_data(i, xindex) .lt. uxmin) cycle
       if (present(uxmax) .and. vis_data(i, xindex) .gt. uxmax) cycle
       lambda = vis_data(i, 1)
       delta_lambda = vis_data(i, 2)
       u = vis_data(i, 3)
       v = vis_data(i, 4)
       vsq = vis_data(i, 5)
       err = vis_data(i, 6)
       if (err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = vis_data(i, xindex) - xzero
          flagged_points(num_flagged, 2) = vsq
          flagged_points(num_flagged, 3) = vsq - err
          flagged_points(num_flagged, 4) = vsq + err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = vis_data(i, xindex) - xzero
          data_points(num_data, 2) = vsq
          data_points(num_data, 3) = vsq + err
          data_points(num_data, 4) = vsq - err
       end if
    end do

    !calculate x range required
    if (present(uxmin)) then
       xmin = uxmin
    else
       xmin = minval(data_points(:num_data, 1))
    end if
    if (present(uxmax)) then
       xmax = uxmax
    else
       xmax = maxval(data_points(:num_data, 1))
    end if

    !calculate model points corresponding to plotted data points
    allocate(model_points(size(vis_data, 1), 2))
    num_model = 0
    do i = 1, size(vis_data, 1)
       lambda = vis_data(i, 1)
       delta_lambda = vis_data(i, 2)
       u = vis_data(i, 3)
       v = vis_data(i, 4)
       if ((vis_data(i, xindex)-xzero) .ge. xmin .and. (vis_data(i, xindex)-xzero) .le. xmax) then
          num_model = num_model + 1
          model_points(num_model, 1) = vis_data(i, xindex) - xzero
          model_points(num_model, 2) = modulus(cmplx_vis(spec, param, lambda, delta_lambda, u, v, mjd))**2.
       end if
    end do

    !calculate y range
    dymin = minval(data_points(:num_data, 4))
    dymax = maxval(data_points(:num_data, 3))
    fymin = minval(flagged_points(:num_flagged, 4))
    fymax = maxval(flagged_points(:num_flagged, 3))
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if (num_data > 0) then
       !exclude flagged extrema by default
       ymin = min(dymin, mymin)
       ymax = max(dymax, mymax)
    else if (num_flagged > 0) then
       ymin = min(fymin, mymin)
       ymax = max(fymax, mymax)
    else
       !nothing to plot
       return
    end if
    !change y range
    if (ymin > 0) then
       ymin = ymin*0.9
    else if (ymin < 0) then !yes, may be -ve
       ymin = ymin*1.2
    end if
    ymax = ymax*1.1

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 7)

    if (allocated(data_points)) deallocate(data_points)
    if (allocated(flagged_points)) deallocate(flagged_points)
    if (allocated(model_points)) deallocate(model_points)

  end subroutine plot_vis

  !============================================================================

  subroutine plot_triple_phase(xindex, spec, param, x_title, y_title, &
       top_title, xzero, uxmin, uxmax, device)

    !plot triple product phase against specified data column
    !xindex gives index into 2nd axis of triple_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    real, intent(in) :: xzero
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, &
         model_points
    real :: xmin, xmax, ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax
    double precision :: lambda, delta_lambda, u1, v1, u2, v2, mjd
    double complex :: vis1, vis2, vis3
    real :: phase, phase_err, model_phase
    integer :: num_data, num_flagged, num_model, i, istat

    !functions
    integer :: pgopen

    ! make up real arrays for unflagged & flagged data points plus model points
    ! columns are x, y, y+delta, y-delta
    allocate(data_points(size(triple_data, 1), 4))
    allocate(flagged_points(size(triple_data, 1), 4))
    allocate(model_points(size(triple_data, 1), 2))
    num_model = 0
    num_data = 0
    num_flagged = 0
    do i = 1, size(triple_data, 1)
       !skip if x value outside plot range
       if (present(uxmin) .and. triple_data(i, xindex) .lt. uxmin) cycle
       if (present(uxmax) .and. triple_data(i, xindex) .gt. uxmax) cycle
       lambda = triple_data(i, 1)
       delta_lambda = triple_data(i, 2)
       phase = modulo(triple_data(i, 9), 360D0)
       if (phase > 180.) phase = phase - 360.
       phase_err = triple_data(i, 10)
       u1 = triple_data(i, 3)
       v1 = triple_data(i, 4)
       u2 = triple_data(i, 5)
       v2 = triple_data(i, 6)
       mjd = triple_data(i, 11)
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
            -(u1+u2), -(v1+v2), mjd)
       num_model = num_model + 1
       model_points(num_model, 1) = triple_data(i, xindex) - xzero
       model_phase = modulo(rad2deg*argument(vis1*vis2*vis3), 360D0)
       if (model_phase > 180.) model_phase = model_phase - 360.
       model_points(num_model, 2) = model_phase
       if (phase_err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = triple_data(i, xindex) - xzero
          flagged_points(num_flagged, 2) = phase
          flagged_points(num_flagged, 3) = phase - phase_err
          flagged_points(num_flagged, 4) = phase + phase_err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = triple_data(i, xindex) - xzero
          data_points(num_data, 2) = phase
          data_points(num_data, 3) = phase + phase_err
          data_points(num_data, 4) = phase - phase_err
       end if
    end do

    !calculate y range required
    dymin = minval(data_points(:num_data, 4))
    dymax = maxval(data_points(:num_data, 3))
    fymin = minval(flagged_points(:num_flagged, 4))
    fymax = maxval(flagged_points(:num_flagged, 3))
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if (num_data > 0) then
       !exclude flagged extrema by default
       ymin = min(dymin, mymin)
       ymax = max(dymax, mymax)
    else if (num_flagged > 0) then
       ymin = min(fymin, mymin)
       ymax = max(fymax, mymax)
    else
       !nothing to plot
       return
    end if
    !change y range
    if (ymax > 180.) then
       ymax = 180.
    else if (ymax < -180.) then
       ymax = -180.
    end if

    !calculate x range required
    if (present(uxmin)) then
       xmin = uxmin
    else
       xmin = minval(data_points(:num_data, 1))
    end if
    if (present(uxmax)) then
       xmax = uxmax
    else
       xmax = maxval(data_points(:num_data, 1))
    end if

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 7)

    if (allocated(data_points)) deallocate(data_points)
    if (allocated(flagged_points)) deallocate(flagged_points)
    if (allocated(model_points)) deallocate(model_points)

  end subroutine plot_triple_phase

  !============================================================================

  subroutine plot_triple_amp(xindex, spec, param, x_title, y_title, &
       top_title, xzero, uxmin, uxmax, device)

    !plot triple product phase against specified data column
    !xindex gives index into 2nd axis of triple_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    real, intent(in) :: xzero
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, &
         model_points
    real :: xmin, xmax, ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax
    double precision :: lambda, delta_lambda, u1, v1, u2, v2, mjd
    double complex :: vis1, vis2, vis3
    real :: amp, amp_err
    integer :: num_data, num_flagged, num_model, i, istat

    !functions
    integer :: pgopen

    ! make up real arrays for unflagged & flagged data points plus model points
    ! columns are x, y, y+delta, y-delta
    allocate(data_points(size(triple_data, 1), 4))
    allocate(flagged_points(size(triple_data, 1), 4))
    allocate(model_points(size(triple_data, 1), 2))
    num_model = 0
    num_data = 0
    num_flagged = 0
    do i = 1, size(triple_data, 1)
       !skip if x value outside plot range
       if (present(uxmin) .and. triple_data(i, xindex) .lt. uxmin) cycle
       if (present(uxmax) .and. triple_data(i, xindex) .gt. uxmax) cycle
       lambda = triple_data(i, 1)
       delta_lambda = triple_data(i, 2)
       amp = triple_data(i, 7)
       amp_err = triple_data(i, 8)
       u1 = triple_data(i, 3)
       v1 = triple_data(i, 4)
       u2 = triple_data(i, 5)
       v2 = triple_data(i, 6)
       mjd = triple_data(i, 11)
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
            -(u1+u2), -(v1+v2), mjd)
       num_model = num_model + 1
       model_points(num_model, 1) = triple_data(i, xindex) - xzero
       model_points(num_model, 2) = modulus(vis1*vis2*vis3)
       if (amp_err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = triple_data(i, xindex) - xzero
          flagged_points(num_flagged, 2) = amp
          flagged_points(num_flagged, 3) = amp - amp_err
          flagged_points(num_flagged, 4) = amp + amp_err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = triple_data(i, xindex) - xzero
          data_points(num_data, 2) = amp
          data_points(num_data, 3) = amp + amp_err
          data_points(num_data, 4) = amp - amp_err
       end if
    end do

    !calculate y range required
    dymin = minval(data_points(:num_data, 4))
    dymax = maxval(data_points(:num_data, 3))
    fymin = minval(flagged_points(:num_flagged, 4))
    fymax = maxval(flagged_points(:num_flagged, 3))
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if (num_data > 0) then
       !exclude flagged extrema by default
       ymin = min(dymin, mymin)
       ymax = max(dymax, mymax)
    else if (num_flagged > 0) then
       ymin = min(fymin, mymin)
       ymax = max(fymax, mymax)
    else
       !nothing to plot
       return
    end if
    !change y range
    if (ymin > 0) then
       ymin = ymin*0.9
    else if (ymin < 0) then !yes, may be -ve
       ymin = ymin*1.2
    end if
    ymax = ymax*1.1

    !calculate x range required
    if (present(uxmin)) then
       xmin = uxmin
    else
       xmin = minval(data_points(:num_data, 1))
    end if
    if (present(uxmax)) then
       xmax = uxmax
    else
       xmax = maxval(data_points(:num_data, 1))
    end if

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 7)

    if (allocated(data_points)) deallocate(data_points)
    if (allocated(flagged_points)) deallocate(flagged_points)
    if (allocated(model_points)) deallocate(model_points)

  end subroutine plot_triple_amp

  !============================================================================

end module Plot
