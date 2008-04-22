!$Id: plot.f90,v 1.24 2008/04/22 12:35:12 jsy1001 Exp $

module Plot
  
  use Maths
  use Visibility
  use Bayes, only : vis_data, triple_data, num_vis, num_triple

  implicit none

  private :: plot_vis_xval, get_plot_vis_data, &
       plot_triple_xval, get_plot_triple_data
  public

  !public subroutines contained:
  !
  !plot_pts - call PGPLOT to plot points; save points to file
  !
  !plot_model_pts - wrapper for plot_pts, used for model points
  !
  !plot_data_pts - call PGPLOT to plot data points; optionally save to file
  !
  !plot_vis_bas - plot squared visibility against projected baseline
  !
  !plot_uv - plot uv coverage
  !
  !plot_vis - plot squared visibility against specified quantity
  !
  !plot_triple_phase - plot closure phase against specified quantity
  !
  !plot_triple_amp - plot triple product amplitude specified quantity

contains

  !============================================================================

  function plot_vis_xval(i, j)

    !returns x-axis value for plot
    !
    !i is index into 1st axis of vis_data
    !If positive, j is index into 2nd axis of vis_data, otherwise:
    !j=-1: return projected baseline length /mega-lambda
    !j=-2: return GMST /h
    !j=-3: return baseline position angle /deg

    !function arguments
    integer, intent(in) :: i, j
    double precision :: plot_vis_xval

    !local variables
    double precision :: u, v, lambda, mjd

    !functions
    double precision :: sla_gmst

    if (j>=1 .and. j<=size(vis_data, 2)) then
       plot_vis_xval = vis_data(i, j)
    else if (j == -1) then
       lambda = vis_data(i, 1)
       u = vis_data(i, 3)
       v = vis_data(i, 4)
       plot_vis_xval = 1000.*sqrt(u**2. + v**2.)/lambda
    else if (j == -2) then
       mjd = vis_data(i, 7)
       !note sla_gmsta would give better precision
       plot_vis_xval = rad2deg/15D0*sla_gmst(mjd) !convert to hours
    else if (j == -3) then
       u = vis_data(i, 3)
       v = vis_data(i, 4)
       plot_vis_xval = rad2deg*atan2(u, v)
    else
       stop 'Illegal value for j in plot_vis_xval()'
    end if

  end function plot_vis_xval

  !============================================================================

  subroutine get_plot_vis_data(xindex, xzero, &
       data_points, num_data, flagged_points, num_flagged, xrange, yrange, &
       uxmin, uxmax)

    !If positive, xindex is index into 2nd axis of vis_data, otherwise:
    !-1: projected baseline length /mega-lambda
    !-2: GMST /h
    !-3: baseline position angle /deg

    !subroutine arguments
    integer, intent(in) :: xindex
    double precision, intent(in), optional :: uxmin, uxmax
    real, intent(out) :: xzero
    real, dimension(:, :), intent(out) :: data_points, flagged_points
    integer, intent(out) :: num_data, num_flagged
    real, dimension(2), intent(out) :: xrange, yrange

    !local variables
    integer :: i
    real :: xmin, x, xspan
    double precision :: lambda, delta_lambda, u, v, vsq, err

    !choose xzero
    xzero = 0.
    if(xindex >= 1 .and. xindex <= size(vis_data, 2)) then
       xmin = minval(vis_data(:, xindex))
       if(abs(xmin) > 100.*(maxval(vis_data(:, xindex)) - xmin)) &
            xzero = xmin
    end if

    !make up data arrays
    num_data = 0
    num_flagged = 0
    do i = 1, num_vis
       x = plot_vis_xval(i, xindex)
       !skip if x value outside plot range
       if (present(uxmin)) then
          if ((x-xzero) < uxmin) cycle
       end if
       if (present(uxmax)) then
          if ((x-xzero) > uxmax) cycle
       end if
       lambda = vis_data(i, 1)
       delta_lambda = vis_data(i, 2)
       u = vis_data(i, 3)
       v = vis_data(i, 4)
       vsq = vis_data(i, 5)
       err = vis_data(i, 6)
       if (err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = x - xzero
          flagged_points(num_flagged, 2) = vsq
          flagged_points(num_flagged, 3) = vsq - err
          flagged_points(num_flagged, 4) = vsq + err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = x - xzero
          data_points(num_data, 2) = vsq
          data_points(num_data, 3) = vsq + err
          data_points(num_data, 4) = vsq - err
       end if
    end do
    if(num_data + num_flagged == 0) return

    !calculate x range required
    if (present(uxmin)) then
       xrange(1) = uxmin
    else
       if (num_data > 0) then
          !exclude flagged extrema by default
          xrange(1) = minval(data_points(:num_data, 1))
       else
          xrange(1) = minval(flagged_points(:num_flagged, 1))
       end if
    end if
    if (present(uxmax)) then
       xrange(2) = uxmax
    else
       if (num_data > 0) then
          xrange(2) = maxval(data_points(:num_data, 1))
       else
          xrange(2) = maxval(flagged_points(:num_flagged, 1))
       end if
    end if
    xspan = xrange(2) - xrange(1)
    xrange(1) = xrange(1) - 0.05*xspan
    xrange(2) = xrange(2) + 0.05*xspan

    !calculate y range
    if (num_data > 0) then
       yrange(1) = minval(data_points(:num_data, 4))
       yrange(2) = maxval(data_points(:num_data, 3))
    else
       yrange(1) = minval(flagged_points(:num_flagged, 4))
       yrange(2) = maxval(flagged_points(:num_flagged, 3))
    end if

    !expand y range
    if (yrange(1) > 0) then
       yrange(1) = yrange(1)*0.9
    else if (yrange(1) < 0) then !yes, may be -ve
       yrange(1) = yrange(1)*1.2
    end if
    yrange(2) = yrange(2)*1.1

  end subroutine get_plot_vis_data

  !============================================================================

  function plot_triple_xval(i, j)

    !returns x-axis value for plot
    !
    !i is index into 1st axis of triple_data
    !If positive, j is index into 2nd axis of triple_data, otherwise:
    !j=-1: return longest projected baseline /mega-lambda
    !j=-2: return GMST /h

    !function arguments
    integer, intent(in) :: i, j
    double precision :: plot_triple_xval

    !local variables
    double precision :: u1, v1, u2, v2, lambda, mjd
    double precision, dimension(3) :: bas

    !functions
    double precision :: sla_gmst

    if (j>=1 .and. j<=size(triple_data, 2)) then
       plot_triple_xval = triple_data(i, j)
    else if (j == -1) then
       u1 = triple_data(i, 3)
       v1 = triple_data(i, 4)
       u2 = triple_data(i, 5)
       v2 = triple_data(i, 6)
       lambda = triple_data(i, 1)
       bas(1) = 1000.*sqrt(u1**2. + v1**2.)/lambda
       bas(2) = 1000.*sqrt(u2**2. + v2**2.)/lambda
       bas(3) = 1000.*sqrt((u1+u2)**2. + (v1+v2)**2.)/lambda
       plot_triple_xval = maxval(bas)
    else if (j == -2) then
       mjd = triple_data(i, 11)
       !note sla_gmsta would give better precision
       plot_triple_xval = rad2deg/15D0*sla_gmst(mjd) !convert to hours
    else
       stop 'Illegal value for j in plot_triple_xval()'
    end if

  end function plot_triple_xval

  !============================================================================

  subroutine get_plot_triple_data(xindex, xzero, get_amp, &
       data_points, num_data, flagged_points, num_flagged, xrange, yrange, &
       uxmin, uxmax)

    !If positive, xindex is index into 2nd axis of triple_data, otherwise:
    !-1: projected baseline length /mega-lambda
    !-2: GMST /h

    !subroutine arguments
    integer, intent(in) :: xindex
    real, intent(out) :: xzero
    logical, intent(in) :: get_amp
    double precision, intent(in), optional :: uxmin, uxmax
    real, dimension(:, :), intent(out) :: data_points, flagged_points
    integer, intent(out) :: num_data, num_flagged
    real, dimension(2), intent(out) :: xrange, yrange

    !local variables
    integer :: i
    real :: xmin, x, y, y_err, xspan

    !choose xzero
    xzero = 0.
    if(xindex >= 1 .and. xindex <= size(triple_data, 2)) then
       xmin = minval(triple_data(:, xindex))
       if(abs(xmin) > 100.*(maxval(triple_data(:, xindex)) - xmin)) &
            xzero = xmin
    end if

    !make up data arrays
    num_data = 0
    num_flagged = 0
    do i = 1, num_triple
       x = plot_triple_xval(i, xindex)
       !skip if x value outside plot range
       if (present(uxmin)) then
          if ((x-xzero) < uxmin) cycle
       end if
       if (present(uxmax)) then
          if ((x-xzero) > uxmax) cycle
       end if
       if (get_amp) then
          y = triple_data(i, 7)
          y_err = triple_data(i, 8)
       else
          y = modulo(triple_data(i, 9), 360D0)
          if (y > 180.) y = y - 360.
          y_err = triple_data(i, 10)
       end if
       if (y_err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = x - xzero
          flagged_points(num_flagged, 2) = y
          flagged_points(num_flagged, 3) = y - y_err
          flagged_points(num_flagged, 4) = y + y_err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = x - xzero
          data_points(num_data, 2) = y
          data_points(num_data, 3) = y + y_err
          data_points(num_data, 4) = y - y_err
       end if
    end do
    if(num_data + num_flagged == 0) return

    !calculate x range required
    if (present(uxmin)) then
       xrange(1) = uxmin
    else
       if (num_data > 0) then
          !exclude flagged extrema by default
          xrange(1) = minval(data_points(:num_data, 1))
       else
          xrange(1) = minval(flagged_points(:num_flagged, 1))
       end if
    end if
    if (present(uxmax)) then
       xrange(2) = uxmax
    else
       if (num_data > 0) then
          xrange(2) = maxval(data_points(:num_data, 1))
       else
          xrange(2) = maxval(flagged_points(:num_flagged, 1))
       end if
    end if
    xspan = xrange(2) - xrange(1)
    xrange(1) = xrange(1) - 0.05*xspan
    xrange(2) = xrange(2) + 0.05*xspan

    !calculate y range
    if (num_data > 0) then
       yrange(1) = minval(data_points(:num_data, 4))
       yrange(2) = maxval(data_points(:num_data, 3))
    else
       yrange(1) = minval(flagged_points(:num_flagged, 4))
       yrange(2) = maxval(flagged_points(:num_flagged, 3))
    end if

    !expand y range
    if(get_amp) then
       if (yrange(1) > 0) then
          yrange(1) = yrange(1)*0.9
       else if (yrange(1) < 0) then !yes, may be -ve
          yrange(1) = yrange(1)*1.2
       end if
       yrange(2) = yrange(2)*1.1
    else
       if (yrange(1) > 180.) then
          yrange(1) = 180.
       else if (yrange(1) < -180.) then
          yrange(1) = -180.
       end if
    end if

  end subroutine get_plot_triple_data

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

  subroutine plot_vis_bas(spec, param, mod_line, xlabel, y_title, top_title, &
       uxmin, uxmax, device)

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    logical, intent(in) :: mod_line !plot continuous line for model?
    character(len=*), intent(in) :: xlabel, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    character(len=128) :: x_title
    ! columns are x, y, y+delta, y-delta:
    real, dimension(num_vis, 4) :: data_points, flagged_points
    real, dimension(:, :), allocatable :: model_points
    integer :: num_data, num_flagged, num_model, i, istat
    double precision :: u, v, lambda, delta_lambda, mjd
    real, dimension(2) :: xrange, yrange
    real :: xzero, bas, mymin, mymax
    integer, parameter :: grid_size = 200

    !functions
    integer :: pgopen

    call get_plot_vis_data(-1, xzero, & 
         data_points, num_data, flagged_points, num_flagged, xrange, yrange, &
         uxmin, uxmax)
    if(num_data + num_flagged == 0) return
    if(abs(xzero) > 1E-6) then
       write (x_title, '(a, a, f9.1)') trim(xlabel), ' - ', xzero
    else
       x_title = xlabel
    end if

    if (mod_line) then
       !XXX can't have time-dependent model?
       !calculate grid of model points
       num_model = grid_size
       allocate(model_points(num_model, 2))
       lambda = vis_data(1, 1)
       delta_lambda = vis_data(1, 2)
       do i = 1, num_model
          u = xrange(1)*lambda/1000. &
               + ((i-1.)/(num_model-1.))*(xrange(2)-xrange(1))*lambda/1000.
          model_points(i, 1) = 1000.*u/lambda
          model_points(i, 2) = modulus( &
               cmplx_vis(spec, param, lambda, delta_lambda, u, 0D0, mjd))**2.
       end do
    else
       !calculate model points corresponding to plotted data points
       allocate(model_points(num_vis, 2))
       num_model = 0
       do i = 1, num_vis
          bas = plot_vis_xval(i, -1)
          if (bas >= xrange(1) .and. bas <= xrange(2)) then
             num_model = num_model + 1
             model_points(num_model, 1) = bas
             lambda = vis_data(i, 1)
             delta_lambda = vis_data(i, 2)
             u = vis_data(i, 3)
             v = vis_data(i, 4)
             mjd = vis_data(i, 7)
             model_points(num_model, 2) = modulus( &
                  cmplx_vis(spec, param, lambda, delta_lambda, u, v, mjd))**2.
          end if
       end do
    end if

    !expand y range so model points visible
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if(mymin < yrange(1)) yrange(1) = mymin
    if(mymax > yrange(2)) yrange(2) = mymax

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xrange(1), xrange(2), yrange(1), yrange(2), 0, 1)
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

    if (allocated(model_points)) deallocate(model_points)

  end subroutine plot_vis_bas

  !============================================================================

  subroutine plot_uv(x_title, y_title, top_title, device)

    !subroutine arguments
    character(len=*), intent(in) :: x_title, y_title, top_title
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(num_vis, 2) :: data_points, flagged_points
    real :: u, v, lambda, delta_lambda, bas, max
    integer :: num_data, num_flagged, i, istat

    !functions
    integer :: pgopen

    ! make up real arrays for unflagged & flagged data points
    num_data = 0
    num_flagged = 0
    max = 0.
    do i = 1, num_vis
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

  end subroutine plot_uv

  !============================================================================

  subroutine plot_vis(xindex, spec, param, &
       xlabel, y_title, top_title, uxmin, uxmax, device)

    !plot squared visibility against specified quanity
    !
    !xindex gives index into 2nd axis of vis_data for independent
    !variable, or indicates derived quantity :
    !-1: projected baseline length /mega-lambda
    !-2: GMST /h
    !-3: baseline position angle /deg

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: xlabel, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    character(len=128) :: x_title
    ! columns are x, y, y+delta, y-delta:
    real, dimension(num_vis, 4) :: data_points, flagged_points
    real, dimension(num_vis, 2) :: model_points
    double precision :: x, lambda, delta_lambda, u, v, mjd
    integer :: num_data, num_flagged, num_model, i, istat
    real, dimension(2) :: xrange, yrange
    real :: xzero, mymin, mymax

    !functions
    integer :: pgopen

    call get_plot_vis_data(xindex, xzero, & 
         data_points, num_data, flagged_points, num_flagged, xrange, yrange, &
         uxmin, uxmax)
    if(num_data + num_flagged == 0) return
    if(abs(xzero) > 1E-6) then
       write (x_title, '(a, a, f9.1)') trim(xlabel), ' - ', xzero
    else
       x_title = xlabel
    end if

    !calculate model points corresponding to plotted data points
    num_model = 0
    do i = 1, num_vis
       x = plot_vis_xval(i, xindex)
       if ((x-xzero) >= xrange(1) .and. (x-xzero) <= xrange(2)) then
          num_model = num_model + 1
          model_points(num_model, 1) = x - xzero
          lambda = vis_data(i, 1)
          delta_lambda = vis_data(i, 2)
          u = vis_data(i, 3)
          v = vis_data(i, 4)
          mjd = vis_data(i, 7)
          model_points(num_model, 2) = modulus( &
               cmplx_vis(spec, param, lambda, delta_lambda, u, v, mjd))**2.
       end if
    end do

    !expand y range so model points visible
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if(mymin < yrange(1)) yrange(1) = mymin
    if(mymax > yrange(2)) yrange(2) = mymax

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xrange(1), xrange(2), yrange(1), yrange(2), 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 13)

  end subroutine plot_vis

  !============================================================================

  subroutine plot_triple_phase(xindex, spec, param, xlabel, y_title, &
       top_title, uxmin, uxmax, device)

    !plot triple product phase against specified data column
    !xindex gives index into 2nd axis of triple_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: xlabel, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    character(len=128) :: x_title
    ! columns are x, y, y+delta, y-delta:
    real, dimension(num_triple, 4) :: data_points, flagged_points
    real, dimension(num_triple, 2) :: model_points
    integer :: num_data, num_flagged, num_model, i, istat
    double precision :: x, lambda, delta_lambda, u1, v1, u2, v2, mjd
    double complex :: vis1, vis2, vis3
    real, dimension(2) :: xrange, yrange
    real :: xzero, model_phase, mymin, mymax

    !functions
    integer :: pgopen

    call get_plot_triple_data(xindex, xzero, .false., & 
         data_points, num_data, flagged_points, num_flagged, xrange, yrange, &
         uxmin, uxmax)
    if(num_data + num_flagged == 0) return
    if(xzero > 1E-6) then
       write (x_title, '(a, a, f9.1)') trim(xlabel), ' - ', xzero
    else
       x_title = xlabel
    end if

    !calculate model points corresponding to plotted data points
    num_model = 0
    do i = 1, num_triple
       x = plot_triple_xval(i, xindex)
       if ((x-xzero) >= xrange(1) .and. (x-xzero) <= xrange(2)) then
          num_model = num_model + 1
          model_points(num_model, 1) = x - xzero
          lambda = triple_data(i, 1)
          delta_lambda = triple_data(i, 2)
          u1 = triple_data(i, 3)
          v1 = triple_data(i, 4)
          u2 = triple_data(i, 5)
          v2 = triple_data(i, 6)
          mjd = triple_data(i, 11)
          vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
          vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
          vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
               -(u1+u2), -(v1+v2), mjd)
          model_phase = modulo(rad2deg*argument(vis1*vis2*vis3), 360D0)
          if (model_phase > 180.) model_phase = model_phase - 360.
          model_points(num_model, 2) = model_phase
       end if
    end do

    !expand y range so model points visible
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if(mymin < yrange(1)) yrange(1) = mymin
    if(mymax > yrange(2)) yrange(2) = mymax

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xrange(1), xrange(2), yrange(1), yrange(2), 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 7)

  end subroutine plot_triple_phase

  !============================================================================

  subroutine plot_triple_amp(xindex, spec, param, xlabel, y_title, &
       top_title, uxmin, uxmax, device)

    !plot triple product phase against specified data column
    !xindex gives index into 2nd axis of triple_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: xlabel, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    character(len=128) :: x_title
    ! columns are x, y, y+delta, y-delta:
    real, dimension(num_triple, 4) :: data_points, flagged_points
    real, dimension(num_triple, 2) :: model_points
    integer :: num_data, num_flagged, num_model, i, istat
    double precision :: x, lambda, delta_lambda, u1, v1, u2, v2, mjd
    double complex :: vis1, vis2, vis3
    real, dimension(2) :: xrange, yrange
    real :: xzero, mymin, mymax

    !functions
    integer :: pgopen

    call get_plot_triple_data(xindex, xzero, .true., & 
         data_points, num_data, flagged_points, num_flagged, xrange, yrange, &
         uxmin, uxmax)
    if(num_data + num_flagged == 0) return
    if(xzero > 1E-6) then
       write (x_title, '(a, a, f9.1)') trim(xlabel), ' - ', xzero
    else
       x_title = xlabel
    end if

    !calculate model points corresponding to plotted data points
    num_model = 0
    do i = 1, num_triple
       x = plot_triple_xval(i, xindex)
       if ((x-xzero) >= xrange(1) .and. (x-xzero) <= xrange(2)) then
          num_model = num_model + 1
          model_points(num_model, 1) = x - xzero
          lambda = triple_data(i, 1)
          delta_lambda = triple_data(i, 2)
          u1 = triple_data(i, 3)
          v1 = triple_data(i, 4)
          u2 = triple_data(i, 5)
          v2 = triple_data(i, 6)
          mjd = triple_data(i, 11)
          vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
          vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
          vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
               -(u1+u2), -(v1+v2), mjd)
          model_points(num_model, 2) = modulus(vis1*vis2*vis3)
       end if
    end do

    !expand y range so model points visible
    mymin = minval(model_points(:num_model, 2))
    mymax = maxval(model_points(:num_model, 2))
    if(mymin < yrange(1)) yrange(1) = mymin
    if(mymax > yrange(2)) yrange(2) = mymax

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xrange(1), xrange(2), yrange(1), yrange(2), 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(2)
    call plot_data_pts(x_title, y_title, num_flagged, flagged_points, &
         9, .false.)
    call pgsci(1)
    call plot_data_pts(x_title, y_title, num_data, data_points, 2, .true.)
    call pgsci(3)
    call plot_model_pts(x_title, y_title, num_model, model_points, 7)

  end subroutine plot_triple_amp

  !============================================================================

end module Plot
