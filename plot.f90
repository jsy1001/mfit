!$Id: plot.f90,v 1.13 2005/05/24 12:53:06 jsy1001 Exp $

module Plot
  
  use Model
  use Fit

  !subroutines contained:
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
  !plot_post - plot 1d cut through -ln(posterior)
  !
  !plot_vis - plot squared visibility against specified data column
  !
  !plot_triple_phase - plot closure phase against specified data column
  !
  !plot_triple_amp - plot triple product amplitude specified data column

  implicit none

contains

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
    double precision :: lambda, delta_lambda, u1, v1, u2, v2
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
       bas(1) = 1000.*sqrt(u1**2. + v1**2.)/lambda
       bas(2) = 1000.*sqrt(u2**2. + v2**2.)/lambda
       bas(3) = 1000.*sqrt((u1+u2)**2. + (v1+v2)**2.)/lambda
       if (present(uxmin) .and. maxval(bas) .lt. uxmin) cycle!out of plot range
       if (present(uxmax) .and. maxval(bas) .gt. uxmax) cycle!out of plot range
       phase = modulo(triple_data(i, 9), 360D0)
       if (phase > 180.) phase = phase - 360.
       phase_err = triple_data(i, 10)
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, -(u1+u2), -(v1+v2))
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
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 9)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(1)
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(3)
    !!call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)
    call pgpt(num_model, model_points(:, 1), model_points(:, 2), 13)

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
    double precision :: lambda, delta_lambda, u1, v1, u2, v2
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
       bas(1) = 1000.*sqrt(u1**2. + v1**2.)/lambda
       bas(2) = 1000.*sqrt(u2**2. + v2**2.)/lambda
       bas(3) = 1000.*sqrt((u1+u2)**2. + (v1+v2)**2.)/lambda
       if (present(uxmin) .and. maxval(bas) .lt. uxmin) cycle!out of plot range
       if (present(uxmax) .and. maxval(bas) .gt. uxmax) cycle!out of plot range
       amp = triple_data(i, 7)
       amp_err = triple_data(i, 8)
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, -(u1+u2), -(v1+v2))
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
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 9)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(1)
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(3)
    call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)

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
    double precision :: u, v, lambda, delta_lambda
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
       !calculate grid of model points
       num_model = grid_size
       allocate(model_points(num_model, 2))
       lambda = vis_data(1, 1)
       delta_lambda = vis_data(1, 2)
       do i = 1, num_model
          u = xmin*lambda/1000. + ((i-1.)/(num_model-1.))*(xmax-xmin)*lambda/1000.
          model_points(i, 1) = 1000.*u/lambda
          model_points(i, 2) = modulus(cmplx_vis(spec, param, lambda, delta_lambda, u, 0D0))**2.
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
          bas = 1000.*sqrt(u**2. + v**2.)/lambda
          if (bas .ge. xmin .and. bas .le. xmax) then
             num_model = num_model + 1
             model_points(num_model, 1) = bas
             model_points(num_model, 2) = modulus(cmplx_vis(spec, param, lambda, delta_lambda, u, v))**2.
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
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 9)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(1)
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(3)
    if (mod_line) then
       call pgline(num_model, model_points(:, 1), model_points(:, 2))
    else
       !!call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)
       call pgpt(num_model, model_points(:, 1), model_points(:, 2), 13)
    end if

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

  end subroutine plot_uv

  !============================================================================

  subroutine plot_post(spec, param, index, x_title, y_title, top_title, &
       uxmin, uxmax, device)

    !index is into x_pos array of variable parameters
    !x_pos, x_info must be initialised on entry to this routine i.e. must have
    !called minimiser() from Fit module

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    integer, intent(in) :: index
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in) :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    double precision, dimension(:), allocatable :: var_param
    double precision, dimension(:, :), allocatable :: all_param
    real, dimension(:, :), allocatable :: post_points
    double precision :: val, lhd, pri
    integer :: nvar, num_points, i, istat
    real :: xmin, xmax, ymin, ymax
    integer, parameter :: grid_size = 200

    !functions
    integer :: pgopen

    !return if nothing to plot
    nvar = size(x_pos, 1)
    if (index < 1 .or. index > nvar) return
    if (uxmax < x_info(index, 2)) return
    if (uxmin > x_info(index, 3)) return

    !adjust x range if necessary - need to keep within x_info limits
    if (uxmin < x_info(index, 2)) then
       xmin = x_info(index, 2)
    else
       xmin = uxmin
    end if
    if (uxmax > x_info(index, 3)) then
       xmax = x_info(index, 3)
    else
       xmax = uxmax
    end if

    !copy model parameters
    allocate(all_param(size(param, 1), 17))
    all_param = param
    allocate(var_param(nvar))
    do i = 1, nvar
       var_param(i) = param(x_pos(i, 1), x_pos(i, 2))
    end do

    !calculate grid of points
    num_points = grid_size
    allocate(post_points(num_points, 2))
    do i = 1, num_points
       val = xmin + ((i-1.)/(num_points-1.))*(xmax-xmin)
       all_param(x_pos(index, 1), x_pos(index, 2)) = val
       var_param(index) = val
       lhd = likelihood(vis_data, triple_data, spec, all_param)
       pri = prior(var_param, x_pos, all_param, model_prior)
       post_points(i, 1) = val
       post_points(i, 2) = lhd + pri
    end do

    !calculate y range
    ymin = minval(post_points(:, 2))
    ymax = maxval(post_points(:, 2))

    !plot posterior
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    end if
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgsci(1)
    call pgline(num_points, post_points(:, 1), post_points(:, 2))

  end subroutine plot_post

  !============================================================================

  subroutine plot_vis(xindex, spec, param, x_title, y_title, top_title, &
       uxmin, uxmax, device)

    !plot squared visibility against specified data column
    !xindex gives index into 2nd axis of vis_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, &
         model_points
    double precision :: u, v, lambda, delta_lambda
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
          flagged_points(num_flagged, 1) = vis_data(i, xindex)
          flagged_points(num_flagged, 2) = vsq
          flagged_points(num_flagged, 3) = vsq - err
          flagged_points(num_flagged, 4) = vsq + err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = vis_data(i, xindex)
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
       if (vis_data(i, xindex) .ge. xmin .and. vis_data(i, xindex) .le. xmax) then
          num_model = num_model + 1
          model_points(num_model, 1) = vis_data(i, xindex)
          model_points(num_model, 2) = modulus(cmplx_vis(spec, param, lambda, delta_lambda, u, v))**2.
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
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 9)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(1)
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(3)
    call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)

  end subroutine plot_vis

  !============================================================================

  subroutine plot_triple_phase(xindex, spec, param, x_title, y_title, &
       top_title, uxmin, uxmax, device)

    !plot triple product phase against specified data column
    !xindex gives index into 2nd axis of triple_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, &
         model_points
    real :: xmin, xmax, ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax
    double precision :: lambda, delta_lambda, u1, v1, u2, v2
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
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, -(u1+u2), -(v1+v2))
       num_model = num_model + 1
       model_points(num_model, 1) = triple_data(i, xindex)
       model_phase = modulo(rad2deg*argument(vis1*vis2*vis3), 360D0)
       if (model_phase > 180.) model_phase = model_phase - 360.
       model_points(num_model, 2) = model_phase
       if (phase_err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = triple_data(i, xindex)
          flagged_points(num_flagged, 2) = phase
          flagged_points(num_flagged, 3) = phase - phase_err
          flagged_points(num_flagged, 4) = phase + phase_err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = triple_data(i, xindex)
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
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 9)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(1)
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(3)
    call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)

  end subroutine plot_triple_phase

  !============================================================================

  subroutine plot_triple_amp(xindex, spec, param, x_title, y_title, &
       top_title, uxmin, uxmax, device)

    !plot triple product phase against specified data column
    !xindex gives index into 2nd axis of triple_data for independent variable

    !subroutine arguments
    integer :: xindex
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    double precision, intent(in), optional :: uxmin, uxmax
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, &
         model_points
    real :: xmin, xmax, ymin, ymax, dymin, dymax, fymin, fymax, mymin, mymax
    double precision :: lambda, delta_lambda, u1, v1, u2, v2
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
       vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1)
       vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2)
       vis3 = cmplx_vis(spec, param, lambda, delta_lambda, -(u1+u2), -(v1+v2))
       num_model = num_model + 1
       model_points(num_model, 1) = triple_data(i, xindex)
       model_points(num_model, 2) = modulus(vis1*vis2*vis3)
       if (amp_err <= 0D0) then
          num_flagged = num_flagged + 1
          flagged_points(num_flagged, 1) = triple_data(i, xindex)
          flagged_points(num_flagged, 2) = amp
          flagged_points(num_flagged, 3) = amp - amp_err
          flagged_points(num_flagged, 4) = amp + amp_err
       else
          num_data = num_data + 1
          data_points(num_data, 1) = triple_data(i, xindex)
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
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 9)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(1)
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(3)
    call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)

  end subroutine plot_triple_amp

  !============================================================================

end module Plot
