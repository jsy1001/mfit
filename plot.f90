module Plot
  
  use Model
  use Fit

  !subroutines contained:
  !
  !plot_vis - plot squared visibility against projected baseline
  !
  !plot_triple - plot closure phase against longest projected baseline
  !in triangle
  !
  !plot_uv - plot uv coverage

  implicit none

contains

  !==============================================================================

  subroutine plot_triple(spec, param, x_title, y_title, top_title, device)

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    character(len=*), intent(in) :: x_title, y_title, top_title
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, model_points
    real :: xmin, xmax, ymin, ymax
    double precision :: lambda, u1, v1, u2, v2
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
    num_model = size(triple_data, 1)
    allocate(model_points(num_model, 2))
    num_data = 0
    num_flagged = 0
    do i = 1, size(triple_data, 1)
       lambda = triple_data(i, 1)
       u1 = triple_data(i, 2)
       v1 = triple_data(i, 3)
       u2 = triple_data(i, 4)
       v2 = triple_data(i, 5)
       phase = modulo(triple_data(i, 8), 360D0)
       if (phase > 180.) phase = phase - 360.
       phase_err = triple_data(i, 9)
       bas(1) = 1000.*sqrt(u1**2. + v1**2.)/lambda
       bas(2) = 1000.*sqrt(u2**2. + v2**2.)/lambda
       bas(3) = 1000.*sqrt((u1+u2)**2. + (v1+v2)**2.)/lambda
       vis1 = cmplx_vis(spec, param, lambda, u1, v1)
       vis2 = cmplx_vis(spec, param, lambda, u2, v2)
       vis3 = cmplx_vis(spec, param, lambda, -(u1+u2), -(v1+v2))
       model_points(i, 1) = maxval(bas)
       model_phase = modulo(rad2deg*argument(vis1*vis2*vis3), 360D0)
       if (model_phase > 180.) model_phase = model_phase - 360.
       model_points(i, 2) = model_phase
       if (phase_err < 0D0) then
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

    !calculate ranges required
    xmin = minval(data_points(:num_data, 1))
    xmax = maxval(data_points(:num_data, 1))
    ymin = minval(data_points(:num_data, 4))
    ymax = maxval(data_points(:num_data, 3))
    
    !change ranges
    xmin = xmin*0.9
    xmax = xmax*1.1
    if (ymax > 180.) then
       ymax = 180.
    else if (ymax < -180.) then
       ymax = -180.
    end if

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    endif
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(2)
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 2)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(3)
    call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)

  end subroutine plot_triple

  !==============================================================================

  subroutine plot_vis(spec, param, symm, x_title, y_title, top_title, device)

    !subroutine arguments
    character(len=128), dimension(:,:), intent(in) :: spec
    double precision, dimension(:,:), intent(in) :: param
    logical, intent(in) :: symm
    character(len=*), intent(in) :: x_title, y_title, top_title
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points, model_points
    double precision :: u, v, lambda
    integer :: num_data, num_flagged, num_model, i, istat
    real :: bas, vsq, err, xmin, xmax, ymin, ymax
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
       u = vis_data(i, 2)
       v = vis_data(i, 3)
       bas = 1000.*sqrt(u**2. + v**2.)/lambda
       vsq = vis_data(i, 4)
       err = vis_data(i, 5)
       if (err < 0D0) then
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

    !calculate ranges required
    xmin = minval(data_points(:num_data, 1))
    xmax = maxval(data_points(:num_data, 1))
    ymin = minval(data_points(:num_data, 4))
    ymax = maxval(data_points(:num_data, 3))
    
    !change ranges
    xmin = 0.0
    xmax = xmax*1.1
    if (ymin > 0) then
       ymin = ymin*0.9
    else if (ymin < 0) then
       ymin = ymin*1.2
    end if
    if (ymax < 1.1) then
       ymax = 1.1
    end if

    if (symm) then !and single wavelength
       !calculate grid of model points
       num_model = grid_size
       allocate(model_points(num_model, 2))
       lambda = vis_data(1, 1)
       do i = 1, num_model
          u = xmin + ((i-1.)/(num_model-1.))*(xmax-xmin)
          model_points(i, 1) = 1000.*u/lambda
          model_points(i, 2) = modulus(cmplx_vis(spec, param, lambda, u, 0D0))**2.
       end do
    else
       !calculate model points corresponding to data points
       num_model = size(vis_data, 1)
       allocate(model_points(num_model, 2))
       do i = 1, num_model
          lambda = vis_data(i, 1)
          u = vis_data(i, 2)
          v = vis_data(i, 3)
          model_points(i, 1) = 1000.*sqrt(u**2. + v**2.)/lambda
          model_points(i, 2) = modulus(cmplx_vis(spec, param, lambda, u, v))**2.
       end do
    endif

    !plot data
    if (present(device)) then
       istat = pgopen(device)
       if (istat <= 0) return
    endif
    call pgsci(1)
    call pgenv(xmin, xmax, ymin, ymax, 0, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 2)
    call pgerry(num_data, data_points(:, 1), &
         data_points(:, 3), data_points(:, 4), 1.0)
    call pgsci(2)
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 2)
    call pgerry(num_flagged, flagged_points(:, 1), &
         flagged_points(:, 3), flagged_points(:, 4), 1.0)
    call pgsci(3)
    if (symm) then
       call pgline(num_model, model_points(:, 1), model_points(:, 2))
    else
       call pgpt(num_model, model_points(:, 1), model_points(:, 2), 7)
    end if

  end subroutine plot_vis

  !==============================================================================

  subroutine plot_uv(x_title, y_title, top_title, device)

    !subroutine arguments
    character(len=*), intent(in) :: x_title, y_title, top_title
    character(len=*), intent(in), optional :: device

    !local variables
    real, dimension(:, :), allocatable :: data_points, flagged_points
    real :: u, v, lambda, bas, max
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
       u = 1000.*vis_data(i, 2)/lambda
       v = 1000.*vis_data(i, 3)/lambda
       bas = sqrt(u**2. + v**2.)
       if (bas > max) max = bas
       if (vis_data(i, 5) < 0D0) then
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
    endif
    call pgsci(1)
    call pgenv(-max, max, -max, max, 1, 1)
    call pglab(trim(x_title), trim(y_title), trim(top_title))
    call pgpt(num_data, data_points(:, 1), data_points(:, 2), 17)
    call pgpt(num_data, -data_points(:, 1), -data_points(:, 2), 22)
    call pgsci(2)
    call pgpt(num_flagged, flagged_points(:, 1), flagged_points(:, 2), 17)
    call pgpt(num_flagged, -flagged_points(:, 1), -flagged_points(:, 2), 22)

  end subroutine plot_uv

  !==============================================================================

end module Plot
