module Plot

!subroutines contained
!plot_vis_data - nominal plot of visibility as a function of physical 
!baseline (u^2 + v^2) [i.e. not rho, the rotated/stretched baseline]

implicit none

contains

!==============================================================================

subroutine plot_vis_data(info, vis_data, model_vis_data, lambda, x_title, &
                         y_title, top_title)

  !subroutine arguments
  character(len=128) :: info, x_title, y_title, top_title
  double precision, dimension(:,:), allocatable :: vis_data
  double precision, dimension(:,:), allocatable :: model_vis_data
  double precision :: lambda
  
  !local variables
  real, dimension(:,:), allocatable :: plot_data, plot_data2
  double precision :: u, v
  integer :: num_points, num_points2, symbol, status
  logical :: errors, points, lines, newplot
  
  !plotting variables
  real :: xmin, xmax, ymin, ymax
  integer :: i, num_points !loop
  
  !calculate x-points (real physical baseline)
  num_points = size(vis_data,1)
  num_points2 = size(model_vis_data,1)
  allocate(plot_data(num_points,2))
  allocate(plot_data2(num_points2,2))
  do i = 1, num_points
     u = vis_data(i,2)
     v = vis_data(i,3)
     plot_data(i,1) = real((sqrt((u**2)+(v**2)))/(lambda*1D-3))
  end do
  do i = 1, num_points2
     u = model_vis_data(i,2)
     v = model_vis_data(i,3)
     plot_data2(i,1) = real((sqrt((u**2)+(v**2)))/(lambda*1D-3))
  end do
  
  !calculate y-points (visiblity)
  plot_data(:,2) = vis_data(:,4)
  plot_data2(:,2) = model_vis_data(:,4)
  
  !calculate ranges required
  xmin = minval(plot_data(:,1))
  xmax = maxval(plot_data(:,1))
  ymin = minval(plot_data(:,2))
  ymax = maxval(plot_data(:,2))
  
  !change ranges
  xmin = 0.0
  xmax = xmax*1.1
  if (ymin > 0) then
     ymin = 0.0
  else if (ymin < 0) then
     ymin = ymin*1.2
  end if
  if (ymax < 1.1) then
     ymax = 1.1
  end if
  
  !initialise plot
  call pgbegin(0, '/xserve', 1, 1)
  
  !set ranges
  call pgenv(xmin, xmax, ymin, ymax, 0, 1)
  
  !label plot
  call pglabel(trim(x_title), trim(y_title), trim(top_title))
  
  !plot data
  symbol = 2
  call pgpoint(num_points, plot_data(:,1), plot_data(:,2), symbol)
  call pgline(1000, plot_data2(:,1), plot_data2(:,2))
  
  !end plotting
  call pgend
  
  return
  
end subroutine plot_vis_data

!==============================================================================

end module Plot
