module Inout

!subroutines contained
!
!read_vis
!read_mapdat
!read_model

implicit none

contains

!==============================================================================

subroutine read_vis(info, file_name, source, max_lines, vis_data, lambda, &
                    calib_error)
    
  !Reads vis file with up to max_lines (excluding blanks)
  !Fills in vis_data array. Note projected baseline sqrt(u**2 + v**2) gets
  !put in u column, with zero in v column
 
  !subroutine arguments
  character(len=128) :: info, file_name, source
  integer :: max_lines
  double precision, dimension(:,:), allocatable :: vis_data
  double precision :: lambda, calib_error

  !local variables
  character(len=32) :: dummy
  integer :: i, data_items
  double precision :: vis, baseline, default_error

  !default absolute error on the visibility^2 data points
  default_error = 0.001D0
  
  !check for zero length filename
  if (file_name == '') goto 90    
  
  !read file header
  open (unit=11, err=91, status='old', action='read', file=file_name)
  read (11, *, err=92, end=92) (source, data_items)
  source = trim(source(1:index(source,'~')-1))
  close (11)
  
  !check valid number of lines
  open (unit=11, action='read', file=file_name)
  do i = 1, max_lines+1
     read (11, *, err=94, end=1) dummy
  end do
  goto 93
1 close (11)
  
  !allocate array size
  allocate(vis_data(data_items,5))
  
  !read vis data properly and close
  !vis_data: lambda, u, v, cal_vis, err
  !where cal_vis is the squared visibility (clearly always positive)
  open (unit=11, action='read', file=file_name)
  read (11, *, err=94) (dummy, dummy) !skip header
  do i = 1, data_items
     read (11, *, err=94) (vis, baseline)
       vis_data(i,1) = lambda
       vis_data(i,2) = abs(baseline)/1D+3
       vis_data(i,3) = 0D0
       vis_data(i,4) = vis**2D0
       vis_data(i,5) = vis_data(i,4) * sqrt( calib_error**2D0 + &
            (default_error/vis_data(i,4))**2D0 )
  end do

  !clean-up and return
200 continue
    
  return

  !error trapping 
90 info = 'blank filename'
  goto 200
91 info = 'cannot open file'
  goto 200
92 info = 'cannot read from file'
  close (11)
  goto 200
93 info = 'file exceeds maximum permitted length'
  close (11)
  goto 200
94 info = 'unknown read problem - possible invalid file format'
  close (11)
  goto 200
  
  end subroutine read_vis

!==============================================================================

subroutine read_nvis(info, file_name, source, max_lines, vis_data, lambda, &
                    calib_error)
    
  !Reads nvis file with up to max_lines (excluding comments)
  !Fills in vis_data array. Note projected baseline sqrt(u**2 + v**2) gets
  !put in u column, with zero in v column
 
  !subroutine arguments
  character(len=128), intent(in) :: file_name
  character(len=128), intent(out) :: info, source
  integer, intent(in) :: max_lines
  double precision, dimension(:,:), allocatable, intent(out) :: vis_data
  double precision, intent(in) :: lambda, calib_error

  !local variables
  character(len=256) :: line
  integer :: i, data_items
  double precision :: vis, err, baseline
  
  !check for zero length filename
  if (file_name == '') then
     info = 'blank filename'
     return
  end if

  !count number of data lines
  data_items = 0
  open (unit=11, action='read', file=file_name)
  do i = 1, max_lines+1
     read (11, '(a)', err=94, end=1) line
     if (line(1:1) /= '#') then
        data_items = data_items + 1
     else if (i == 1) then
        source = trim(line(2:)) !by convention, 1st comment is source name
     end if
  end do
  info = 'file exceeds maximum permitted length'
  close (11)
  return

1 close (11)
  
  !allocate array size
  allocate(vis_data(data_items,5))
  
  !read vis data properly and close
  !vis_data: lambda, u, v, cal_vis, err
  !where cal_vis is the squared visibility
  i = 1
  open (unit=11, action='read', file=file_name)
  do
     read (11, '(a)', err=92, end=2) line
     if (line(1:1) /= '#') then
        read (line, *, err=94) (vis, err, baseline)
        vis_data(i,1) = lambda
        vis_data(i,2) = abs(baseline)/1D+3
        vis_data(i,3) = 0D0
        vis_data(i,4) = vis**2D0
        !calculate abs error on squared visibility; err is abs error
        !on (assumed real) visibility
        vis_data(i,5) = vis_data(i,4) * sqrt( calib_error**2D0 + &
             (2D0*err/vis_data(i,4))**2D0 )
        i = i + 1
        if (i > max_lines) then
           info = 'file exceeds maximum permitted length'
           close (11)
           return
        end if
    end if
  end do

2 close(11)
  return

  !error trapping 
92 info = 'cannot read from file'
  close (11)
  return
94 info = 'unknown read problem - possible invalid file format'
  close (11)
  return
  
  end subroutine read_nvis

!==============================================================================

subroutine read_mapdat(info, file_name, max_lines, vis_data, triple_data, &
                       wavebands, calib_error)

  !Reads mapdat file with up to max_lines (excluding blanks)
  !
  !Fills in vis_data and triple_data arrays:
  !
  !vis_data: lambda, u, v, vis, err
  !          vis is squared visibility amplitude (may be -ve for data points)
  !triple_data: lambda, u1, v1, u2, v2, amp, err, cp, err
  !
  !wavebands array is list of different wavelength values encountered
  
  !subroutine arguments
  character(len=128) :: info, file_name
  integer :: max_lines
  double precision :: calib_error
  double precision, dimension(:,:), allocatable :: vis_data
  double precision, dimension(:,:), allocatable :: triple_data
  double precision, dimension(:), allocatable :: wavebands

  !local variables
  character(len=32) :: dummy
  integer :: i, j, i1, i2, lines, num
  double precision :: vis, vis_err, amp, amp_err, cp, cp_err, swap
  double precision, dimension(:), allocatable :: lambdas
  
  !check for zero length filename
  if (file_name == '') goto 90
  
  !read and count vis/triple occurences and number of lines
  open (unit=11, err=91, status='old', action='read', file=file_name)
  i1 = 0
  i2 = 0
  do i = 1, max_lines+1
     read (11, *, err=92, end=1) dummy
     if (dummy == 'vis') then 
        i1 = i1 + 1
     else if (dummy == 'triple') then 
        i2 = i2 + 1
     end if
  end do
  goto 93
1 lines = i-1
  close (11)
 
  if ((i1 == 0) .and. (i2 == 0)) goto 94
 
  !allocate vis and triple data array sizes
  allocate(vis_data(i1,5))
  allocate(triple_data(i2,9))
  
  !read data properly and close
  open (unit=11, action='read', file=file_name)
  open (unit=12, action='read', file=file_name)
  i1 = 0 !counters for vis and triple data item
  i2 = 0
  do i = 1, lines
     read (11, *) dummy

     if (dummy == 'vis') then

        i1 = i1 + 1
        read (12,*,err=95) (dummy, dummy, dummy, vis_data(i1,1), dummy, &
             dummy, dummy, dummy, vis_data(i1,2:3), vis, vis_err)

        !V^2 positive stored as sqroot(V^2) in file
        !V^2 negative stored as -sqroot(-V^2) in file
        vis_data(i1,4) = vis**2D0
        if (vis<0D0) vis_data(i1,4) = -vis_data(i1,4)

        !"vis" error in file is the modulus of the fractional error
        !(in the squared visibility) divided by 2. Need to convert to
        !absolute error
        vis_data(i1,5) = vis_data(i1,4) * sqrt(calib_error**2D0 + &
                         (2D0*vis_err)**2D0 )

        if (vis_err<0D0) vis_data(i1,5) = -vis_data(i1,5)

     else if (dummy == 'triple') then

        i2 = i2 + 1
        read (12,*,err=95) (dummy, dummy, dummy, dummy, triple_data(i2,1), &
             dummy, dummy, dummy, dummy, triple_data(i2,2:5), &
             amp, amp_err, cp, cp_err)

        !amplitude in file is cube root of amplitude
        triple_data(i2,6) = amp**3D0

        !error in amplitude is 1/3 of fractional error in the amplitude
        triple_data(i2,7) = abs(3D0*amp_err*triple_data(i2,6))
        !preserve sign of tp error - this is flagging info
        if (amp_err<0D0) triple_data(i2,7) = -triple_data(i2,7)

        triple_data(i2,8) = cp
        triple_data(i2,9) = cp_err

     else
        read (12,*,err=95) dummy
     end if
  end do
  close (11)
  close (12)

  !make array of wavebands detected:
  !count number of different wavelengths and create vector lambdas
  !containing only one occurence of each wavelength
  allocate(lambdas(size(vis_data,1)))
  lambdas = vis_data(:,1)
  num = 0
  do i = 1, size(lambdas,1)
     if (lambdas(i)/=-1) then
        num = num + 1
        do j = i+1, size(lambdas,1)
           if (lambdas(j)==lambdas(i)) lambdas(j) = -1
        end do
     end if
  end do
  !make list of wavelengths found
  allocate(wavebands(num))
  j = 0
  do i = 1, size(lambdas,1)
     if (lambdas(i)/=-1) then
        j = j + 1
        wavebands(j) = lambdas(i)
     end if
  end do
  !sort into ascending order
  do i = num, 2, -1
     do j = 1, i-1
        if (wavebands(j)>wavebands(j+1)) then
           swap = wavebands(j) 
           wavebands(j) = wavebands(j+1)
           wavebands(j+1) = swap
        end if
     end do
  end do
  
  !clean-up and return
  if (allocated(lambdas)) deallocate(lambdas)
200 continue

  return
  
  !error trapping 
90 info = 'blank filename'
  goto 200
91 info = 'cannot open file'
  goto 200
92 info = 'cannot read from file'
  close (11)
  goto 200
93 info = 'file exceeds maximum permitted length'
  close (11)
  goto 200
94 info = 'no vis or triple data found'
  goto 200
95 info = 'unknown read problem - possible invalid file format'
  close (11)
  goto 200
  
end subroutine read_mapdat

!==============================================================================

subroutine read_model(info, file_name, source, max_lines, &
                      model_spec, model_param, model_prior, limits)
    
  !Reads model files and puts data into model_spec and model_param 
  !
  !model_spec holds 3 items per component: 1 name, 2 shape and 3 LD type
  !model_param holds 17 items: 1 LD order, 2-7 r, theta, brightness, major 
  !axis, phi, epsilon, 8-17 LD coeffs 1-10.
  !model_param holds priors for all 17 items in model_param (LD order prior
  !forced to zero to prevent order variation).
  
  !subroutine arguments
  character(len=128) :: info, file_name, source
  integer :: max_lines 
  character(len=128), dimension(:,:), allocatable :: model_spec
  double precision, dimension(:,:), allocatable :: model_param
  double precision, dimension(:,:), allocatable :: model_prior
  double precision, dimension(17,2) :: limits
  
  !local variables
  character(len=32) :: dummy, source1, source2
  integer :: i, j, k, lines, comps, order
  
  !check for zero length filename
  if (file_name == '') goto 90    
  
  !read and count number of components and number of lines
  open (unit=11, err=91, status='old', action='read', file=file_name)
  comps = 0
  do i = 1, max_lines+1
     read (11, *, err=92, end=1) dummy
     if (dummy == 'component') comps = comps+1
  end do
  goto 93
1 lines = i-1
  close (11)
  
  !check if legal number of components
  if ((comps<1).or.(comps>10)) goto 97

  !allocate model data arrays
  allocate(model_spec(comps,3))
  allocate(model_param(comps,17))
  allocate(model_prior(comps,17))
  model_param = 0D0
  model_prior = 0D0

  !Limits array stores the lower and upper acceptable limits
  !on the parameter values, it is passed into the minimising routines
  !later to ensure parameter values stay legal. Model values not within
  !these limits trigger error on reading the model
  limits(1:7,1) = dble((/0,0,-720,-100,0,-720,0/))
  limits(8:17,1) = -100D0
  limits(1:7,2) = dble((/10,100,720,100,1000,720,10/))
  limits(8:17,2) = 100D0

  !read source name
  open (unit=11, action='read', file=file_name)
  open (unit=12, action='read', file=file_name)
  do i = 1, lines
     read (11, *, err=94) dummy
     if (dummy == 'source') then
        read (12, *, err=94) (dummy, source1, source2)
        source = trim(source1) // ' ' // trim(source2)
     else
        read (12, *, err=94)
     end if
  end do
  close (11)
  close (12)
  
  !read spec and param info for components
  open (unit=11, action='read', file=file_name)
  j = 0 !component counter
  do i = 1, lines
     read (11,*,err=94,end=2) dummy
     if (dummy == 'component') then
        j = j+1
        
        !read component name
        read (11,*,err=95,end=95) (dummy, model_spec(j,1))
        if (dummy /= 'name') goto 96
        
        !read shape type
        read (11,*,err=95,end=95) (dummy, model_spec(j,2))
        if (dummy /= 'shape_type') goto 96
        
        !read LD type (set as uniform for point case)
        select case (trim(model_spec(j,2)))
           case ('point')
              read (11,*,err=95) dummy
              model_spec(j,3) = 'uniform'
           case default
              read (11,*,err=95,end=95) (dummy, model_spec(j,3))
        end select
        if (dummy /= 'ld_type') goto 96
        
        !Read LD order (preset for some LD type cases)
        !Non-integer will cause error in read statement
        select case (trim(model_spec(j,3)))
           case ('uniform','gaussian')
              read (11,*,err=95,end=95) dummy
              order = 0
           case ('square-root')
              read (11,*,err=95,end=95) dummy
              order = 2
           case ('hestroffer')
              read (11,*,err=95,end=95) dummy
              order = 1
           case ('taylor','gauss-hermite') 
              read (11,*,err=95,end=95) (dummy, order)
           case default
              goto 98
           end select
        if (dummy /= 'ld_order') goto 96
        model_param(j,1) = dble(order)
        
        !read position r and theta
        read (11,*,err=95,end=95) (dummy, model_param(j,2:3))
        if (dummy /= 'position') goto 96
        
        !read position prior, must be non-negative
        read (11,*,err=95,end=95) (dummy, model_prior(j,2:3))
        if (dummy /= 'position_prior') goto 96
        
        !read flux
        read (11,*,err=95,end=95) (dummy, model_param(j,4))
        if (dummy /= 'flux') goto 96
        
        !read flux prior, must be non-negative
        read (11,*,err=95,end=95) (dummy, model_prior(j,4))
        if (dummy /= 'flux_prior') goto 96  
        
        !Read shape parameters a, phi, epsilon depending on
        !shape type. Also read priors
        !Shape a and epsilon must be non-negative
        !All priors must be non-negative
        model_param(j,5) = 0D0
        model_param(j,6) = 0D0
        model_param(j,7) = 1D0
        select case (trim(model_spec(j,2)))
           case ('point')
              read (11,*,err=95,end=95) dummy
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) dummy
           case ('disc')
              read (11,*,err=95,end=95) (dummy, model_param(j,5))
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) (dummy, model_prior(j,5))
           case ('ellipse')
              read (11,*,err=95,end=95) (dummy, model_param(j,5:7))
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) (dummy, model_prior(j,5:7))
           case default
              goto 99
        end select
        if (dummy /= 'shape_param_prior') goto 96

        !read LD parameters
        select case (order)
           case (0)
              read (11,*,err=95,end=95) dummy
              if (dummy /= 'ld_param') goto 96
              read (11,*,err=95,end=95) dummy
           case default
              read (11,*,err=95,end=95) (dummy, model_param(j,8:7+order))
              if (dummy /= 'ld_param') goto 96
              read (11,*,err=95,end=95) (dummy, model_prior(j,8:7+order))
        end select

        !check to ensure everything inside limits
        do k=1, size(limits,1)
           if (model_param(j,k)<limits(k,1)) goto 100
           if (model_param(j,k)>limits(k,2)) goto 100
        end do
        !ensure alpha>0 for hestroffer case
        select case(trim(model_spec(j,3)))
           case('hestroffer')
              if (.not.(model_param(j,8) > 0D0)) goto 100
        end select

        if (dummy /= 'ld_param_prior') goto 96
        
     end if
  end do
  
2 continue 

  close (11)

  !clean-up and return
200 continue
  return
  
  !error trapping 
90 info = 'blank filename'
  goto 200
91 info = 'cannot open file'
  goto 200
92 info = 'cannot read from file'
  close (11)
  goto 200
93 info = 'file exceeds maximum permitted length'
  close (11)
  goto 200
94 info = 'unknown read problem - possible invalid file format'
  close (11)
  goto 200
95 info = 'unknown read problem whilst reading component data'
  close (11)
  goto 200
96 info = 'invalid keyword or keyword order'
  close (11)
  goto 200
97 info = 'illegal number of components present'
  goto 200
98 info = 'invalid limb darkening type'
  close (11)
  goto 200
99 info = 'invalid shape type'
  close (11)
  goto 200
100 info = 'illegal parameter value detected'
  close (11)
  goto 200
  
end subroutine read_model

!==============================================================================

end module Inout








