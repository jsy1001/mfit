!$Id: inout.f90,v 1.4 2002/11/18 10:28:58 jsy1001 Exp $

module Inout

!subroutines contained
!
!read_vis
!read_nvis
!read_mapdat

implicit none

contains

!==============================================================================

subroutine read_vis(info, file_name, source, max_lines, vis_data, lambda, &
                    calib_error)
    
  !Reads vis file with up to max_lines (excluding blanks)
  !(Re)Allocates and fills in vis_data array. Note projected
  !baseline sqrt(u**2 + v**2) gets
  !put in u column, with zero in v column
  !Adds calib_error, assumed to be constant fractional error in mod V
 
  !subroutine arguments
  character(len=128), intent(out) :: info, source
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: max_lines
  double precision, dimension(:,:), allocatable, intent(out) :: vis_data
  double precision, intent(in) :: calib_error
  double precision, intent(out) :: lambda

  !local variables
  character(len=32) :: dummy
  integer :: i, data_items
  double precision :: vis, baseline, default_error, frac_error

  !default absolute error on the visibility data points
  !(for consistency with analyse)
  default_error = 0.001D0
  
  !check for zero length filename
  if (file_name == '') then
     info = 'blank filename'
     return
  end if
  
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
  info = 'file exceeds maximum permitted length'
  close (11)
  return
1 close (11)
  
  !allocate array size
  if(allocated(vis_data)) deallocate(vis_data)
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
       !calculate fractional error on mod V
       frac_error = sqrt( calib_error**2D0 + (default_error/vis)**2D0 )
       !hence calculate abs error on squared visibility
       vis_data(i,5) = 2D0*frac_error*vis_data(i,4)
  end do

  close(11)
  return

  !error trapping 
91 info = 'cannot open file'
  return
92 info = 'cannot read from file'
  close (11)
  return
94 info = 'unknown read problem - possible invalid file format'
  close (11)
  return
  
  end subroutine read_vis

!==============================================================================

subroutine read_nvis(info, file_name, source, max_lines, vis_data, lambda, &
                    calib_error)
    
  !Reads nvis file with up to max_lines (excluding comments)
  !(Re)Allocates and fills in vis_data array. Note projected
  !baseline sqrt(u**2 + v**2) gets
  !put in u column, with zero in v column
  !Adds calib_error, assumed to be constant fractional error in mod V
 
  !subroutine arguments
  character(len=*), intent(in) :: file_name
  character(len=128), intent(out) :: info, source
  integer, intent(in) :: max_lines
  double precision, dimension(:,:), allocatable, intent(out) :: vis_data
  double precision, intent(in) :: lambda, calib_error

  !local variables
  character(len=256) :: line
  integer :: i, data_items
  double precision :: vis, err, baseline, frac_error
  
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
     if (line(1:1) /= '#' .and. len_trim(line) > 0) then
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
  if(allocated(vis_data)) deallocate(vis_data)
  allocate(vis_data(data_items,5))
  
  !read vis data properly and close
  !vis_data: lambda, u, v, cal_vis, err
  !where cal_vis is the squared visibility
  i = 1
  open (unit=11, action='read', file=file_name)
  do
     read (11, '(a)', err=92, end=2) line
     if (line(1:1) /= '#' .and. len_trim(line) > 0) then
        read (line, *, err=94) (vis, err, baseline)
        vis_data(i,1) = lambda
        vis_data(i,2) = abs(baseline)/1D+3
        vis_data(i,3) = 0D0
        vis_data(i,4) = vis**2D0
        !calculate fractional error on mod V
        frac_error = sqrt( calib_error**2D0 + (err/vis)**2D0 )
        !hence calculate abs error on squared visibility
        vis_data(i,5) = 2D0*frac_error*vis_data(i,4)
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

subroutine read_mapdat(info, file_name, source, max_lines, &
     vis_data, triple_data, wavebands, calib_error)

  !Reads mapdat file with up to max_lines (excluding blanks)
  !
  !(Re)Allocates and fills in vis_data and triple_data arrays:
  !
  !vis_data: lambda, u, v, vis, err
  !          vis is squared visibility amplitude (may be -ve for data points)
  !triple_data: lambda, u1, v1, u2, v2, amp, err, cp, err
  !
  !Adds calib_error, assumed to be constant fractional error in mod V
  !
  !wavebands array (allocated here) is list of different wavelength values
  !encountered
  
  !subroutine arguments
  character(len=128), intent(out) :: info, source
  character(len=*), intent(in) :: file_name
  integer, intent(in) :: max_lines
  double precision, intent(in) :: calib_error
  double precision, dimension(:,:), allocatable, intent(out) :: vis_data
  double precision, dimension(:,:), allocatable, intent(out) :: triple_data
  double precision, dimension(:), allocatable, intent(out) :: wavebands

  !local variables
  character(len=32) :: dummy, source1, source2
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
  if(allocated(vis_data)) deallocate(vis_data)
  if(allocated(triple_data)) deallocate(triple_data)
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
        !(in the squared visibility) divided by 2. Need to add calibration error
        !and convert to absolute error
        vis_data(i1,5) = vis_data(i1,4) * sqrt((2D0*calib_error)**2D0 + &
                         (2D0*vis_err)**2D0 )

        if (vis_err<0D0) vis_data(i1,5) = -vis_data(i1,5)

     else if (dummy == 'triple') then

        i2 = i2 + 1
        read (12,*,err=95) (dummy, dummy, dummy, dummy, triple_data(i2,1), &
             dummy, dummy, dummy, dummy, triple_data(i2,2:5), &
             amp, amp_err, cp, cp_err)

        !reverse signs of u1, v1, u2, v2 as Caltech baseline convention is
        !opposite to the one used here
        triple_data(i2,2:5) = -triple_data(i2,2:5)

        !amplitude in file is cube root of amplitude
        triple_data(i2,6) = amp**3D0

        !error in amplitude is 1/3 of fractional error in the amplitude
        triple_data(i2,7) = abs(3D0*amp_err*triple_data(i2,6))
        !preserve sign of tp error - this is flagging info
        if (amp_err<0D0) triple_data(i2,7) = -triple_data(i2,7)

        triple_data(i2,8) = cp
        triple_data(i2,9) = cp_err

     else if (dummy == 'source') then

        read (12,*,err=95) dummy, source1, source2
        source = trim(source1)//' '//trim(source2)

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
  if(allocated(wavebands)) deallocate(wavebands)
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

end module Inout








