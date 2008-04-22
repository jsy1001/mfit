!$Id: inout.f90,v 1.23 2008/04/22 11:42:07 jsy1001 Exp $

module Inout

  implicit none

  private

  !public subroutines contained:
  !
  !read_vis
  !read_nvis
  !read_calib
  !read_mapdat
  !read_oi_fits
  !yesno

  public :: read_vis, read_nvis, read_calib, read_mapdat, read_oi_fits, yesno

  !private subroutines contained:
  !read_oi_wavelength
  !read_oi_vis2
  !read_oi_t3
  !DateToMjd

  !public module variables contained:

  public :: release

  character(len=5), parameter :: release = '1.5.1' !! Package release number

contains
  
  !============================================================================
  
  subroutine read_vis(info, file_name, source, max_lines, &
       vis_data, num_vis, waveband, calib_error)

    !Reads vis file with up to max_lines (excluding blanks)
    !(Re)Allocates and fills in vis_data array. Note projected
    !baseline sqrt(u**2 + v**2) gets
    !put in u column, with zero in v column
    !Adds calib_error, assumed to be constant fractional error in mod V

    !subroutine arguments
    character(len=*), intent(out) :: info, source
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: max_lines
    integer, intent(out) :: num_vis
    double precision, dimension(:,:), allocatable, intent(out) :: vis_data
    double precision, intent(in) :: calib_error
    double precision, dimension(2), intent(in) :: waveband

    !local variables
    character(len=32) :: dummy
    integer :: i
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
    read (11, *, err=92, end=92) source, num_vis
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
1   close (11)

    !allocate array size
    if(allocated(vis_data)) deallocate(vis_data)
    allocate(vis_data(num_vis,7))

    !read vis data properly and close
    !vis_data: lambda, delta_lambda, u, v, (vis^2), err
    open (unit=11, action='read', file=file_name)
    read (11, *, err=94) dummy, dummy !skip header
    do i = 1, num_vis
       read (11, *, err=94) vis, baseline
       vis_data(i,1:2) = waveband
       vis_data(i,3) = abs(baseline)/1D+3
       vis_data(i,4) = 0D0
       vis_data(i,5) = vis**2D0
       !calculate fractional error on mod V
       frac_error = sqrt( calib_error**2D0 + (default_error/vis)**2D0 )
       !hence calculate abs error on squared visibility
       vis_data(i,6) = 2D0*frac_error*vis_data(i,5)
       vis_data(i,7) = -1D0 !MJD unknown
    end do

    close(11)
    return

    !error trapping 
91  info = 'cannot open file'
    return
92  info = 'cannot read from file'
    close (11)
    return
94  info = 'unknown read problem - possible invalid file format'
    close (11)
    return

  end subroutine read_vis

  !============================================================================

  subroutine read_nvis(info, file_name, source, max_lines, &
       vis_data, num_vis, waveband, calib_error)

    !Reads nvis file with up to max_lines (excluding comments)
    !(Re)Allocates and fills in vis_data array. Note projected
    !baseline sqrt(u**2 + v**2) gets
    !put in u column, with zero in v column
    !Adds calib_error, assumed to be constant fractional error in mod V

    !subroutine arguments
    character(len=*), intent(in) :: file_name
    character(len=*), intent(out) :: info, source
    integer, intent(in) :: max_lines
    integer, intent(out) :: num_vis
    double precision, dimension(:,:), allocatable, intent(out) :: vis_data
    double precision, dimension(2), intent(in) :: waveband
    double precision, intent(in) :: calib_error

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
    open (unit=11, status='old', action='read', err=92, file=file_name)
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

1   close (11)

    !allocate array size
    if(allocated(vis_data)) deallocate(vis_data)
    allocate(vis_data(data_items,7))

    !read vis data properly and close
    !vis_data: lambda, delta_lambda, u, v, (vis^2), err
    i = 1
    open (unit=11, status='old', action='read', err=92, file=file_name)
    do
       read (11, '(a)', err=92, end=2) line
       if (line(1:1) /= '#' .and. len_trim(line) > 0) then
          read (line, *, err=94) vis, err, baseline
          vis_data(i,1:2) = waveband
          vis_data(i,3) = abs(baseline)/1D+3
          vis_data(i,4) = 0D0
          vis_data(i,5) = vis**2D0
          if (vis /= 0) then
             !calculate fractional error on mod V
             frac_error = sqrt( calib_error**2D0 + (err/vis)**2D0 )
             !hence calculate abs error on squared visibility
             vis_data(i,6) = 2D0*frac_error*vis_data(i,5)
          else
             vis_data(i,6) = 2D0*err
          end if
          vis_data(i,7) = -1D0 !MJD unknown
          i = i + 1
          if (i > max_lines) then
             info = 'file exceeds maximum permitted length'
             close (11)
             return
          end if
       end if
    end do

2   close(11)
    return

    !error trapping 
92  info = 'cannot read from file'
    close (11)
    return
94  info = 'unknown read problem - possible invalid file format'
    close (11)
    return

  end subroutine read_nvis

  !============================================================================

  subroutine read_calib(info, file_name, source, max_lines, vis_data, num_vis,&
       waveband, wavebands, calib_error)

    !Reads ouput produced by the MSC V2 calibration applications
    !wbCalib (wide band) and nbCalib (narrow band)
    !
    !Determines which application output it is by looking at item 7 in table
    !if integer (no '.') it is NChannels and nbCalib data,
    !if floating point then it is wavelength and wbCalib data

    character(len=*), intent(in) :: file_name
    character(len=*), intent(out) :: info, source
    integer, intent(in) :: max_lines
    integer, intent(out) :: num_vis
    double precision, dimension(:,:), allocatable, intent(out) :: vis_data
    double precision, intent (in) :: waveband
    double precision, dimension(:,:), allocatable, intent(out) :: wavebands
    double precision, intent(in) :: calib_error

    !local variables
    character(len=1024) :: line, dummy
    integer :: i, j, k, loc, channels, bands
    logical :: newband
    double precision, dimension (:), allocatable :: data
    double precision :: frac_error

    !check for zero length filename
    if (file_name == '') then
       info = 'blank filename'
       return
    end if

    !count number of data lines
    !if nbCalib data increment data lines by NChannels
    num_vis = 0
    open (unit=11, status='old', action='read', err=92, file=file_name)
    do i = 1, max_lines+1
       read (11, '(a)', err=94, end=1) line
       line = trim (adjustl (line))
       if (line(1:1) /= '#' .and. len_trim(line) > 0) then
          !remove '/' character (confuses read statement)
          do
             loc = index (line, '/')
             if (loc == 0) exit
             line (loc:loc) = '.'
          end do
          !determine if wbCalib or nbCalib and read data from line
          read (line, *, err=94) source, dummy, dummy, dummy, dummy, dummy, dummy
          if (index (dummy, '.') /= 0) then
             num_vis = num_vis + 1
          else
             read (dummy, *, err=94) channels
             num_vis = num_vis + channels
          end if
       end if
    end do
    info = 'file exceeds maximum permitted length'
    close (11)
    return

1   close (11)

    !allocate array size
    if(allocated(vis_data)) deallocate(vis_data)
    allocate(vis_data(num_vis,7))

    !read *Calib data and close
    i = 1
    bands = 0
    open (unit=11, status='old', action='read', err=92, file=file_name)
    do
       !read line and check if comment or data
       read (11, '(a)', err=92, end=2) line
       line = trim (adjustl (line))
       if (line(1:1) /= '#' .and. len_trim(line) > 0) then
          !remove '/' character (confuses read statement)
          do
             loc = index (line, '/')
             if (loc == 0) exit
             line (loc:loc) = '.'
          end do

          !determine if wbCalib or nbCalib and read data from line
          read (line, *, err=94) dummy, dummy, dummy, dummy, dummy, dummy, dummy
          if (index (dummy, '.') /= 0) then
             !we have wbCalib data
             channels = 1
             allocate(data (7))
             read (line, *, err=94) dummy, dummy, dummy, dummy, dummy, dummy, &
                  data(1:7), dummy, vis_data(i,3:4)
          else
             !we have nbCalib data
             read (dummy, *, err=94) channels
             allocate(data (channels*7))
             read (line, *, err=94) dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
                  data(1:channels*7), dummy, vis_data(i,3:4)
          end if

          !conform data
          do j = 0, channels-1
             vis_data(i+j,1) = data(j*7+1) * 1D+3
             newband = .true.
             do k = 1, i+j-1
                if (vis_data(k,1) == vis_data(i+j,1)) newband = .false.
             end do
             if (newband) bands = bands + 1
             vis_data(i+j,2) = waveband
             vis_data(i+j,3:4) = vis_data(i,3:4)
             vis_data(i+j,5) = data(j*7+2)
             !calculate fractional error on mod V
             frac_error = sqrt( calib_error**2D0 + (0.5*data(j*7+3)/sqrt(abs(vis_data(i+j,5))))**2D0 )
             !hence calculate abs error on squared visibility
             vis_data(i+j,6) = 2D0*frac_error*vis_data(i+j,5)
             vis_data(i+j,7) = -1D0 !MJD unknown
          end do

          !deallocate arrays
          deallocate(data)

          !increment i
          i = i + channels
       end if
    end do

2   close(11)

    !set waveband data
    allocate(wavebands(bands,2))
    do i = 1, bands
       do j = 1, num_vis
          newband = .true.
          do k = 1, i
             if (vis_data(j,1) == wavebands(k,1)) newband = .false.
          end do
          if (newband) then
             wavebands(i, 1) = vis_data (j, 1)
             exit
          end if
       end do
       wavebands(i, 2) = waveband
    end do

    return

    !error trapping 
92  info = 'cannot read from file'
    close (11)
    return
94  info = 'unknown read problem - possible invalid file format'
    close (11)
    return

  end subroutine read_calib

  !============================================================================

  subroutine read_mapdat(info, file_name, source, max_lines, &
       vis_data, num_vis, triple_data, num_triple, wavebands, calib_error)

    !Reads mapdat file with up to max_lines (excluding blanks)
    !
    !(Re)Allocates and fills in vis_data and triple_data arrays:
    !
    !vis_data: lambda, delta_lambda, u, v, (vis^2), err
    !      (vis^2) is squared visibility amplitude (may be -ve for data points)
    !triple_data: lambda, delta_lambda, u1, v1, u2, v2, amp, err, cp, err
    !
    !Adds calib_error, assumed to be constant fractional error in mod V
    !
    !wavebands array (allocated here) gives different (wavelength, bandwidth)
    !pairs encountered

    !subroutine arguments
    character(len=*), intent(out) :: info, source
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: max_lines
    integer, intent(out) :: num_vis, num_triple
    double precision, intent(in) :: calib_error
    double precision, dimension(:,:), allocatable, intent(out) :: vis_data
    double precision, dimension(:,:), allocatable, intent(out) :: triple_data
    double precision, dimension(:,:), allocatable, intent(out) :: wavebands

    !local variables
    character(len=32) :: dummy, source1, source2
    integer :: i, j, i1, i2, lines, num, year, month, day, hour, min
    real :: sec
    double precision :: vis, vis_err, amp, amp_err, cp, cp_err, mjd
    double precision, dimension(2) :: swap
    double precision, dimension(:,:), allocatable :: all_wb

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
1   lines = i-1
    close (11)

    if ((i1 == 0) .and. (i2 == 0)) goto 94

    !allocate vis and triple data array sizes
    if(allocated(vis_data)) deallocate(vis_data)
    if(allocated(triple_data)) deallocate(triple_data)
    allocate(vis_data(i1,7))
    allocate(triple_data(i2,11))
    num_vis = i1
    num_triple = i2

    !read data properly and close
    open (unit=11, action='read', file=file_name)
    open (unit=12, action='read', file=file_name)
    i1 = 0 !counters for vis and triple data item
    i2 = 0
    do i = 1, lines
       read (11, *) dummy

       if (dummy == 'group_date') then

          read (12,*,err=95) dummy, year, month, day

       else if (dummy == 'vis') then

          i1 = i1 + 1
          read (12,*,err=95) dummy, dummy, dummy, vis_data(i1,1:2), &
               hour, min, sec, vis_data(i1,3:4), vis, vis_err

          !reverse signs of u, v, to correct for inconsistent sign in FT
          !doesn't actually matter for V^2 data
          vis_data(i1,3:4) = -vis_data(i1,3:4)

          !V^2 positive stored as sqroot(V^2) in file
          !V^2 negative stored as -sqroot(-V^2) in file
          vis_data(i1,5) = vis**2D0
          if (vis<0D0) vis_data(i1,5) = -vis_data(i1,5)

          !"vis" error in file is the modulus of the fractional error
          !(in the squared visibility) divided by 2. Need to add calibration error
          !and convert to absolute error
          vis_data(i1,6) = abs(vis_data(i1,5)) * sqrt((2D0*calib_error)**2D0 + &
               (2D0*vis_err)**2D0 )
          mjd = DateToMjd(year, month, day) + &
               ((((hour*60) + min)*60) + sec)/86400.
          vis_data(i1,7) = mjd

          if (vis_err<0D0) vis_data(i1,6) = -vis_data(i1,6)

       else if (dummy == 'triple') then

          i2 = i2 + 1
          read (12,*,err=95) dummy, dummy, dummy, dummy, triple_data(i2,1:2), &
               hour, min, sec, triple_data(i2,3:6), &
               amp, amp_err, cp, cp_err

          !reverse signs of u1, v1, u2, v2 to correct for inconsistent sign in FT
          triple_data(i2,3:6) = -triple_data(i2,3:6)

          !amplitude in file is cube root of amplitude
          triple_data(i2,7) = amp**3D0

          !error in amplitude is 1/3 of fractional error in the amplitude
          triple_data(i2,8) = abs(3D0*amp_err*triple_data(i2,7))
          !preserve sign of tp error - this is flagging info
          if (amp_err<0D0) triple_data(i2,8) = -triple_data(i2,8)

          triple_data(i2,9) = cp
          triple_data(i2,10) = cp_err
          mjd = DateToMjd(year, month, day) + &
               ((((hour*60) + min)*60) + sec)/86400.
          triple_data(i2,11) = mjd

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
    !collate all (wavelength, bandwidth) pairs
    allocate(all_wb(num_vis+num_triple, 2))
    all_wb(:num_vis, :) = vis_data(:, 1:2)
    all_wb(num_vis+1:, :) = triple_data(:, 1:2)
    !mark duplicates
    num = 0
    do i = 1, size(all_wb,1)
       if (all_wb(i, 1) /= -1) then
          num = num + 1
          do j = i+1, size(all_wb,1)
             if (all_wb(j, 1) == all_wb(i, 1) &
                  .and. all_wb(j, 2) == all_wb(i, 2)) all_wb(j, 1) = -1
          end do
       end if
    end do
    !make list of unique wavebands found
    if(allocated(wavebands)) deallocate(wavebands)
    allocate(wavebands(num, 2))
    j = 0
    do i = 1, size(all_wb,1)
       if (all_wb(i, 1) /= -1) then
          j = j + 1
          wavebands(j, :) = all_wb(i, :)
       end if
    end do
    !sort into ascending order of wavelengths (ignore bw)
    do i = num, 2, -1
       do j = 1, i-1
          if (wavebands(j, 1) > wavebands(j+1, 1)) then
             swap = wavebands(j, :) 
             wavebands(j, :) = wavebands(j+1, :)
             wavebands(j+1, :) = swap
          end if
       end do
    end do

    !clean-up and return
    if (allocated(all_wb)) deallocate(all_wb)
200 continue

    return

    !error trapping 
90  info = 'blank filename'
    goto 200
91  info = 'cannot open file'
    goto 200
92  info = 'cannot read from file'
    close (11)
    goto 200
93  info = 'file exceeds maximum permitted length'
    close (11)
    goto 200
94  info = 'no vis or triple data found'
    goto 200
95  info = 'unknown read problem - possible invalid file format'
    close (11)
    goto 200

  end subroutine read_mapdat

  !============================================================================

  subroutine read_oi_fits(info, file_name, user_target_id, source, &
       vis_data, num_vis, triple_data, num_triple, wavebands, calib_error)

    !Reads OIFITS format
    !
    !(Re)Allocates and fills in vis_data and triple_data arrays:
    !
    !vis_data: lambda, u, v, (vis^2), err
    !      (vis^2) is squared visibility amplitude (may be -ve for data points)
    !triple_data: lambda, u1, v1, u2, v2, amp, err, cp, err
    !
    !wavebands array (allocated here) gives different (wavelength, bandwidth)
    !pairs encountered.
    !Reads rows with TARGET_ID == user_target_id,
    !unless user_target_id is -1, when 1st TARGET_ID in OI_TARGET table is used
    !
    !Deals with arbitrary table positions in file, and references to
    !OI_WAVELENGTH tables via INSNAME. Code to actually read the data is quite
    !complicated because returned data arrays are not a good match for the
    !structure of the tables in the file
    !
    !Adds calib_error, assumed to be constant fractional error in mod V
    !
    !On error, returns message in info (otherwise empty string)

    !subroutine arguments
    character(len=*), intent(out) :: info, source
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: user_target_id
    integer, intent(out) :: num_vis, num_triple
    double precision, intent(in) :: calib_error
    double precision, dimension(:, :), allocatable, intent(out) :: vis_data
    double precision, dimension(:, :), allocatable, intent(out) :: triple_data
    double precision, dimension(:, :), allocatable, intent(out) :: wavebands

    !local variables
    integer, parameter :: maxhdu = 40 !must exceed no. of tables in file
    integer, dimension(maxhdu) :: vis2_hdus, t3_hdus, wl_hdus
    integer, dimension(maxhdu) :: vis2_rows, t3_rows, vis2_nwave, t3_nwave
    character(len=80), dimension(maxhdu) :: vis2_insname, t3_insname, wl_insname
    integer :: status, unit, blocksize, hdutype, nhdu, maxwave, ntarget
    integer :: nvis2, nt3, nwl, iwl, irow, ntot, nwave, iwave, itab, i, j
    integer :: colnum, naxis, count, num_wb
    integer :: target_id, sel_target_id
    integer, dimension(2) :: vis2sta
    integer, dimension(3) :: naxes, t3sta
    character(len=80) :: keyval, comment
    real, dimension(:, :), allocatable :: tab_wb, all_wb
    double precision, dimension(2) :: swap
    double precision, dimension(:), allocatable :: vis2data, vis2err
    double precision, dimension(:), allocatable :: t3amp, t3amperr
    double precision, dimension(:), allocatable :: t3phi, t3phierr
    logical, dimension(:), allocatable :: flag
    integer, dimension(:), allocatable :: all_ids
    character (len=16), dimension(:), allocatable :: all_names
    double precision :: mjd, ucoord, vcoord, u2coord, v2coord
    double precision :: frac_error, tot_error

    !Open FITS file
    status = 0
    call ftgiou(unit, status) !get unused unit number
    call ftopen(unit, file_name, 0, blocksize, status)
    if (status /= 0) then
       print *, 'FITSIO Error opening file:'
       call ftgerr(status, info)
       return
    end if

    !Get name & target_id of chosen target
    !call ftmnhd(unit, 2, 'OI_TARGET', 0, status) !move to OI_TARGET
    !can't link to ftmnhd for some reason, so find OI_TARGET by brute force
    nhdu = 1
    hduloop0: do
       call ftmrhd(unit, 1, hdutype, status)
       if (status /= 0) exit hduloop0 !most likely EOF
       nhdu = nhdu + 1
       call ftgkys(unit, 'EXTNAME', keyval, comment, status)
       if (trim(keyval) == 'OI_TARGET') then
          call ftgkyj(unit, 'NAXIS2', ntarget, comment, status)
          allocate(all_ids(ntarget))
          allocate(all_names(ntarget))
          call read_oi_target(unit, ntarget, all_ids, all_names, status)
          if (user_target_id == -1) then
             sel_target_id = all_ids(1)
          else
             sel_target_id = user_target_id
          end if
          print *, 'Read OI_TARGET table containing:'
          do irow = 1, ntarget
             print '(i4, 2x, a16)', all_ids(irow), all_names(irow)
             if (all_ids(irow) == sel_target_id) then
                print *, '^^^^^^^^^^^^^^^^^^^^^^'
                source = all_names(irow)
             end if
          end do
          deallocate(all_ids)
          deallocate(all_names)
          exit hduloop0
       end if
    end do hduloop0
    if (status /= 0 .and. status /= 107) then
       !return if error wasn't EOF
       print *, 'FITSIO Error scanning for OI_TARGET:'
       call ftgerr(status, info)
       return
    end if

    !Pass 1 - Scan HDUs
    !Want positions and sizes of data tables, and positions of
    !OI_WAVELENGTH tables. Also get INSNAME values for cross-referencing
    nvis2 = 0
    nt3 = 0
    nwl = 0
    maxwave = 0 !will be upper limit on number of distinct wavebands in file
    nhdu = 1
    call ftmahd(unit, nhdu, hdutype, status) !go back to start
    hduloop: do
       call ftmrhd(unit, 1, hdutype, status)
       if (status /= 0) exit hduloop !most likely EOF
       nhdu = nhdu + 1
       if(nhdu > maxhdu) stop 'Too many HDUs. Recompile to change limit'
       call ftgkys(unit, 'EXTNAME', keyval, comment, status)
       if (trim(keyval) == 'OI_VIS2') then
          !Get table size & INSNAME for OI_VIS2
          nvis2 = nvis2 + 1
          vis2_hdus(nvis2) = nhdu !remember position
          call ftgkyj(unit, 'NAXIS2', ntot, comment, status)
          call ftgcno(unit, .false., 'VIS2DATA', colnum, status)
          call ftgtdm(unit, colnum, 3, naxis, naxes, status)
          vis2_nwave(nvis2) = naxes(1)
          maxwave = maxwave + vis2_nwave(nvis2)
          !read all rows to find no. matching sel_target_id
          vis2_rows(nvis2) = 0
          allocate(vis2data(vis2_nwave(nvis2)))
          allocate(vis2err(vis2_nwave(nvis2)))
          allocate(flag(vis2_nwave(nvis2)))
          do irow = 1, ntot
             call read_oi_vis2(unit, irow, vis2_nwave(nvis2), target_id, mjd, &
                  vis2data, vis2err, ucoord, vcoord, vis2sta, flag, status)
             if (target_id == sel_target_id) vis2_rows(nvis2) = vis2_rows(nvis2) + 1
          end do
          call ftgkys(unit, 'INSNAME', vis2_insname(nvis2), comment, status)
          print *, nhdu, ': OI_VIS2 table with', vis2_rows(nvis2), &
               'matching rows, INSNAME=', trim(vis2_insname(nvis2))
          deallocate(vis2data)
          deallocate(vis2err)
          deallocate(flag)
       else if (trim(keyval) == 'OI_T3') then
          !Get table size & INSNAME for OI_T3
          nt3 = nt3 + 1
          t3_hdus(nt3) = nhdu !remember position
          call ftgkyj(unit, 'NAXIS2', ntot, comment, status)
          call ftgcno(unit, .false., 'T3AMP', colnum, status)
          call ftgtdm(unit, colnum, 3, naxis, naxes, status)
          t3_nwave(nt3) = naxes(1)
          maxwave = maxwave + t3_nwave(nt3)
          !read all rows to find no. matching sel_target_id
          t3_rows(nt3) = 0
          allocate(t3amp(t3_nwave(nt3)))
          allocate(t3amperr(t3_nwave(nt3)))
          allocate(t3phi(t3_nwave(nt3)))
          allocate(t3phierr(t3_nwave(nt3)))
          allocate(flag(t3_nwave(nt3)))
          do irow = 1, ntot
             call read_oi_t3(unit, irow, t3_nwave(nt3), target_id, mjd, &
                  t3amp, t3amperr, t3phi, t3phierr, &
                  ucoord, vcoord, u2coord, v2coord, t3sta, flag, status)
             if (target_id == sel_target_id) t3_rows(nt3) = t3_rows(nt3) + 1
          end do
          call ftgkys(unit, 'INSNAME', t3_insname(nt3), comment, status)
          print *, nhdu, ': OI_T3 table with', t3_rows(nt3), &
               'matching rows, INSNAME=', trim(t3_insname(nt3))
          deallocate(t3amp)
          deallocate(t3amperr)
          deallocate(t3phi)
          deallocate(t3phierr)
          deallocate(flag)
       else if (trim(keyval) == 'OI_WAVELENGTH') then
          !Get INSNAME for OI_WAVELENGTH
          nwl = nwl + 1
          wl_hdus(nwl) = nhdu !remember position
          call ftgkys(unit, 'INSNAME', wl_insname(nwl), comment, status)
          print *, nhdu, ': OI_WAVELENGTH table with INSNAME=', &
               trim(wl_insname(nwl))
       else
          print *, nhdu, ': Skipping ', trim(keyval), ' table'
       end if
    end do hduloop
    if (status /= 107) then
       !return if error wasn't EOF
       print *, 'FITSIO Error in Pass 1:'
       call ftgerr(status, info)
       return
    end if
    status = 0

    !Allocate storage
    if(allocated(vis_data)) deallocate(vis_data)
    if(allocated(triple_data)) deallocate(triple_data)
    num_vis = 0
    do itab = 1, nvis2
       num_vis = num_vis + vis2_rows(itab)*vis2_nwave(itab)
    end do
    num_triple = 0
    do itab = 1, nt3
       num_triple = num_triple + t3_rows(itab)*t3_nwave(itab)
    end do
    allocate(vis_data(num_vis, 7))
    allocate(triple_data(num_triple, 11))
    allocate(all_wb(maxwave, 2))
    num_wb = 0

    !Pass 2 - read the data
    !Read OI_VIS2 tables in turn
    count = 1
    do itab = 1, nvis2
       nwave = vis2_nwave(itab)
       do iwl = 1, nwl
          if (trim(wl_insname(iwl)) == trim(vis2_insname(itab))) then
             !Read 1st matching OI_WAVELENGTH table
             allocate(tab_wb(nwave, 2))
             call ftmahd(unit, wl_hdus(iwl), hdutype, status)
             call read_oi_wavelength(unit, nwave, tab_wb(:, 1), tab_wb(:, 2), &
                  status)
             tab_wb = tab_wb/1e-9 !store wavelength/bw in nm
             exit
          end if
       end do
       if (.not. allocated(tab_wb)) then
          write(info,*) 'No OI_WAVELENGTH with INSNAME=', &
               trim(vis2_insname(itab))
          return
       end if
       !Read OI_VIS2 table
       allocate(vis2data(nwave))
       allocate(vis2err(nwave))
       allocate(flag(nwave))
       call ftmahd(unit, vis2_hdus(itab), hdutype, status)
       !read 1 row at a time to save memory
       do irow = 1, vis2_rows(itab)
          call read_oi_vis2(unit, irow, nwave, target_id, mjd, &
               vis2data, vis2err, ucoord, vcoord, vis2sta, flag, status)
          !map each binary table row to nwave rows of vis_data array
          wave1: do iwave = 1, nwave
             vis_data(count, 1:2) = tab_wb(iwave, :)
             !reverse signs of u, v, to correct for inconsistent sign in FT
             !doesn't actually matter for V^2 data
             vis_data(count, 3) = -ucoord
             vis_data(count, 4) = -vcoord
             vis_data(count, 5) = vis2data(iwave) !may be -ve
             if (vis2data(iwave) /= 0D0) then
                frac_error = abs(vis2err(iwave)/vis2data(iwave))
                tot_error = abs(vis2data(iwave)) &
                     * sqrt((2D0*calib_error)**2D0 + frac_error**2D0 )
             else
                tot_error = vis2err(iwave)
             end if
             if (flag(iwave)) then
                if (tot_error == 0D0) then
                   vis_data(count, 6) = -1.0D0
                else
                   vis_data(count, 6) = -tot_error
                end if
             else
                vis_data(count, 6) = tot_error
             end if
             vis_data(count, 7) = mjd
             count = count + 1
             do j = 1, num_wb
                if (tab_wb(iwave, 1) .eq. all_wb(j, 1) &
                     .and. tab_wb(iwave, 2) .eq. all_wb(j, 2)) cycle wave1
             end do
             !add new waveband to list
             num_wb = num_wb + 1
             all_wb(num_wb, :) = tab_wb(iwave, :)
          end do wave1
       end do
       deallocate(tab_wb)
       deallocate(vis2data)
       deallocate(vis2err)
       deallocate(flag)
    end do

    !Read OI_T3 tables in turn
    count = 1
    do itab = 1, nt3
       nwave = t3_nwave(itab)
       do iwl = 1, nwl
          if (trim(wl_insname(iwl)) == trim(t3_insname(itab))) then
             !Read 1st matching OI_WAVELENGTH table
             allocate(tab_wb(nwave, 2))
             call ftmahd(unit, wl_hdus(iwl), hdutype, status)
             call read_oi_wavelength(unit, nwave, tab_wb(:, 1), tab_wb(:, 2), &
                  status)
             tab_wb = tab_wb/1e-9 !store wavelength/bw in nm
             exit
          end if
       end do
       if (.not. allocated(tab_wb)) then
          write(info,*) 'No OI_WAVELENGTH with INSNAME=', &
               trim(t3_insname(itab))
          return
       end if
       !Read OI_T3 table
       allocate(t3amp(nwave))
       allocate(t3amperr(nwave))
       allocate(t3phi(nwave))
       allocate(t3phierr(nwave))
       allocate(flag(nwave))
       call ftmahd(unit, t3_hdus(itab), hdutype, status)
       do irow = 1, t3_rows(itab)
          call read_oi_t3(unit, irow, nwave, target_id, mjd, &
               t3amp, t3amperr, t3phi, t3phierr, &
               ucoord, vcoord, u2coord, v2coord, t3sta, flag, status)
          !map each binary table row to nwave rows of triple_data array
          wave2: do iwave = 1, nwave
             triple_data(count, 1:2) = tab_wb(iwave, :)
             !reverse signs of u1, v1, u2, v2 to correct for inconsistent sign in FT
             triple_data(count, 3) = -ucoord
             triple_data(count, 4) = -vcoord
             triple_data(count, 5) = -u2coord
             triple_data(count, 6) = -v2coord
             triple_data(count, 7) = t3amp(iwave)
             triple_data(count, 9) = t3phi(iwave)
             !read_oi_t3 returns -ve t3amperr if NULL t3amp
             if (flag(iwave)) then
                !flag amplitude and phase
                if (t3amperr(iwave) /= 0) then
                   triple_data(count, 8) = -abs(t3amperr(iwave))
                else
                   triple_data(count, 8) = -1.0D0
                end if
                if (t3phierr(iwave) /= 0) then
                   triple_data(count, 10)= -t3phierr(iwave)
                else
                   triple_data(count, 10) = -99D0
                end if
             else
                triple_data(count, 8) = t3amperr(iwave)
                triple_data(count, 10) = t3phierr(iwave)
             end if
             triple_data(count, 11) = mjd
             count = count + 1
             do j = 1, num_wb
                if (tab_wb(iwave, 1) .eq. all_wb(j, 1) &
                     .and. tab_wb(iwave, 2) .eq. all_wb(j, 2)) cycle wave2
             end do
             !add new waveband to list
             num_wb = num_wb + 1
             all_wb(num_wb, :) = tab_wb(iwave, :)
          end do wave2
       end do
       deallocate(tab_wb)
       deallocate(t3amp)
       deallocate(t3amperr)
       deallocate(t3phi)
       deallocate(t3phierr)
       deallocate(flag)
    end do

    allocate(wavebands(num_wb, 2))
    wavebands = all_wb(:num_wb, :)
    deallocate(all_wb)
    !sort into ascending order of wavelengths (ignore bw)
    do i = num_wb, 2, -1
       do j = 1, i-1
          if (wavebands(j, 1) > wavebands(j+1, 1)) then
             swap = wavebands(j, :) 
             wavebands(j, :) = wavebands(j+1, :)
             wavebands(j+1, :) = swap
          end if
       end do
    end do
    call ftclos(unit, status)
    call ftfiou(unit, status) !free unit number

  end subroutine read_oi_fits

  !============================================================================

  subroutine read_oi_target(unit, ntarget, target_id, target_name, status)

    !Read OI_TARGET table: TARGET_ID and TARGET columns
    !OI_TARGET must be current HDU
    !To read all rows, set ntarget equal to value of NAXIS2

    integer, intent(in) :: unit
    integer, intent(in) :: ntarget
    integer, dimension(ntarget), intent(out) :: target_id
    character (len=16), dimension(ntarget), intent(out) :: target_name
    integer, intent(inout) :: status

    !local variables
    integer colnum, nullint
    character (len=1) :: nullstr
    logical anyf

    call ftgcno(unit, .false., 'TARGET_ID', colnum, status)
    call ftgcvj(unit, colnum, 1, 1, ntarget, nullint, &
         target_id, anyf, status)
    call ftgcno(unit, .false., 'TARGET', colnum, status)
    call ftgcvs(unit, colnum, 1, 1, ntarget, nullstr, &
         target_name, anyf, status)

  end subroutine read_oi_target

  !============================================================================

  subroutine read_oi_wavelength(unit, nwave, eff_wave, eff_band, status)

    !Read columns of OI_WAVELENGTH FITS binary table using fitsio calls
    !OI_WAVELENGTH must be current HDU

    integer, intent(in) :: unit, nwave
    real, dimension(nwave), intent(out) :: eff_wave, eff_band
    integer, intent(inout) :: status

    !local variables
    integer colnum
    real nullreal
    logical anyf

    call ftgcno(unit, .false., 'EFF_WAVE', colnum, status)
    call ftgcve(unit, colnum, 1, 1, nwave, nullreal, &
         eff_wave, anyf, status)
    call ftgcno(unit, .false., 'EFF_BAND', colnum, status)
    call ftgcve(unit, colnum, 1, 1, nwave, nullreal, &
         eff_band, anyf, status)

  end subroutine read_oi_wavelength

  !============================================================================

  subroutine read_oi_vis2(unit, row, nwave, target_id, mjd, &
       vis2data, vis2err, ucoord, vcoord, sta_index, flag, status)

    !Read specified row of OI_VIS2 FITS binary table using fitsio calls
    !OI_VIS2 must be current HDU

    integer, intent(in) :: unit, row, nwave
    double precision, dimension(nwave), intent(out) :: vis2data, vis2err
    double precision, intent(out) :: mjd, ucoord, vcoord
    integer, intent(out) :: target_id
    integer, dimension(2), intent(out) :: sta_index
    logical, dimension(nwave), intent(out) :: flag
    integer, intent(inout) :: status

    !local variables
    integer, dimension(1) :: target_id_array
    double precision, dimension(1) :: mjd_array, uv_array
    integer colnum, nullint
    double precision nullval
    logical anyf

    call ftgcno(unit, .false., 'TARGET_ID', colnum, status)
    call ftgcvj(unit, colnum, row, 1, 1, nullint, &
         target_id_array, anyf, status)
    target_id = target_id_array(1)
    call ftgcno(unit, .false., 'MJD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         mjd_array, anyf, status)
    mjd = mjd_array(1)
    call ftgcno(unit, .false., 'VIS2DATA', colnum, status)
    call ftgcvd(unit, colnum, row, 1, nwave, nullval, &
         vis2data, anyf, status)
    call ftgcno(unit, .false., 'VIS2ERR', colnum, status)
    call ftgcvd(unit, colnum, row, 1, nwave, nullval, &
         vis2err, anyf, status)
    call ftgcno(unit, .false., 'UCOORD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         uv_array, anyf, status)
    ucoord = uv_array(1)
    call ftgcno(unit, .false., 'VCOORD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         uv_array, anyf, status)
    vcoord = uv_array(1)
    call ftgcno(unit, .false., 'STA_INDEX', colnum, status)
    call ftgcvj(unit, colnum, row, 1, 2, nullint, &
         sta_index, anyf, status)
    call ftgcno(unit, .false., 'FLAG', colnum, status)
    call ftgcl(unit, colnum, row, 1, nwave, flag, status)

  end subroutine read_oi_vis2

  !============================================================================

  subroutine read_oi_t3(unit, row, nwave, target_id, mjd, &
       t3amp, t3amperr, t3phi, t3phierr, &
       u1coord, v1coord, u2coord, v2coord, sta_index, flag, status)

    !Read specified row of OI_T3 FITS binary table using fitsio calls
    !OI_T3 must be current HDU

    integer, intent(in) :: unit, row, nwave
    double precision, dimension(nwave), intent(out) :: t3amp, t3amperr
    double precision, dimension(nwave), intent(out) :: t3phi, t3phierr
    double precision, intent(out) :: mjd, u1coord, v1coord, u2coord, v2coord
    integer, intent(out) :: target_id
    integer, dimension(3), intent(out) :: sta_index
    logical, dimension(nwave), intent(out) :: flag
    integer, intent(inout) :: status

    !local variables
    integer, dimension(1) :: target_id_array
    double precision, dimension(1) :: mjd_array, uv_array
    integer colnum, iwave, nullint
    double precision nullval
    logical anyf, anynullamp
    logical, dimension(:), allocatable :: isnull

    call ftgcno(unit, .false., 'TARGET_ID', colnum, status)
    call ftgcvj(unit, colnum, row, 1, 1, nullint, &
         target_id_array, anyf, status)
    target_id = target_id_array(1)
    call ftgcno(unit, .false., 'MJD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         mjd_array, anyf, status)
    mjd = mjd_array(1)
    allocate(isnull(nwave))
    call ftgcno(unit, .false., 'T3AMP', colnum, status)
    call ftgcfd(unit, colnum, row, 1, nwave, &
         t3amp, isnull, anynullamp, status)
    call ftgcno(unit, .false., 'T3AMPERR', colnum, status)
    call ftgcvd(unit, colnum, row, 1, nwave, nullval, &
         t3amperr, anyf, status)
    ! make errors for NULL T3AMP values -ve, put dummy values in for T3AMP
    if (anynullamp) then
       do iwave = 1, nwave
          if (isnull(iwave)) then
             t3amperr(iwave) = -t3amperr(iwave) !zero error => ignored also
             t3amp(iwave) = 1.0D0
          end if
       end do
    end if
    call ftgcno(unit, .false., 'T3PHI', colnum, status)
    call ftgcvd(unit, colnum, row, 1, nwave, nullval, &
         t3phi, anyf, status)
    call ftgcno(unit, .false., 'T3PHIERR', colnum, status)
    call ftgcvd(unit, colnum, row, 1, nwave, nullval, &
         t3phierr, anyf, status)
    call ftgcno(unit, .false., 'U1COORD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         uv_array, anyf, status)
    u1coord = uv_array(1)
    call ftgcno(unit, .false., 'V1COORD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         uv_array, anyf, status)
    v1coord = uv_array(1)
    call ftgcno(unit, .false., 'U2COORD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         uv_array, anyf, status)
    u2coord = uv_array(1)
    call ftgcno(unit, .false., 'V2COORD', colnum, status)
    call ftgcvd(unit, colnum, row, 1, 1, nullval, &
         uv_array, anyf, status)
    v2coord = uv_array(1)
    call ftgcno(unit, .false., 'STA_INDEX', colnum, status)
    call ftgcvj(unit, colnum, row, 1, 3, nullint, &
        sta_index, anyf, status)
    call ftgcno(unit, .false., 'FLAG', colnum, status)
    call ftgcl(unit, colnum, row, 1, nwave, flag, status)
    deallocate(isnull)

  end subroutine read_oi_t3

  !============================================================================

  function DateToMjd(year, month, day)

    !Return Modified Julian Day at 0h UT on specified date.
    !
    !Two-digit years are not supported: 99 means the year 99 AD
    !month, day start at 1
    !
    !Code translated from routine julday in Numerical Recipes in C, Section 1.1

    !function arguments
    integer, intent(in) :: year, month, day
    double precision :: DateToMjd

    !local variables
    integer :: ja, jy, jm
    integer*8 :: jul
    !Gregorian calendar adopted 15 Oct 1582:
    integer*8, parameter :: greg = 15 + 31*(10+12*1582)

    if (year == 0) then
       print *, 'DateToMjd: There is no year zero.'
       return
    end if
    if (year < 0) then
       jy = year + 1
    else
       jy = year
    end if
    if (month > 2) then
       jm = month + 1
    else
       jy = jy - 1
       jm = month + 13
    end if
    jul = floor(365.25*jy) + floor(30.6001*jm) + day + 1720995
    if (day + 31*(month+12*year) >= greg) then
       ja = int(0.01*year)
       jul = jul + 2 - ja + int(0.25*ja)
    end if
    !jul is Julian Day for noon on given date
    DateToMjd = jul - 2400001

  end function DateToMjd

  !============================================================================

  function yesno(prompt, default)

    !Get yes/no response
    character(len=*), intent(in) :: prompt, default
    logical yesno
    character(len=128) :: response

    print *, trim(prompt), ' [', default, '] ?'
    read (*, '(a)') response
    if (len_trim(response) == 0) then
       response = default
       print '(a)', trim(response)
    end if
    if (response(1:1) == 'y' .or. response(1:1) == 'Y') then
       yesno = .true.
    else
       yesno = .false.
    end if

  end function yesno

  !============================================================================

end module Inout








