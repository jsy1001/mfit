!$Id: clnest.f90,v 1.3 2009/07/16 13:27:11 jsy1001 Exp $

program Main

  use f2kcli
  use Inout
  use Model
  use Bayes
  use Wrap
  use nestwrapper
  use ReadMC

  implicit none
  
  !============================================================================
  !variables
  
  !parameters
  !pi, deg/rad conversions all picked up from Maths module
  integer, parameter :: max_lines = 5000      !max. lines in data file
  integer, parameter :: width = 78            !for spacer lines
  double precision, parameter :: sig = 0.1D0  !waveband must match to sig nm

  !other local variables
  double precision, allocatable :: wavebands(:,:)
  double precision :: wb(2), wl(2)
  character(len=width) :: spacer_line
  character(len=128) :: switch, arg, info, file_name, ext, source
  character(len=128) :: cvs_rev, revision
  integer :: narg, iarg, i, n, user_target_id
  integer :: useful_vis, useful_amp, useful_cp, mode, degfreedom
  double precision :: chisqrd, normchisqrd
  double precision :: calib_error
  logical :: mmodal, force_symm
  type(allparam) :: mean_param
  
  !posterior samples
  integer nsamp
  character(len=128) :: sampFilename
  double precision, allocatable :: sol(:), err(:)
  integer, allocatable :: var_pos(:,:)
  character(len=model_desc_len), allocatable :: var_desc(:)

  !----------------------------------------------------------------------------
  !Formatting stuff
  spacer_line = repeat('-', width)

  !----------------------------------------------------------------------------
  !Introduction

  cvs_rev = '$Revision: 1.3 $'
  revision = cvs_rev(scan(cvs_rev, ':')+2:scan(cvs_rev, '$', .true.)-1)
  print *,' '
  print *,spacer_line
  print *,' '
  print *,'  clnest - command-line OI nested sampler'
  print *,'  package release ',release
  print *,'  [clnest revision ',trim(revision),']'
  print *,' '
  print *,spacer_line

  !----------------------------------------------------------------------------
  !Parse command-line arguments

  narg = command_argument_count()
  iarg = 1
  !defaults
  wb = (/-1.0D0, -1.0D0/)
  wl = (/-1.0D0, -1.0D0/)
  user_target_id = -1 !default to 1st target in OI_TARGET
  calib_error = 0D0
  mmodal = .false.
  mode = -1
  do
     if (iarg == narg - 1) exit
     if (iarg > narg - 1) then
        call print_usage()
        stop 'Not enough arguments'
     end if
     call get_command_argument(iarg, switch)
     select case(switch)
        case('-c', '--calerr')
           call get_command_argument(iarg+1, arg)
           read(arg, *) calib_error
           iarg = iarg + 1
        case('-w', '--waveband')
           call get_command_argument(iarg+1, arg)
           read(arg, *) wb(1)
           call get_command_argument(iarg+2, arg)
           read(arg, *) wb(2)
           iarg = iarg + 2
        case('-r', '--waverange')
           call get_command_argument(iarg+1, arg)
           read(arg, *) wl(1)
           call get_command_argument(iarg+2, arg)
           read(arg, *) wl(2)
           iarg = iarg + 2
        case('-t', '--target_id')
           call get_command_argument(iarg+1, arg)
           read(arg, *) user_target_id
           iarg = iarg + 1
        case('-m', '--multimodal')
           mmodal = .true.
        case('-n', '--mode')
           call get_command_argument(iarg+1, arg)
           read(arg, *) mode
           iarg = iarg + 1
        case default
           print *, 'Ignoring invalid command-line argument: ', arg
        end select
        iarg = iarg + 1
  end do

  !----------------------------------------------------------------------------
  !Read vis/triple data
  !
  !End up with vis_data and triple_data arrays. 
  !
  !vis_data holds the visibility data points of lambda, u, v, vis, vis_err, 
  !where vis is the squared visibility amplitude and vis_err is its abs error
  !
  !triple_data holds the triple product data points of lambda, u1,v1,u2,v2,
  ! amp, amp_err, cp, cp_err
  !where amp is the triple product amplitude (NOT its square) and amp_err is
  !the absolute error. cp is the closure phase (deg) and cp_err is the absolute
  !error.
  !
  !Negative or zero absolute errors are stored if the mapdat/vis file contains
  !negative or zero entries for its errors (defined differently). These points
  !are retained in the file for potential plotting purposes but are completely 
  !ignored in fitting (i.e. ignored in likelihood calculation and goodness of
  !fit calculations).

  call get_command_argument(narg-1, file_name)
  ext = trim(file_name(scan(file_name,'.',.true.)+1:len(file_name)))

  !call appropriate reading routine for the file type
  info = ''

  if (ext == 'mapdat') then
     !Can be non-centrosymmetric (complex model visibilities)
     !.mapdat can have visibility and triple product measurements
     !Can contain different waveband observations, these are all read into
     !data arrays in case multi-waveband functionality is added later

     !read_mapdat allocates vis_data, triple_data, and wavebands
     call read_mapdat(info, file_name, source, max_lines, &
          vis_data, num_vis, triple_data, num_triple, wavebands, calib_error)
     force_symm = .false.

  else if (ext == 'vis') then
     !Model is forced to be centrosymmetric (real visibilities)
     !.vis only has visibility measurements (no errors) and
     !projected baselines sqrt(u**2 + v**2), stored in u column

     if (wb(1) .eq. -1.0D0) &
          stop 'need observing wavelength and bandwidth (nm) for vis-format data'
     print '(1x, a, 2f8.2, a)', 'Assuming waveband for data is', wb, ' nm'
     !read_vis allocates vis_data
     call read_vis(info, file_name, source, max_lines, vis_data, num_vis, &
          wb, calib_error)
     allocate(triple_data(0, 0))
     num_triple = 0
     allocate(wavebands(1, 2))
     wavebands(1, :) = wb
     force_symm = .true.

  else if (ext == 'nvis') then
     !Model is forced to be centrosymmetric (real visibilities)
     !.nvis only has visibility measurements & errors, plus
     !projected baselines sqrt(u**2 + v**2), stored in u column
     if (wb(1) .eq. -1.0D0) &
          stop 'need observing wavelength and bandwidth (nm) for nvis-format data'
     print '(1x, a, 2f8.2, a)', 'Assuming waveband for data is', wb, ' nm'
     !read_nvis allocates vis_data
     call read_nvis(info, file_name, source, max_lines, vis_data, num_vis, &
          wb, calib_error)
     allocate(triple_data(0, 0))
     num_triple = 0
     allocate(wavebands(1, 2))
     wavebands(1, :) = wb
     force_symm = .true.

  else if (ext(len_trim(ext)-3:len_trim(ext)) == 'fits') then
     !read_oi_fits allocates vis_data, triple_data, and wavebands
     call read_oi_fits(info, file_name, user_target_id, source, &
          vis_data, num_vis, triple_data, num_triple, wavebands, calib_error)
     force_symm = .false.

  else if (ext(len_trim(ext)-4:len_trim(ext)) == 'calib') then
     !.*Calib has no waveband data so must be supplied by wavelength range
     if (wb(2) .eq. -1.0D0) &
          stop 'need waveband for wbCalib format'
     print '(1x, a, f8.2, a)', 'Assuming waveband for data is', wb(2), ' nm'
     call read_calib(info, file_name, source, max_lines, vis_data, num_vis, &
          wb(2), wavebands, calib_error)
     allocate(triple_data(0, 0))
     force_symm = .false.

  else
     info = 'file type "'//trim(ext)//'" not handled'

  end if
  if (info /= '') then
     print *, trim(info)
     stop
  end if

  !Filter here to reduce vis and triple data to chosen waveband(s)
  if (size(wavebands,1) > 1) then
     print *,' '
     print *,'multiple waveband data found with the following wavelengths'
     do i = 1, size(wavebands,1)
        print 59, 'waveband', i, 'wavelength (nm)', real(wavebands(i, 1)), &
             'bandwidth (nm)', real(wavebands(i, 2))
59      format(a, 1x, i3, 1x, a, 1x, f8.2, 1x, a, 1x, f7.2)
     end do

     if (wb(1) /= -1.0D0) then
        call filt_by_wb(vis_data, wb, sig, num_vis)
        call filt_by_wb(triple_data, wb, sig, num_triple)
        num_wb = 1
        allocate(sel_wavebands(num_wb, 2))
        sel_wavebands(1, :) = wb
        print *, ' '
        print '(1x, a, 1x, 2f8.2)', 'Using specified waveband:', wb
     else if (wl(1) /= -1.0D0) then
        call filt_by_wl(vis_data, wl(1), wl(2), num_vis)
        call filt_by_wl(triple_data, wl(1), wl(2), num_triple)
        allocate(sel_wavebands(size(wavebands, 1), 2))
        sel_wavebands = wavebands
        call filt_by_wl(sel_wavebands, wl(1), wl(2), num_wb)
        print *, ' '
        print '(1x, a, 1x, 2f8.2)', 'Using specified wavelength range:', wl
     else
        !use all wavebands present
        num_wb = size(wavebands, 1)
        allocate(sel_wavebands(num_wb, 2))
        sel_wavebands = wavebands
     end if

  else
     allocate(sel_wavebands(1, 2))
     sel_wavebands = wavebands
  end if

  !Count number of useful data points (i.e. with +ve errors, points with
  !-ve or zero errors are ignored in fitting but are held in data arrays
  !for plotting purposes so mustn't count them)
  useful_vis = 0
  useful_amp = 0
  useful_cp = 0
  do i = 1, num_vis
     if (vis_data(i,6) > 0D0) useful_vis = useful_vis + 1
  end do
  do i = 1, num_triple
     if (triple_data(i,8) > 0D0) useful_amp = useful_amp + 1
     if (triple_data(i,10) > 0D0) useful_cp = useful_cp + 1
  end do
  if ((useful_vis + useful_amp + useful_cp) == 0) stop 'No unflagged data'

  print *,' '
  print *, trim(ext), ' file visibility data for ', trim(source),':'

  print *, ' '
  print *, '    MJD   wave   band     baseline coords     sqrd      abs'
  print *, '        length  width         u         v      vis    error'
  print *, '          (nm)   (nm)       (m)       (m)'  
  do i = 1, num_vis
     write(*,60) vis_data(i,7), vis_data(i,:6)
  end do
60 format(f8.2, f7.1, 1x, f6.1, 1x, f9.4, 1x, f9.4, 1x, f8.5, 1x, f8.5)

  print *, ' '
  print *, trim(ext),' file triple product data for ',trim(source),':'
  print *, ' '
  print '(1x, 2a)', &
       '    MJD   wave  band                baseline triangle coords', &
       '                      triple product'
  print '(1x, 2a)', &
       '        length width        u1        v1        u2        v2', &
       '  amplitude     error   phase    err'
  print *, '          (nm)  (nm)       (m)       (m)       (m)       (m)'
  do i = 1, num_triple
     write(*,61) triple_data(i,11), triple_data(i,:10)
  end do
61 format(f8.2, f7.1, 1x, f5.1, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, &
        e10.3, 1x, e9.3, 1x, f7.2, 1x, f6.2, 1x)

  print *, ' '
  print *, useful_vis, 'unflagged visibility data points out of', num_vis
  print *, useful_amp, 'unflagged triple product amplitudes out of', num_triple
  print *, useful_cp, 'unflagged closure phases out of', num_triple
  print *, ' '
  print *, 'total of', useful_vis+useful_amp+useful_cp, &
       'unflagged data points out of', num_vis+2*num_triple
  print *, ' '
  print *, spacer_line

  !-------------------------------------------------------------------------
  !Read model
  !
  !Model file describing model as per documentation is read in.
  !Checks are made to ensure parameters lie within acceptable ranges
  !(defined by the model_limits array), the same limits are used later in
  !the fitting routines (refer to them).
  !Free parameters should be supplied with (non zero and positive)
  !prior widths which for clnest *ONLY* are interpreted as the
  !half-width of a uniform prior, centred on the specified parameter value
  !(which is otherwise unused).

  call get_command_argument(narg, file_name)
  print *, 'reading model...'
  info = ''
  call read_model(info, file_name, sel_wavebands)
  if (info /= '') then
     print *, trim(info)
     stop
  end if

  print *, '...done'
  print *, ' '
  print *, 'model details:'
  call print_model()
  print *, ' '
  print *, spacer_line

  !-------------------------------------------------------------------------
  !run nested sampler

  print *, ' '
  print *, 'running sampler...'

  !check model freedoms
  if (.not. model_valid(info, force_symm)) then
     print *, trim(info)
     stop
  end if

  !run nest sampler
  call model_nvar(n)
  call nest_Sample(n, mmodal)

  !display results
  if(mode > 0) then
     sampFilename = trim(nest_root)//'post_separate.dat'
  else
     sampFilename = trim(nest_root)//'.txt'
  end if
  call mc_count_lines(sampFilename, mode, nsamp)
  print '(1x, a, i6, a, a)', &
       'Reading', nsamp, ' samples from ', sampFilename
  print *, 'Mean results are in ', trim(nest_root)//'stats.dat'
  allocate(sol(n), err(n))
  allocate(var_pos(n,2), var_desc(n))
  call mc_get_params(sampFilename, nest_root, nsamp, n, mode, sol, err)
  call model_getvar(n, var_pos, var_desc)
  print '(1x, a, 49x, a)', '                   ', '  mean     standard'
  print '(1x, a, 49x, a)', 'num  parameter name', ' value    deviation'
  do i = 1, n
     write(*,62) i, var_desc(i), sol(i), err(i)
  end do
62 format(' (', i2, ') ', A55, 1x, f13.6, 1x, f12.6) 

  !Display chi-squared
  call allparam_init(mean_param, model_param, model_limits, n, var_pos)
  call allparam_setvar(mean_param, sol)
  chisqrd = 2d0*likelihood(vis_data, triple_data, model_spec, mean_param%param)
  call allparam_free(mean_param)
  degfreedom = useful_vis + useful_amp + useful_cp - n
  print *, ' '
  print *, '           chi squared =',real(chisqrd) 
  print *, '    degrees of freedom =',degfreedom
  if (degfreedom > 0) then
     normchisqrd = chisqrd/degfreedom
     print *, 'chi sqrd / deg freedom =',real(normchisqrd)
  end if
  

  !-------------------------------------------------------------------------
  !Deallocate storage
  call free_model() !model_*

  if (allocated(vis_data)) deallocate(vis_data)
  if (allocated(triple_data)) deallocate(triple_data)
  if (allocated(wavebands)) deallocate(wavebands)
  if (allocated(sel_wavebands)) deallocate(sel_wavebands)

  if (allocated(sol)) deallocate(sol)
  if (allocated(err)) deallocate(err)
  if (allocated(var_pos)) deallocate(var_pos)
  if (allocated(var_desc)) deallocate(var_desc)

contains

!==============================================================================

  subroutine filt_by_wb(data, wb, sig, new_dim1)

    !Filter 2d double precision data array by 1st two columns
    !Rows are kept if 1st two columns match wb to precision sig
    !Reallocates data array
    double precision, allocatable, intent(inout) :: data(:,:)
    double precision, intent(in) :: wb(2)
    double precision, intent(in) :: sig
    integer, intent(out) :: new_dim1

    !local variables
    integer :: dim2
    logical, allocatable :: mask(:)
    double precision, allocatable :: filt_data(:,:)

    dim2 = size(data,2)
    allocate(mask(size(data,1)))
    mask = (data(:, 1) >= wb(1)-sig .and. data(:, 1) <= wb(1)+sig &
         .and. data(:, 2) >= wb(2)-sig .and. data(:, 2) <= wb(2)+sig)
    new_dim1 = count(mask)
    allocate(filt_data(new_dim1, dim2))
    filt_data = reshape(pack(data, spread(mask, 2, dim2)), (/new_dim1, dim2/))
    deallocate(data)
    allocate(data(new_dim1, dim2))
    data = filt_data
    new_dim1 = size(data,1)
    deallocate(filt_data)
    deallocate(mask)

  end subroutine filt_by_wb

!==============================================================================

  subroutine filt_by_wl(data, wlmin, wlmax, new_dim1)

    !Filter 2d double precision data array by 1st column
    !Rows are kept if 1st column between wlmin and wlmax
    !Reallocates data array
    double precision, allocatable, intent(inout) :: data(:,:)
    double precision, intent(in) :: wlmin, wlmax
    integer, intent(out) :: new_dim1

    !local variables
    integer :: dim2
    logical, allocatable :: mask(:)
    double precision, allocatable :: filt_data(:, :)

    dim2 = size(data,2)
    allocate(mask(size(data,1)))
    mask = (data(:, 1) >= wlmin .and. data(:, 1) <= wlmax)
    new_dim1 = count(mask)
    allocate(filt_data(new_dim1, dim2))
    filt_data = reshape(pack(data, spread(mask, 2, dim2)), (/new_dim1, dim2/))
    deallocate(data)
    allocate(data(new_dim1, dim2))
    data = filt_data
    deallocate(filt_data)
    deallocate(mask)

  end subroutine filt_by_wl

!==============================================================================

  subroutine print_usage()

    print *, 'Usage:'
    print *, ' '
    print *, 'clnest [options] datafile modelfile'
    print *, ' '
    print *, 'Options (may appear in any order):'
    print *, ' '
    print *, '-w|--waveband CWL BW    select this waveband;'
    print *, '                        must specify wb this way for .(n)vis data'
    print *, '-r|--waverange WL1 WL2  select all wavebands in this range'
    print *, '-t|--target_id ID       select this target (default is 1st in OI_TARGET table)'
    print *, '-c|--calerr FRAC        add calibration error (frac. err. in system visib.)'
    print *, '-m|--multimodal         do multimodal sampling'
    print *, '-n|--mode               process only samples for this mode'

  end subroutine print_usage

end program Main
