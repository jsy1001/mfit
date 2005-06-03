!$Id: clfit.f90,v 1.17 2005/06/03 09:47:57 jsy1001 Exp $

program Main

  use Inout
  use Plot
  use Visibility
  use Fit
  use Model
  use f2kcli

  implicit none
  
  !============================================================================
  !variables
  
  !parameters
  !pi, deg/rad conversions all picked up from Maths module
  integer, parameter :: max_lines = 5000      !max. lines in data file
  integer, parameter :: width = 78            !for spacer lines
  double precision, parameter :: sig = 0.1D0  !waveband must match to sig nm

  !arrays for fit results
  character(len=55), dimension(:), allocatable :: desc
  double precision, dimension(:, :), allocatable :: sol, hes, cov, cor

  !other local variables
  double precision, dimension(:, :), allocatable :: wavebands
  double precision, dimension(2) :: wb, wl
  character(len=width) :: spacer_line
  character(len=128) :: switch, arg, device, info, file_name, ext, source
  character(len=128) :: xlabel, top_title, cvs_rev, revision
  character(len=8) :: sel_plot
  integer :: narg, iarg, i, j, n, length, flag, index, user_target_id
  integer :: degfreedom, useful_vis, useful_amp, useful_cp
  double precision :: nlposterior, nlevidence, chisqrd, normchisqrd
  double precision :: calib_error, uxmin, uxmax, x0, sigma
  real :: xzero
  logical :: force_symm, nofit, zoom, mod_line

  integer :: pgopen, istat

  !----------------------------------------------------------------------------
  !Formatting stuff
  spacer_line = repeat('-', width)

  !----------------------------------------------------------------------------
  !Introduction

  cvs_rev = '$Revision: 1.17 $'
  revision = cvs_rev(scan(cvs_rev, ':')+2:scan(cvs_rev, '$', .true.)-1)
  print *,' '
  print *,spacer_line
  print *,' '
  print *,'  clfit - command-line mfit'
  print *,'  package release ',release
  print *,'  [clfit revision ',trim(revision),']'
  print *,' '
  print *,spacer_line

  !----------------------------------------------------------------------------
  !Parse command-line arguments

  narg = command_argument_count()
  iarg = 1
  !defaults
  sel_plot = '' !no plot
  device = '/xserv'
  nofit = .false.
  zoom = .false.
  wb = (/-1.0D0, -1.0D0/)
  wl = (/-1.0D0, -1.0D0/)
  user_target_id = -1 !default to 1st target in OI_TARGET
  calib_error = 0D0
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
        case('-d', '--device')
           call get_command_argument(iarg+1, device)
           iarg = iarg + 1
        case('-n', '--nofit')
           nofit = .true.
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
        case('-p', '--plot')
           call get_command_argument(iarg+1, sel_plot)
           iarg = iarg + 1
           if (sel_plot == 'post') then
              call get_command_argument(iarg+1, arg)
              read(arg, *) index
              iarg = iarg + 1
           end if
        case('-z', '--zoomplot')
           zoom = .true.
           call get_command_argument(iarg+1, sel_plot)
           if (sel_plot == 'post') then
              call get_command_argument(iarg+2, arg)
              read(arg, *) index
              iarg = iarg + 1
           end if
           call get_command_argument(iarg+2, arg)
           read(arg, *) uxmin
           call get_command_argument(iarg+3, arg)
           read(arg, *) uxmax
           iarg = iarg + 3
        case('-t', '--target_id')
           call get_command_argument(iarg+1, arg)
           read(arg, *) user_target_id
           iarg = iarg + 1
        case default
           print *, 'Ignoring invalid command-line argument: ', arg
        end select
        iarg = iarg + 1
  end do
  ! If plot specified, open plot device
  if (sel_plot /= '') istat = pgopen(device)

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
          vis_data, triple_data, wavebands, calib_error)
     force_symm = .false.

  else if (ext == 'vis') then
     !Model is forced to be centrosymmetric (real visibilities)
     !.vis only has visibility measurements (no errors) and
     !projected baselines sqrt(u**2 + v**2), stored in u column

     if (wb(1) .eq. -1.0D0) &
          stop 'need observing wavelength and bandwidth (nm) for vis-format data'
     print '(1x, a, 2f8.2, a)', 'Assuming waveband for data is', wb, ' nm'
     !read_vis allocates vis_data
     call read_vis(info, file_name, source, max_lines, vis_data, wb, &
          calib_error)
     allocate(triple_data(0, 0))
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
     call read_nvis(info, file_name, source, max_lines, vis_data, wb, &
          calib_error)
     allocate(triple_data(0, 0))
     allocate(wavebands(1, 2))
     wavebands(1, :) = wb
     force_symm = .true.

  else if (ext(len_trim(ext)-3:len_trim(ext)) == 'fits') then
     !read_oi_fits allocates vis_data, triple_data, and wavebands
     call read_oi_fits(info, file_name, user_target_id, source, &
          vis_data, triple_data, wavebands, calib_error)
     force_symm = .false.

  else
     info = 'file type "'//trim(ext)//'" not handled'

  end if
  if (info /= '') then
     print *, info
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
        call filt_by_wb(vis_data, wb, sig)
        call filt_by_wb(triple_data, wb, sig)
        allocate(sel_wavebands(1, 2))
        sel_wavebands(1, :) = wb
        print *, ' '
        print '(1x, a, 1x, 2f8.2)', 'Using specified waveband:', wb
     else if (wl(1) /= -1.0D0) then
        call filt_by_wl(vis_data, wl(1), wl(2))
        call filt_by_wl(triple_data, wl(1), wl(2))
        allocate(sel_wavebands(size(wavebands, 1), 2))
        sel_wavebands = wavebands
        call filt_by_wl(sel_wavebands, wl(1), wl(2))
        print *, ' '
        print '(1x, a, 1x, 2f8.2)', 'Using specified wavelength range:', wl
     else
        !use all wavebands present
        allocate(sel_wavebands(size(wavebands, 1), 2))
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
  do i = 1, size(vis_data,1)
     if (vis_data(i,6) > 0D0) useful_vis = useful_vis + 1
  end do
  do i = 1, size(triple_data,1)
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
  do i = 1, size(vis_data,1)
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
  do i = 1, size(triple_data,1)
     write(*,61) triple_data(i,11), triple_data(i,:10)
  end do
61 format(f8.2, f7.1, 1x, f5.1, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, &
        e10.3, 1x, e9.3, 1x, f7.2, 1x, f6.2, 1x)

  print *, ' '
  print *, useful_vis, 'unflagged visibility data points out of', &
       size(vis_data,1)
  print *, useful_amp, 'unflagged triple product amplitudes out of', &
       size(triple_data,1)
  print *, useful_cp, 'unflagged closure phases out of', size(triple_data,1)
  print *, ' '
  print *, 'total of', useful_vis+useful_amp+useful_cp, &
       'unflagged data points out of', size(vis_data,1)+(2*size(triple_data,1))
  print *, ' '
  print *, spacer_line

  if (sel_plot == 'uv' .and. ext /= 'vis' .and. ext /= 'nvis' &
       .and. useful_vis > 0) &
       call plot_uv('u /M\gl', 'v /M\gl', trim(source))

  !-------------------------------------------------------------------------
  !Read model data
  !
  !Model file describing model as per documentation is read in.
  !Checks are made to ensure parameters lie within acceptable ranges
  !(defined by the model_limits array), the same limits are used later in
  !the fitting routines (refer to them).
  !Free parameters should be supplied with (non zero and positive) prior
  !widths which are the 1-sigma width of the gaussian prior
  !distributions as per DB thesis chapter 2.
  !Checks are not made here as to the legality of the freedom in the model -
  !refer to the fit module.

  call get_command_argument(narg, file_name)
  print *, 'reading model...'
  info = ''
  call read_model(info, file_name, sel_wavebands)
  if (info /= '') then
     print *, info
     stop
  end if

  print *, '...done'
  print *, ' '
  print *, 'model details:'
  call print_model()
  print *, ' '
  print *, spacer_line

  if (nofit) then
     !-------------------------------------------------------------------------
     ! report chi-squared, plot initial model
     call gof(model_spec, model_param, chisqrd)
     !go through model arrays and count number of variable parameters
     n = 0
     do i = 1, size(model_param,1)
        do j = 1, size(model_param,2)
           if (model_prior(i,j) /= 0D0) n=n+1
        end do
     end do
     degfreedom = useful_vis + useful_amp + useful_cp - n
     print *, ' '
     print *, '           chi squared =',real(chisqrd) 
     print *, '    degrees of freedom =',degfreedom
     if (degfreedom > 0) then
        normchisqrd = chisqrd/degfreedom
        print *, 'chi sqrd / deg freedom =',real(normchisqrd)
     end if
     mod_line = (symm .and. size(sel_wavebands, 1) == 1)
     top_title = trim(source)//' - initial model: '//trim(model_name)
     if (sel_plot == 'vis2' .and. useful_vis > 0) then
        if (zoom) then 
           call plot_vis_bas(model_spec, model_param, mod_line, &
                'Baseline /M\gl', 'Squared visibility', top_title, &
                uxmin, uxmax)
        else
           call plot_vis_bas(model_spec, model_param, mod_line, &
                'Baseline /M\gl', 'Squared visibility', top_title)
        end if
     end if
     if (sel_plot == 'vis2wl' .and. useful_vis > 0) then
        if (zoom) then 
           call plot_vis(1, model_spec, model_param, &
                'Wavelength /nm', 'Squared visibility', top_title, &
                0., uxmin, uxmax)
        else
           call plot_vis(1, model_spec, model_param, &
                'Wavelength /nm', 'Squared visibility', top_title, 0.)
        end if
     end if
     if (sel_plot == 'vis2mjd' .and. useful_vis > 0) then
        xzero = floor(minval(vis_data(:,7)))
        write (xlabel, *) 'Modified Julian Day -', xzero
        if (zoom) then 
           call plot_vis(7, model_spec, model_param, &
                xlabel, 'Squared visibility', top_title, xzero, uxmin, uxmax)
        else
           call plot_vis(7, model_spec, model_param, &
                xlabel, 'Squared visibility', top_title, xzero)
        end if
     end if
     if (sel_plot == 't3amp' .and. useful_amp > 0) then
        if (zoom) then 
           call plot_triple_amp_bas(model_spec, model_param, &
                'Longest baseline /M\gl', 'Triple amplitude', top_title, &
                uxmin, uxmax)
        else
           call plot_triple_amp_bas(model_spec, model_param, &
                'Longest baseline /M\gl', 'Triple amplitude', top_title)
        end if
     end if
     if (sel_plot == 't3ampwl' .and. useful_amp > 0) then
        if (zoom) then 
           call plot_triple_amp(1, model_spec, model_param, &
                'Wavelength /nm', &
                'Triple amplitude', top_title, 0., uxmin, uxmax)
        else
           call plot_triple_amp(1, model_spec, model_param, &
                'Wavelength /nm', &
                'Triple amplitude', top_title, 0.)
        end if
     end if
     if (sel_plot == 't3ampmjd' .and. useful_amp > 0) then
        xzero = floor(minval(triple_data(:,11)))
        write (xlabel, *) 'Modified Julian Day -', xzero
        if (zoom) then 
           call plot_triple_amp(11, model_spec, model_param, &
                xlabel, 'Triple amplitude', top_title, xzero, uxmin, uxmax)
        else
           call plot_triple_amp(11, model_spec, model_param, &
                xlabel, 'Triple amplitude', top_title, xzero)
        end if
     end if
     if (sel_plot == 't3phi' .and. useful_cp > 0) then
        if (zoom) then 
           call plot_triple_phase_bas(model_spec, model_param, &
                'Longest baseline /M\gl', 'Closure phase /'//char(176), &
                top_title, uxmin, uxmax)
        else
           call plot_triple_phase_bas(model_spec, model_param, &
                'Longest baseline /M\gl', 'Closure phase /'//char(176), &
                top_title)
        end if
     end if
     if (sel_plot == 't3phiwl' .and. useful_cp > 0) then
        if (zoom) then 
           call plot_triple_phase(1, model_spec, model_param, &
                'Wavelength /nm', &
                'Closure phase /'//char(176), top_title, 0., uxmin, uxmax)
        else
           call plot_triple_phase(1, model_spec, model_param, &
                'Wavelength /nm', &
                'Closure phase /'//char(176), top_title, 0.)
        end if
     end if
     if (sel_plot == 't3phimjd' .and. useful_cp > 0) then
        xzero = floor(minval(triple_data(:,11)))
        write (xlabel, *) 'Modified Julian Day -', xzero
        if (zoom) then 
           call plot_triple_phase(11, model_spec, model_param, &
                xlabel, 'Closure phase /'//char(176), top_title, xzero, &
                uxmin, uxmax)
        else
           call plot_triple_phase(11, model_spec, model_param, &
                xlabel, 'Closure phase /'//char(176), top_title, xzero)
        end if
     end if
     !cannot do 'post' plot if --nofit
  else
     !-------------------------------------------------------------------------
     !fit model to the data by minimising negative log posterior
     !
     !Refer to fit module. the fit solution is returned along with various 
     !diagnostic quantities.
     print *, ' '
     print *, 'fitting model by minimising negative log posterior...'

     info = ''
     ! minimiser allocates fit_param, x, x_pos, sol, desc, hes, cov, cor
     call minimiser(info, force_symm, sol, flag, desc, &
          hes, cov, cor, chisqrd, nlposterior, nlevidence)
     if (info /= '') then
        print *,'*****'
        print *,trim(info)
        print *, '*****'
     end if

     if (flag > -1 .and. flag < 4) then
        degfreedom = useful_vis + useful_amp + useful_cp - size(sol,1)

        print *, ' '
        print *, 'negative log posterior =',real(nlposterior)
        print *, 'negative log evidence  =',real(nlevidence)
        print *, ' '
        print *, '           chi squared =',real(chisqrd) 
        print *, '    degrees of freedom =',degfreedom
        if (degfreedom > 0) then
           normchisqrd = chisqrd/degfreedom
           print *, 'chi sqrd / deg freedom =',real(normchisqrd)
        end if

        print *, ' '
        print *, 'solution details:'
        print '(1x, a, 49x, a)', '                   ', 'fitted      hessian'
        print '(1x, a, 49x, a)', 'num  parameter name', ' value        error'
        do i = 1, size(sol,1)
           write(*,62) i, desc(i), sol(i,:)
        end do
62      format(' (', i2, ') ', A55, 1x, f13.6, 1x, f12.6) 

        ! Display alternative error bars for fitted parameters
        ! Estimates assuming data errors scaled so chi sqrd/deg freedom = unity
        print *, ' '
        print '(1x, a, 49x, a)', 'SCALED DATA ERRORS ->', 'fitted      hessian'
        print '(1x, a, 49x, a)', '  num  parameter name', ' value        error'
        do i = 1, size(sol,1)
           write(*,63) i, desc(i), sol(i,1), sqrt(normchisqrd)*sol(i,2)
        end do
63      format('  *(', i2, ') ', A55, 1x, f13.6, 1x, f12.6, '*') 

        length = size(hes,1)
        print *, ' '
        print *, 'hessian matrix'
        do i = 1, length
           write(*,64) hes(i,:)
        end do

        print *, ' '
        print *, 'covariance matrix'
        do i = 1, length
           write(*,64) cov(i,:)
        end do

        print *, ' '
        print *, 'correlation matrix'
        do i = 1, length
           write(*,64) cor(i,:)
        end do

64      format(1x, 10(e11.4,1x))

        !----------------------------------------------------------------------
        !plot
        mod_line = (symm .and. size(sel_wavebands, 1) == 1)
        top_title = trim(source)//' - final model: '//trim(model_name)
        if (sel_plot == 'vis2' .and. useful_vis > 0) then
           if (zoom) then 
              call plot_vis_bas(model_spec, fit_param, mod_line, &
                   'Baseline /M\gl', 'Squared visibility', top_title, &
                   uxmin, uxmax)
           else
              call plot_vis_bas(model_spec, fit_param, mod_line, &
                   'Baseline /M\gl', 'Squared visibility', top_title)
           end if
        end if
        if (sel_plot == 'vis2wl' .and. useful_vis > 0) then
           if (zoom) then 
              call plot_vis(1, model_spec, fit_param, &
                   'Wavelength /nm', 'Squared visibility', top_title, &
                   0., uxmin, uxmax)
           else
              call plot_vis(1, model_spec, fit_param, &
                   'Wavelength /nm', 'Squared visibility', top_title, 0.)
           end if
        end if
        if (sel_plot == 'vis2mjd' .and. useful_vis > 0) then
           xzero = floor(minval(vis_data(:,7)))
           write (xlabel, *) 'Modified Julian Day -', xzero
           if (zoom) then 
              call plot_vis(7, model_spec, fit_param, &
                   xlabel, 'Squared visibility', top_title, xzero, &
                   uxmin, uxmax)
           else
              call plot_vis(7, model_spec, fit_param, &
                   xlabel, 'Squared visibility', top_title, xzero)
           end if
        end if
        if (sel_plot == 't3amp' .and. useful_amp > 0) then
           if (zoom) then 
              call plot_triple_amp_bas(model_spec, fit_param, &
                   'Longest baseline /M\gl', 'Triple amplitude', &
                   top_title, uxmin, uxmax)
           else
              call plot_triple_amp_bas(model_spec, fit_param, &
                   'Longest baseline /M\gl', 'Triple amplitude', &
                   top_title)
           end if
        end if
        if (sel_plot == 't3ampwl' .and. useful_amp > 0) then
           if (zoom) then 
              call plot_triple_amp(1, model_spec, fit_param, &
                   'Wavelength /nm', &
                   'Triple amplitude', top_title, 0., uxmin, uxmax)
           else
              call plot_triple_amp(1, model_spec, fit_param, &
                   'Wavelength /nm', &
                   'Triple amplitude', top_title, 0.)
           end if
        end if
        if (sel_plot == 't3ampmjd' .and. useful_amp > 0) then
           xzero = floor(minval(triple_data(:,11)))
           write (xlabel, *) 'Modified Julian Day -', xzero
           if (zoom) then 
              call plot_triple_amp(11, model_spec, fit_param, &
                   xlabel, 'Triple amplitude', top_title, xzero, uxmin, uxmax)
           else
              call plot_triple_amp(11, model_spec, fit_param, &
                   xlabel, 'Triple amplitude', top_title, xzero)
           end if
        end if
        if (sel_plot == 't3phi' .and. useful_cp > 0) then
           if (zoom) then 
              call plot_triple_phase_bas(model_spec, fit_param, &
                   'Longest baseline /M\gl', &
                   'Closure phase /'//char(176), top_title, uxmin, uxmax)
           else
              call plot_triple_phase_bas(model_spec, fit_param, &
                   'Longest baseline /M\gl', &
                   'Closure phase /'//char(176), top_title)
           end if
        end if
        if (sel_plot == 't3phiwl' .and. useful_cp > 0) then
           if (zoom) then 
              call plot_triple_phase(1, model_spec, fit_param, &
                   'Wavelength /nm', &
                   'Closure phase /'//char(176), top_title, 0., uxmin, uxmax)
           else
              call plot_triple_phase(1, model_spec, fit_param, &
                   'Wavelength /nm', &
                   'Closure phase /'//char(176), top_title, 0.)
           end if
        end if
        if (sel_plot == 't3phimjd' .and. useful_cp > 0) then
           xzero = floor(minval(triple_data(:,11)))
           write (xlabel, *) 'Modified Julian Day -', xzero
           if (zoom) then 
              call plot_triple_phase(11, model_spec, fit_param, &
                   xlabel, 'Closure phase /'//char(176), top_title, xzero, &
                   uxmin, uxmax)
           else
              call plot_triple_phase(11, model_spec, fit_param, &
                   xlabel, 'Closure phase /'//char(176), top_title, xzero)
           end if
        end if
        if (sel_plot == 'post') then
           if (.not. zoom) then
              x0 = sol(index, 1)
              if (sol(index, 2) == 0D0) then
                 !was problem calculating errors from hessian
                 sigma = model_prior(x_pos(index, 1), x_pos(index, 2))
              else
                 sigma = sol(index, 2)
              end if
              uxmin = x0 - 3*sigma
              uxmax = x0 + 3*sigma
           end if
           call plot_post(model_spec, fit_param, index, &
                desc(index), '-ln(post)', top_title, uxmin, uxmax)
        end if
     end if
  end if

  print *,' '
  print *,spacer_line

  !-------------------------------------------------------------------------
  !Deallocate model/fitting storage
  call free_model() !model_*
  call free_fit() !fit_param, x_pos, x_info
  if (allocated(desc)) deallocate(desc)
  if (allocated(sol)) deallocate(sol)
  if (allocated(hes)) deallocate(hes)
  if (allocated(cov)) deallocate(cov)
  if (allocated(cor)) deallocate(cor)


  !----------------------------------------------------------------------------
  if (sel_plot /= '') call pgend !close graphics

  !Deallocate data storage
  if (allocated(vis_data)) deallocate(vis_data)
  if (allocated(triple_data)) deallocate(triple_data)
  if (allocated(wavebands)) deallocate(wavebands)
  if (allocated(sel_wavebands)) deallocate(sel_wavebands)

contains

!==============================================================================

  subroutine filt_by_wb(data, wb, sig)

    !Filter 2d double precision data array by 1st two columns
    !Rows are kept if 1st two columns match wb to precision sig
    !Reallocates data array
    double precision, dimension(:, :), allocatable :: data
    double precision, dimension(2) :: wb
    double precision :: sig
    integer dim2

    !local variables
    integer nfilt
    logical, dimension(:), allocatable :: mask
    double precision, dimension(:, :), allocatable :: filt_data

    dim2 = size(data,2)
    allocate(mask(size(data,1)))
    mask = (data(:, 1) >= wb(1)-sig .and. data(:, 1) <= wb(1)+sig &
         .and. data(:, 2) >= wb(2)-sig .and. data(:, 2) <= wb(2)+sig)
    nfilt = count(mask)
    allocate(filt_data(nfilt, dim2))
    filt_data = reshape(pack(data, spread(mask, 2, dim2)), (/nfilt, dim2/))
    deallocate(data)
    allocate(data(nfilt, dim2))
    data = filt_data
    deallocate(filt_data)
    deallocate(mask)

  end subroutine filt_by_wb

!==============================================================================

  subroutine filt_by_wl(data, wlmin, wlmax)

    !Filter 2d double precision data array by 1st column
    !Rows are kept if 1st column between wlmin and wlmax
    !Reallocates data array
    double precision, dimension(:, :), allocatable :: data
    double precision :: wlmin, wlmax
    integer dim2

    !local variables
    integer nfilt
    logical, dimension(:), allocatable :: mask
    double precision, dimension(:, :), allocatable :: filt_data

    dim2 = size(data,2)
    allocate(mask(size(data,1)))
    mask = (data(:, 1) >= wlmin .and. data(:, 1) <= wlmax)
    nfilt = count(mask)
    allocate(filt_data(nfilt, dim2))
    filt_data = reshape(pack(data, spread(mask, 2, dim2)), (/nfilt, dim2/))
    deallocate(data)
    allocate(data(nfilt, dim2))
    data = filt_data
    deallocate(filt_data)
    deallocate(mask)

  end subroutine filt_by_wl

!==============================================================================

  subroutine print_usage()

    print *, 'Usage:'
    print *, ' '
    print *, 'clfit [options] datafile modelfile'
    print *, ' '
    print *, 'Options (may appear in any order):'
    print *, ' '
    print *, '-n|--nofit              just report initial chi^2 and make plot (if specified)'
    print *, '                        for initial model'
    print *, '-w|--waveband CWL BW    select this waveband;'
    print *, '                        must specify wb this way for .(n)vis data'
    print *, '-r|--waverange WL1 WL2  select all wavebands in this range'
    print *, '-t|--target_id ID       select this target (default is 1st in OI_TARGET table)'
    print *, '-c|--calerr FRAC        add calibration error (frac. err. in system visib.)'
    print *, '-p|--plot      uv|post N|vis2|t3amp|t3phi|vis2wl|t3ampwl|t3phiwl|vis2mjd|t3ampmjd|t3phimjd'
    print *, '                        make specified plot'
    print *, '-z|--zoomplot  uv|post N|vis2|t3amp|t3phi|vis2wl|t3ampwl|t3phiwl|vis2mjd|t3ampmjd|t3phimjd XMIN XMAX'
    print *, '                        plot with specified x-axis range'
    print *, '-d|--device DEV         PGPLOT device to use'
    print *, ' '

  end subroutine print_usage

end program Main
