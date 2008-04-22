!$Id: clfit.f90,v 1.24 2008/04/22 09:51:19 jsy1001 Exp $

program Main

  use f2kcli
  use Inout
  use Model
  use Fit
  use Wrap
  use Plot
  use PostPlot

  implicit none
  
  !============================================================================
  !variables
  
  !parameters
  !pi, deg/rad conversions all picked up from Maths module
  integer, parameter :: max_lines = 5000      !max. lines in data file
  integer, parameter :: width = 78            !for spacer lines
  double precision, parameter :: sig = 0.1D0  !waveband must match to sig nm

  !local variables to do with fit
  type(allparam) :: allpar
  integer, allocatable :: var_pos(:,:)
  character(len=model_desc_len), allocatable :: var_desc(:)
  double precision, allocatable :: sol(:), err(:)
  double precision, allocatable :: hes(:,:), cov(:,:), cor(:,:)

  !other local variables
  double precision, allocatable :: wavebands(:,:)
  double precision, allocatable :: errguess(:)
  double precision :: wb(2), wl(2), alt_err(2)
  character(len=width) :: spacer_line
  character(len=128) :: switch, arg, device, info, file_name, ext, source
  character(len=128) :: xlabel, ylabel, top_title, cvs_rev, revision
  character(len=8) :: sel_plot
  integer :: narg, iarg, i, n, user_target_id, margerr_var
  integer :: indx(2)
  integer :: degfreedom, useful_vis, useful_amp, useful_cp, xindex
  double precision :: nlposterior, nlevidence, nlnpost, chisqrd, normchisqrd
  double precision :: calib_error, uxmin, uxmax, uymin, uymax
  double precision :: x0, y0, xsig, ysig
  logical :: force_symm, nofit, zoom, mod_line, marg
  logical :: fit_ok, found_min, hes_valid
  logical :: pl_vis, pl_amp, pl_cp

  integer :: pgopen, istat

  !----------------------------------------------------------------------------
  !Formatting stuff
  spacer_line = repeat('-', width)

  !----------------------------------------------------------------------------
  !Introduction

  cvs_rev = '$Revision: 1.24 $'
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
  margerr_var = -1
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
           if (sel_plot(:4) == 'post' .or. sel_plot(:5) == 'mpost') then
              call get_command_argument(iarg+1, arg)
              read(arg, *) indx(1)
              iarg = iarg + 1
           end if
           if (sel_plot == 'post2d' .or. sel_plot == 'mpost2d') then
              call get_command_argument(iarg+1, arg)
              read(arg, *) indx(2)
              iarg = iarg + 1
           end if
        case('-z', '--zoomplot')
           zoom = .true.
           call get_command_argument(iarg+1, sel_plot)
           if (sel_plot(:4) == 'post' .or. sel_plot(:5) == 'mpost') then
              call get_command_argument(iarg+2, arg)
              read(arg, *) indx(1)
              iarg = iarg + 1
           end if
           if (sel_plot == 'post2d' .or. sel_plot == 'mpost2d') then
              call get_command_argument(iarg+2, arg)
              read(arg, *) indx(2)
              iarg = iarg + 1
           end if
           call get_command_argument(iarg+2, arg)
           read(arg, *) uxmin
           call get_command_argument(iarg+3, arg)
           read(arg, *) uxmax
           if (sel_plot == 'post2d' .or. sel_plot == 'mpost2d') then
              call get_command_argument(iarg+4, arg)
              read(arg, *) uymin
              call get_command_argument(iarg+5, arg)
              read(arg, *) uymax
              iarg = iarg + 2
           end if
           iarg = iarg + 3
        case('-t', '--target_id')
           call get_command_argument(iarg+1, arg)
           read(arg, *) user_target_id
           iarg = iarg + 1
        case('-m', '--margerr')
           call get_command_argument(iarg+1, arg)
           read(arg, *) margerr_var
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

  else if (ext(len_trim(ext)-4:len_trim(ext)) == 'calib') then
     !.*Calib has no waveband data so must be supplied by wavelength range
     if (wb(2) .eq. -1.0D0) &
          stop 'need waveband for wbCalib format'
     print '(1x, a, f8.2, a)', 'Assuming waveband for data is', wb(2), ' nm'
     call read_calib(info, file_name, source, max_lines, vis_data, &
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

  if (nofit) then
     !-------------------------------------------------------------------------
     ! report reduced chi-squared, plot initial model
     call gof(model_spec, model_param, chisqrd)
     call model_nvar(n)
     degfreedom = useful_vis + useful_amp + useful_cp - n
     print *, ' '
     print *, '           chi squared =',real(chisqrd) 
     print *, '    degrees of freedom =',degfreedom
     if (degfreedom > 0) then
        normchisqrd = chisqrd/degfreedom
        print *, 'chi sqrd / deg freedom =',real(normchisqrd)
     end if
     mod_line = (symm_model .and. size(sel_wavebands, 1) == 1)
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
        !cannot do '(m)post(2d)' plot if --nofit
     else
        call parse_plot_arg(sel_plot, pl_vis, pl_amp, pl_cp, &
             xindex, xlabel, ylabel)
        if (pl_vis .and. useful_vis > 0) then
           if (zoom) then 
              call plot_vis(xindex, model_spec, model_param, &
                   xlabel, ylabel, top_title, uxmin, uxmax)
           else
              call plot_vis(xindex, model_spec, model_param, &
                   xlabel, ylabel, top_title)
           end if
        else if (pl_amp .and. useful_amp > 0) then
           if (zoom) then 
              call plot_triple_amp(xindex, model_spec, model_param, &
                   xlabel, ylabel, top_title, uxmin, uxmax)
           else
              call plot_triple_amp(xindex, model_spec, model_param, &
                   xlabel, ylabel, top_title)
           end if
        else if (pl_cp .and. useful_cp > 0) then
           if (zoom) then 
              call plot_triple_phase(xindex, model_spec, model_param, &
                   xlabel, ylabel, top_title, uxmin, uxmax)
           else
              call plot_triple_phase(xindex, model_spec, model_param, &
                   xlabel, ylabel, top_title)
           end if
        end if
     end if
  else
     !-------------------------------------------------------------------------
     !fit model to the data by minimising negative log posterior

     print *, ' '
     print *, 'fitting model by minimising negative log posterior...'

     !check model freedoms
     if (.not. model_valid(info, force_symm)) then
        print *, trim(info)
        stop
     end if

     !allocate and initialise arrays related to variable model parameters
     call model_nvar(n)
     degfreedom = useful_vis + useful_amp + useful_cp - n
     allocate(var_pos(n,2), var_desc(n))
     call model_getvar(n, var_pos, var_desc)

     !allocate arrays for fit results
     allocate(sol(n), err(n), hes(n,n), cov(n,n), cor(n,n))

     call allparam_init(allpar, model_param, model_limits, n, var_pos)

     !calculate initial goodness-of-fit
     call gof(model_spec, allpar%param, chisqrd)
     print *,' '
     print *,'initial chi squared =',real(chisqrd) 

     !minimise
     call minimiser(allpar, sol, chisqrd, nlposterior, fit_ok, info)
     if (info /= '') then
        print *,'*****'
        print *,trim(info)
        print *, '*****'
     end if

     if (fit_ok) then
        if (.not. is_minimum(n, sol)) &
             print *, 'NOT A LOCAL MINIMUM'
        call err_est(n, sol, &
             found_min, hes_valid, hes, cov, cor, err, nlevidence, info)
        if (info /= '') then
           print *,'*****'
           print *,trim(info)
           print *, '*****'
        end if
     end if

     if (found_min) then
        print *, ' '
        print *, 'negative log posterior =',real(nlposterior)
        if (hes_valid) print *, 'negative log evidence  =',real(nlevidence)
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
        do i = 1, n
           write(*,62) i, var_desc(i), sol(i), err(i)
        end do
62      format(' (', i2, ') ', A55, 1x, f13.6, 1x, f12.6) 

        if (hes_valid) then
           ! Display alternative error bars for fitted parameters: estimates
           ! assuming data errors scaled so chi sqrd/deg freedom = unity
           print *, ' '
           print '(1x, a, 49x, a)', 'SCALED DATA ERRORS ->', &
                'fitted      hessian'
           print '(1x, a, 49x, a)', '  num  parameter name', &
                ' value        error'
           do i = 1, n
              write(*,63) i, var_desc(i), sol(i), sqrt(normchisqrd)*err(i)
           end do
63         format('  *(', i2, ') ', A55, 1x, f13.6, 1x, f12.6, '*') 

           print *, ' '
           print *, 'hessian matrix'
           do i = 1, n
              write(*,64) hes(i,:)
           end do

           print *, ' '
           print *, 'covariance matrix'
           do i = 1, n
              write(*,64) cov(i,:)
           end do

           print *, ' '
           print *, 'correlation matrix'
           do i = 1, n
              write(*,64) cor(i,:)
           end do

64         format(1x, 10(e11.4,1x))
        end if

        !----------------------------------------------------------------------
        !plot
        call allparam_setvar(allpar, sol)
        mod_line = (symm_model .and. size(sel_wavebands, 1) == 1)
        top_title = trim(source)//' - final model: '//trim(model_name)
        if (sel_plot == 'post' .or. sel_plot == 'mpost') then
           if (indx(1) > n) stop 'Invalid parameter number'
           if (.not. zoom) then
              x0 = sol(indx(1))
              if (.not. hes_valid) then
                 !was problem calculating errors from hessian
                 xsig = model_prior(var_pos(indx(1), 1), var_pos(indx(1), 2))
              else
                 xsig = err(indx(1))
              end if
              uxmin = x0 - 3*xsig
              uxmax = x0 + 3*xsig
           end if
           if (sel_plot == 'mpost') then
              marg = .true.
              if (hes_valid) then
                 nlnpost = nlevidence
              else
                 nlnpost = nlposterior
              end if
              ylabel = '-ln(MARG. postprob)'
           else
              marg = .false.
              if (hes_valid) then
                 nlnpost = nlposterior - 0.5D0*log(2D0*pi) + &
                      0.5D0*log(hes(indx(1),indx(1)))
              else
                 nlnpost = nlposterior !just want a plot
              end if
              ylabel = '-ln(postprob)'
           end if
           call plot_post1d(marg, nlnpost, allpar, indx(1), &
                var_desc(indx(1)), ylabel, top_title, uxmin, uxmax)
        else if (sel_plot == 'post2d' .or. sel_plot == 'mpost2d') then
           if (indx(1) > n .or. indx(2) > n) &
                stop 'Invalid parameter number'
           if (.not. zoom) then
              x0 = sol(indx(1))
              if (.not. hes_valid) then
                 !was problem calculating errors from hessian
                 xsig = model_prior(var_pos(indx(1), 1), var_pos(indx(1), 2))
              else
                 xsig = err(indx(1))
              end if
              uxmin = x0 - 3*xsig
              uxmax = x0 + 3*xsig
              y0 = sol(indx(2))
              if (.not. hes_valid) then
                 !was problem calculating errors from hessian
                 ysig = model_prior(var_pos(indx(2), 1), var_pos(indx(2), 2))
              else
                 ysig = err(indx(2))
              end if
              uymin = y0 - 3*ysig
              uymax = y0 + 3*ysig
           end if
           if (sel_plot == 'mpost2d') then
              marg = .true.
              top_title = '-ln(MARG. postprob) - '//trim(top_title)
           else
              marg = .false.
              top_title = '-ln(postprob) - '//trim(top_title)
           end if
           call plot_post2d(marg, nlposterior, allpar, indx, &
                var_desc(indx(1)), var_desc(indx(2)), top_title, &
                uxmin, uxmax, uymin, uymax)
        else if (sel_plot == 'vis2' .and. useful_vis > 0) then
           if (zoom) then 
              call plot_vis_bas(model_spec, allpar%param, mod_line, &
                   'Baseline /M\gl', 'Squared visibility', top_title, &
                   uxmin, uxmax)
           else
              call plot_vis_bas(model_spec, allpar%param, mod_line, &
                   'Baseline /M\gl', 'Squared visibility', top_title)
           end if
        else
           call parse_plot_arg(sel_plot, pl_vis, pl_amp, pl_cp, &
                xindex, xlabel, ylabel)
           if (pl_vis .and. useful_vis > 0) then
              if (zoom) then 
                 call plot_vis(xindex, model_spec, allpar%param, &
                      xlabel, ylabel, top_title, uxmin, uxmax)
              else
                 call plot_vis(xindex, model_spec, allpar%param, &
                      xlabel, ylabel, top_title)
              end if
           else if (pl_amp .and. useful_amp > 0) then
              if (zoom) then 
                 call plot_triple_amp(xindex, model_spec, allpar%param, &
                      xlabel, ylabel, top_title, uxmin, uxmax)
              else
                 call plot_triple_amp(xindex, model_spec, allpar%param, &
                      xlabel, ylabel, top_title)
              end if
           else if (pl_cp .and. useful_cp > 0) then
              if (zoom) then 
                 call plot_triple_phase(xindex, model_spec, allpar%param, &
                      xlabel, ylabel, top_title, uxmin, uxmax)
              else
                 call plot_triple_phase(xindex, model_spec, allpar%param, &
                      xlabel, ylabel, top_title)
              end if
           end if
        end if

        !----------------------------------------------------------------------
        if (margerr_var >= 1 .and. margerr_var <= n) then
           !check specified error bar
           print *, 'checking error bar',margerr_var,' by marginalisation...'
           allocate(errguess(n))
           errguess = err
           do i = 1, n
              if (err(i) < sqrt(1D0/hes(i,i))) &
                   errguess(i) = sqrt(1D0/hes(i,i))
           end do
           call marg_err(allpar, errguess, margerr_var, alt_err)
           print '(1x, a, 49x, a)', '                   ', 'fitted'
           print '(1x, a, 49x, a)', 'num  parameter name', ' value'
           write(*,77) margerr_var, var_desc(margerr_var), sol(margerr_var), &
                ' +', alt_err(1), ' -', alt_err(2)
77         format(' (', i2, ') ', a55, 1x, f13.6, a, f12.6, a, f12.6) 
           deallocate(errguess)
        end if
     end if
  end if

  print *,' '
  print *,spacer_line

  !-------------------------------------------------------------------------
  !Deallocate model/fitting storage
  call allparam_free(allpar)
  call free_model() !model_*
  call free_fit()
  if (allocated(var_pos)) deallocate(var_pos)
  if (allocated(var_desc)) deallocate(var_desc)
  if (allocated(sol)) deallocate(sol)
  if (allocated(err)) deallocate(err)
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

  subroutine parse_plot_arg(plot_id, plot_v2, plot_t3amp, plot_t3phi, &
       xindex, xlabel, ylabel)

    character(len=*), intent(in) :: plot_id
    logical, intent(out) :: plot_v2, plot_t3amp, plot_t3phi
    integer, intent(out) :: xindex
    character(len=*), intent(out) :: xlabel, ylabel

    plot_v2 = .false.
    plot_t3amp = .false.
    plot_t3phi = .false.

    select case (plot_id)

    case ('vis2')
       plot_v2 = .true.
       xindex = -1
       xlabel = 'Baseline /M\gl'
    case ('t3amp')
       plot_t3amp = .true.
       xindex = -1
       xlabel = 'Longest baseline /M\gl'
    case ('t3phi')
       plot_t3phi = .true.
       xindex = -1
       xlabel = 'Longest baseline /M\gl'
    case ('vis2-wl')
       plot_v2 = .true.
       xindex = 1
       xlabel = 'Wavelength /nm'
    case ('t3amp-wl')
       plot_t3amp = .true.
       xindex = 1
       xlabel = 'Wavelength /nm'
    case ('t3phi-wl')
       plot_t3phi = .true.
       xindex = 1
       xlabel = 'Wavelength /nm'
    case ('vis2-mjd')
       plot_v2 = .true.
       xindex = 7
       xlabel = 'Modified Julian Day'
    case ('t3amp-mjd')
       plot_t3amp = .true.
       xindex = 11
       xlabel = 'Modified Julian Day'
    case ('t3phi-mjd')
       plot_t3phi = .true.
       xindex = 11
       xlabel = 'Modified Julian Day'
    case ('vis2-st')
       plot_v2 = .true.
       xindex = -2
       xlabel = 'GMST /h'
    case ('t3amp-st')
       plot_t3amp = .true.
       xindex = -2
       xlabel = 'GMST /h'
    case ('t3phi-st')
       plot_t3phi = .true.
       xindex = -2
       xlabel = 'GMST /h'
    case ('vis2-pa')
       plot_v2 = .true.
       xindex = -3
       xlabel = 'Baseline PA /deg'
    case default
       return
    end select

    if (plot_v2) ylabel = 'Squared visibility'
    if (plot_t3amp) ylabel = 'Triple amplitude'
    if (plot_t3phi) ylabel = 'Closure phase /'//char(176)

  end subroutine parse_plot_arg
 
!==============================================================================

  subroutine filt_by_wb(data, wb, sig)

    !Filter 2d double precision data array by 1st two columns
    !Rows are kept if 1st two columns match wb to precision sig
    !Reallocates data array
    double precision, allocatable :: data(:,:)
    double precision :: wb(2)
    double precision :: sig
    integer :: dim2

    !local variables
    integer :: nfilt
    logical, allocatable :: mask(:)
    double precision, allocatable :: filt_data(:,:)

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
    double precision, allocatable :: data(:,:)
    double precision :: wlmin, wlmax
    integer :: dim2

    !local variables
    integer :: nfilt
    logical, allocatable :: mask(:)
    double precision, allocatable :: filt_data(:, :)

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
    print *, '-p|--plot      uv|post N|mpost N|post2d M N|mpost2d M N|vis2|t3amp|t3phi'
    print *, '                |vis2-wl|t3amp-wl|t3phi-wl|vis2-mjd|t3amp-mjd|t3phi-mjd'
    print *, '                |vis2-st|t3amp-st|t3phi-st|vis2-pa'
    print *, '                        make specified plot'
    print *, '-z|--zoomplot  uv|post N|mpost N|post2d M N|mpost2d M N|vis2|t3amp|t3phi'
    print *, '                |vis2-wl|t3amp-wl|t3phi-wl|vis2-mjd|t3amp-mjd|t3phi-mjd'
    print *, '                |vis2-st|t3amp-st|t3phi-st|vis2-pa XMIN XMAX'
    print *, '                        plot with specified x-axis range'
    print *, '-d|--device DEV         PGPLOT device to use'
    print *, '-m|--margerr N          alternate error bar N by brute-force marginalisation'
    print *, ' '

  end subroutine print_usage

end program Main
