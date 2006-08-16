!$Id: main.f90,v 1.24 2006/08/16 09:20:22 jsy1001 Exp $

program Main

  use Inout
  use Plot
  use Visibility
  use Fit
  use Model

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
  double precision, dimension(2) :: wb, wave
  character(len=width) :: spacer_line
  character(len=128) :: info, file_name, ext, source
  character(len=128) :: xlabel, ylabel, top_title, xrange, wavestr
  integer :: i, length, flag
  integer :: degfreedom, useful_vis, useful_amp, useful_cp
  integer, dimension(2) :: indx
  double precision :: nlposterior, nlevidence, nlnpost, chisqrd, normchisqrd
  double precision :: calib_error, uxmin, uxmax, uymin, uymax, temp
  double precision :: x0, y0, xsig, ysig
  logical :: force_symm, mod_line, marg

  integer :: pgopen, istat

  !----------------------------------------------------------------------------
  !Formatting stuff
  spacer_line = repeat('-', width)

  !----------------------------------------------------------------------------
  !Introduction

  print *,' '
  print *,spacer_line
  print *,' '
  print *,'  mfit model fitting program'
  print *,'  package release ',release
  print *,' '
  print *,spacer_line

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

  read_data: do
     print *,' '
     print *,'enter data filename in current directory'
     read (*, '(a)') file_name
     ext = trim(file_name(scan(file_name,'.',.true.)+1:len(file_name)))

     if (ext(len_trim(ext)-3:len_trim(ext)) == 'fits') &
          print *, 'OIFITS-format data often includes calibration error already'
     do
        print *, ' '
        print *, 'enter calibration error to add (fractional error in system visibility)'
        read *, calib_error
        if (calib_error >= 0D0) exit
        print *, 'must specify a positive calibration error'
     end do

     if (ext == 'vis' .or. ext == 'nvis') then
        print *, 'enter observing wavelength and bandwidth (nm) for (n)vis-format data'
        read *, wb
        allocate(wavebands(1, 2))
        wavebands(1, :) = wb
     end if

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

        !read_vis allocates vis_data
        call read_vis(info, file_name, source, max_lines, vis_data, wb, &
             calib_error)
        allocate(triple_data(0, 0))
        force_symm = .true.

     else if (ext == 'nvis') then
        !Model is forced to be centrosymmetric (real visibilities)
        !.nvis only has visibility measurements & errors, plus
        !projected baselines sqrt(u**2 + v**2), stored in u column

        !read_nvis allocates vis_data
        call read_nvis(info, file_name, source, max_lines, vis_data, wb, &
             calib_error)
        allocate(triple_data(0, 0))
        force_symm = .true.

     else if (ext(len_trim(ext)-3:len_trim(ext)) == 'fits') then
        !read_oi_fits allocates vis_data, triple_data, and wavebands
        call read_oi_fits(info, file_name, -1, source, &
             vis_data, triple_data, wavebands, calib_error)
        force_symm = .false.

     else
        info = 'file type "'//trim(ext)//'" not handled'

     end if
     if (info == '') exit
     print *, trim(info)

  end do read_data

  !Filter here to reduce vis and triple data to chosen wavebands
  if (size(wavebands,1) > 1) then
     print *,' '
     print *,'multiple waveband data found with the following wavelengths'
     do i = 1, size(wavebands,1)
        print 57, 'waveband', i, 'wavelength (nm)', real(wavebands(i, 1)), &
             'bandwidth (nm)', real(wavebands(i, 2))
57      format(a, 1x, i3, 1x, a, 1x, f8.2, 1x, a, 1x, f7.2)
     end do
     choose_wave: do
        print *, ' '
        print *, 'enter waveband number or wavelength range to be used in fitting'
        wave = (/-1.0D0, -1.0D0/)
        read (*, '(a)') wavestr
        read (wavestr, *, end=58) wave(1), wave(2)
        !two values - wavelength range
        if (wave(2) .lt. wave(1)) then
           temp = wave(2)
           wave(2) = wave(1)
           wave(1) = temp
        end if
        call filt_by_wl(vis_data, wave(1), wave(2))
        call filt_by_wl(triple_data, wave(1), wave(2))
        allocate(sel_wavebands(size(wavebands,1), 2))
        sel_wavebands = wavebands
        call filt_by_wl(sel_wavebands, wave(1), wave(2))
        exit choose_wave
58      if (wave(1) .eq. -1.0D0) then
           !nothing entered
           cycle choose_wave
        end if
        !one value - waveband number
        if (wave(1) < 1 .or. wave(1) > size(wavebands, 1)) then
           print *, 'Waveband number should be in range 1 to ', size(wavebands, 1)
           cycle choose_wave
        end if
        call filt_by_wb(vis_data, wavebands(int(wave(1)), :), sig)
        call filt_by_wb(triple_data, wavebands(int(wave(1)), :), sig)
        allocate(sel_wavebands(1, 2))
        sel_wavebands(1, :) = wavebands(int(wave(1)), :)
        exit choose_wave
     end do choose_wave

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
  print *, '    MJD   wave   band   baseline coords       sqrd      abs'
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
  if (ext /= 'vis' .and. ext /= 'nvis' .and. useful_vis > 0) then
     call plot_uv('u /M\gl', 'v /M\gl', trim(source), '?')
  else
     istat = pgopen('?')
  end if

  ! allow repeated model fits to same data
  model_loop: do

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

     do
        print *, ' '
        print *, 'enter model filename in current directory (or [return] to exit)'
        read (*, '(a)') file_name
        if (file_name == '') exit model_loop
        print *, ' '
        print *, 'reading model...'
        info = ''
        call read_model(info, file_name, sel_wavebands)
        if (info == '') exit
        print *, trim(info)
     end do

     print *, '...done'
     print *, ' '
     print *, 'model details:'
     call print_model()
     print *, ' '
     print *, spacer_line

     !-------------------------------------------------------------------------
     ! plot initial model
     mod_line = (symm .and. size(sel_wavebands, 1) == 1)
     top_title = trim(source)//' - initial model: '//trim(model_name)
     if (useful_vis > 0) &
          call plot_vis_bas(model_spec, model_param, mod_line, &
          'Baseline /M\gl', 'Squared visibility', top_title)
     if (useful_amp > 0) &
          call plot_triple_amp_bas(model_spec, model_param, &
          'Longest baseline /M\gl', 'Triple amplitude', top_title)
     if (useful_cp > 0) &
          call plot_triple_phase_bas(model_spec, model_param, &
          'Longest baseline /M\gl', &
          'Closure phase /'//char(176), top_title)

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

        print *,' '
        print *,spacer_line

        !----------------------------------------------------------------------
        !plot

        mod_line = (symm .and. size(sel_wavebands, 1) == 1)
        top_title = trim(source)//' - final model: '//trim(model_name)
        if (useful_vis > 0) then
           call plot_vis_bas(model_spec, fit_param, mod_line, 'Baseline /M\gl', &
                'Squared Visibility', top_title)
           print *, 'enter x-axis range for replot ([return] to skip)'
           read (*, '(a)') xrange
           if (len_trim(xrange) .gt. 0) then
              read (xrange, *) uxmin, uxmax
              call plot_vis_bas(model_spec, fit_param, mod_line, 'Baseline /M\gl', &
                   'Squared Visibility', top_title, uxmin, uxmax)
           end if
        end if
        if (useful_amp > 0) then 
           call plot_triple_amp_bas(model_spec, fit_param, &
                'Longest baseline /M\gl', 'Triple amplitude', top_title)
           print *, 'enter x-axis range for replot ([return] to skip)'
           read (*, '(a)') xrange
           if (len_trim(xrange) .gt. 0) then
              read (xrange, *) uxmin, uxmax
              call plot_triple_amp_bas(model_spec, fit_param, &
                   'Longest baseline /M\gl', 'Triple amplitude', top_title, &
                   uxmin, uxmax)
             end if
        end if
        if (useful_cp > 0) then
           call plot_triple_phase_bas(model_spec, fit_param, &
                'Longest baseline /M\gl', 'Closure phase /'//char(176), &
                top_title)
           print *, 'enter x-axis range for replot ([return] to skip)'
           read (*, '(a)') xrange
           if (len_trim(xrange) .gt. 0) then
              read (xrange, *) uxmin, uxmax
              call plot_triple_phase_bas(model_spec, fit_param, &
                   'Longest baseline /M\gl', 'Closure phase /'//char(176), &
                   top_title, uxmin, uxmax)
             end if
        end if
        if (yesno('Plot 1d cut through -ln(postprob)', 'no')) then
           do
              print *, 'enter variable number for x axis'
              read *, indx(1)
              if (indx(1) <= length) exit
           end do
           x0 = sol(indx(1), 1)
           if (sol(indx(1), 2) == 0D0) then
              !was problem calculating errors from hessian
              xsig = model_prior(x_pos(indx(1), 1), x_pos(indx(1), 2))
           else
              xsig = sol(indx(1), 2)
           end if
           uxmin = x0 - 3*xsig
           uxmax = x0 + 3*xsig
           if (yesno('Marginalise over other variables', 'no')) then
              marg = .true.
              nlnpost = nlevidence
              ylabel = '-ln(MARG. postprob)'
           else
              marg = .false.
              nlnpost = nlposterior - 0.5D0*log(2D0*pi) + &
                   0.5D0*log(hes(indx(1),indx(1)))
              ylabel = '-ln(postprob)'
           end if
           call plot_post(marg, nlnpost, indx(1), desc(indx(1)), &
                ylabel, top_title, uxmin, uxmax)
           print *, 'enter x-axis range for replot ([return] to skip)'
           read (*, '(a)') xrange
           if (len_trim(xrange) .gt. 0) then
              read (xrange, *) uxmin, uxmax
              call plot_post(marg, nlnpost, indx(1), desc(indx(1)), &
                   ylabel, top_title, uxmin, uxmax)
           end if
        end if
        if (yesno('Plot 2d cut through -ln(postprob)', 'no')) then
           do
              print *, 'enter variable numbers for x and y axes'
              read *, indx
              if (indx(1) <= length .and. indx(2) <= length) exit
           end do
           x0 = sol(indx(1), 1)
           if (sol(indx(1), 2) == 0D0) then
              !was problem calculating errors from hessian
              xsig = model_prior(x_pos(indx(1), 1), x_pos(indx(1), 2))
           else
              xsig = sol(indx(1), 2)
           end if
           uxmin = x0 - 3*xsig
           uxmax = x0 + 3*xsig
           y0 = sol(indx(2), 1)
           if (sol(indx(2), 2) == 0D0) then
              !was problem calculating errors from hessian
              ysig = model_prior(x_pos(indx(2), 1), x_pos(indx(2), 2))
           else
              ysig = sol(indx(2), 2)
           end if
           uymin = y0 - 3*ysig
           uymax = y0 + 3*ysig
           if (yesno('Marginalise over other variables', 'no')) then
              marg = .true.
              top_title = trim(top_title)//': -ln(MARG. postprob)'
           else
              marg = .false.
              top_title = trim(top_title)//': -ln(postprob)'
           end if
           call plot_post2d(marg, nlposterior, indx, desc(indx(1)), &
                desc(indx(2)), top_title, uxmin, uxmax, uymin, uymax)
           print *, 'enter x-axis range for replot ([return] to skip)'
           read (*, '(a)') xrange
           if (len_trim(xrange) .gt. 0) then
              print *, 'enter y-axis range for replot'
              read *, uymin, uymax
              call plot_post2d(marg, nlposterior, indx, desc(indx(1)), &
                   desc(indx(2)), top_title, uxmin, uxmax, uymin, uymax)
           end if
        end if
     end if

     !-------------------------------------------------------------------------
     !Deallocate model/fitting storage
     call free_model() !model_*
     call free_fit() !fit_param, x_pos, x_bound
     if (allocated(desc)) deallocate(desc)
     if (allocated(sol)) deallocate(sol)
     if (allocated(hes)) deallocate(hes)
     if (allocated(cov)) deallocate(cov)
     if (allocated(cor)) deallocate(cor)

  end do model_loop

  !----------------------------------------------------------------------------
  call pgend !close graphics

  !Deallocate data storage
  if (allocated(vis_data)) deallocate(vis_data)
  if (allocated(triple_data)) deallocate(triple_data)
  if (allocated(wavebands)) deallocate(wavebands)

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

end program Main
