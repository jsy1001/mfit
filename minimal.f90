program Main

  use f2kcli
  use Inout, only: read_oi_fits
  use Bayes, only: vis_data, triple_data, num_vis, num_triple
  use Model, only: read_model, model_valid, model_nvar, print_model, free_model
  use Fit, only: minimiser_wrap, err_est, is_minimum, free_fit

  implicit none
  
  !============================================================================
  !variables
  
  !parameters
  integer, parameter :: max_lines = 5000      !max. lines in data file
  integer, parameter :: width = 78            !for spacer lines

  !local variables to do with fit
  double precision, allocatable :: sol(:), err(:)
  double precision, allocatable :: hes(:,:), cov(:,:), cor(:,:)

  !other local variables
  double precision, allocatable :: wavebands(:,:)
  character(len=width) :: spacer_line
  character(len=128) :: switch, arg, info, file_name, ext, source
  integer :: narg, iarg, i, n, user_target_id
  integer :: degfreedom, useful_vis, useful_amp, useful_cp
  double precision :: nlposterior, nlevidence, chisqrd, normchisqrd
  double precision :: calib_error
  logical :: fit_ok, found_min, hes_valid

  !----------------------------------------------------------------------------
  !Formatting stuff
  spacer_line = repeat('-', width)

  !----------------------------------------------------------------------------
  !Parse command-line arguments

  narg = command_argument_count()
  iarg = 1
  !defaults
  user_target_id = -1 !default to 1st target in OI_TARGET
  calib_error = 0D0
  do
     if (iarg == narg - 1) exit
     if (iarg > narg - 1) then
        stop 'Not enough arguments'
     end if
     call get_command_argument(iarg, switch)
     select case(switch)
        case('-c', '--calerr')
           call get_command_argument(iarg+1, arg)
           read(arg, *) calib_error
           iarg = iarg + 1
        case('-t', '--target_id')
           call get_command_argument(iarg+1, arg)
           read(arg, *) user_target_id
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
  if (ext(len_trim(ext)-3:len_trim(ext)) == 'fits') then
     !read_oi_fits allocates vis_data, triple_data, and wavebands
     call read_oi_fits(info, file_name, user_target_id, source, calib_error)
  else
     info = 'file type "'//trim(ext)//'" not handled'
  end if
  if (info /= '') then
     print *, trim(info)
     stop
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
  call read_model(info, file_name, wavebands)
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
  !fit model to the data by minimising negative log posterior
  
  print *, ' '
  print *, 'fitting model by minimising negative log posterior...'

  !check model freedoms
  if (.not. model_valid(info, .false.)) then
     print *, trim(info)
     stop
  end if

  !allocate and initialise arrays related to variable model parameters
  call model_nvar(n)
  degfreedom = useful_vis + useful_amp + useful_cp - n
  !allocate(var_pos(n,2), var_desc(n))
  !call model_getvar(n, var_pos, var_desc)

  !allocate arrays for fit results
  allocate(sol(n), err(n), hes(n,n), cov(n,n), cor(n,n))

  !minimise
  call minimiser_wrap(n, sol, chisqrd, nlposterior, fit_ok, info)
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
!      do i = 1, n
!         write(*,62) i, var_desc(i), sol(i), err(i)
!      end do
! 62   format(' (', i2, ') ', A55, 1x, f13.6, 1x, f12.6) 
      do i = 1, n
         write(*,62) i, sol(i), err(i)
      end do
 62   format(' (', i2, ') ', f13.6, 1x, f12.6) 

     if (hes_valid) then
        ! Display alternative error bars for fitted parameters: estimates
        ! assuming data errors scaled so chi sqrd/deg freedom = unity
        print *, ' '
        print '(1x, a, 49x, a)', 'SCALED DATA ERRORS ->', &
             'fitted      hessian'
        print '(1x, a, 49x, a)', '  num  parameter name', &
             ' value        error'
!         do i = 1, n
!            write(*,63) i, var_desc(i), sol(i), sqrt(normchisqrd)*err(i)
!         end do
! 63      format('  *(', i2, ') ', A55, 1x, f13.6, 1x, f12.6, '*') 
         do i = 1, n
            write(*,63) i, sol(i), sqrt(normchisqrd)*err(i)
         end do
 63      format('  *(', i2, ') ', f13.6, 1x, f12.6, '*') 

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

64      format(1x, 10(e11.4,1x))
     end if
  end if

  print *,' '
  print *,spacer_line

  !-------------------------------------------------------------------------
  !Deallocate model/fitting storage
  call free_model() !model_*
  call free_fit()
  !if (allocated(var_pos)) deallocate(var_pos)
  !if (allocated(var_desc)) deallocate(var_desc)
  if (allocated(sol)) deallocate(sol)
  if (allocated(err)) deallocate(err)
  if (allocated(hes)) deallocate(hes)
  if (allocated(cov)) deallocate(cov)
  if (allocated(cor)) deallocate(cor)

  !----------------------------------------------------------------------------

  !Deallocate data storage
  if (allocated(vis_data)) deallocate(vis_data)
  if (allocated(triple_data)) deallocate(triple_data)
  if (allocated(wavebands)) deallocate(wavebands)

end program Main
