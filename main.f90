!$Id: main.f90,v 1.6 2002/11/29 10:40:05 jsy1001 Exp $

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
  integer, parameter :: max_lines = 1000 ! max. lines in data file

  !arrays for manipulating data
  double precision, dimension(:,:), allocatable :: swap
  double precision, dimension(:), allocatable :: wavebands
  logical, dimension(:), allocatable :: keep

  !arrays for fit results
  character(len=35), dimension(:), allocatable :: desc
  double precision, dimension(:,:), allocatable :: sol, hes, cov, cor

  !arrays for plotting
  double precision, dimension(:,:), allocatable :: model_vis_data

  !variables
  integer, parameter :: width=80 !for spacer lines
  character(len=width) :: spacer_line
  character(len=128) :: info, file_name, ext, source, top_title
  integer :: i, j, length, flag, degfreedom
  integer :: degfreedom, useful_vis, useful_amp, useful_cp
  double precision :: nlposterior, chisqrd, normchisqrd
  double precision :: version, lambda, u1, v1, u2, v2, step, calib_error, cp
  double complex :: W1, W2, W3, T
  logical :: force_symm

  external myhandler
  integer :: ieeer, ieee_handler, myhandler

  !----------------------------------------------------------------------------
  !Set up Floating Point Exception handler
  ieeer = ieee_handler('set','common',myhandler)

  !----------------------------------------------------------------------------
  !Formatting stuff

  spacer_line = repeat('-',width)

  !----------------------------------------------------------------------------
  !Introduction

  version = 1.10D0
  print *,' '
  print *,spacer_line
  print *,' '
  print *,'  mfit model fitting program'
  print *,'  version ',version
  print *,' '
  print *,spacer_line

  !----------------------------------------------------------------------------
  !Read vis/triple data
  !
  !End up with vis_data and triple_data arrays. 
  !
  !vis_data holds the visibility data points of lambda, u, v, vis, vis_err, 
  !where vis is the sqaured visibility amplitude and vis_err is its abs error
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
     read *,file_name
     ext = trim(file_name(index(file_name,'.')+1:len(file_name)))

     do
        print *,' '
        print *,'enter calibration error (fractional error in system visibility)'
        read *,calib_error
        if (calib_error>=0D0) exit
        print *,'must specify a positive calibration error'
     end do

     if (ext == 'vis' .or. ext == 'nvis') then
        print *,'enter observing wavelength (nm) for (n)vis-format data'
        read *,lambda
     end if

     print *,' '
     print *,'reading ',trim(ext),' file data...'

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

        !Filter here to reduce vis and triple data to single waveband, selected
        !by user, this will not be necessary if multi-waveband fitting is
        !implemented
        !** hive this off
        if (size(wavebands,1)>1) then
           print *,' '
           print *,'multiple waveband data found with the following wavelengths'
           do i = 1, size(wavebands,1)
              print *,'waveband',i,'wavelength (nm)',real(wavebands(i))
           end do
           do
              print *,' '
              print *,'enter waveband number to be used in fitting'
              read *,j
              if (j>=1.and.j<=size(wavebands,1)) then
                 exit
              else
                 print *,'illegal selection'
              end if
           end do
           lambda = wavebands(j)      
           !do the filtering on vis data
           length = size(vis_data,1)
           if (length>0) then
              allocate(keep(length))
              keep = .false.
              do i = 1, length
                 if (vis_data(i,1)==lambda) keep(i) = .true.
              end do
              allocate(swap(count(keep,1),size(vis_data, 2)))
              j = 0
              do i = 1, length
                 if (keep(i)) then
                    j = j + 1
                    swap(j,:) = vis_data(i,:)
                 end if
              end do
              deallocate(vis_data)
              allocate(vis_data(count(keep,1),size(vis_data, 2)))
              vis_data = swap
              deallocate(swap)
              deallocate(keep)
           end if
           !do the filtering on triple data
           length = size(triple_data,1)
           if (length>1) then
              allocate(keep(length))
              keep = .false.
              do i = 1, length
                 if (triple_data(i,1)==lambda) keep(i) = .true.
              end do
              allocate(swap(count(keep,1),9))
              j = 0
              do i = 1, length
                 if (keep(i)) then
                    j = j + 1
                    swap(j,:) = triple_data(i,:)
                 end if
              end do
              deallocate(triple_data)
              allocate(triple_data(count(keep,1),9))
              triple_data = swap
              deallocate(swap)
              deallocate(keep)
           end if
        end if

     else if (ext == 'vis') then
        !Model is forced to be centrosymmetric (real visibilities)
        !.vis only has visibility measurements (no errors) and
        !projected baselines sqrt(u**2 + v**2), stored in u column

        !read_vis allocates vis_data
        call read_vis(info, file_name, source, max_lines, vis_data, lambda, &
             calib_error)
        allocate(triple_data(0,0))
        force_symm = .true.

     else if (ext == 'nvis') then
        !Model is forced to be centrosymmetric (real visibilities)
        !.nvis only has visibility measurements & errors, plus
        !projected baselines sqrt(u**2 + v**2), stored in u column

        !read_nvis allocates vis_data
        call read_nvis(info, file_name, source, max_lines, vis_data, lambda, &
             calib_error)
        allocate(triple_data(0,0))
        force_symm = .true.

     else
        info = 'file type not handled'

     end if
     if (info == '') exit
     print *,info

  end do read_data

  !Count number of useful data points (i.e. with +ve errors, points with
  !-ve or zero errors are ignored in fitting but are held in data arrays
  !for plotting purposes so mustn't count them)
  useful_vis = 0
  useful_amp = 0
  useful_cp = 0
  do i = 1, size(vis_data,1)
     if (vis_data(i,5)>0D0) useful_vis = useful_vis + 1
  end do
  do i = 1, size(triple_data,1)
     if (triple_data(i,7)>0D0) useful_amp = useful_amp + 1
     if (triple_data(i,9)>0D0) useful_cp = useful_cp + 1
  end do

  print *,'...done'
  print *,' '
  print *,trim(ext),' file visibility data for ',trim(source),':'

  print *,' '
  print *,'  wave   baseline coords     sqrd      abs'
  print *,'length        u        v      vis    error'
  print *,'  (nm)      (m)      (m)'  
  do i = 1, size(vis_data,1)
     write(*,60) vis_data(i,:)
  end do
60 format(f7.1, 1x, f8.4, 1x, f8.4, 1x, f8.5, 1x, f8.5)

  print *,' '
  print *,trim(ext),' file triple product data for ',trim(source),':'
  print *,' '
  print *,'  wave            baseline triangle coords', &
       '                        triple product'
  print *,'length       u1       v1       u2       v2', &
       '  amplitude      error   phase   error'
  print *,'  (nm)      (m)      (m)      (m)      (m)'
  do i = 1, size(triple_data,1)
     write(*,61) triple_data(i,:)
  end do
61 format(f7.1, 1x, f8.4, 1x, f8.4, 1x, f8.4, 1x, f8.4, 1x, &
       e10.3, 1x, e10.3, 1x, f7.2, 1x, f7.2, 1x)

  print *,' '
  print *,useful_vis,'unflagged visibility data points out of', size(vis_data,1)
  print *,useful_amp,'unflagged triple product amplitudes out of', &
       size(triple_data,1)
  print *,useful_cp,'unflagged closure phases out of',size(triple_data,1)
  print *,' '
  print *,'total of',useful_vis+useful_amp+useful_cp, &
       'unflagged data points out of',size(vis_data,1)+(2*size(triple_data,1))
  print *,' '
  print *,spacer_line
  call plot_uv('u /M\gl', 'v /M\gl', trim(source), '?')

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
        print *,' '
        print *,'enter model filename in current directory (or [return] to exit)'
        read (*,'(a)') file_name
        if (file_name == '') exit model_loop
        print *,' '
        print *,'reading model...'
        info = ''
        call read_model(info, file_name)
        if (info == '') exit
        print *, info
     end do

     print *,'...done'
     print *,' '
     print *,'model details:'
     call print_model()
     print *,' '
     print *,spacer_line

     !-------------------------------------------------------------------------
     ! plot initial model
     top_title = trim(source)//' - initial model: '//trim(model_name)
     if (useful_vis > 0) &
          call plot_vis(model_spec, model_param, symm, &
          'Baseline /M\gl', 'Squared visibility', top_title)
     if (useful_cp > 0) &
          call plot_triple(model_spec, model_param, 'Longest baseline /M\gl',&
          'Closure phase /'//char(176), top_title)

     !-------------------------------------------------------------------------
     !fit model to the data by minimising negative log posterior
     !
     !Refer to fit module. the fit solution is returned along with various 
     !diagnostic quantities.

     print *,' '
     print *,'fitting model by minimising negative log posterior...'

     info = ''
     ! minimiser allocates fit_param, x, x_pos, sol, desc, hes, cov, cor
     call minimiser(info, symm, sol, flag, desc, &
          hes, cov, cor, chisqrd, nlposterior)
     if (info /= '') then
        print *,'*****'
        print *,trim(info)
        print *, '*****'
     end if

     if (flag < 4) then
        degfreedom = useful_vis + useful_amp + useful_cp - size(sol,1)
        normchisqrd = chisqrd/degfreedom

        print *,' '
        print *,'negative log posterior =',real(nlposterior)
        print *,' '
        print *,'           chi squared =',real(chisqrd) 
        print *,'    degrees of freedom =',degfreedom
        print *,'chi sqrd / deg freedom =',real(normchisqrd)

        print *,' '
        print *,'solution details:'
        print *,'                                                ', &
             'fitted      hessian'
        print *,'num  parameter name                             ', &
             ' value        error'
        do i = 1, size(sol,1)
           write(*,62) (i, desc(i), sol(i,:))
        end do
62      format(' (', i2, ') ', A35, 1x, f13.6, 1x, f12.6) 

        length = size(hes,1)
        print *,' '
        print *,'hessian matrix'
        do i = 1, length
           write(*,64) (hes(i,:))
        end do

        print *,' '
        print *,'covariance matrix'
        do i = 1, length
           write(*,64) (cov(i,:))
        end do

        print *,' '
        print *,'correlation matrix'
        do i = 1, length
           write(*,64) (cor(i,:))
        end do

64      format(1x,10(e11.4,1x))

        print *,' '
        print *,spacer_line

        !----------------------------------------------------------------------
        !plot

        top_title = trim(source)//' - final model: '//trim(model_name)
        if (useful_vis > 0) &
             call plot_vis(model_spec, fit_param, symm, 'Baseline /M\gl', &
             'Squared Visibility', top_title)
        if (useful_cp > 0) &
             call plot_triple(model_spec, fit_param, 'Longest baseline /M\gl',&
             'Closure phase /'//char(176), top_title)
        call pgend
     end if

     !-------------------------------------------------------------------------
     !Deallocate model/fitting storage
     call free_model() !model_*
     call free_fit() !fit_param, x_pos, x_info
     if (allocated(desc)) deallocate(desc)
     if (allocated(sol)) deallocate(sol)
     if (allocated(hes)) deallocate(hes)
     if (allocated(cov)) deallocate(cov)
     if (allocated(cor)) deallocate(cor)

  end do model_loop

  !----------------------------------------------------------------------------
  !Deallocate data storage
  if (allocated(vis_data)) deallocate(vis_data)
  if (allocated(triple_data)) deallocate(triple_data)
  if (allocated(wavebands)) deallocate(wavebands)

end program Main


!==============================================================================
integer function myhandler(sig,code,context) 
  integer sig, code, context(5)
  myhandler = 0 !avoids compiler warning
  call abort()
end function myhandler

