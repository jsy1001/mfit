!$Id: model.f90,v 1.20 2013/02/18 13:18:16 jsy1001 Exp $

module Model

  use Maths
  use Search

  implicit none

  include 'fftw3.f'

  private

  !public subroutines contained
  !
  !read_model
  !model_valid
  !model_nvar
  !model_getvar
  !print_model
  !write_model
  !free_model

  public :: read_model, model_valid, model_nvar, model_getvar
  public :: print_model, write_model, free_model


  !public module variables contained:
  public :: model_spec, model_pos_relto, model_wldep, nwave, model_wb
  public :: model_param, model_prior, model_limits, max_order
  public :: model_desc, model_desc_len, model_name, symm_model
  public :: nxsiz, clv_mdiam, clv_mvis, clv_mbase

  !! Specifies nature of model cpts -
  !! holds 3 items per component: 1 name, 2 shape and 3 LD type
  character(len=128), allocatable :: model_spec(:,:)

  !! If model_pos_relto(icpt) is a valid model component, r & theta for
  !! component icpt are treated as multipliers for r & theta for the specified
  !! (other) component
  integer, allocatable :: model_pos_relto(:)

  !! model_wldep(i) is unity if position, flux, shape_param,
  !! ld_param respectively is fn(waveband), else zero
  integer :: model_wldep(4)
  !! Model has nwave values for each wavelength-dependent parameter
  integer :: nwave

  !! Wavebands for wavelength-dependent model parameters
  !! (if supplied to read_model) and CLV
  real, allocatable :: model_wb(:, :)

  !! Initial model parameters (both variable and non-variable)
  !! for each cpt (1st axis): 2nd axis is
  !! LD order
  !! {r, theta(t0), t0, d(theta)/dt} possibly repeated for each waveband
  !! {brightness}                    possibly repeated for each waveband
  !! {major axis, phi, epsilon}      possibly repeated for each waveband
  !!  {LD coeffs 1-max_order}         possibly repeated for each waveband
  !! EXCEPT for 'ld_type two-layer', which has
  !!  {LD coeffs 1-5}                 coeff 5 possibly repeated for each waveband
  double precision, allocatable :: model_param(:,:)
  integer max_order !at least one component has this many LD coeffs

  !! Gaussian prior widths for parameters; prior centres are in model_param
  double precision, allocatable :: model_prior(:,:)

  !! Limits on legal values for parameters (both variable and non-variable)
  double precision, allocatable :: model_limits(:,:,:)

  integer, parameter :: model_desc_len = 55

  !! Descriptions of model parameters
  character(len=model_desc_len), allocatable :: model_desc(:, :)

  !! Name from model file (2 words max)
  character(len=128) :: model_name

  !! Centrosymmetric model?
  logical :: symm_model

  !numerical CLV data
  real, allocatable :: clv_rad(:)
  real, allocatable :: clv_inten(:, :)
  integer :: nclv, clvcomp
  !model visibilities calculated from numerical CLV
  integer, parameter :: nxsiz = 4096, modelsize = 120
  real :: clv_mbase(nxsiz+1)
  real, allocatable :: clv_mvis(:, :)
  real :: clv_mdiam

contains

  !============================================================================

  !! Read consecutive lines which are blank or start with specified character
  !! Backspace on reading first non-blank, non-comment line
  function skip_lines(iunit, mark, eof)

    !subroutine arguments
    integer, intent(in) :: iunit
    character (len=1), intent(in) :: mark
    logical, intent(out) :: eof
    integer skip_lines

    !local variables
    !!character(len=256) :: line
    character(len=80) :: line

    skip_lines = 0
    eof = .false.
    do
       read(iunit,'(a)', end=1, err=2) line
       if(len_trim(line) > 0 .and. line(1:1) /= mark) then
          backspace(iunit)
          return
       end if
       skip_lines = skip_lines + 1
       print *, line
    end do

1   eof = .true.
    return

    ! I/O error, rely on calling routine to re-detect
2   return

  end function skip_lines

  !============================================================================

  !! Read model files
  !! (Re-)Allocates and assigns to module variables
  subroutine read_model(info, file_name, wavebands)

    !subroutine arguments
    character(len=*), intent(out) :: info !! Error message
    character(len=*), intent(in) :: file_name !! Filename to read
    !! Wavebands being used in data
    double precision, intent(in), optional :: wavebands(:,:) 

    !local variables
    logical eof
    character(len=1) :: comment
    character(len=2) :: numbers(10)
    character(len=32) :: keyw, cpt
    character(len=32), allocatable :: wb(:)
    character(len=128) :: clv_filename, shape_type, ld_type
    character(len=256) :: line
    integer :: i, j, k, n, line_no, iwave, loc, relto
    integer :: comps, order, pos, npar, ipar, nread, nn

    !check for zero length filename
    if (file_name == '') then
       info = 'blank filename'
       return
    end if

    !read and count number of components and number of lines
    !also determine which parameters are wavelength-dependent (model_wldep),
    !how many wavebands values are supplied for (nwave),
    !and the maximum ld_order (max_order)
    open (unit=11, err=91, status='old', action='read', file=file_name)
    comment = '!'
    comps = 0
    model_wldep = 0
    nwave = -1
    max_order = 0
    line_no = 1
    pass1:    do
       line_no = line_no + skip_lines(11, comment, eof)
       if(eof) exit pass1
       read (11, '(a)', err=92) line
       line_no = line_no + 1
       if (len_trim(line) == 0) cycle !blank line
       read (line, *) keyw !read keyword & qualifiers
       if (keyw == 'component') comps = comps+1
       if (keyw(:8) == 'position' .and. index(keyw, '#fofwave') /= 0) then
          model_wldep(1) = 1
          ! qualifiers could be any combination of #fofwave, #rotate, #reltoN
          if (index(keyw, '#rotate') /= 0) then
             !{r, theta(t0), t0, d(theta)/dt} poss. repeated for each waveband
             nn = (countsym(line)-1)/4
          else
             !file contains {r, theta} poss. repeated for each waveband
             nn = (countsym(line)-1)/2
          end if
          if (nwave == -1) then
             nwave = nn
          else if (nn /= nwave) then
             write(info, *) 'line ', line_no, &
                  ': expecting values for ', nwave, ' wavebands, got ', nn
             close(11)
             return
          end if
       else if (keyw == 'flux#fofwave') then
          model_wldep(2) = 1
          nn = countsym(line)-1
          if (nwave == -1) then
             nwave = nn
          else if (nn /= nwave) then
             write(info, *) 'line ', line_no, &
                  ': expecting values for ', nwave, ' wavebands, got ', nn
             close(11)
             return
          end if
       else if (keyw == 'shape_type') then
          read (line, *, end=95) keyw, shape_type
       else if (keyw == 'shape_param#fofwave') then
          select case (shape_type)
          case ('point')
             print *, 'cannot have wavelength-dependent point: no shape parameters'
          case ('disc')
             model_wldep(3) = 1
             nn = countsym(line)-1
          case ('ellipse')
             model_wldep(3) = 1
             nn = (countsym(line)-1)/3
          case default
             info = 'invalid shape type'
             close (11)
             return
          end select
          if (nwave == -1) then
             nwave = nn
          else if (nn /= nwave) then
             write(info, *) 'line ', line_no, &
                  ': expecting values for ', nwave, ' wavebands, got ', nn
             close(11)
             return
          end if
       else if (keyw == 'ld_type') then
          read (line, *, end=95) keyw, ld_type
          !Read LD order (preset for some LD type cases)
          !Non-integer will cause error in read statement
          select case (trim(ld_type))
          case ('uniform','gaussian','thin-ring')
             read (11,*,err=95,end=95) keyw
             order = 0
          case ('square-root')
             read (11,*,err=95,end=95) keyw
             order = 2
          case ('hestroffer')
             read (11,*,err=95,end=95) keyw
             order = 1
          case ('taylor','gauss-hermite') 
             read (11,*,err=95,end=95) keyw, order
          case ('two-layer')
             read (11,*,err=95,end=95) keyw
             order = 5
          case default
             if (ld_type(1:1) == '<') then
                read (11,*,err=95,end=95) keyw
                order = 0
             else
                !not <filename>
                info = 'invalid limb darkening type'
                close (11)
                return
             end if
          end select
          line_no = line_no + 1
          if (order > max_order) max_order = order
       else if (keyw == 'ld_param#fofwave') then
          model_wldep(4) = 1
          if (trim(ld_type) == 'two-layer') then
             nn = countsym(line)-5
          else
             nn = (countsym(line)-1)/order
          end if
          if (nwave == -1) then
             nwave = nn
          else if (nn /= nwave) then
             write(info, *) 'line ', line_no, &
                  ': expecting values for ', nwave, ' wavebands, got ', nn
             close(11)
             return
          end if
       end if
    end do pass1

    close(11)

    !check if legal number of components
    if ((comps<1).or.(comps>10)) then
       info = 'illegal number of components present'
       return
    end if

    !check wavebands
    if (nwave /= -1) then
       !wavelength-dependent model
       if (present(wavebands)) then
          if (nwave < size(wavebands, 1)) then
             write(info, '(i3, a, i3, a)') size(wavebands, 1), &
                  ' waveband(s) in data, but only ', nwave, &
                  ' in wavelength-dependent model'
             return
          else if (nwave > size(wavebands, 1)) then
             nwave = size(wavebands, 1)
             print *, 'Reading first ', nwave, ' waveband(s) from model'
          endif
          !wavebands not present if read_model called from modelplot
       end if
    else
       !could still have wavelength-dependent CLV
       if (present(wavebands)) then
          nwave = size(wavebands, 1)
       else
          nwave = 1
       end if
    end if

    !allocate model data arrays
    !each of position, flux, shape_param, ld_param can be either
    !wavelength-independent (default) or wavelength-dependent, FOR ALL
    !components
    !shape_type, ld_type, ld_order of a component are the same for all wavebands
    call free_model()
    allocate(model_spec(comps,3))
    if (present(wavebands)) then
       allocate(model_wb(size(wavebands, 1), 2))
       model_wb = wavebands
    end if
    npar = 1 !LD order
    npar = npar + 4*(1 + model_wldep(1)*(nwave-1)) !r, theta(t0), t0, d(theta)/dt
    npar = npar + 1*(1 + model_wldep(2)*(nwave-1)) !flux
    npar = npar + 3*(1 + model_wldep(3)*(nwave-1)) !major axis, phi, epsilon
    npar = npar + max_order*(1 + model_wldep(4)*(nwave-1)) !LD coeffs
    allocate(model_param(comps,npar))
    allocate(model_prior(comps,npar))
    allocate(model_limits(comps,npar,2))
    allocate(model_desc(comps,npar))
    allocate(model_pos_relto(comps))
    model_param = 0D0
    model_prior = 0D0
    clvcomp = 0 !no model component yet has numerical CLV

    !assign strings describing model wavebands
    allocate(wb(nwave))
    do iwave = 1, nwave
       if (allocated(model_wb)) then
          write(wb(iwave), '(a, f7.1, a, f6.1, a)') ' [', model_wb(iwave,1), &
               '/', model_wb(iwave,2), ']'
       else
          write(wb(iwave), '(a, i3, a)') '[', iwave, ']'
       end if
    end do

    !Limits array stores the lower and upper acceptable limits
    !on the parameter values, it is passed into the minimising routines
    !later to ensure parameter values stay legal. Model values not within
    !these limits trigger error on reading the model
    !Remainder filled in as we read model file
    numbers = (/' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10'/)
    do j = 1, comps
       model_limits(j,1,1:2) = dble((/0,10/))
       cpt = 'cpt '// trim(adjustl(numbers(j))) // ', '
       model_desc(j, 1) = trim(cpt) // ' LD order'
    end do

    !read spec and param info for components
    open (unit=11, action='read', file=file_name)
    j = 0 !component counter
    model_name = 'model' !default if "source" keyword missing
    pass2: do
       read (11,*,err=94,end=12) keyw
       if (keyw == 'source') then
          backspace(11)
          read (11,'(a)',err=95,end=95) line
          line = line(index(line, 'source')+6:) !strip 'source'
          pos = verify(line, ' '//achar(9))
          model_name = line(pos:) !strip leading spaces and tabs
       end if
       if (keyw == 'component') then
          j = j+1
          cpt = 'cpt '// trim(adjustl(numbers(j))) // ', '

          !read component name
          n = skip_lines(11, comment, eof)
          read (11,*,err=95,end=95) keyw, model_spec(j,1)
          if (keyw /= 'name') goto 96

          !read shape type
          n = skip_lines(11, comment, eof)
          read (11,*,err=95,end=95) keyw, model_spec(j,2)
          if (keyw /= 'shape_type') goto 96

          !read LD type (set as uniform for point case)
          select case (trim(model_spec(j,2)))
          case ('point')
             read (11,*,err=95) keyw
             model_spec(j,3) = 'uniform'
          case default
             n = skip_lines(11, comment, eof)
             read (11,'(a)',err=95,end=95) line
             read (line,*,err=95,end=95) keyw
             !from last non-trailing space or tab character:
             pos = scan(trim(line), ' '//achar(9), back=.true.) + 1
             model_spec(j,3) = line(pos:)
          end select
          if (keyw /= 'ld_type') goto 96

          !Read LD order (preset for some LD type cases)
          !Non-integer will cause error in read statement
          n = skip_lines(11, comment, eof)
          select case (trim(model_spec(j,3)))
          case ('uniform','gaussian','thin-ring')
             read (11,*,err=95,end=95) keyw
             order = 0
          case ('square-root')
             read (11,*,err=95,end=95) keyw
             order = 2
          case ('hestroffer')
             read (11,*,err=95,end=95) keyw
             order = 1
          case ('taylor','gauss-hermite') 
             read (11,*,err=95,end=95) keyw, order
          case ('two-layer')
             read (11,*,err=95,end=95) keyw
             order = 5
          case default
             if (model_spec(j,3)(1:1) == '<') then
                if (clvcomp /= 0) then
                   info = 'only one numerical CLV component allowed'
                   close (11)
                   return
                end if
                clvcomp = j
                clv_filename = model_spec(j,3)(2:(len_trim(model_spec(j,3))-1))
                if (allocated(model_wb)) then
                   call read_clv(info, clv_filename, model_wb)
                else
                   call read_clv(info, clv_filename)
                end if
                if (info /= '') then
                   close (11)
                   return
                end if
                read (11,*,err=95,end=95) keyw
                order = 0
             else
                !not <filename>
                info = 'invalid limb darkening type'
                close (11)
                return
             end if
          end select
          if (keyw /= 'ld_order') goto 96
          model_param(j,1) = dble(order)

          !read position r and theta etc.
          n = skip_lines(11, comment, eof)
          read (11,*,err=95,end=95) keyw
          if (keyw(:8) /= 'position') goto 96
          loc = index(keyw, '#relto')
          if (loc /= 0) then
             read (keyw(loc+6:), *) model_pos_relto(j)
          else
             model_pos_relto(j) = 0
          end if
          relto = model_pos_relto(j)
          !loop over wavebands, assign parameter limits and descriptions
          do k = 1, 1+model_wldep(1)*(nwave-1)
             model_limits(j,2+(k-1)*4,1:2) = dble((/0, 500/))
             model_limits(j,3+(k-1)*4,1:2) = dble((/-720, 720/))
             model_limits(j,4+(k-1)*4,1:2) = dble((/10000, 100000/))
             model_limits(j,5+(k-1)*4,1:2) = dble((/-1000, 1000/))
             if (relto >= 1 .and. relto <= comps) then
                !allow negative radius multipliers
                model_limits(j,2+(k-1)*4,1:2) = dble((/-500, 500/))
             end if
             !defaults for t0 and d(theta)/dt
             model_param(j, 4+(k-1)*4) = 10000D0
             model_param(j, 5+(k-1)*4) = 0D0
             !descriptions
             if (model_wldep(1) == 1) then
                if (relto >= 1 .and. relto <= comps) then
                   model_desc(j, 2+(k-1)*4) = trim(cpt) // trim(wb(k)) &
                        // ' position radius as multiple of cpt ' &
                        // trim(adjustl(numbers(relto))) // ' radius'
                   model_desc(j, 3+(k-1)*4) = trim(cpt) // trim(wb(k)) &
                        // ' position angle(t0) as multiple of cpt ' &
                        // trim(adjustl(numbers(relto))) // ' angle'
                else
                   model_desc(j, 2+(k-1)*4) = trim(cpt) // trim(wb(k)) &
                        // ' position radius (mas)'
                   model_desc(j, 3+(k-1)*4) = trim(cpt) // trim(wb(k)) &
                        // ' position angle(t0) (deg)'
                end if
                model_desc(j, 4+(k-1)*4) = trim(cpt) // trim(wb(k)) &
                     // ' t0 (MJD)'
                model_desc(j, 5+(k-1)*4) = trim(cpt) // trim(wb(k)) &
                     // ' d(PA)/dt (deg/day)'
             else
                if (relto >= 1 .and. relto <= comps) then
                   model_desc(j, 2+(k-1)*4) = trim(cpt) &
                        // ' position radius as multiple of cpt ' &
                        // trim(adjustl(numbers(relto))) // ' radius'
                   model_desc(j, 3+(k-1)*4) = trim(cpt) &
                        // ' position angle(t0) as multiple of cpt ' &
                        // trim(adjustl(numbers(relto))) // ' angle'
                else
                   model_desc(j, 2+(k-1)*4) = trim(cpt) // ' position radius (mas)'
                   model_desc(j, 3+(k-1)*4) = trim(cpt) // ' position angle(t0) (deg)'
                end if
                model_desc(j, 4+(k-1)*4) = trim(cpt) // ' t0 (MJD)'
                model_desc(j, 5+(k-1)*4) = trim(cpt) // ' d(PA)/dt (deg/day)'
             end if
          end do
          ! qualifiers could be any combination of #fofwave, #rotate, #reltoN
          if (index(keyw, '#rotate') /= 0) then
             !{r, theta(t0), t0, d(theta)/dt} poss. repeated for each waveband
             backspace(11)
             nread = 4*(1 + model_wldep(1)*(nwave-1))
             read (11,*,err=95,end=95) keyw, model_param(j, 2:1+nread)
             n = skip_lines(11, comment, eof)
             read (11,*,err=95,end=95) keyw, model_prior(j, 2:1+nread)
             if (keyw /= 'position_prior') goto 96
          else
             !file contains {r, theta} poss. repeated for each waveband
             backspace(11)
             read (11,*,err=95,end=95) keyw, &
                  (model_param(j, 2+(k-1)*4:3+(k-1)*4), &
                  k=1,(1+model_wldep(1)*(nwave-1)))
             n = skip_lines(11, comment, eof)
             read (11,*,err=95,end=95) keyw, &
                  (model_prior(j, 2+(k-1)*4:3+(k-1)*4), &
                  k=1,(1+model_wldep(1)*(nwave-1)))
             if (keyw /= 'position_prior') goto 96
          end if

          !read flux
          ipar = 6+4*model_wldep(1)*(nwave-1)
          n = skip_lines(11, comment, eof)
          read (11,*,err=95,end=95) keyw, &
               model_param(j, ipar:ipar+model_wldep(2)*(nwave-1))
          if (keyw(:4) /= 'flux') goto 96
          do k = 1, 1+model_wldep(2)*(nwave-1)
             model_limits(j,ipar+(k-1),1:2) = dble((/-100, 100/))
             if (model_wldep(2) == 1) then
                model_desc(j, ipar+(k-1)) = trim(cpt) // trim(wb(k)) // ' flux (arb units)'
             else
                model_desc(j, ipar+(k-1)) = trim(cpt) // ' flux (arb units)'
             end if
          end do

          !read flux prior, must be non-negative
          n = skip_lines(11, comment, eof)
          read (11,*,err=95,end=95) keyw, &
               model_prior(j, ipar:ipar+model_wldep(2)*(nwave-1))
          if (keyw /= 'flux_prior') goto 96  

          !Read shape parameters a, phi, epsilon depending on
          !shape type. Also read priors
          !Shape a and epsilon must be non-negative
          !All priors must be non-negative
          ipar = 7 + (4*model_wldep(1) + model_wldep(2))*(nwave-1)
          do k=1, 1+model_wldep(3)*(nwave-1)
             model_param(j,ipar+(k-1)*3) = 0D0
             model_param(j,ipar+(k-1)*3+1) = 0D0
             model_param(j,ipar+(k-1)*3+2) = 1D0
             model_limits(j,ipar+(k-1)*3,1:2) = dble((/0, 1000/))
             model_limits(j,ipar+(k-1)*3+1,1:2) = dble((/-720, 720/))
             model_limits(j,ipar+(k-1)*3+2,1:2) = dble((/0, 10/))
             if (model_wldep(3) == 1) then
                model_desc(j,ipar+(k-1)*3) = trim(cpt) // trim(wb(k)) // ' major axis (mas)'
                model_desc(j,ipar+(k-1)*3+1) = trim(cpt) // trim(wb(k)) // ' orientation (deg)'
                model_desc(j,ipar+(k-1)*3+2) = trim(cpt) // trim(wb(k)) // ' ellipticity'
             else
                model_desc(j,ipar+(k-1)*3) = trim(cpt) // ' major axis (mas)'
                model_desc(j,ipar+(k-1)*3+1) = trim(cpt) // ' orientation (deg)'
                model_desc(j,ipar+(k-1)*3+2) = trim(cpt) // ' ellipticity'
             end if
          end do
          n = skip_lines(11, comment, eof)
          select case (trim(model_spec(j,2)))
          case ('point')
             read (11,*,err=95,end=95) keyw
             if (keyw(:11) /= 'shape_param') goto 96
             n = skip_lines(11, comment, eof)
             read (11,*,err=95,end=95) keyw
          case ('disc')
             read (11,*,err=95,end=95) keyw, (model_param(j,ipar+(k-1)*3), &
                  k=1,(1+model_wldep(3)*(nwave-1)))
             if (keyw(:11) /= 'shape_param') goto 96
             n = skip_lines(11, comment, eof)
             read (11,*,err=95,end=95) keyw, (model_prior(j,ipar+(k-1)*3), &
                  k=1,(1+model_wldep(3)*(nwave-1)))
          case ('ellipse')
             nread = 3*(1 + model_wldep(3)*(nwave-1))
             read (11,*,err=95,end=95) keyw, model_param(j,ipar:ipar+nread-1)
             if (keyw(:11) /= 'shape_param') goto 96
             n = skip_lines(11, comment, eof)
             read (11,*,err=95,end=95) keyw, model_prior(j,ipar:ipar+nread-1)
          case default
             info = 'invalid shape type'
             close (11)
             return
          end select
          if (keyw /= 'shape_param_prior') goto 96

          !read LD parameters
          n = skip_lines(11, comment, eof)
          ipar = 10 + (4*model_wldep(1) + model_wldep(2) &
               + 3*model_wldep(3))*(nwave-1)
          select case (order)
          case (0)
             read (11,*,err=95,end=95) keyw
             if (keyw(:8) /= 'ld_param') goto 96
             n = skip_lines(11, comment, eof)
             read (11,*,err=95,end=95) keyw
          case default
             if (trim(model_spec(j,3)) == 'two-layer') then
                nread = 5 + model_wldep(4)*(nwave-1)
             else
                nread = order*(1 + model_wldep(4)*(nwave-1))
             end if
             read (11,*,err=95,end=95) keyw, &
                  model_param(j, ipar:ipar+nread-1)
             if (keyw(:8) /= 'ld_param') goto 96
             n = skip_lines(11, comment, eof)
             read (11,*,err=95,end=95) keyw, &
                  model_prior(j, ipar:ipar+nread-1)
          end select
          if (trim(model_spec(j,3)) == 'two-layer') then
             !different limits for LD parameters if 2-layer component
             model_limits(j,ipar:ipar+2,1) = 0D0
             model_limits(j,ipar:ipar+2,2) = 1D5
             model_limits(j,ipar+3,1) = 1D0
             model_limits(j,ipar+3,2) = 100D0
             do k = 1, 1+model_wldep(4)*(nwave-1)
                model_limits(j,ipar+3+k,1) = 0D0
                model_limits(j,ipar+3+k,2) = 100D0
             end do
          else
             model_limits(j, ipar: &
                  ipar+max_order*(1 + model_wldep(4)*(nwave-1))-1,1) = -100D0
             model_limits(j, ipar: &
                  ipar+max_order*(1 + model_wldep(4)*(nwave-1))-1,2) = 100D0
          end if
          !LD param descriptions
          if (trim(model_spec(j,3)) == 'two-layer') then
             model_desc(j,ipar) = trim(cpt) // ' NStep '
             model_desc(j,ipar+1) = trim(cpt) // ' Temp1 /K '
             model_desc(j,ipar+2) = trim(cpt) // ' Temp2 /K '
             model_desc(j,ipar+3) = trim(cpt) // ' R2/R1 '
             if (model_wldep(4) == 1) then
                do iwave = 1, nwave
                   model_desc(j,ipar+3+iwave) = trim(cpt) &
                        // trim(wb(iwave)) // ' Optical depth '
                end do
             else
                model_desc(j,ipar+4) = trim(cpt) // ' Optical depth '
             end if
          else
             do k = 1, int(model_param(j, 1)) !order for this cpt
                if (model_wldep(4) == 1) then
                   do iwave = 1, nwave
                      model_desc(j, ipar+max_order*(iwave-1)+k-1) = trim(cpt) &
                           // trim(wb(iwave)) // ' LD parameter ' // adjustl(numbers(k))
                   end do
                else
                   model_desc(j, ipar+k-1) = trim(cpt) // ' LD parameter ' &
                        // adjustl(numbers(k))
                end if
             end do
          end if

          !amend model_limits:
          !alpha(1) in hestroffer model must be > 0
          if (trim(model_spec(j,3)) == 'hestroffer') then
             if (model_wldep(4) == 1) then
                do iwave = 1, nwave
                   k = 10+(4*model_wldep(1)+model_wldep(2) &
                        +3*model_wldep(3))*(nwave-1) &
                        + max_order*model_wldep(4)*(iwave-1)
                   model_limits(j,k,1) = 0D0
                end do
             else
                k = 10+(4*model_wldep(1)+model_wldep(2) &
                     +3*model_wldep(3))*(nwave-1)
                model_limits(j,k,1) = 0D0
             end if
          end if

          !check to ensure everything inside limits
          do k=1, size(model_limits,2)
             if (model_param(j,k) < model_limits(j,k,1)) goto 100
             if (model_param(j,k) > model_limits(j,k,2)) goto 100
          end do

          !if numerical CLV, calculate visibilities for initial guess diameter
          !wavelength-dep CLV with wavelength-dep diameter doesn't make sense
          !- this is trapped in minimiser()
          if (model_spec(j,3)(1:1) == '<') &
               call calcvis_clv(model_param(j, &
               7+(4*model_wldep(1)+model_wldep(2))*(nwave-1)))

          if (keyw /= 'ld_param_prior') goto 96

       end if
    end do pass2

12  continue 
    close (11)

    !is this a centrosymmetric model?
    symm_model = .true.
    do i = 1, size(model_param, 1)
       !need all epsilon fixed at unity
       ipar = 9 + (4*model_wldep(1) + model_wldep(2))*(nwave-1)
       do j = 1, 1+model_wldep(3)*(nwave-1)
          if (.not.((model_param(i,ipar+(j-1)*3)==1D0) &
               .and.(model_prior(i,ipar+(j-1)*3)==0D0))) symm_model = .false.
       end do
       !and all r fixed at zero
       do j = 1, 1+model_wldep(1)*(nwave-1)
          if (.not.((model_param(i,2+(j-1)*4)==0D0) &
               .and.(model_prior(i,2+(j-1)*4)==0D0))) symm_model = .false.
       end do
    end do

    deallocate(wb) !free local storage
    return

    !error trapping 
91  info = 'cannot open file '//trim(file_name)
    return
92  info = 'cannot read from file '//trim(file_name)
    close(11)
    return
94  info = 'unknown read problem - possible invalid file format'
    close (11)
    return
95  info = 'unknown read problem whilst reading component data'
    close (11)
    return
96  info = 'invalid keyword or keyword order'
    close (11)
    return
100 info = 'illegal parameter value detected'
    close (11)
    return

  end subroutine read_model

  !============================================================================

  !! Are variable parameters sensible?
  function model_valid(info, force_symm)

    logical :: model_valid

    !function arguments
    character(len=*), intent(out) :: info !! Error message if model invalid
    logical, intent(in) :: force_symm !! True if expecting a symmetric model

    !local variables
    integer :: nvar, i, iwave, ipar

    info = ''
    model_valid = .true.

    !must have at least 1 freedom to minimise with
    call model_nvar(nvar)
    if (nvar == 0) model_valid = .false.

    !for single component model cannot vary r/theta (change in position
    !only constant phase offset) or B (flux is normalised anyway) 
    if (size(model_param,1) == 1) then
       do iwave = 1, nwave
          if (model_prior(1, 2+4*model_wldep(1)*(iwave-1)) /= 0D0) &
               model_valid = .false.
          if (model_prior(1, 3+4*model_wldep(1)*(iwave-1)) /= 0D0) &
               model_valid = .false.
          if (model_prior(1, &
               6+4*model_wldep(1)*(nwave-1)+model_wldep(2)*(iwave-1)) /= 0D0) &
               model_valid = .false.
       end do
    end if

    !for any component cannot vary theta if r is zero and not free to vary
    !(position angle has no effect if position radius fixed at zero)
    do i = 1, size(model_param,1)
       do iwave = 1, nwave
          ipar = 2+4*model_wldep(1)*(iwave-1)
          if ((model_param(i,ipar)==0D0) .and. (model_prior(i,ipar)==0D0) &
               .and. (model_prior(i,ipar+1)/=0D0)) model_valid = .false.
       end do
    end do

    !for any component cannot vary phi if epsilon is unity and not free to vary
    !(orientation is meaningless if ellipse reduced to circular disc)
    do i = 1, size(model_param,1)
       do iwave = 1, nwave
          ipar = 7 + (4*model_wldep(1)+model_wldep(2))*(nwave-1) &
               + model_wldep(3)*(iwave-1)
          if ((model_param(i,ipar+2)==1D0) .and. (model_prior(i,ipar+2)==0D0) &
               .and. (model_prior(i,ipar+1)/=0D0)) model_valid = .false.
       end do
    end do

    !doesn't make sense to have wavelength-dependent numerical CLV with
    !wavelength-dependent diameter
    if (allocated(clv_mvis)) then
       if (size(clv_mvis, 2) > 1 .and. model_wldep(3) == 1) &
            model_valid = .false.
    end if

    if (.not. model_valid) &
         info = 'illegal freedom(s) in model'

    !if centrosymmetric model is forced then must have 
    !eccentricity epsilon fixed to be unity and position radius fixed at zero
    if (force_symm .and. .not. symm_model) then
       model_valid = .false.
       info = 'for vis/nvis data must have guaranteed symmetric model'
    end if

  end function model_valid

  !============================================================================

  !! Get number of variable parameters in model
  subroutine model_nvar(nvar)

    !subroutine arguments
    integer, intent(out) :: nvar

    !local variables
    integer :: i, j

    nvar = 0
    do i = 1, size(model_param,1)
       do j = 1, size(model_param,2)
          if (model_prior(i,j) /= 0D0)  nvar = nvar + 1
       end do
    end do

  end subroutine model_nvar

  !============================================================================

  !! Extract locations and descriptions for variable model parameters
  subroutine model_getvar(nvar, pos, desc)

    !subroutine arguments
     !! Number of variables, e.g. from model_nvar
    integer, intent(in) :: nvar
    !! pos(i,:) gives indices of i'th variable parameter in model_param
    !! and model_prior arrays
    integer, intent(out) :: pos(nvar,2)
    !! Variable descriptions
    character(len=model_desc_len), intent(out) :: desc(nvar)

    !local variables
    integer :: i, j, ivar

    ivar = 0
    do i = 1, size(model_param,1)
       do j = 1, size(model_param,2)
          if (model_prior(i,j) /= 0D0) then
             ivar = ivar + 1
             pos(ivar,1) = i
             pos(ivar,2) = j
             desc(ivar) = model_desc(i,j)
          end if
       end do
    end do

  end subroutine model_getvar

  !============================================================================

  !! Print model details
  subroutine print_model()

    !local variables
    integer :: i, ipar, jpar, iwave, order
    character(len=32), allocatable :: wb(:)
    character(len=512) :: line1, line2, ch1, ch2

    allocate(wb(nwave))
    do iwave = 1, nwave
       if (allocated(model_wb)) then
          write(wb(iwave), '(a, f7.1, a, f6.1, a)') ' [', model_wb(iwave,1), &
               '/', model_wb(iwave,2), ']'
       else
          write(wb(iwave), '(a, i3, a)') '[', iwave, ']'
       end if
    end do

    do i = 1, size(model_param,1)
       print *,' '
       print *,'component:',i
       print *,'name     : ',trim(model_spec(i,1))
       order = int(model_param(i,1))
       if (model_spec(i,3)(1:1) == '<') then
          print *,'LD type  : ',trim(model_spec(i,3)), &
               ' :',nclv,'points'
       else
          print *,'LD type  : ',trim(model_spec(i,3)), &
               ' of order',order
       end if
       print *,'shape    : ',trim(model_spec(i,2))
       !position
       ipar = 6 + 4*model_wldep(1)*(nwave-1)
       jpar = 7 + (4*model_wldep(1) + model_wldep(2))*(nwave-1)
       print 10,'position :',real(model_param(i,2:ipar-1))
       print 10,'    prior:',real(model_prior(i,2:ipar-1))
       !flux
       line1 = 'flux     :'
       line2 = '    prior:'
       do iwave = 1, nwave
          if (model_wldep(2) == 1) then
             line1 = trim(line1) // trim(wb(iwave)) //  ':'
             line2 = trim(line2) // trim(wb(iwave)) //  ':'
          end if
          write (ch1, '(f10.3)') real(model_param(i,ipar+iwave-1))
          write (ch2, '(f10.3)') real(model_prior(i,ipar+iwave-1))
          line1 = trim(line1) // trim(ch1)
          line2 = trim(line2) // trim(ch2)
          if (model_wldep(2) == 0) exit
       end do
       print '(1x, a)', trim(line1)
       print '(1x, a)', trim(line2)
       !shape
       ipar = jpar
       line1 = 'shape par:'
       line2 = '    prior:'
       do iwave = 1, nwave
          if (model_wldep(3) == 1) then
             line1 = trim(line1) // trim(wb(iwave)) //  ':'
             line2 = trim(line2) // trim(wb(iwave)) //  ':'
          end if
          write (ch1, '(3f10.3)') &
               real(model_param(i,ipar+3*(iwave-1):ipar+3*(iwave-1)+2))
          write (ch2, '(3f10.3)') &
               real(model_prior(i,ipar+3*(iwave-1):ipar+3*(iwave-1)+2))
          line1 = trim(line1) // trim(ch1)
          line2 = trim(line2) // trim(ch2)
          if (model_wldep(3) == 0) exit
       end do
       print '(1x, a)', trim(line1)
       print '(1x, a)', trim(line2)
       !ld params
       if (order >= 1) then
          if (trim(model_spec(i,3)) == 'two-layer') then
             ipar = 10 + (4*model_wldep(1) + model_wldep(2) &
                  + 3*model_wldep(3))*(nwave-1)
             write (line1, '(a, 4f9.3)') 'LD params:', &
                  real(model_param(i,ipar:ipar+3))
             write (line2, '(a, 4f9.3)') '    prior:', &
                  real(model_prior(i,ipar:ipar+3))
             do iwave = 1, nwave
                if (model_wldep(4) == 1) then
                   line1 = trim(line1) // trim(wb(iwave)) //  ':'
                   line2 = trim(line2) // trim(wb(iwave)) //  ':'
                end if
                write (ch1, '(f7.3)') real(model_param(i,ipar+3+iwave))
                write (ch2, '(f7.3)') real(model_prior(i,ipar+3+iwave))
                line1 = trim(line1) // trim(ch1)
                line2 = trim(line2) // trim(ch2)
                if (model_wldep(4) == 0) exit
             end do
          else
             ipar = 10 + (4*model_wldep(1) + model_wldep(2) &
                  + 3*model_wldep(3))*(nwave-1)
             line1 = 'LD params:'
             line2 = '    prior:'
             do iwave = 1, nwave
                if (model_wldep(4) == 1) then
                   line1 = trim(line1) // trim(wb(iwave)) //  ':'
                   line2 = trim(line2) // trim(wb(iwave)) //  ':'
                end if
                write (ch1, '(10f7.3)') real(model_param(i,ipar+max_order*(iwave-1):ipar+max_order*(iwave-1)+order-1))
                write (ch2, '(10f7.3)') real(model_prior(i,ipar+max_order*(iwave-1):ipar+max_order*(iwave-1)+order-1))
                line1 = trim(line1) // trim(ch1)
                line2 = trim(line2) // trim(ch2)
                if (model_wldep(4) == 0) exit
             end do
          end if
          print '(1x, a)', trim(line1)
          print '(1x, a)', trim(line2)
       end if
10     format (1x, a, 12f10.3)
    end do
    deallocate(wb)

  end subroutine print_model

  !============================================================================

  !! Write model file containing given parameter values
  subroutine write_model(param)

    !subroutine arguments
    double precision :: param(:,:)

    print *, 'write_model not implemented'

  end subroutine write_model

  !============================================================================

  !! Deallocate model storage
  subroutine free_model()

    if (allocated(model_spec)) deallocate(model_spec)
    if (allocated(model_wb)) deallocate(model_wb)
    if (allocated(model_param)) deallocate(model_param)
    if (allocated(model_prior)) deallocate(model_prior)
    if (allocated(model_limits)) deallocate(model_limits)
    if (allocated(model_desc)) deallocate(model_desc)
    if (allocated(model_pos_relto)) deallocate(model_pos_relto)
    call free_clv()

  end subroutine free_model
  
  !============================================================================
  
  !! Read numerical CLV from file and calculate model visibilities
  subroutine read_clv(info, file_name, wavebands)

    !subroutine arguments
    character(len=*), intent(out) :: info !! Error message
    character(len=*), intent(in) :: file_name !! Filename to read
    !model can be wavelength-dependent:
    real, intent(in), optional :: wavebands(:,:)

    !local variables
    integer, parameter :: iunit = 12, iunit2 = 13
    integer iclv, iwl, nline, iline, nwl, iwb, istart, iend
    character(len=256) :: line, ext
    real dummy1, dummy2
    real, allocatable :: wl(:)
    real, allocatable :: inten(:,:)

    !decide what sort of CLV file we have
    ext = trim(file_name(scan(file_name,'.',.true.)+1:len(file_name)))

    if (ext == 'clv') then

       !wavelength-independent CLV; two columns giving r, I(r)

       !count number of data lines, allocate storage
       open (unit=iunit, file=file_name, status='old', action='read', err=91)
       nclv = 0
       do while(.true.)
          read(iunit,'(a)',end=20) line
          !don't count comments or blank lines
          if (line(1:1) /= '#' .and. len_trim(line) > 0) nclv = nclv + 1
       end do
20     close (iunit)
       if (nclv == 0) then
          info = 'file contains no data'
          return
       end if
       allocate(clv_rad(nclv))
       allocate(clv_inten(nclv, 1))
       allocate(clv_mvis(nxsiz+1, 1))

       !read two columns of data
       open (unit=iunit, file=file_name, action='read', err=91)
       iclv = 1
       do while(.true.)
          read(iunit,'(a)',end=30) line
          if (line(1:1) /= '#' .and. len_trim(line) > 0) then
             read(line,*) clv_rad(iclv), clv_inten(iclv, 1)
             iclv = iclv + 1
          end if
       end do
30     close (iunit)

    else

       !wavelength-dependent CLV
       !format as in Kurucz models supplied to Chris Haniff by Bill Tango!
       if (.not. present(wavebands)) then
          info = 'observing waveband(s) unspecified for wavelength-dependent numerical CLV'
          return
       end if

       open (unit=iunit, file=file_name, status='old', action='read', err=91)

       !1st non-comment line gives no. of mu values
       do while(.true.)
          read(iunit,'(a)',end=90) line
          if (line(1:1) /= '#' .and. len_trim(line) > 0) exit
       end do
       read (line,*) nclv
       nclv = nclv + 1 !extra point for intensity beyond mu=0
       allocate(clv_rad(nclv))
       allocate(clv_inten(nclv, size(wavebands, 1)))
       allocate(clv_mvis(nxsiz+1, size(wavebands, 1)))
       read (iunit,*) clv_rad(2:nclv)
       !convert from mu to r
       do iclv = 2, nclv
          clv_rad(iclv) = sqrt(1. - clv_rad(iclv)**2)
       end do
       clv_rad(1) = 1.0001

       !subsequent lines give Temp(K) logg wavelength(nm) I(mu=0) ... I(mu=1)
       !count number of data lines, allocate local storage
       nline = 0 !number of reads to EOF
       nwl = 0
       do while(.true.)
          read(iunit,'(a)',end=40) line
          nline = nline + 1
          !don't count comments or blank lines
          if (line(1:1) /= '#' .and. len_trim(line) > 0) nwl = nwl + 1
       end do
40     do iline = 1, nline
          backspace(iunit)
       end do
       allocate(wl(nwl))
       allocate(inten(nclv, nwl))

       !read CLV data (typically on fine wavelength grid)
       iwl = 0
       do while(.true.)
          read(iunit,'(a)',end=50) line
          if (line(1:1) /= '#' .and. len_trim(line) > 0) then
             iwl = iwl + 1
             inten(1, iwl) = 0.0
             read(line,*) dummy1, dummy2, wl(iwl), inten(2:nclv, iwl)
          end if
       end do
50     close (iunit)

       !integrate to get required (top-hat) bandpasses
       do iwb = 1, size(wavebands, 1)
          istart = locate(wl, wavebands(iwb, 1)-0.5*wavebands(iwb, 2))
          iend = locate(wl, wavebands(iwb, 1)+0.5*wavebands(iwb, 2))
          clv_inten(:, iwb) = 0.
          do iwl = istart, iend
             clv_inten(:, iwb) = clv_inten(:, iwb) + inten(:, iwl)
          end do
          clv_inten(:, iwb) = clv_inten(:, iwb) / (iend-istart+1)
       end do

    end if

    return

    !error trapping 
90  info = 'file contains no data'
    return
91  info = 'cannot open file '//trim(file_name)
    return

  end subroutine read_clv

  !============================================================================
  
  !! Calculate visibilities from numerical CLV
  subroutine calcvis_clv(model_diam)

    !subroutine arguments
    double precision, intent(in) :: model_diam !! Diameter of model star /mas

    !local variables
    integer threshold, ixcen, iycen, ix, iy, i, j, lookup, sign
    !integer plan
    integer*8 :: plan
    integer icalc, ncalc
    real flux, radius, max_rad, frac, delta, factor
    double precision :: map2d(nxsiz+1,nxsiz+1)
    double precision :: map1d(2*nxsiz+1), map1d_rft(2*nxsiz+1)

    !limits
    max_rad = maxval(clv_rad)
    threshold = nint(max_rad*modelsize) + 2
    ncalc = size(clv_inten, 2) !either unity or no. of wavebands

    if (ncalc > 20) then
       call dfftw_plan_r2r_1d(plan, 2*nxsiz+1, map1d, map1d_rft, &
            FFTW_R2HC, FFTW_MEASURE)
    else
       call dfftw_plan_r2r_1d(plan, 2*nxsiz+1, map1d, map1d_rft, &
            FFTW_R2HC, FFTW_ESTIMATE)
    end if

    do icalc = 1, ncalc

       !do a brute force Hankel transform:
       !first, make a 2d image of the top left quadrant of the disk
       map2d = 0.
       flux = 0.
       ixcen = nxsiz + 1
       iycen = nxsiz + 1
       do i = 1, nxsiz + 1
          ix = i - ixcen
          if (abs(ix) <= threshold) then
             do j = 1, nxsiz + 1
                iy = j - iycen
                if (abs(iy) <= threshold) then
                   !radius in units of photospheric radius
                   !(modelsize is photospheric radius in pixels)
                   radius = sqrt(real(ix*ix + iy*iy))/real(modelsize)
                   if (radius <= max_rad) then
                      !use linear interpolation here
                      lookup = locate(clv_rad, radius)
                      if (lookup == 0) then
                         map2d(i,j) = clv_inten(1, icalc)
                      else if (lookup == nclv) then
                         map2d(i,j) = clv_inten(nclv, icalc)
                      else
                         !CAH had this wrong
                         if (clv_rad(lookup+1) == clv_rad(lookup)) then
                            map2d(i,j) = clv_inten(lookup, icalc)
                         else
                            frac = (radius - clv_rad(lookup)) / &
                                 (clv_rad(lookup+1) - clv_rad(lookup))
                            delta = clv_inten(lookup+1, icalc) &
                                 - clv_inten(lookup, icalc)
                            map2d(i,j) = clv_inten(lookup, icalc) + (frac*delta)
                         end if
                      end if
                      flux = flux + map2d(i,j)
                   end if
                end if
             end do
          end if
       end do

       !now project to 1 dimension, and reflect
       do i = 1, nxsiz + 1
          map1d(i) = 0D0
          ix = i - ixcen
          if (abs(ix) <= threshold) then
             do j = 1, nxsiz + 1
                map1d(i) = map1d(i) + map2d(i,j)
             end do
             map1d(i) = map1d(i)/flux
             map1d(2*nxsiz + 2 - i) = map1d(i)
          end if
          write(30, *) i, map1d(i)
       end do

       !and take the fourier transform of this real (symmetric) array
       call dfftw_execute_r2r(plan, map1d, map1d_rft)

       !normalise to zero freq. value and fudge the sign
       sign = 1
       do i = 1, nxsiz + 1
          clv_mvis(i, icalc) = sign*map1d_rft(i)/map1d_rft(1)
          sign = -sign
       end do

    end do

    !scale baselines (which are in Mega-lambda) so that the model
    !corresponds to a star with photospheric diameter model_diam mas
    clv_mdiam = model_diam
    factor = 1.0e-6*modelsize/(clv_mdiam*mas2rad)/nxsiz
    do i = 1, nxsiz + 1
       clv_mbase(i) = real(i-1)*factor
       write(20, *) clv_mbase(i), clv_mvis(i, 1)
    end do

    call dfftw_destroy_plan(plan)

  end subroutine calcvis_clv

  !============================================================================

  !! Deallocate numerical clv storage
  subroutine free_clv()

    if (allocated(clv_rad)) deallocate(clv_rad)
    if (allocated(clv_inten)) deallocate(clv_inten)
    if (allocated(clv_mvis)) deallocate(clv_mvis)

  end subroutine free_clv

  !============================================================================

  !! Count number of symbols in line
  function countsym(line)

    integer :: countsym

    !function arguments
    character(len=*), intent(in) :: line !! Line to parse

    !local variables
    integer istart, iend, lline
    character (len=2) :: sep='  ' !valid separators

    sep(2:2) = achar(9) !tab
    lline = len(line)
    countsym = 0
    iend = 0

    do while (iend < lline)
       istart = verify(line((iend+1):lline),sep)+iend !start of next symbol
       if (istart == iend) exit
       countsym = countsym + 1
       iend = scan(line((istart+1):lline),sep)+istart-1 !end of symbol
       if (iend < istart) iend = lline !scan returned 0 (no match)
    end do

  end function countsym

  !============================================================================

end module Model
