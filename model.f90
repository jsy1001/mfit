!$Id: model.f90,v 1.8 2005/05/24 12:53:06 jsy1001 Exp $

module Model

!callable subroutines contained
!
!read_model
!print_model
!write_model
!free_model
!read_clv
!calcvis_clv
!free_clv

use Maths

implicit none

include 'fftw_f77.i'

!module variables contained:

!specifies nature of model cpts -
!holds 3 items per component: 1 name, 2 shape and 3 LD type
character(len=128), dimension(:,:), allocatable :: model_spec

!model_wldep(i) is unity if position, flux, shape_param, ld_param respectively
!is fn(waveband), else zero
integer, dimension(4) :: model_wldep
!model has nwave values for each wavelength-dependent parameter
integer nwave

!wavebands for wavelength-dependent model parameters (if supplied to read_model)
!and CLV
real, dimension(:, :), allocatable :: model_wb

!initial model parameters (both variable and non-variable)
!for each cpt (1st axis): 2nd axis is
! LD order
! {r, theta(t0), t0, d(theta)/dt} possibly repeated for each waveband
! {brightness}                    possibly repeated for each waveband
! {major axis, phi, epsilon}      possibly repeated for each waveband
! {LD coeffs 1-max_order}         possibly repeated for each waveband
double precision, dimension(:,:), allocatable :: model_param
integer max_order !at least one component has this many LD coeffs

!Gaussian prior widths for parameters; prior centres are in model_param
double precision, dimension(:,:), allocatable :: model_prior

!limits on legal values for parameters (both variable and non-variable)
double precision, dimension(:,:), allocatable :: model_limits

!descriptions of model parameters
character(len=55), dimension(:, :), allocatable :: model_desc

!name from model file (2 words max)
character(len=128) :: model_name

!centrosymmetric model?
logical symm

!numerical CLV data
real, dimension(:), allocatable :: clv_rad
real, dimension(:, :), allocatable :: clv_inten
integer nclv, clvcomp
!model visibilities calculated from numerical CLV
integer, parameter :: nxsiz = 4096, modelsize = 120
real, dimension(nxsiz+1) :: clv_mbase
real, dimension(:, :), allocatable :: clv_mvis
real clv_mdiam

contains

!==============================================================================

subroutine read_model(info, file_name, wavebands)
    
  !Reads model files
  !(Re-)Allocates and assigns to module variables
  
  !subroutine arguments
  character(len=128), intent(out) :: info
  character(len=*), intent(in) :: file_name
  double precision, dimension(:, :), intent(in), optional :: wavebands
  
  !local variables
  integer, parameter :: max_lines = 1000
  character(len=2), dimension(10) :: numbers
  character(len=32) :: dummy, source1, source2, cpt
  character(len=32), dimension(:), allocatable :: wb
  character(len=128) :: clv_filename, shape_type, ld_type
  character(len=256) :: line
  integer :: i, j, k, iwave
  integer :: nlines, comps, order, pos, npar, ipar, nread, nn
  
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
  comps = 0
  model_wldep = 0
  nwave = -1
  max_order = 0
  do i = 1, max_lines+1
     read (11, '(a)', err=92, end=18) line
     if (len_trim(line) == 0) cycle
     read (line, *) dummy !read keyword & qualifiers
     if (dummy == 'component') comps = comps+1
     if (dummy(1:16) == 'position#fofwave' &
          .or. dummy == 'position#rotate#fofwave') then
        model_wldep(1) = 1
        ! qualifiers could be "#fofwave#rotate" or "#rotate#fofwave"
        if (dummy(9:15) == '#rotate' .or. dummy(17:23) == '#rotate') then
           !{r, theta(t0), t0, d(theta)/dt} poss. repeated for each waveband
           nn = (countsym(line)-1)/4
        else
           !file contains {r, theta} poss. repeated for each waveband
           nn = (countsym(line)-1)/2
        end if
        if (nwave == -1) then
           nwave = nn
        else if (nn /= nwave) then
           write(info, *) 'line ', i, ': expecting values for ', nwave, &
                 ' wavebands, got ', nn
           close(11)
           return
        end if
     else if (dummy == 'flux#fofwave') then
        model_wldep(2) = 1
        nn = countsym(line)-1
        if (nwave == -1) then
           nwave = nn
        else if (nn /= nwave) then
           write(info, *) 'line ', i, ': expecting values for ', nwave, &
                 ' wavebands, got ', nn
           close(11)
           return
        end if
     else if (dummy == 'shape_type') then
        read (line, *, end=95) dummy, shape_type
     else if (dummy == 'shape_param#fofwave') then
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
           write(info, *) 'line ', i, ': expecting values for ', nwave, &
                 ' wavebands, got ', nn
           close(11)
           return
        end if
     else if (dummy == 'ld_type') then
        read (line, *, end=95) dummy, ld_type
        !Read LD order (preset for some LD type cases)
        !Non-integer will cause error in read statement
        select case (trim(ld_type))
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
              read (11,*,err=95,end=95) dummy, order
           case default
              if (ld_type(1:1) == '<') then
                 read (11,*,err=95,end=95) dummy
                 order = 0
              else
                 !not <filename>
                 info = 'invalid limb darkening type'
                 close (11)
                 return
              end if
        end select
        if (order > max_order) max_order = order
     else if (dummy == 'ld_param#fofwave') then
        model_wldep(4) = 1
        nn = (countsym(line)-1)/order
        if (nwave == -1) then
           nwave = nn
        else if (nn /= nwave) then
           write(info, *) 'line ', i, ': expecting values for ', nwave, &
                 ' wavebands, got ', nn
           close(11)
           return
        end if
     end if
  end do
  info = 'file exceeds maximum permitted length'
  close (11)
  return

18 nlines = i-1
  close (11)
  
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
  allocate(model_limits(npar,2))
  allocate(model_desc(comps,npar))
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
  model_limits(1,1:2) = dble((/0,10/))
  numbers = (/' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10'/)
  do j = 1, comps
     cpt = 'component '// trim(adjustl(numbers(j))) // ', '
     model_desc(j, 1) = trim(cpt) // ' LD order'
  end do

  !read source name
  open (unit=11, action='read', file=file_name)
  open (unit=12, action='read', file=file_name)
  do i = 1, nlines
     read (11, *, err=94) dummy
     if (dummy == 'source') then
        read (12, *, err=94) dummy, source1, source2
        model_name = trim(source1) // ' ' // trim(source2)
        exit
     else
        read (12, *, err=94) dummy
     end if
  end do
  close (11)
  close (12)
  
  !read spec and param info for components
  open (unit=11, action='read', file=file_name)
  j = 0 !component counter
  do i = 1, nlines
     read (11,*,err=94,end=2) dummy
     if (dummy == 'component') then
        j = j+1
        cpt = 'component '// trim(adjustl(numbers(j))) // ', '
        
        !read component name
        read (11,*,err=95,end=95) dummy, model_spec(j,1)
        if (dummy /= 'name') goto 96
        
        !read shape type
        read (11,*,err=95,end=95) dummy, model_spec(j,2)
        if (dummy /= 'shape_type') goto 96
        
        !read LD type (set as uniform for point case)
        select case (trim(model_spec(j,2)))
           case ('point')
              read (11,*,err=95) dummy
              model_spec(j,3) = 'uniform'
           case default
              read (11,'(a)',err=95,end=95) line
              read (line,*,err=95,end=95) dummy
              pos = scan(trim(line), ' '//achar(9), back=.true.) + 1
              model_spec(j,3) = line(pos:)
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
              read (11,*,err=95,end=95) dummy, order
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
                 read (11,*,err=95,end=95) dummy
                 order = 0
              else
                 !not <filename>
                 info = 'invalid limb darkening type'
                 close (11)
                 return
              end if
           end select
        if (dummy /= 'ld_order') goto 96
        model_param(j,1) = dble(order)
        
        !read position r and theta etc.
        read (11,*,err=95,end=95) dummy
        if (dummy(:8) /= 'position') goto 96
        do k = 1, 1+model_wldep(1)*(nwave-1)
           if (j == 1) then
              model_limits(2+(k-1)*4,1:2) = dble((/0, 500/))
              model_limits(3+(k-1)*4,1:2) = dble((/-720, 720/))
              model_limits(4+(k-1)*4,1:2) = dble((/10000, 100000/))
              model_limits(5+(k-1)*4,1:2) = dble((/-1000, 1000/))
           end if
           !defaults for t0 and d(theta)/dt
           model_param(j, 4+(k-1)*4) = 10000D0
           model_param(j, 5+(k-1)*4) = 0D0
           !descriptions
           if (model_wldep(1) == 1) then
              model_desc(j, 2+(k-1)*4) = trim(cpt) // trim(wb(k)) // ' position radius (mas)'
              model_desc(j, 3+(k-1)*4) = trim(cpt) // trim(wb(k)) // ' position angle(t0) (deg)'
              model_desc(j, 4+(k-1)*4) = trim(cpt) // trim(wb(k)) // ' t0 (MJD)'
              model_desc(j, 5+(k-1)*4) = trim(cpt) // trim(wb(k)) // ' d(PA)/dt (deg/day)'
           else
              model_desc(j, 2+(k-1)*4) = trim(cpt) // ' position radius (mas)'
              model_desc(j, 3+(k-1)*4) = trim(cpt) // ' position angle(t0) (deg)'
              model_desc(j, 4+(k-1)*4) = trim(cpt) // ' t0 (MJD)'
              model_desc(j, 5+(k-1)*4) = trim(cpt) // ' d(PA)/dt (deg/day)'
           end if
        end do
        ! qualifiers could be "#fofwave#rotate" or "#rotate#fofwave"
        if (dummy(9:15) == '#rotate' .or. dummy(17:23) == '#rotate') then
           !{r, theta(t0), t0, d(theta)/dt} poss. repeated for each waveband
           backspace(11)
           nread = 4*(1 + model_wldep(1)*(nwave-1))
           read (11,*,err=95,end=95) dummy, model_param(j, 2:1+nread)
           read (11,*,err=95,end=95) dummy, model_prior(j, 2:1+nread)
           if (dummy /= 'position_prior') goto 96
        else
           !file contains {r, theta} poss. repeated for each waveband
           backspace(11)
           read (11,*,err=95,end=95) dummy, &
                (model_param(j, 2+(k-1)*4:3+(k-1)*4), &
                k=1,(1+model_wldep(1)*(nwave-1)))
           read (11,*,err=95,end=95) dummy, &
                (model_prior(j, 2+(k-1)*4:3+(k-1)*4), &
                k=1,(1+model_wldep(1)*(nwave-1)))
           if (dummy /= 'position_prior') goto 96
        end if
        
        !read flux
        ipar = 6+4*model_wldep(1)*(nwave-1)
        read (11,*,err=95,end=95) dummy, &
             model_param(j, ipar:ipar+model_wldep(2)*(nwave-1))
        if (dummy(:4) /= 'flux') goto 96
        do k = 1, 1+model_wldep(2)*(nwave-1)
           if (j == 1) &
                model_limits(ipar+(k-1),1:2) = dble((/-100, 100/))
           if (model_wldep(2) == 1) then
              model_desc(j, ipar+(k-1)) = trim(cpt) // trim(wb(k)) // ' flux (arb units)'
           else
              model_desc(j, ipar+(k-1)) = trim(cpt) // ' flux (arb units)'
           end if
        end do
        
        !read flux prior, must be non-negative
        read (11,*,err=95,end=95) dummy, &
             model_prior(j, ipar:ipar+model_wldep(2)*(nwave-1))
        if (dummy /= 'flux_prior') goto 96  
        
        !Read shape parameters a, phi, epsilon depending on
        !shape type. Also read priors
        !Shape a and epsilon must be non-negative
        !All priors must be non-negative
        ipar = 7 + (4*model_wldep(1) + model_wldep(2))*(nwave-1)
        do k=1, 1+model_wldep(3)*(nwave-1)
           model_param(j,ipar+(k-1)*3) = 0D0
           model_param(j,ipar+(k-1)*3+1) = 0D0
           model_param(j,ipar+(k-1)*3+2) = 1D0
           if (j == 1) then
              model_limits(ipar+(k-1)*3,1:2) = dble((/0, 1000/))
              model_limits(ipar+(k-1)*3+1,1:2) = dble((/-720, 720/))
              model_limits(ipar+(k-1)*3+2,1:2) = dble((/0, 10/))
           end if
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
        select case (trim(model_spec(j,2)))
           case ('point')
              read (11,*,err=95,end=95) dummy
              if (dummy(:11) /= 'shape_param') goto 96
              read (11,*,err=95,end=95) dummy
           case ('disc')
              read (11,*,err=95,end=95) dummy, (model_param(j,ipar+(k-1)*3), &
                   k=1,(1+model_wldep(3)*(nwave-1)))
              if (dummy(:11) /= 'shape_param') goto 96
              read (11,*,err=95,end=95) dummy, (model_prior(j,ipar+(k-1)*3), &
                   k=1,(1+model_wldep(3)*(nwave-1)))
           case ('ellipse')
              nread = 3*(1 + model_wldep(3)*(nwave-1))
              read (11,*,err=95,end=95) dummy, model_param(j,ipar:ipar+nread-1)
              if (dummy(:11) /= 'shape_param') goto 96
              read (11,*,err=95,end=95) dummy, model_prior(j,ipar:ipar+nread-1)
           case default
              info = 'invalid shape type'
              close (11)
              return
        end select
        if (dummy /= 'shape_param_prior') goto 96

        !read LD parameters
        ipar = 10 + (4*model_wldep(1) + model_wldep(2) &
             + 3*model_wldep(3))*(nwave-1)
        select case (order)
           case (0)
              read (11,*,err=95,end=95) dummy
              if (dummy(:8) /= 'ld_param') goto 96
              read (11,*,err=95,end=95) dummy
           case default
              nread = order*(1 + model_wldep(4)*(nwave-1))
              read (11,*,err=95,end=95) dummy, &
                   model_param(j, ipar:ipar+nread-1)
              if (dummy(:8) /= 'ld_param') goto 96
              read (11,*,err=95,end=95) dummy, &
                   model_prior(j, ipar:ipar+nread-1)
        end select
        if (j == 1) then
           model_limits(ipar: &
                ipar+max_order*(1 + model_wldep(4)*(nwave-1))-1,1) = -100D0
           model_limits(ipar: &
                ipar+max_order*(1 + model_wldep(4)*(nwave-1))-1,2) = 100D0
        end if
        !LD param descriptions
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

        !check to ensure everything inside limits
        do k=1, size(model_limits,1)
           if (model_param(j,k) < model_limits(k,1)) goto 100
           if (model_param(j,k) > model_limits(k,2)) goto 100
        end do
        !ensure alpha>0 for hestroffer case
        if (trim(model_spec(j,3)) == 'hestroffer') then
           if (model_wldep(4) == 1) then
              do iwave = 1, nwave
                 if (.not.(model_param(j, ipar+(iwave-1)*max_order) > 0D0)) &
                      goto 100
              end do
           else
              if (.not.(model_param(j, ipar) > 0D0)) goto 100
           end if
        end if
        !if numerical CLV, calculate visibilities for initial guess diameter
        !wavelength-dep CLV with wavelength-dep diameter doesn't make sense
        !- this is trapped in minimiser()
        if (model_spec(j,3)(1:1) == '<') &
             call calcvis_clv(model_param(j, &
             7+(4*model_wldep(1)+model_wldep(2))*(nwave-1)))

        if (dummy /= 'ld_param_prior') goto 96
        
     end if
  end do
  
2 continue 
  close (11)

  !is this a centrosymmetric model?
  symm = .true.
  do i = 1, size(model_param, 1)
     !need all epsilon fixed at unity
     if (.not.((model_param(i,7)==1D0).and.(model_prior(i,7)==0D0))) &
          symm = .false.
     !and all r fixed at zero
     if (.not.((model_param(i,2)==0D0).and.(model_prior(i,2)==0D0))) &
          symm = .false.
  end do

  deallocate(wb) !free local storage
  return
  
  !error trapping 
91 info = 'cannot open file '//trim(file_name)
  return
92 info = 'cannot read from file '//trim(file_name)
  close(11)
  return
94 info = 'unknown read problem - possible invalid file format'
  close (11)
  return
95 info = 'unknown read problem whilst reading component data'
  close (11)
  return
96 info = 'invalid keyword or keyword order'
  close (11)
  return
100 info = 'illegal parameter value detected'
  close (11)
  return
  
end subroutine read_model

!==============================================================================

subroutine print_model()

  !Print model details

  integer :: i, ipar, jpar, iwave, order
  character(len=32), dimension(:), allocatable :: wb
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
        print '(1x, a)', trim(line1)
        print '(1x, a)', trim(line2)
     end if
10   format (1x, a, 12f10.3)
  end do
  deallocate(wb)

end subroutine print_model

!==============================================================================

subroutine write_model(param)

  !Write model file containing given parameter values

  !subroutine arguments
  double precision, dimension(:,:) :: param

  print *, 'write_model not implemented'

end subroutine write_model

!==============================================================================

subroutine free_model()

  !Deallocate model storage
  if (allocated(model_spec)) deallocate(model_spec)
  if (allocated(model_wb)) deallocate(model_wb)
  if (allocated(model_param)) deallocate(model_param)
  if (allocated(model_prior)) deallocate(model_prior)
  if (allocated(model_limits)) deallocate(model_limits)
  if (allocated(model_desc)) deallocate(model_desc)
  call free_clv()

end subroutine free_model

!==============================================================================

subroutine read_clv(info, file_name, wavebands)

  !Read numerical CLV from file and calculate model visibilities

  !subroutine arguments
  character(len=128), intent(out) :: info
  character(len=*), intent(in) :: file_name
  !model can be wavelength-dependent:
  real, dimension(:, :), intent(in), optional :: wavebands

  !local variables
  integer, parameter :: iunit = 12, iunit2 = 13
  integer iclv, iwl, nwl, iwb, nwb, start, end
  character(len=256) :: line, ext
  real dummy1, dummy2
  real, dimension(:), allocatable :: wl
  real, dimension(:, :), allocatable :: inten

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
20   close (iunit)
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
30   close (iunit)

  else

     !wavelength-dependent CLV
     !format as in Kurucz models supplied to Chris Haniff by Bill Tango!
     if (.not. present(wavebands)) then
        info = 'observing waveband(s) unspecified for wavelength-dependent numerical CLV'
        return
     end if

     open (unit=iunit, file=file_name, status='old', action='read', err=91)
     open (unit=iunit2, file=file_name, action='read')

     !1st non-comment line gives no. of mu values
     do while(.true.)
        read(iunit,'(a)',end=90) line
        read(iunit2,'(a)') line
        if (line(1:1) /= '#' .and. len_trim(line) > 0) exit
     end do
     read (line,*) nclv
     nclv = nclv + 1 !extra point for intensity beyond mu=0
     allocate(clv_rad(nclv))
     allocate(clv_inten(nclv, size(wavebands, 1)))
     allocate(clv_mvis(nxsiz+1, size(wavebands, 1)))
     read (iunit,'(a)') line
     read (iunit2,*) clv_rad(2:nclv)
     !convert from mu to r
     do iclv = 2, nclv
        clv_rad(iclv) = sqrt(1. - clv_rad(iclv)**2)
     end do
     clv_rad(1) = 1.0001

     !subsequent lines give Temp(K) logg wavelength(nm) I(mu=0) ... I(mu=1)
     !count number of data lines, allocate local storage
     nwl = 0
     do while(.true.)
        read(iunit,'(a)',end=40) line
        !don't count comments or blank lines
        if (line(1:1) /= '#' .and. len_trim(line) > 0) nwl = nwl + 1
     end do
40   close (iunit)
     allocate(wl(nwl))
     allocate(inten(nclv, nwl))

     !read CLV data (typically on fine wavelength grid)
     iwl = 0
     do while(.true.)
        read(iunit2,'(a)',end=50) line
        if (line(1:1) /= '#' .and. len_trim(line) > 0) then
           iwl = iwl + 1
           inten(1, iwl) = 0.0
           read(line,*) dummy1, dummy2, wl(iwl), inten(2:nclv, iwl)
        end if
     end do
50   close (iunit2)

     !integrate to get required (top-hat) bandpasses
     do iwb = 1, size(wavebands, 1)
        start = locate(wl, wavebands(iwb, 1)-0.5*wavebands(iwb, 2))
        end = locate(wl, wavebands(iwb, 1)+0.5*wavebands(iwb, 2))
        clv_inten(:, iwb) = 0.
        do iwl = start, end
           clv_inten(:, iwb) = clv_inten(:, iwb) + inten(:, iwl)
        end do
        clv_inten(:, iwb) = clv_inten(:, iwb) / (end-start+1)
     end do

  end if

  return

  !error trapping 
90 info = 'file contains no data'
  return
91 info = 'cannot open file '//trim(file_name)
  return

end subroutine read_clv

!==============================================================================

subroutine calcvis_clv(model_diam)

  !subroutine arguments
  double precision, intent(in) :: model_diam

  !local variables
  integer threshold, ixcen, iycen, ix, iy, i, j, lookup, sign, plan
  integer icalc, ncalc
  real flux, radius, max_rad, frac, delta, factor
  double precision, dimension(:,:), allocatable :: map2d
  double precision, dimension(:), allocatable :: map1d, map1d_rft

  !allocate storage
  allocate(map2d(nxsiz+1, nxsiz+1))
  allocate(map1d(2*nxsiz+1))
  allocate(map1d_rft(2*nxsiz+1))

  !limits
  max_rad = maxval(clv_rad)
  threshold = nint(max_rad*modelsize) + 2
  ncalc = size(clv_inten, 2) !either unity or no. of wavebands

  if (ncalc > 20) then
     call rfftw_f77_create_plan(plan, 2*nxsiz+1, FFTW_REAL_TO_COMPLEX, &
          FFTW_MEASURE)
  else
     call rfftw_f77_create_plan(plan, 2*nxsiz+1, FFTW_REAL_TO_COMPLEX, &
          FFTW_ESTIMATE)
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
     call rfftw_f77_one(plan, map1d, map1d_rft)

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

  call rfftw_f77_destroy_plan(plan)
  !free storage
  deallocate(map2d)
  deallocate(map1d)
  deallocate(map1d_rft)

end subroutine calcvis_clv

!==============================================================================

subroutine free_clv()

  !Deallocate numerical clv storage
  if (allocated(clv_rad)) deallocate(clv_rad)
  if (allocated(clv_inten)) deallocate(clv_inten)
  if (allocated(clv_mvis)) deallocate(clv_mvis)

end subroutine free_clv

!==============================================================================
function countsym(line)

  ! Count number of symbols in line

  !function arguments
  character(len=*), intent(in) :: line
  integer countsym

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

!==============================================================================

end module Model
