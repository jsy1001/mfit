!$Id: model.f90,v 1.6 2005/01/06 18:45:12 jsy1001 Exp $

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

!initial model parameters (both variable and non-variable)
!17 items per cpt: 1 LD order, 2-7 r, theta, brightness, major axis, phi,
!epsilon, 8-17 LD coeffs 1-10
double precision, dimension(:,:), allocatable :: model_param

!Gaussian prior widths for parameters; prior centres are in model_param
double precision, dimension(:,:), allocatable :: model_prior

!limits on legal values for parameters (both variable and non-variable)
double precision, dimension(17,2) :: model_limits

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
real, dimension(:, :), allocatable :: clv_mvis, clv_wb
real clv_mdiam

contains

!==============================================================================

subroutine read_model(info, file_name, wavebands)
    
  !Reads model files
  !(Re-)Allocates and assigns to module variables
  
  !subroutine arguments
  character(len=128), intent(out) :: info
  character(len=*), intent(in) :: file_name
  double precision, dimension(:, :), intent(in) :: wavebands
  
  !local variables
  integer, parameter :: max_lines = 1000
  character(len=32) :: dummy, source1, source2
  character(len=128) :: clv_filename
  character(len=256) :: line
  integer :: i, j, k, nlines, comps, order, pos
  
  !check for zero length filename
  if (file_name == '') then
     info = 'blank filename'
     return
  end if
  
  !read and count number of components and number of lines
  open (unit=11, err=91, status='old', action='read', file=file_name)
  comps = 0
  do i = 1, max_lines+1
     read (11, *, err=92, end=1) dummy
     if (dummy == 'component') comps = comps+1
  end do
  info = 'file exceeds maximum permitted length'
  close (11)
  return

1 nlines = i-1
  close (11)
  
  !check if legal number of components
  if ((comps<1).or.(comps>10)) then
     info = 'illegal number of components present'
     return
  end if

  !allocate model data arrays
  call free_model()
  allocate(model_spec(comps,3))
  allocate(model_param(comps,17))
  allocate(model_prior(comps,17))
  model_param = 0D0
  model_prior = 0D0
  clvcomp = 0 !no model component yet has numerical CLV

  !Limits array stores the lower and upper acceptable limits
  !on the parameter values, it is passed into the minimising routines
  !later to ensure parameter values stay legal. Model values not within
  !these limits trigger error on reading the model
  model_limits(1:7,1) = dble((/0,0,-720,-100,0,-720,0/))
  model_limits(8:17,1) = -100D0
  model_limits(1:7,2) = dble((/10,100,720,100,1000,720,10/))
  model_limits(8:17,2) = 100D0

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
                 call read_clv(info, clv_filename, wavebands)
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
        
        !read position r and theta
        read (11,*,err=95,end=95) dummy, model_param(j,2:3)
        if (dummy /= 'position') goto 96
        
        !read position prior, must be non-negative
        read (11,*,err=95,end=95) dummy, model_prior(j,2:3)
        if (dummy /= 'position_prior') goto 96
        
        !read flux
        read (11,*,err=95,end=95) dummy, model_param(j,4)
        if (dummy /= 'flux') goto 96
        
        !read flux prior, must be non-negative
        read (11,*,err=95,end=95) dummy, model_prior(j,4)
        if (dummy /= 'flux_prior') goto 96  
        
        !Read shape parameters a, phi, epsilon depending on
        !shape type. Also read priors
        !Shape a and epsilon must be non-negative
        !All priors must be non-negative
        model_param(j,5) = 0D0
        model_param(j,6) = 0D0
        model_param(j,7) = 1D0
        select case (trim(model_spec(j,2)))
           case ('point')
              read (11,*,err=95,end=95) dummy
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) dummy
           case ('disc')
              read (11,*,err=95,end=95) dummy, model_param(j,5)
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) dummy, model_prior(j,5)
           case ('ellipse')
              read (11,*,err=95,end=95) dummy, model_param(j,5:7)
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) dummy, model_prior(j,5:7)
           case default
              info = 'invalid shape type'
              close (11)
              return
        end select
        if (dummy /= 'shape_param_prior') goto 96

        !read LD parameters
        select case (order)
           case (0)
              read (11,*,err=95,end=95) dummy
              if (dummy /= 'ld_param') goto 96
              read (11,*,err=95,end=95) dummy
           case default
              read (11,*,err=95,end=95) dummy, model_param(j,8:7+order)
              if (dummy /= 'ld_param') goto 96
              read (11,*,err=95,end=95) dummy, model_prior(j,8:7+order)
        end select

        !check to ensure everything inside limits
        do k=1, size(model_limits,1)
           if (model_param(j,k)<model_limits(k,1)) goto 100
           if (model_param(j,k)>model_limits(k,2)) goto 100
        end do
        !ensure alpha>0 for hestroffer case
        select case(trim(model_spec(j,3)))
           case('hestroffer')
              if (.not.(model_param(j,8) > 0D0)) goto 100
        end select
        !if numerical CLV, calculate visibilities for initial guess diameter
        if (model_spec(j,3)(1:1) == '<') &
             call calcvis_clv(model_param(j,5))

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

  integer :: i

  do i = 1, size(model_param,1)
     print *,' '
     print *,'component:',i
     print *,'name     : ',trim(model_spec(i,1))
     if (model_spec(i,3)(1:1) == '<') then
        print *,'LD type  : ',trim(model_spec(i,3)), &
             ' :',nclv,'points'
     else
        print *,'LD type  : ',trim(model_spec(i,3)), &
             ' of order',int(model_param(i,1))
     end if
     print *,'shape    : ',trim(model_spec(i,2))
     print 10,'position :',real(model_param(i,2:3))
     print 10,'    prior:',real(model_prior(i,2:3))
     print 10,'flux     :',real(model_param(i,4))
     print 10,'    prior:',real(model_prior(i,4))
     print 10,'shape par:',real(model_param(i,5:7))
     print 10,'    prior:',real(model_prior(i,5:7))
     print 20,'LD params:',real(model_param(i,8:17))
     print 20,'    prior:',real(model_prior(i,8:17))
10   format (1x, a, 3f10.3)
20   format (1x, a, 10f7.3)
  end do

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
  if (allocated(model_param)) deallocate(model_param)
  if (allocated(model_prior)) deallocate(model_prior)
  call free_clv()

end subroutine free_model

!==============================================================================

subroutine read_clv(info, file_name, wavebands)

  !Read numerical CLV from file and calculate model visibilities

  !subroutine arguments
  character(len=128), intent(out) :: info
  character(len=*), intent(in) :: file_name
  double precision, dimension(:, :), intent(in) :: wavebands !model can be wavelength-dependent

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
     !clv_wb not required

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
     allocate(clv_wb(size(wavebands, 1), 2))
     clv_wb = wavebands
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
     do iwb = 1, size(clv_wb, 1)
        start = locate(wl, clv_wb(iwb, 1)-0.5*clv_wb(iwb, 2))
        end = locate(wl, clv_wb(iwb, 1)+0.5*clv_wb(iwb, 2))
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
  if (allocated(clv_wb)) deallocate(clv_wb)

end subroutine free_clv

!==============================================================================

end module Model
