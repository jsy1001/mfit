module Model

!callable subroutines contained
!
!read_model
!print_model
!write_model
!free_model

implicit none

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

contains

!==============================================================================

subroutine read_model(info, file_name)
    
  !Reads model files
  !(Re-)Allocates and assigns to module variables
  
  !subroutine arguments
  character(len=128), intent(out) :: info
  character(len=*), intent(in) :: file_name
  
  !local variables
  integer, parameter :: max_lines = 1000
  character(len=32) :: dummy, source1, source2
  integer :: i, j, k, lines, comps, order
  
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

1 lines = i-1
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
  do i = 1, lines
     read (11, *, err=94) dummy
     if (dummy == 'source') then
        read (12, *, err=94) (dummy, source1, source2)
        model_name = trim(source1) // ' ' // trim(source2)
     else
        read (12, *, err=94)
     end if
  end do
  close (11)
  close (12)
  
  !read spec and param info for components
  open (unit=11, action='read', file=file_name)
  j = 0 !component counter
  do i = 1, lines
     read (11,*,err=94,end=2) dummy
     if (dummy == 'component') then
        j = j+1
        
        !read component name
        read (11,*,err=95,end=95) (dummy, model_spec(j,1))
        if (dummy /= 'name') goto 96
        
        !read shape type
        read (11,*,err=95,end=95) (dummy, model_spec(j,2))
        if (dummy /= 'shape_type') goto 96
        
        !read LD type (set as uniform for point case)
        select case (trim(model_spec(j,2)))
           case ('point')
              read (11,*,err=95) dummy
              model_spec(j,3) = 'uniform'
           case default
              read (11,*,err=95,end=95) (dummy, model_spec(j,3))
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
              read (11,*,err=95,end=95) (dummy, order)
           case default
              info = 'invalid limb darkening type'
              close (11)
              return
           end select
        if (dummy /= 'ld_order') goto 96
        model_param(j,1) = dble(order)
        
        !read position r and theta
        read (11,*,err=95,end=95) (dummy, model_param(j,2:3))
        if (dummy /= 'position') goto 96
        
        !read position prior, must be non-negative
        read (11,*,err=95,end=95) (dummy, model_prior(j,2:3))
        if (dummy /= 'position_prior') goto 96
        
        !read flux
        read (11,*,err=95,end=95) (dummy, model_param(j,4))
        if (dummy /= 'flux') goto 96
        
        !read flux prior, must be non-negative
        read (11,*,err=95,end=95) (dummy, model_prior(j,4))
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
              read (11,*,err=95,end=95) (dummy, model_param(j,5))
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) (dummy, model_prior(j,5))
           case ('ellipse')
              read (11,*,err=95,end=95) (dummy, model_param(j,5:7))
              if (dummy /= 'shape_param') goto 96
              read (11,*,err=95,end=95) (dummy, model_prior(j,5:7))
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
              read (11,*,err=95,end=95) (dummy, model_param(j,8:7+order))
              if (dummy /= 'ld_param') goto 96
              read (11,*,err=95,end=95) (dummy, model_prior(j,8:7+order))
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
91 info = 'cannot open file'
  return
92 info = 'cannot read from file'
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
     print *,'LD type  : ',trim(model_spec(i,3)), &
          ' of order',int(model_param(i,1))
     print *,'shape    : ',trim(model_spec(i,2))
     print *,'position :',real(model_param(i,2:3))
     print *,'    prior:',real(model_prior(i,2:3))
     print *,'flux     :',real(model_param(i,4))
     print *,'    prior:',real(model_prior(i,4))
     print *,'shape par:',real(model_param(i,5:7))
     print *,'    prior:',real(model_prior(i,5:7))
     print *,'LD params:',real(model_param(i,8:17))
     print *,'    prior:',real(model_prior(i,8:17))
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

end subroutine free_model

!==============================================================================

end module Model
