! Copyright (C) 2003-2018 John Young, Matthew Worsley
!
! This file is part of mfit.
!
! Mfit is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/ .

module Wrap

  implicit none

  public

  type allparam

     logical :: done_init = .false.
     double precision, pointer :: param(:,:)
     double precision, pointer :: svar(:)
     double precision, pointer :: limits(:,:,:)
     integer :: nvar
     integer, pointer :: var_pos(:,:)
     double precision, pointer :: var_scale(:)
     double precision, pointer :: var_offset(:)

  end type allparam

contains

  !============================================================================

  !! Initialise allparam instance
  subroutine allparam_init(this, param, limits, nvar, var_pos, var_scale, var_offset)

    !subroutine arguments
    type(allparam), intent(out) :: this
    double precision, intent(in) :: param(:,:)
    double precision, intent(in) :: limits(:,:,:)
    integer, intent(in) :: nvar
    integer, intent(in) :: var_pos(nvar,2)
    double precision, intent(in), optional :: var_scale(nvar)
    double precision, intent(in), optional :: var_offset(nvar)

    this%done_init = .true.
    allocate(this%param(size(param,1), size(param,2)))
    this%param = param
    allocate(this%svar(nvar))
    this%svar = 0D0  !incorrect until first call to allparam_setvar()
    allocate(this%limits(size(limits,1), size(limits,2), 2))
    this%limits = limits
    this%nvar = nvar
    allocate(this%var_pos(nvar,2))
    this%var_pos = var_pos
    allocate(this%var_scale(nvar))
    if (present(var_scale)) then
       this%var_scale = var_scale
    else
       this%var_scale = 1D0
    end if
    allocate(this%var_offset(nvar))
    if (present(var_offset)) then
       this%var_offset = var_offset
    else
       this%var_offset = 0D0
    end if

  end subroutine allparam_init

  !============================================================================

  !! Copy from one instance to another, perhaps fixing subset of the variables
  subroutine allparam_copy(src, dest, fix_var)

    !subroutine arguments
    type(allparam), intent(in) :: src
    type(allparam), intent(out) :: dest
    logical, intent(in), optional :: fix_var(:)

    !local variables
    integer isrc, idest

    if (present(fix_var)) then
       if (size(fix_var,1) < src%nvar) &
            stop 'invalid fix_var dimension in allparam_copy'
    end if
    if (.not. src%done_init) &
         stop 'src not initialised in allparam_copy'

    !if necessary, deallocate dest storage
    if (dest%done_init) call allparam_free(dest)

    dest%done_init = .true.

    !copy src%param -> dest%param
    allocate(dest%param(size(src%param,1), size(src%param,2)))
    dest%param = src%param

    !copy src%svar-> dest%svar
    allocate(dest%svar(size(src%svar,1)))
    dest%svar = src%svar

    !copy src%limits -> dest%limits
    allocate(dest%limits(size(src%limits,1), size(src%limits,2), 2))
    dest%limits = src%limits

    !assign dest%nvar, dest%var_pos, dest%var_scale and dest%var_offset...
    if (present(fix_var)) then
       !...pass 1 - count variables
       dest%nvar = 0
       do isrc = 1, src%nvar
          if (.not. fix_var(isrc)) dest%nvar = dest%nvar + 1
       end do
       allocate(dest%var_pos(dest%nvar,2))
       allocate(dest%var_scale(dest%nvar))
       allocate(dest%var_offset(dest%nvar))
       !...pass 2 - copy variable positions and scaling
       idest = 0
       do isrc = 1, src%nvar
          if (.not. fix_var(isrc)) then
             idest = idest + 1
             dest%var_pos(idest,:) = src%var_pos(isrc,:)
             dest%var_scale(idest) = src%var_scale(isrc)
             dest%var_offset(idest) = src%var_offset(isrc)
          end if
       end do
    else
       !all variables are still variable!
       dest%nvar = src%nvar
       allocate(dest%var_pos(dest%nvar,2))
       dest%var_pos = src%var_pos
       allocate(dest%var_scale(dest%nvar))
       dest%var_scale = src%var_scale
       allocate(dest%var_offset(dest%nvar))
       dest%var_offset = src%var_offset
    end if

  end subroutine allparam_copy

  !==========================================================================

  !! Set scaling for future calls to allparam_inlimit, _setvar, and _setone
  subroutine allparam_setscale(this, var_scale)

    !subroutine arguments
    type(allparam), intent(inout) :: this
    double precision, intent(in) :: var_scale(:)

    if (.not. this%done_init) stop 'allparam instance not initialised'
    this%var_scale = var_scale

  end subroutine allparam_setscale

  !==========================================================================

  !! Set offset for future calls to allparam_inlimit, _setvar, and _setone
  subroutine allparam_setoffset(this, var_offset)

    !subroutine arguments
    type(allparam), intent(inout) :: this
    double precision, intent(in) :: var_offset(:)

    if (.not. this%done_init) stop 'allparam instance not initialised'
    this%var_offset = var_offset

  end subroutine allparam_setoffset

  !==========================================================================

  !! Disable variable scaling for future calls to allparam_inlimit,
  !! _setvar, and _setone
  subroutine allparam_setnoscale(this)

    !subroutine arguments
    type(allparam), intent(inout) :: this

    if (.not. this%done_init) stop 'allparam instance not initialised'
    this%var_scale = 1D0

  end subroutine allparam_setnoscale

  !==========================================================================

  !! Disable variable offsetting for future calls to allparam_inlimit,
  !! _setvar, and _setone
  subroutine allparam_setnooffset(this)

    !subroutine arguments
    type(allparam), intent(inout) :: this

    if (.not. this%done_init) stop 'allparam instance not initialised'
    this%var_offset = 0D0

  end subroutine allparam_setnooffset

  !==========================================================================

  !! Print values of allparam components
  subroutine allparam_print(this)

    !subroutine arguments
    type(allparam), intent(in) :: this

    if (.not. this%done_init) stop 'allparam instance not initialised'
    print *, 'param:', this%param
    print *, 'svar:', this%svar
    print *, 'nvar:', this%nvar
    print *, 'var_pos:', this%var_pos
    print *, 'var_scale:', this%var_scale
    print *, 'var_offset:', this%var_offset

  end subroutine allparam_print

  !==========================================================================

  !! Return true if all supplied (scaled) variable values within limits
  function allparam_inlimit(this, var) result(ok)

    logical :: ok

    !subroutine arguments
    type(allparam), intent(in) :: this
    double precision, intent(in) :: var(:)

    !local variables
    integer :: ivar, indx1, indx2
    double precision :: scale, offset

    if (.not. this%done_init) stop 'allparam instance not initialised'
    ok = .false.
    do ivar = 1, this%nvar
       indx1 = this%var_pos(ivar,1)
       indx2 = this%var_pos(ivar,2)
       scale = this%var_scale(ivar)
       offset = this%var_offset(ivar)
       if ((var(ivar)*scale + offset) < this%limits(indx1, indx2, 1)) return
       if ((var(ivar)*scale + offset) > this%limits(indx1, indx2, 2)) return
    end do
    ok = .true.

  end function allparam_inlimit

  !============================================================================

  !! Update all variable parameters from supplied (scaled) values
  !!
  !! Doesn't check against limits
  !! Want to optimise this, as its needed for each posterior evaluation
  subroutine allparam_setvar(this, var)

    !subroutine arguments
    type(allparam), intent(inout) :: this
    double precision, intent(in) :: var(:)

    !local variables
    integer :: ivar

    if (.not. this%done_init) stop 'allparam instance not initialised'
    do ivar = 1, this%nvar
       this%svar(ivar) = var(ivar)*this%var_scale(ivar) + this%var_offset(ivar)
       this%param(this%var_pos(ivar,1), this%var_pos(ivar,2)) = this%svar(ivar)
    end do

  end subroutine allparam_setvar

  !============================================================================

  !! Update single variable parameter from supplied (scaled) value
  !! Doesn't check against limits
  subroutine allparam_setone(this, ivar, val)

    !subroutine arguments
    type(allparam), intent(inout) :: this
    integer, intent(in) :: ivar
    double precision, intent(in) :: val

    if (.not. this%done_init) stop 'allparam instance not initialised'
    this%svar(ivar) = val*this%var_scale(ivar) + this%var_offset(ivar)
    this%param(this%var_pos(ivar,1), this%var_pos(ivar,2)) = this%svar(ivar)

  end subroutine allparam_setone

  !============================================================================

  !! Free storage for components of allparam instance
  subroutine allparam_free(this)

    type(allparam), intent(inout) :: this

    if (.not. this%done_init) return !harmless
    this%done_init = .false.
    deallocate(this%param)
    deallocate(this%svar)
    deallocate(this%limits)
    deallocate(this%var_pos)
    deallocate(this%var_scale)
    deallocate(this%var_offset)

  end subroutine allparam_free

  !============================================================================

end module Wrap
