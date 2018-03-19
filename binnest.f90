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

program Main

  use f2kcli
  use nestwrapper
  use ReadMC

  implicit none
  
  !============================================================================
  !variables
  integer, parameter :: width = 78            !for spacer lines
  character(len=width) :: spacer_line
  integer narg, iarg, mode, n, nsamp, i
  character(len=128) :: switch, arg, sampFilename
  double precision, allocatable :: sol(:), err(:)

  !----------------------------------------------------------------------------
  !Formatting stuff
  spacer_line = repeat('-', width)

  !----------------------------------------------------------------------------
  !Introduction

  print *,' '
  print *,spacer_line
  print *,' '
  print *,'  Bin nested samples'
  print *,' '
  print *,spacer_line

  !----------------------------------------------------------------------------
  !Parse command-line arguments

  narg = command_argument_count()
  iarg = 1
  !defaults
  mode = -1
  do
     if (iarg == narg + 1) exit
     if (iarg > narg - 1) then
        call print_usage()
        stop 'Not enough arguments'
     end if
     call get_command_argument(iarg, switch)
     select case(switch)
        case('-n', '--mode')
           call get_command_argument(iarg+1, arg)
           read(arg, *) mode
           iarg = iarg + 1
        case default
           print *, 'Ignoring invalid command-line argument: ', arg
        end select
        iarg = iarg + 1
  end do

  !display results
  if(mode > 0) then
     sampFilename = trim(nest_root)//'post_separate.dat'
  else
     sampFilename = trim(nest_root)//'.txt'
  end if
  call mc_count_lines(sampFilename, mode, nsamp)
  call mc_count_vars(sampFilename, n)
  print '(1x, a, i6, a, i3, a, a)', &
       'Reading', nsamp, ' samples for', n, ' variables from ', sampFilename
  print *, 'Mean results are in ', trim(nest_root)//'stats.dat'
  allocate(sol(n), err(n))
  call mc_get_params(sampFilename, nest_root, nsamp, n, mode, sol, err)
  print '(1x, a, 10x, a)', '   ', '  mean     standard'
  print '(1x, a, 10x, a)', 'num', ' value    deviation'
  do i = 1, n
     write(*,62) i, sol(i), err(i)
  end do
62 format(' (', i2, ') ', 1x, f13.6, 1x, f12.6) 

  !-------------------------------------------------------------------------
  !Deallocate storage
  if (allocated(sol)) deallocate(sol)
  if (allocated(err)) deallocate(err)

contains

!==============================================================================

  subroutine print_usage()

    print *, 'Usage:'
    print *, ' '
    print *, 'binnest [options]'
    print *, ' '
    print *, 'Options (may appear in any order):'
    print *, ' '
    print *, '-n|--mode               process only samples for this mode'

  end subroutine print_usage

end program Main
