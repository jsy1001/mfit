!$Id: binnest.f90,v 1.1 2009/07/14 16:35:21 jsy1001 Exp $

program Main

  use f2kcli
  use nestwrapper
  use ReadMC

  implicit none
  
  !============================================================================
  !variables
  integer, parameter :: width = 78            !for spacer lines
  character(len=width) :: spacer_line
  integer narg, iarg, mode, n, nsamp
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
