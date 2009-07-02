!$Id: readmc.f90,v 1.1 2009/07/02 13:48:56 jsy1001 Exp $

module ReadMC

  implicit none

  private

  public :: mc_count_lines, mc_get_params

contains
  
  !! Count number of lines in specified text file
  subroutine mc_count_lines(filename, nline)
 
    !! Pathname of text file
    character(len=*), intent(in) :: filename
    !! On exit, number of lines in text file
    integer, intent(out) :: nline

    integer, parameter :: iunit = 12
    character(len=1) :: dummy

    nline = 0
    open(unit=iunit, action='read', err=2, file=filename)
    do
       read(iunit, *, end=2) dummy
       nline = nline + 1
    end do
2   close(iunit)

  end subroutine mc_count_lines


  !! Calculate parameter means from [root].txt file written by MultiNest
  !! Writes marginalised posterior pdf for each parameter
  subroutine mc_get_params(filename, nest_root, nsamp, nvar, mean, stddev)

    !! Pathname of text file
    character(len=*), intent(in) :: filename
    !! Prefix for posterior pdf files
    character(len=*), intent(in) :: nest_root
    !! Number of posterior samples in text file (e.g from mc_count_lines)
    integer, intent(in) :: nsamp
    !! Number of variable parameters
    integer, intent(in) :: nvar
    !! On exit, mean parameter values
    double precision, intent(out) :: mean(nvar)
    !! On exit, standard deviations of parameters
    double precision, intent(out) :: stddev(nvar)

    integer :: i, ivar
    integer :: ivar2(2)
    double precision :: prob(nsamp), loglike(nsamp)
    double precision :: var(nvar, nsamp)
    double precision :: sumsq(nvar)
    double precision :: sumwt
    double precision :: range(2)
    double precision :: range2(2, 2)
    character(len=20) :: suffix

    call mc_read_txt(filename, nsamp, nvar, prob, loglike, var)

    !calculate parameter means and standard deviations
    sumwt = 0d0
    mean = 0d0
    sumsq = 0d0
    do i = 1, nsamp
       sumwt = sumwt + prob(i)
       mean = mean + prob(i)*var(:, i)
       sumsq = sumsq + prob(i)*var(:, i)*var(:, i)
    end do
    mean = mean/sumwt
    stddev = sqrt(sumsq/sumwt - mean*mean)

    !write marginalised 1d pdf for each param
    do ivar = 1, nvar
       range(1) = mean(ivar) - 3d0*stddev(ivar)
       range(2) = mean(ivar) + 3d0*stddev(ivar)
       write(suffix, '(a, i1, a)') 'param_', ivar, '.dat'
       call mc_marg_1d(trim(nest_root)//suffix, & 
            nsamp, nvar, ivar, range, prob, loglike, var)
    end do

    !write marginalised 2d pdf for each pair of adjacent params
    do ivar = 1, nvar-1
       ivar2(1) = ivar
       ivar2(2) = ivar+1
       range2(:, 1) = mean(ivar2(:)) - 3d0*stddev(ivar2(:))
       range2(:, 2) = mean(ivar2(:)) + 3d0*stddev(ivar2(:))
       write(suffix, '(a, i1, a, i1, a)') 'param_', ivar2(1), '_', ivar2(2), &
            '.dat'
       call mc_marg_2d(trim(nest_root)//suffix, & 
            nsamp, nvar, ivar2, range2, prob, loglike, var)
    end do

  end subroutine mc_get_params

  
  !! Read [root].txt file written by MultiNest
  subroutine mc_read_txt(filename, nsamp, nvar, prob, loglike, var)

    !! Pathname of text file
    character(len=*), intent(in) :: filename
    !! Number of posterior samples in text file (e.g. from mc_count_lines)
    integer, intent(in) :: nsamp
    !! Number of variable parameters to read sample values for
    integer, intent(in) :: nvar
    !! On exit, sample probability (sample probability is the sample
    !! prior mass multiplied by its likelihood & normalized by the evidence)
    double precision, intent(out) :: prob(nsamp)
    !! On exit, -2*loglikelihood
    double precision, intent(out) :: loglike(nsamp)
    !! On exit, parameter values
    double precision, intent(out) :: var(nvar, nsamp)

    integer, parameter :: iunit = 12
    integer :: i

    open(unit=iunit, action='read', file=filename)
    do i = 1, nsamp
       read(iunit, *) prob(i), loglike(i), var(:, i)
    end do
    close(iunit)

  end subroutine mc_read_txt


  !! Write 1d marginalised posterior pdf for each parameter
  subroutine mc_marg_1d(outfname, nsamp, nvar, &
       ivar, range, prob, loglike, var)

    !! Prefix for posterior pdf files
    character(len=*), intent(in) :: outfname
    integer, intent(in) :: nsamp, nvar, ivar
    double precision, intent(in) :: range(2)
    double precision, intent(in) :: prob(nsamp), loglike(nsamp)
    double precision, intent(in) :: var(nvar, nsamp)

    integer, parameter :: ounit = 11
    integer, parameter :: nbin = 10
    integer :: ibin, i
    double precision :: binw, binl, like, wt, mlike

    open(unit=ounit, file=outfname, action='write')
    binw = (range(2)-range(1))/nbin
    do ibin = 1, nbin
       binl = range(1)+(ibin-1)*binw
       mlike = 1d-9
       do i = 1, nsamp
          if(var(ivar, i) > binl .and. var(ivar, i) <= (binl+binw)) then
             like = exp(loglike(i))
             !weight by (prior mass)/evidence
             wt = prob(i)/like
             mlike = mlike + wt*like  
          end if
       end do
       write(ounit, '(2f12.6)') binl+0.5d0*binw, log(mlike)
    end do
    close(ounit)

  end subroutine mc_marg_1d


  !! Write 2d marginalised posterior pdf for each parameter
  subroutine mc_marg_2d(outfname, nsamp, nvar, &
       ivar, range, prob, loglike, var)

    character(len=*), intent(in) :: outfname
    integer, intent(in) :: nsamp, nvar
    integer, intent(in) :: ivar(2)
    double precision, intent(in) :: range(2, 2)
    double precision, intent(in) :: prob(nsamp), loglike(nsamp)
    double precision, intent(in) :: var(nvar, nsamp)

    integer, parameter :: ounit = 11
    integer, parameter :: nbin = 10
    integer :: ibin, jbin, i
    double precision :: like, wt, mlike
    double precision :: binl(2), binw(2)

    open(unit=ounit, file=outfname, action='write')
    binw(:) = (range(:, 2)-range(:, 1))/nbin
    do ibin = 1, nbin
       binl(1) = range(1, 1)+(ibin-1)*binw(1)
       do jbin = 1, nbin
          binl(2) = range(2, 1)+(jbin-1)*binw(2)
          mlike = 1d-9
          do i = 1, nsamp
             if(var(ivar(1), i) > binl(1) &
                  .and. var(ivar(1), i) <= (binl(1)+binw(1)) &
                  .and. var(ivar(2), i) > binl(2) &
                  .and. var(ivar(2), i) <= (binl(2)+binw(2))) then
                like = exp(loglike(i))
                !weight by (prior mass)/evidence
                wt = prob(i)/like
                mlike = mlike + wt*like  
             end if
          end do
          write(ounit, '(3f12.6)') binl+0.5d0*binw, log(mlike)
       end do
       write(ounit, *) !for gnuplot
    end do
    close(ounit)

  end subroutine mc_marg_2d

end module ReadMC
  
