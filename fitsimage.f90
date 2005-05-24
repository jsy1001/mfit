!$Id: fitsimage.f90,v 1.1 2005/05/24 12:53:06 jsy1001 Exp $

module fitsimage

  implicit none

contains

  subroutine writefits2d(outfile, image, nx, ny, status)

    ! Routine to write out a 2d fits image file
    !
    ! Based on Fortran 77 code by DFB/CAH/JSY
    ! Needs to be linked with fitsio library 
    character(len=*), intent(in) :: outfile
    integer, intent(in) :: nx, ny
    integer, intent(inout) :: status
    real, dimension(:,:), intent(in) :: image

    integer iunit, bitpix, naxis, pcount, gcount
    integer, dimension(2) :: naxes
    integer group
    logical simple, extend
    character(len=30) :: errtxt

    if (status /= 0) return

    !open the new FITS file
    iunit=77
    call ftinit(iunit, outfile, 2880, status)

    !set up some variables
    simple=.true.
    bitpix=-32
    naxis=2
    naxes(1)=nx
    naxes(2)=ny
    pcount=0
    gcount=0
    extend=.false.

    !write the required primary array keywords
    call ftphpr(iunit,simple,bitpix,naxis,naxes,pcount,gcount, &
         extend,status)
    group=1

    !define primary array structure
    call ftrdef(iunit,status)

    !write the primary array of data
    call ftp2de(iunit,group,size(image,1),nx,ny,image,status)

    !close the file and check all is well
    call ftclos(iunit,status)

    if (status <= 0)then
       print *,'*** Wrote fits file successfully ***'
    else
       call ftgerr(status,errtxt)
       print *,'*** ERROR - fits writer did not run successfully ***'
       print *,'status =',status,': ',errtxt
    end if

  end subroutine writefits2d

end module fitsimage
