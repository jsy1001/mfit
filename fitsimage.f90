!$Id: fitsimage.f90,v 1.2 2010/04/08 11:06:02 jsy1001 Exp $

module fitsimage

  implicit none

contains

  subroutine writefits2d(outfile, image, nx, ny, pixscal, status)

    ! Routine to write out a 2d fits image file
    !
    ! Based on Fortran 77 code by DFB/CAH/JSY
    ! Needs to be linked with fitsio library

    !subroutine arguments
    !! FITS filename
    character(len=*), intent(in) :: outfile
    !! Dimensions of image to write
    integer, intent(in) :: nx, ny
    !! Status variable
    integer, intent(inout) :: status
    !! Image pixel array
    real, dimension(:,:), intent(in) :: image
    !! Pixellation [milliarcsec/pixel]
    real, intent(in) :: pixscal

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

    !write PIXEL keyword as required by vis_sim
    call ftpkyf(iunit,'PIXEL',pixscal,4, &
         '[mas/pixel] Pixellation',status)

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
