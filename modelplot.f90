!$Id: modelplot.f90,v 1.2 2005/06/28 16:13:27 jsy1001 Exp $

! Plot sky brightness distribution for mfit-format model
!
! Based on modelplot (for CITVLB-format model), last modified 29/04/99

program Modelplot

  use Model
  use Maths
  use Fitsimage
  use Inout !just for revision no.

  implicit none

  !parameters
  integer unit, maxcmp, npts_max, nlevs_max, max_wb
  parameter(unit=11, npts_max=1500, nlevs_max=8, max_wb=10)
  real scale
  parameter(scale=1000.)

  !local variables
  integer status, istat, i, j, ipar, cmp, npts, nlevs
  real sum, plotxylow, plotxyhigh, xx, yy, disp, peak
  real, dimension(npts_max,npts_max) :: plotxy
  character(len=128) :: modfil, info, outbase, outfil, title
  character(len=32) :: cvs_rev, revision, tag

  !model parameters
  integer iwb
  real radius, theta, x, y, flux
  real axis, phi, epsilon, alpha

  !plot parameters
  logical dogray, invgray, writefits, docont, markcen
  integer symbol
  real xmax
  real, dimension(nlevs_max) :: plevs, levs
  real, dimension(6) :: tr
  real fg, bg
  data plevs/0.5,1.0,2.0,4.0,8.0,16.0,32.0,64.0/

  !functions
  integer pgopen

  cvs_rev = '$Revision: 1.2 $'
  revision = cvs_rev(scan(cvs_rev, ':')+2:scan(cvs_rev, '$', .true.)-1)
  print *,' '
  print *,'  mplot - greyscale/contour plot of mfit-format model'
  print *,'  package release ',release
  print *,'  [mplot revision ',trim(revision),']'
  print *,' '

  !read the model
  do
     print *, 'Enter model filename in current directory:'
     read (*, '(a)') modfil
     print *, ' '
     print *, 'reading model...'
     info = ''
     call read_model(info, modfil)
     !not necessary to normalise fluxes of model components,
     !as greyscale and contours are normalised to max. flux in image
     if (info == '') exit !success
     print *, trim(info)
  end do
  call print_model

  !check for unsupported model features
  do cmp = 1, size(model_spec, 1)
     
     do iwb = 1, nwave
        ipar = 7 + (4*model_wldep(1)+model_wldep(2))*(nwave-1) &
             + 3*model_wldep(3)*(iwb-1)
        epsilon = model_param(cmp, ipar+2)
        if (epsilon /= 1D0) &
             stop '*** elliptical model component not supported'
     end do
     
     select case (trim(model_spec(cmp,3)))
        
     case ('taylor')
     case ('square-root')
     case ('gauss-hermite')
        print *, '*** ld_type ', model_spec(cmp,3), ' not supported'
        stop
     case default
        if (model_spec(cmp,3)(1:1) == '<') &
             stop '*** numerical limb-darkened disc not supported'
     end select

  end do

  !setup for plotting
  print *, ' '
  print *, 'Enter max. radius to plot (mas):'
  read (*, *) xmax
  istat = pgopen('?')
  call pgsubp(nwave, 1)
  call pgpage
  call pgsch(1.5)
  call pgvstand
  call pgwnad(xmax, -xmax, -xmax, xmax)

  !get user input
  print *, 'Enter no. grid points /axis:'
  read (*, *) npts
  if (npts > npts_max) then
     npts = npts_max
     print *, 'Too many points. using ', npts_max
  endif
  dogray = yesno('Plot grayscale', 'yes')
  if (dogray) &
       invgray = yesno('Invert grayscale (black=bright)', 'no')
  writefits = yesno('Write FITS image(s)', 'no')
  if (writefits) then
     print *, 'Enter FITS base filename:'
     read (*, '(a)') outbase
  end if
  docont = yesno('Plot contours', 'yes')
  if (docont) then
     print *, 'Default contours in percent:'
     print *,  plevs
     if (yesno('Change contour levels', 'no')) then
        do
           print *, 'Enter no. of contours (max. ', nlevs_max, ')'
           read (*, *) nlevs
           if (nlevs <= nlevs_max) exit
        end do
        print *, 'Enter', nlevs, ' contours in percent:'
        read (*, *) plevs(:nlevs)
     else
        nlevs = nlevs_max
     end if
  end if
  markcen = yesno('Mark component centres', 'yes')
  print *, 'Enter title for plot:'
  read (*, '(a)') title

  !loop over wavebands in model
  do iwb = 1, nwave

     !build up model brightness distribution in 2d array
     print *, 'Building up model...'

     !mapping from array elements to sky coordinates
     tr(2) = 2.*xmax/real(npts-1)
     tr(6) = 2.*xmax/real(npts-1)
     tr(1) = -xmax - tr(2)
     tr(4) = -xmax - tr(6)
     tr(3) = 0.
     tr(5) = 0.

     plotxyhigh = 0.
     plotxylow = 0.
     sum = 0.

     !loop over 2d array
     do i = 1, npts

        do j = 1, npts

           plotxy(i,j) = 0.

           !centre of current pixel /mas
           xx = tr(1) + i*tr(2) + j*tr(3)
           yy = tr(4) + i*tr(5) + j*tr(6)

           !loop over model components
           do cmp = 1, size(model_spec, 1)
           
              !get parameters for this cpt
              radius = model_param(cmp, 2+4*model_wldep(1)*(iwb-1))
              theta = deg2rad*model_param(cmp, 3+4*model_wldep(1)*(iwb-1))
              flux = model_param(cmp, &
                   6+4*model_wldep(1)*(nwave-1)+model_wldep(2)*(iwb-1))
              ipar = 7 + (4*model_wldep(1)+model_wldep(2))*(nwave-1) &
                   + 3*model_wldep(3)*(iwb-1)
              axis = model_param(cmp, ipar)
              phi = deg2rad*model_param(cmp, ipar+1)
              epsilon = model_param(cmp, ipar+2)
              x = radius*sin(theta)
              y = radius*cos(theta)

              !displacement of component from centre of current pixel
              disp = sqrt((x-xx)**2 + (y-yy)**2)

              select case (trim(model_spec(cmp,3)))
                 
              case ('gaussian')
                 ! - Gaussian
                 peak = scale*flux*4*log(2.)/(pi*axis**2)
                 plotxy(i,j) = plotxy(i,j) &
                      + peak*exp(-4*log(2.)*(disp/axis)**2)
              case ('uniform')
                 ! - uniform disk
                 peak = scale*flux*4/(pi*axis**2)
                 if(disp < axis/2) &
                      plotxy(i,j) = plotxy(i,j) + peak
              case ('hestroffer')
                 ! - Hestroffer limb-darkened disk
                 ipar = 10 + (4*model_wldep(1) + model_wldep(2) &
                      + 3*model_wldep(3))*(nwave-1) &
                      + max_order*model_wldep(4)*(iwb-1)
                 alpha = model_param(cmp, ipar)
                 peak = scale*flux*(4+2*alpha)/(pi*axis**2)
                 plotxy(i,j) = plotxy(i,j) &
                      + peak*inten_hest(alpha, axis, disp)
              case default
                 ! - other types (unsupported for now)
                 ! should be trapped above
                 print *,'*** Ignoring model component ', cmp, &
                      ': ld_type ', model_spec(cmp,3), ' not supported.'
              end select

           end do

           sum = sum + plotxy(i,j)*tr(2)*tr(6)
           if (plotxy(i,j) < plotxylow) plotxylow = plotxy(i,j)
           if (plotxy(i,j) > plotxyhigh) plotxyhigh = plotxy(i,j)

        end do

     end do

     print *, 'Total, max, min: ', sum/scale, plotxyhigh/scale, &
          plotxylow/scale

     !write out image, if desired
     if (writefits) then
        write (tag, '(i1)') iwb
        outfil = trim(outbase)//'_'//trim(tag)//'.fits'
        status = 0
        call writefits2d(outfil, plotxy, npts, npts, status)
     endif

     !do greyscale/contour plot
     if (dogray) then
        if (invgray) then
           fg = 0.
           bg = 2.*plotxyhigh
        else
           fg = plotxyhigh
           bg = 0.
        endif
        call pggray(plotxy, npts_max, npts_max, 1, npts, 1, npts, fg, bg, tr)
     end if
      
     if (docont) then
        do i = 1, nlevs
           levs(i) = plevs(i)/100.*plotxyhigh
        end do
        call pgsci(2)
        call pgcont(plotxy, npts_max, npts_max, 1, npts, 1, &
             npts, levs, nlevs, tr) 
        call pgsci(1)
     endif

     !mark the centre of all components
     if (markcen) then
        call pgsci(3)
        do cmp = 1, size(model_spec, 1)
           select case (trim(model_spec(cmp,3)))
           case ('gaussian')
              if (cmp == 1) then
                 symbol = 5
              else
                 symbol = 3
              end if
           case default
              symbol = 5
           end select
           radius = model_param(cmp, 2+4*model_wldep(1)*(iwb-1))
           theta = deg2rad*model_param(cmp, 3+4*model_wldep(1)*(iwb-1))
           call pgpt1(radius*sin(theta), radius*cos(theta), symbol)
        end do
     end if

     !finally, draw axes and labels
     call pgsci(1)
     call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
     call pglabel('Relative RA (mas)', 'Relative declination (mas)', &
          trim(title))

     if (iwb /= nwave) call pgpage

  end do

  call pgend

contains
      
  function inten_hest(alpha, diam, r)

    !Return intensity for Hestroffer (1997) A&A 327, 199
    !limb-darkened disk model, given:
    !
    !alpha   limb-darkening parameter
    !diam    zero-intensity diameter /mas
    !r       radius /mas
    real, intent(in) :: alpha, diam, r
    real mu
    real inten_hest
      
    if (r > diam/2.0) then
       inten_hest = 0.0
       return
    endif
    mu = sqrt(1.0 - (2.0*r/diam)**2)
    inten_hest = mu**alpha

  end function inten_hest

end program Modelplot
