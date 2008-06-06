! $Id: bayes.f90,v 1.3 2008/06/06 13:02:22 jsy1001 Exp $

module Bayes

  use Maths
  use Visibility

  implicit none

  private

  !public subroutines/functions contained:

  public :: likelihood, prior, gof


  !public module variables contained:

  public :: vis_data, triple_data, sel_wavebands
  public :: num_vis, num_triple, num_wb

  !! Data arrays
  double precision, allocatable :: vis_data(:,:), triple_data(:,:)
  integer :: num_vis, num_triple
  
  !! Wavebands included in data (for wavelength-dependent models).
  !! Not used in this module, but should go with data
  double precision, allocatable :: sel_wavebands(:,:)
  integer :: num_wb

contains

  !============================================================================

  !! Compute the negative log likelihood of data given model
  function likelihood(vis_data, triple_data, model_spec, param)

    double precision :: likelihood

    !function arguments
    !! Visibility data (as module variable)
    double precision, intent(in) :: vis_data(:,:)
    !! Triple product data (as module variable)
    double precision, intent(in) :: triple_data(:,:)
    !! As model_spec: name, shape, LD type of each model cpt
    character(len=128), intent(in) :: model_spec(:,:)
    !! As model_param: all (variable and non-variable) model parameters
    double precision, intent(in) :: param(:,:)

    !local variables
    integer :: i
    double precision :: lambda, delta_lambda, u, v, data_vis, data_vis_err
    double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
    double precision :: data_phase_err, model_vis, model_amp, model_phase, mjd
    double complex :: vis, vis1, vis2, vis3

    likelihood = 0D0

    !sum over the visibility data points
    do i = 1, num_vis
       
       !extract points
       lambda = vis_data(i,1)
       delta_lambda = vis_data(i,2)
       u = vis_data(i,3)
       v = vis_data(i,4)
       data_vis = vis_data(i,5)
       data_vis_err = vis_data(i,6)
       mjd = vis_data(i,7)

       if (data_vis_err>0D0) then
          
          !compute model-predicted visibility amplitude squared
          vis = cmplx_vis(model_spec, param, lambda, delta_lambda, u, v, mjd)
          model_vis = (modulus(vis)**2D0)
      
          !compute contribution to likelihood
          likelihood = likelihood + &
               (((data_vis-model_vis)**2D0)/(2D0*(data_vis_err**2D0)))

       end if

    end do

    !sum over triple product amplitude and phase data points
    do i = 1, num_triple

       !extract points
       lambda = triple_data(i,1)
       delta_lambda = triple_data(i,2)
       u1 = triple_data(i,3)
       v1 = triple_data(i,4)
       u2 = triple_data(i,5)
       v2 = triple_data(i,6)
       data_amp = triple_data(i,7)
       data_amp_err = triple_data(i,8)
       data_phase = triple_data(i,9)
       data_phase_err = triple_data(i,10)
       mjd = triple_data(i,11)

       if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
          !compute model-predicted triple product
          vis1 = cmplx_vis(model_spec, param, lambda, delta_lambda, &
               u1, v1, mjd)
          vis2 = cmplx_vis(model_spec, param, lambda, delta_lambda, &
               u2, v2, mjd)
          vis3 = cmplx_vis(model_spec, param, lambda, delta_lambda, &
               -(u1+u2), -(v1+v2), mjd)
          vis = vis1*vis2*vis3
          model_amp = modulus(vis)
          model_phase = rad2deg*argument(vis)
       end if

       !compute likelihood contributions
       !nb phase calculations modulo 360 degrees
       if (data_amp_err>0D0) then
          likelihood = likelihood + &
               (((data_amp-model_amp)**2D0)/(2D0*(data_amp_err**2D0)))
       end if
       if (data_phase_err>0D0) then
          likelihood = likelihood + &
               (((modx(data_phase,model_phase))**2D0) &
                /(2D0*(data_phase_err**2D0)))
       end if
    end do

  end function likelihood

  !===========================================================================

  !! Return negative log prior probability
  function prior(x_pos, param, model_param, model_prior)

    double precision :: prior

    !function arguments
    integer, intent(in) :: x_pos(:,:)
    double precision, intent(in) :: param(:,:)
    double precision, intent(in) :: model_param(:,:)
    double precision, intent(in) :: model_prior(:,:)

    !local variables
    integer :: i
    double precision :: value, value0, err

    prior = 0D0
    
    !loop over the variable parameters
    do i = 1, size(x_pos,1)

       value0 = model_param(x_pos(i,1),x_pos(i,2)) !parameter original value
       err = model_prior(x_pos(i,1),x_pos(i,2)) !parameter prior width
       value = param(x_pos(i,1),x_pos(i,2)) !parameter value as it is varied
       
       !calculate normalised contribution to the prior
       prior = prior + (((value-value0)**2D0)/(2D0*(err**2D0))) + log(err*sqrt(2D0*pi))

    end do

  end function prior

  !==========================================================================

  !! Calculate goodness-of-fit
  !!
  !! Computes the chisqrd value of a model fit and data set:
  !!  chisqrd = sum[i] { (data[i]-theory[i])^2/data_error[i]^2 }
  !!
  !! chisqrd is thus almost identical to likelihood
  !! could save code here and significantly combine this function
  !! with the likelihood one but maybe less risky to keep them separated...
  subroutine gof(spec, param, chisqrd)

    !subroutine arguments
    !! As model_spec: name, shape, LD type of each model cpt
    character(len=128), intent(in) :: spec(:,:)
    !! As model_param: all (variable and non-variable) model parameters
    double precision, intent(in) :: param(:,:)
    !! Calculated chi-squared
    double precision, intent(out) :: chisqrd

    !local variables
    integer :: i
    double precision :: lambda, delta_lambda, u, v, data_vis, data_vis_err
    double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
    double precision :: data_phase_err, model_vis, model_amp, model_phase, mjd
    double complex :: vis, vis1, vis2, vis3

    chisqrd = 0D0

    !sum over the visibility data points
    do i = 1, num_vis
       
       !extract points
       lambda = vis_data(i,1)
       delta_lambda = vis_data(i,2)
       u = vis_data(i,3)
       v = vis_data(i,4)
       data_vis = vis_data(i,5)
       data_vis_err = vis_data(i,6)
       mjd = vis_data(i,7)

       !ignore if -ve or zero error:
       if (data_vis_err>0D0) then

          !compute model-predicted visibility amplitude squared
          vis = cmplx_vis(spec, param, lambda, delta_lambda, u, v, mjd)
          model_vis = (modulus(vis))**2D0

          !compute contribution to chisqrd
          chisqrd = chisqrd + &
               (((data_vis-model_vis)/data_vis_err)**2D0)

       end if

    end do

    !sum over triple product amplitude and phase data points
    do i = 1, num_triple

       !extract points
       lambda = triple_data(i,1)
       delta_lambda = triple_data(i,2)
       u1 = triple_data(i,3)
       v1 = triple_data(i,4)
       u2 = triple_data(i,5)
       v2 = triple_data(i,6)
       data_amp = triple_data(i,7)
       data_amp_err = triple_data(i,8)
       data_phase = triple_data(i,9)
       data_phase_err = triple_data(i,10)
       mjd = triple_data(i,11)

       if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
          !compute model-predicted triple amplitude
          vis1 = cmplx_vis(spec, param, lambda, delta_lambda, u1, v1, mjd)
          vis2 = cmplx_vis(spec, param, lambda, delta_lambda, u2, v2, mjd)
          vis3 = cmplx_vis(spec, param, lambda, delta_lambda, &
               -(u1+u2), -(v1+v2), mjd)
          vis = vis1*vis2*vis3
          model_amp = modulus(vis)
          model_phase = rad2deg*argument(vis)
       end if

       !compute rmsdev and chisqrd contributions
       !nb phase calculations modulo 360 degrees
       if (data_amp_err>0D0) then
          chisqrd = chisqrd + &
               (((data_amp-model_amp)/(data_amp_err))**2D0)
       end if
       if (data_phase_err>0D0) then
          chisqrd = chisqrd + &
               (modx(data_phase,model_phase)/(data_phase_err))**2D0

       end if

    end do

  end subroutine gof

  !==========================================================================

end module Bayes
