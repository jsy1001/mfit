module nestwrapper

  use Nested  !this module built as part of MultiNest
  use Wrap
  use Model
  use Bayes

  implicit none

  private

  public :: nest_Sample, nest_root

  type(allparam), save :: nest_param


  ! Parameters for Nested Sampler

  !whether to do multimodal sampling
  logical nest_mmodal 

  !sample with constant efficiency
  logical nest_ceff
  parameter(nest_ceff=.false.)

  !max no. of live points
  integer nest_nlive
  parameter(nest_nlive=500)
       
  !tot no. of parameters, should be nvar in most cases but if you need to
  !store some additional parameters with the actual parameters then
  !you need to pass them through the likelihood routine
  integer nest_nPar 

  !seed for nested sampler, -ve means take it from sys clock
  integer nest_rseed 
  parameter(nest_rseed=-1)

  !evidence tolerance factor
  real*8 nest_tol 
  parameter(nest_tol=0.5)

  !enlargement factor reduction parameter
  real*8 nest_efr
  parameter(nest_efr=0.5d0)

  !root for saving posterior files
  character*100 nest_root
  parameter(nest_root='chains/2-')

  !no. of iterations after which the ouput files should be updated
  integer nest_updInt
  parameter(nest_updInt=100)

  !null evidence (set it to very high negative no. if null evidence is unknown)
  real*8 nest_Ztol
  parameter(nest_Ztol=-1.d90)

  !max modes expected, for memory allocation
  integer nest_maxModes 
  parameter(nest_maxModes=10)

  !no. of parameters to cluster (for mode detection)
  integer nest_nClsPar

  !whether to resume from a previous run
  logical nest_resume
  parameter(nest_resume=.false.)
	
  integer nest_maxDim
  parameter(nest_maxDim=20)
  !parameters to wrap around (0 is F & non-zero T)
  integer nest_pWrap(nest_maxDim)

  !feedback on the sampling progress?
  logical nest_fb 
  parameter(nest_fb=.true.)

contains

  !-----*-----------------------------------------------------------------

  subroutine nest_Sample(nvar, mmodal)

    integer, intent(in) :: nvar
    logical, intent(in) :: mmodal

    integer var_pos(nvar,2)
    character(len=model_desc_len) var_desc(nvar)
    double precision var_scale(nvar), var_offset(nvar)
    double precision value0, width
    integer i, context

    if(nvar > nest_maxDim)  stop 'Too many variables. Need to increase nest_maxDim and recompile'

    !initialise nest_param
    call model_getvar(nvar, var_pos, var_desc)
    call allparam_init(nest_param, model_param, model_limits, nvar, var_pos)

    !set scale and offset from unit hypercube
    do i = 1, nvar
       value0 = model_param(var_pos(i,1), var_pos(i,2)) !parameter value
       width = model_prior(var_pos(i,1), var_pos(i,2)) !parameter prior width
       var_scale(i) = 2d0*width
       var_offset(i) = value0 - width
    end do
    call allparam_setscale(nest_param, var_scale)
    call allparam_setoffset(nest_param, var_offset)

    nest_nPar = nvar
    nest_nClsPar = nvar
    nest_mmodal = mmodal
    nest_pWrap = 1D0  !enable wrapping
    call nestRun(nest_mmodal, nest_ceff, nest_nlive, nest_tol, nest_efr, &
         nvar, nest_nPar, nest_nClsPar, nest_maxModes, nest_updInt, &
         nest_Ztol, nest_root, nest_rseed, nest_pWrap, &
         nest_fb, nest_resume, getLogLike, context)

    call allparam_free(nest_param)

  end subroutine nest_Sample

  !-----*-----------------------------------------------------------------

  ! Wrapper around Likelihood Function
  ! Cube(1:n_dim) has nonphysical parameters
  ! scale Cube(1:n_dim) & return the scaled parameters in Cube(1:n_dim) &
  ! additional parameters in Cube(n_dim+1:nPar)
  ! return the log-likelihood in lnew
  subroutine getLogLike(Cube, n_dim, nPar, lnew, context)

    integer, intent(in) :: n_dim, nPar, context
    real*8, intent(inout) :: Cube(nPar)
    real*8, intent(out) :: lnew

    !call your loglike function here   
    call allparam_setvar(nest_param, Cube(1:n_dim))
    Cube(1:n_dim) = nest_param%svar  !scaled (physical) parameters
    lnew = -likelihood(vis_data, triple_data, model_spec, nest_param%param)

  end subroutine getLogLike

  !-----*-----------------------------------------------------------------

end module nestwrapper
