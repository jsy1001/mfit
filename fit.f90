module Fit

!callable subroutines contained
!
!minimiser - wrapper to the minimising routine. this subroutine preps
!            the arrays into true parameter space and fixed model space
!            before calling the actual minimising alogrithm
!
!local subroutines contained
!
!simplex - downhill simplex method minimisation based on numerical 
!          recipies amoeba code
!extrap - extrapolation of new simplex vertex position based on
!         numerical recipies amotry code
!posterior - returns negative log posterior of model and data
!likelihood - returns negative log likelihood of data given model
!prior - returns negative log prior of current model 

use Maths
use Visibility

implicit none

contains

!==============================================================================

subroutine minimiser(info, model_spec, model_param, model_prior, limits, &
                     fit_param, vis_data, triple_data, symm, sol, flag, var, &
                     hes, cov, cor, chisqrd, sumsqr, nlposterior)

!returns the fit_param array of best fitting model parameters to the data
!flag holds sucess state of the minimisation:
!=0 minimisation exceeded number of iterations - bounds provided only
!=1 minimsation suceeded but with -ve on covariance diagonals
!=2 minimisation suceeded with fully valid covariance and correlation matrices

!symm is a logical that forces centrosymmetric only models

!subroutine arguments
character(len=128) :: info
character(len=128), dimension(:,:), allocatable :: model_spec
double precision, dimension(:,:), allocatable :: model_param, model_prior 
double precision, dimension(:,:), allocatable :: fit_param
double precision, dimension(:,:), allocatable :: vis_data, triple_data
double precision, dimension(17,2) :: limits
character(len=35), dimension(:), allocatable :: var
double precision, dimension(:,:), allocatable :: sol, hes, cov, cor
double precision :: nlposterior, nlevidence, chisqrd, sumsqr
integer :: flag
logical :: symm

!local variables
integer :: i, j, k, n
character(len=128) :: name, comp
character(len=2), dimension(10) :: numbers
double precision :: P, Pl, Pu, Pi, Pj, deltai, deltaj, eta, diff
double precision, dimension(:), allocatable ::  x, temp_x
double precision, dimension(:,:), allocatable :: x_info
integer, dimension(:,:), allocatable :: x_pos
logical :: illegal, ok

!go through model arrays and count number of variable parameters
n = 0
do i = 1, size(model_param,1)
   do j = 1, 17
      if (model_prior(i,j) /= 0D0) n=n+1
   end do
end do

!allocate variable array and read in the variables
!x is a vector holding the values of the model parameters that are allowed
!to vary, x_pos holds the positions of the values in the model_param array
!that defines the model.
allocate(x(n))
allocate(temp_x(n))
allocate(x_pos(n,2))
k = 0
do i = 1, size(model_param,1)
   do j = 1, 17
      if (model_prior(i,j) /= 0D0) then
         k = k+1
         x_pos(k,1) = i
         x_pos(k,2) = j
         x(k) = model_param(i,j)
      end if
   end do
end do

!assign names to the variables listed in the sol solution array
allocate(var(n))
var = ''
numbers = (/' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10'/)
do i = 1, n
   comp = adjustl(numbers(x_pos(i,1)))
   select case(x_pos(i,2))
      case (2)
         name = 'position radius (mas)'
      case (3)
         name = 'position angle (deg)'
      case (4)
         name = 'flux (arb units)'
      case (5)
         name = 'major axis (mas)'
      case (6)
         name = 'orientation (deg)'
      case (7)
         name = 'ellipticity'
      case default
         name = 'LD parameter ' // adjustl(numbers(x_pos(i,2)-7))
   end select
   var(i) = 'component ' // trim(comp) // ', ' // trim(name)
end do

!check for illegal types of minimisation:
illegal = .false.
!must have at least 1 freedom to minimise with
if (n==0) illegal = .true.
!for single component model cannot vary r/theta (change in position
!only contstant phase offset) or B (flux is normalised anyway) 
if (size(model_param,1)==1) then
   do i = 1, size(x_pos,1)
      if (x_pos(i,2)==2) illegal = .true.
      if (x_pos(i,2)==3) illegal = .true.
      if (x_pos(i,2)==4) illegal = .true.
   end do
end if
!for any component cannot vary theta if r is zero and not free to vary
!(position angle has no effect if position radius fixed at zero)
do i = 1, size(model_param,1)
   if ((model_param(i,2)==0D0) .and. (model_prior(i,2)==0D0) &
      .and. (model_prior(i,3)/=0D0)) illegal = .true.
end do
!for any component cannot vary phi if epsilon is unity and not free to vary
!(orientation is meaningless is ellipse reduced to circular disc)
do i = 1, size(model_param,1)
   if ((model_param(i,7)==1D0) .and. (model_prior(i,7)==0D0) &
        .and. (model_prior(i,6)/=0D0)) illegal = .true.
end do
if (illegal) goto 90

!if centrosymmetric model is forced then maust have 
!eccentricity epsilon fixed to be unity and position radius fixed at zero
if (symm) then
   do i = 1, size(model_param,1)
      if (.not.((model_param(i,7)==1D0).and.(model_prior(i,7)==0D0))) &
           goto 92
      if (.not.((model_param(i,2)==0D0).and.(model_prior(i,2)==0D0))) &
           goto 92
   end do
end if

!allocate fit_param array
!fit_param is identical to model_param, it holds the definition of the model
!with the difference that fit_param is designed to vary as the model is
!adjusted to fit. x_pos is used to update the variable parameters in 
!the array with their values from x
if (allocated(fit_param)) deallocate(fit_param)
allocate(fit_param(size(model_param,1),size(model_param,2)))
fit_param = model_param

!compile x_info array
!1st column is of typical length scales for variables and set equal to the
!prior width. 2nd and 3rd columns hold lower and upper bounds on the
!value of the parameter in question - used to prevent negative/other illegal
!values causing problems (esp in visibility calculations)
allocate(x_info(n,3))
do i = 1, n
   x_info(i,1) = model_prior(x_pos(i,1),x_pos(i,2))
end do
!preset limits to values from the limits array used in model input
do i = 1, n
   x_info(i,2) = limits(x_pos(i,2),1)
   x_info(i,3) = limits(x_pos(i,2),2)
   !alpha(1) in hestroffer model must be > 0
   if (x_pos(i,2)==8) then
      if (trim(model_spec(x_pos(i,1),3))=='hestroffer') x_info(i,2) = 0D0
   end if
end do

!allocate var, sol, hes, cov, cor arrays
allocate(sol(n,6))
allocate(hes(n,n))
allocate(cov(n,n))
allocate(cor(n,n))
sol = 0D0
hes = 0D0
cov = 0D0
cor = 0D0

print *, 'calculating initial goodness-of-fit...'
!calculate goodness-of-fit statistic chi squared
call gof(chisqrd, sumsqr, model_spec, model_param, vis_data, triple_data)
print *,' '
print *,'initial sum of sqrd deviations =',real(sumsqr)
print *,' '
print *,'           initial chi squared =',real(chisqrd) 
nlposterior = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, model_param, x, x_pos, x_info, n)
print *,'negative log posterior =',real(nlposterior)

!call minimising algorithm
info = ''
print *,'running simplex minimiser...'
call simplex(info, model_spec, model_param, model_prior, &
             vis_data, triple_data, fit_param, x, x_pos, x_info, sol, ok, n)
if (info/='') goto 91
!now have solution x and x_info which provides minimiser est err, and
!min/max bounds on solution

!update fit_param
do i = 1, n
   fit_param(x_pos(i,1),x_pos(i,2)) = sol(i,1)
end do

!negative logarithm of posterior probability at solution
x = sol(:,1)
nlposterior = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, fit_param, x, x_pos, x_info, n)

flag = 2
if (.not.(ok)) flag = 0

print *,'calculating hessian...'
!find error in x position more rigourously by considering hessian
!hessian matrix Aij has components (d2/(dyi*dyj))P({y}) where yi and yj
!are parameters i and j, P({y}) is the neg log posterior for model with
!parameters {y} = y1, y2, ... yi, ... yj, ... yN
!estimate each element Aij via taylor 2nd order symmetric numeric 
!differentiation method, for off-diagonals this is
![P(i+di,j+dj)+P(i-di,j-dj)-P(i+di,j-dj)-P(i-di,j+dj)]/[4di*dj], for diagonals
![P(i+di)+P(i-di)-2P(i)]/[di^2]
!eta is scaling factor in selecting increment in numeric differentiation,
!increment is eta * upper limit on parameter
eta = 1D-4
!find hessian elements (nb Aij=Aji symmetric)
do i = 1, n
   deltai = eta * x_info(i,3)
   do j = 1, i
      deltaj = eta * x_info(j,3)
      temp_x = sol(:,1)
      if (j == i) then
         !P(i)
         P = posterior(model_spec, model_param, model_prior, vis_data, &
                       triple_data, fit_param, temp_x, x_pos, x_info, n)
         !P(i-di)
         temp_x(i) = sol(i,1)-deltai
         Pl = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, fit_param, temp_x, x_pos, x_info, n)
         !P(i+di)
         temp_x(i) = sol(i,1)+deltai
         Pu = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, fit_param, temp_x, x_pos, x_info, n) 
         !hessian by numeric method
         diff = Pu + Pl - (2D0*P)
         if (diff/=0D0) hes(i,i) = diff/(deltai**2)
      else
         !P(i+di,j+dj)
         temp_x(i) = sol(i,1)+deltai
         temp_x(j) = sol(j,1)+deltaj
         Pu = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, fit_param, temp_x, x_pos, x_info, n)
         !P(i+di,j-dj)
         temp_x(j) = sol(j,1)-deltaj
         Pi = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, fit_param, temp_x, x_pos, x_info, n)   
         !P(i-di,j-dj)
         temp_x(i) = sol(i,1)-deltai
         Pl = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, fit_param, temp_x, x_pos, x_info, n)
         !P(i-di,j+dj)
         temp_x(j) = sol(j,1)+deltaj
         Pj = posterior(model_spec, model_param, model_prior, vis_data, &
                        triple_data, fit_param, temp_x, x_pos, x_info, n)  
         !hessian by numeric method
         diff = Pu + Pl - Pi - Pj
         if (diff/=0D0) hes(i,j) = diff/(4D0*deltai*deltaj)
         hes(j,i) = hes(i,j) !since hessian is symmetric
      end if
   end do
end do

print *,'calculating covariance...'
!calculate covariance matrix is (inverse of the hessian)
cov = hes
call inv_rsis(cov) !inversion of real symmetric indeterminate sqaure matrix

print *,'calculating correlation...'
!calculate correlation matrix
!correlation defined as: corr(x,y) = cov(x,y) / sqroot(var(x)*var(y))
!noting that var(x) = cov(x,x)
do i = 1, n
   do j = 1, i
      if (cov(i,i)*cov(j,j)>0D0) then
         cor(i,j) = cov(i,j)/sqrt(cov(i,i)*cov(j,j))
      else
         flag = 1
      end if
      cor(j,i) = cor(i,j)
   end do
end do

print *, 'estimating errors on parameters...'
!calculate estimated error on solution
!ideal case:
!error on a single parameter x is taken as the sqroot(cov(x,x)) plus the
!minima location error (rms deviation of amoeba vertices) - the error on
!the parameter cannot be lower than this value
!negative cov(x,x) case: estimated error reported as zero (as occurs for
!minimsation failure case too)
do i = 1, n
   if (cov(i,i)>0) then
      sol(i,5) = sqrt(cov(i,i))
      sol(i,6) = sol(i,5) + sol(i,4)
   end if
end do

print *, 'calculating final goodness-of-fit...'
!calculate goodness-of-fit statistic chi squared
call gof(chisqrd, sumsqr, model_spec, fit_param, vis_data, triple_data)

!clean-up and return
200 continue

!deallocate local storage
if (allocated(x)) deallocate(x)
if (allocated(x_info)) deallocate(x_info)
if (allocated(x_pos)) deallocate(x_pos)
if (allocated(temp_x)) deallocate(temp_x)

return

!error trapping 
90 info = 'illegal freedom(s) in model'
goto 200
91 info = 'minimisation routine error: ' // info
goto 200
92 info = 'for vis/nvis data must have gauranteed symmetric model'
goto 200

end subroutine minimiser

!==============================================================================

subroutine simplex(info, model_spec, model_param, model_prior, vis_data, &
                   triple_data, fit_param, x, x_pos, x_info, sol, ok, ndim)

!based on amoeba code (numerical recipies) for downhill simplex algorithm
!in multidimensions (uses amotry code also in numerical recipies)

!subroutine arguments
character(len=128) :: info
character(len=128), dimension(:,:), allocatable :: model_spec
double precision, dimension(:,:), allocatable :: model_param, model_prior 
double precision, dimension(:,:), allocatable :: fit_param, x_info
double precision, dimension(:,:), allocatable :: vis_data, triple_data
double precision, dimension(:,:), allocatable :: sol
double precision, dimension(:), allocatable :: x
logical :: ok
integer, dimension(:,:), allocatable :: x_pos
integer :: ndim

!local variables
integer :: s
integer :: iter, nmax, itmax
double precision :: zeta, eta, ftol
double precision, dimension(:,:), allocatable :: p
double precision, dimension(:), allocatable :: y
integer i, ihi, ilo, inhi, j, m, n
double precision :: rtol, sum, swap, ysave, ytry, min, max
double precision, dimension(:), allocatable :: psum

!routine parameter setting
!nmax limit on dimensions is a carryover from original code
!setting equal to ndim effectively puts onus of ndim into minimiser hands
!??
nmax = ndim
!set fractional tolerance to machine precision, zeta scale factor is
!to prevent rounding errors causing problems (refer to numerical recipies)
zeta = 1D6
ftol = zeta*machine_precision()
print *, 'ftol=', ftol
!initial bounding of minimum accomplished by bounding moving to lower
!and upper limits on each variable. the limit is chosen by moving a distance
!of prior * eta away from the model value. 
!normal set to 5D0 for 5 standard deviation initial bounding
eta = 5D0
!itmax limit of maximum number of function evaluations
itmax = 200000

!allocations
allocate(p(1:ndim+1,1:ndim))
allocate(y(1:ndim+1))
allocate(psum(nmax))

!form initial simplex
!rows of p are the vertices of the starting simplex
!1st vertex is the model point
! ** NO IT WASN'T??
p(1,:) = x - eta*x_info(:,1)

!check and correct for range violations
do s = 1, ndim
   if (p(1,s) < x_info(s,2)) p(1,s) = x_info(s,2)
   if (p(1,s) > x_info(s,3)) p(1,s) = x_info(s,3)
end do
!remaining vertices are model points offset by eta*prior in each 
!variable direction in parameter space
do s = 2, ndim+1
   p(s,:) = x
   p(s,s-1) = p(s,s-1) + eta*x_info(s-1,1)
   !check and correct for range violations
   if (p(s,s-1) < x_info(s-1,2)) p(s,s-1) = x_info(s-1,2)
   if (p(s,s-1) > x_info(s-1,3)) p(s,s-1) = x_info(s-1,3)
end do

!form initial vector with elements equal to function value at
!each of the simplex vertices
do s = 1, ndim+1
   x = p(s,:)
   print *, x
   y(s) = posterior(model_spec, model_param, model_prior, vis_data, &
                    triple_data, fit_param, x, x_pos, x_info, ndim)
end do

!amoeba code
!enter here when starting or have just overall contracted
1 continue
do n = 1, ndim
   sum=0D0 !recompute psum
   do m = 1, ndim+1
      sum = sum + p(m,n)
   end do
   psum(n) = sum
end do

!enter here when have just changed a single point
2 continue
ilo=1
if (y(1) > y(2)) then !determine which point is the highest (worst),
   ihi=1              !next-highest and lowest (best) 
   inhi=2
else
   ihi=2
   inhi=1
end if

do i=1, ndim+1 !by looping over points in the simplex
   if (y(i) < y(ilo)) ilo=i
   if (y(i) > y(ihi)) then
      inhi=ihi
      ihi=i
   else if (y(i) > y(inhi)) then
      if (i /= ihi) inhi=i
   end if
end do
print *, '--', iter, y(:), ihi, inhi, ilo

!compute the fractional range from highest to lowest and return if satisfactory
rtol = 2D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))

if (rtol < ftol) then !if returning, put best point and value in slot 1
   swap = y(1)
   y(1) = y(ilo)
   y(ilo) = swap
   do n=1, ndim
      swap = p(1,n)
      p(1,n) = p(ilo,n)
      p(ilo,n) = swap
   end do
   ok = .true.
   goto 3 !termination
end if

if (iter > itmax) then
   ok = .false.
   goto 3 !termination
end if   

iter = iter + 2

!begin new iteration, first extrapolate by a factor -1 through the face of
!the simplex across from the high point, i.e. reflect the simplex from the
!high point
ytry = extrap(p, y, psum, ndim, ihi, -1D0, nmax, model_spec, model_param, &
              model_prior, vis_data, triple_data, fit_param, x_pos, x_info)
! ** NR has <= here
if (ytry < y(ilo)) then !gives better result than the best point so
                        !try additonal extrapolation by a factor 2
   ytry = extrap(p, y, psum, ndim, ihi, 2D0, nmax, model_spec, model_param, &
                 model_prior, vis_data, triple_data, fit_param, x_pos, x_info)
! ** NR has >= y(inhi) here
else if (ytry > y(ihi)) then !the reflected point is worse than the second
                             !highest so look for intermediate lower point
                             !i.e. do one dimensional contraction
   ysave = y(ihi)
   ytry = extrap(p, y, psum, ndim, ihi, 0.5D0, nmax, model_spec, model_param, &
                 model_prior, vis_data, triple_data, fit_param, x_pos, x_info)
! ** NR has >= here
   if (ytry > ysave) then !cant seem to get rid of that high point, better
                          !contract around the lowest (best) point
      do i=1, ndim+1
         if (i /= ilo) then
            do j=1, ndim
               psum(j) = 0.5D0*(p(i,j)+p(ilo,j))
               p(i,j) = psum(j)
            end do
            y(i) = posterior(model_spec, model_param, model_prior, vis_data, &
                             triple_data, fit_param, psum, x_pos, x_info, ndim)
         end if
      end do
      iter = iter + ndim !keep track of function evaluations
      goto 1 !go back for the test of doneness and the next iteration
   end if
else
   iter = iter - 1 !correct the evaluation count
end if
goto 2

3 continue

!write position to solution matrix (row 1 of p has best point)
sol(:,1) = p(1,:)

!return x_info with info on final solution, - and + tolerances
!all obtained by looking over simplex vertices
do n = 1, ndim !loop over variables
   sum = 0D0
   min = p(1,n)
   max = p(1,n)
   do m = 1, ndim+1 !loop over vertices
      if (p(m,n)<min) min=p(m,n)
      if (p(m,n)>max) max=p(m,n)
      sum = sum + ((p(m,n)-p(1,n))**2D0)
   end do
   sol(n,2) = min
   sol(n,3) = max
   sol(n,4) = sqrt(sum/(dble(ndim))) !mean sqrd deviation 
end do

!clean-up and return
200 continue

!deallocations
if (allocated(p)) deallocate(p)
if (allocated(y)) deallocate(y)
if (allocated(psum)) deallocate(psum)

return

!error trapping
90 info = ''
goto 200

end subroutine simplex

!==============================================================================

subroutine gof(chisqrd, sumsqr, model_spec, fit_param, vis_data, triple_data)

!goodness of fit information calculation
!computes the chisqrd value of a model fit and data set
!chisqrd = sum[i] { (data[i]-theory[i])^2/data_error[i]^2 }
!chisqrd is thus almost identical to likelihood
!could save code here and significantly combine this function
!with the likeihood one but maybe less risky to keep them separated...
!could in future add code here to look up critical chi sqaured values
!for num of data points and num of free parameters

!subroutine arguments
double precision :: chisqrd, sumsqr
character(len=128), dimension(:,:), allocatable :: model_spec
double precision, dimension(:,:), allocatable :: fit_param
double precision, dimension(:,:), allocatable :: vis_data, triple_data

!local variables
integer :: i
double precision :: lambda, u, v, data_vis, data_vis_err, model_vis
double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
double precision :: data_phase_err, model_amp, model_phase
double complex :: vis, vis1, vis2, vis3

chisqrd = 0D0
sumsqr = 0D0

!sum over the visibility data points
do i = 1, size(vis_data,1)
   
   !extract points
   lambda = vis_data(i,1)
   u = vis_data(i,2)
   v = vis_data(i,3)
   data_vis = vis_data(i,4)
   data_vis_err = vis_data(i,5)

   !ignore if -ve or zero error:
   if (data_vis_err>0D0) then
      
      !compute model-predicted visibility amplitude squared
      vis = cmplx_vis(model_spec, fit_param, lambda, u, v)
      model_vis = (modulus(vis))**2D0
      
      !compute contribution to rmsdev and chisqrd
      sumsqr = sumsqr + ((data_vis-model_vis)**2D0)
      chisqrd = chisqrd + &
                (((data_vis-model_vis)/data_vis_err)**2D0)

      print '(a, 3f7.4, 2f9.2)','vis',data_vis,model_vis,data_vis_err, &
           ((data_vis-model_vis)/data_vis_err)**2D0, chisqrd

   end if      

end do

!sum over triple product amplitude and phase data points
do i = 1, size(triple_data,1)

   !extract points
   lambda = triple_data(i,1)
   u1 = triple_data(i,2)
   v1 = triple_data(i,3)
   u2 = triple_data(i,4)
   v2 = triple_data(i,5)
   data_amp = triple_data(i,6)
   data_amp_err = triple_data(i,7)
   data_phase = triple_data(i,8)
   data_phase_err = triple_data(i,9)

   if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
      !compute model-predicted triple amplitude
      vis1 = cmplx_vis(model_spec, fit_param, lambda, u1, v1)
      vis2 = cmplx_vis(model_spec, fit_param, lambda, u2, v2)
      vis3 = cmplx_vis(model_spec, fit_param, lambda, -(u1+u2), -(v1+v2))
      vis = vis1*vis2*vis3
      model_amp = modulus(vis)
      model_phase = rad2deg*argument(vis)
   end if

   !compute rmsdev and chisqrd contributions
   !nb phase calculations modulo 360 degrees
   if (data_amp_err>0D0) then
      sumsqr = sumsqr + ((data_amp-model_amp)**2D0)
      chisqrd = chisqrd + &
                (((data_amp-model_amp)/(data_amp_err))**2D0)
   end if
   if (data_phase_err>0D0) then
      sumsqr = sumsqr + ((modx(data_phase,model_phase))**2D0)
      chisqrd = chisqrd + &
           (modx(data_phase,model_phase)/data_phase_err)**2D0

      print '(a, 3f8.2, 2f9.2)','t_phase',data_phase,model_phase,data_phase_err, &
           (modx(data_phase,model_phase)/data_phase_err)**2D0, chisqrd
   end if 


end do

end subroutine gof

!==============================================================================

function extrap(p, y, psum, ndim, ihi, fac, nmax, model_spec, model_param, &
                model_prior, vis_data, triple_data, fit_param, x_pos, x_info)

!called from amoeba code in subroutine simplex
!identical to amotry code given in numerical recipies
!function extrapolates by a factor fac through the face of the simplex across
!from the high point, tries it, and replaces the high point if the new point
!is better

!function arguments
integer :: ihi, ndim, nmax
double precision :: extrap, fac
double precision, dimension(:,:), allocatable :: p
double precision, dimension(:), allocatable :: psum, y
character(len=128), dimension(:,:), allocatable :: model_spec
double precision, dimension(:,:), allocatable :: model_param, model_prior 
double precision, dimension(:,:), allocatable :: fit_param, x_info
double precision, dimension(:,:), allocatable :: vis_data, triple_data
integer, dimension(:,:), allocatable :: x_pos

!local variables
integer :: j
double precision :: fac1, fac2, ytry
double precision, dimension(:), allocatable :: ptry

!allocations
allocate(ptry(nmax))

!amotry start
fac1 = (1D0-fac)/ndim
fac2 = fac1-fac
do j = 1, ndim
   ptry(j) = (psum(j)*fac1) - (p(ihi,j)*fac2)
end do

!evaluate the function at the highest point, it its better than the highest
!then replace the highest
ytry = posterior(model_spec, model_param, model_prior, vis_data, &
                 triple_data, fit_param, ptry, x_pos, x_info, ndim)
if (ytry < y(ihi)) then
   y(ihi) = ytry
   do j = 1, ndim
      psum(j) = psum(j) - p(ihi,j) + ptry(j)
      p(ihi,j) = ptry(j)
   end do
end if

extrap = ytry

!deallocate
deallocate(ptry)

return

end function extrap

!==============================================================================

function posterior(model_spec, model_param, model_prior, vis_data, &
                   triple_data, fit_param, x, x_pos, x_info, n)

!subroutine arguments
double precision :: posterior
character(len=128), dimension(:,:), allocatable :: model_spec
double precision, dimension(:,:), allocatable :: model_param, model_prior 
double precision, dimension(:,:), allocatable :: fit_param, x_info
double precision, dimension(:,:), allocatable :: vis_data, triple_data
double precision, dimension(:), allocatable :: x
integer, dimension(:,:), allocatable :: x_pos
integer :: n

!local variables
double precision :: lhd, pri
integer :: i
logical :: violation

!check for range violations in current position
!if position is out of range then return very large energy 
violation = .false.
do i = 1, n !over variable parameters
   if (x(i) < x_info(i,2)) then
      violation = .true.
      print *,'out of bounds in posterior()',x(i),x_info(i,2),x_info(i,3)
   end if
   if (x(i) > x_info(i,3)) then
      violation = .true.
      print *,'out of bounds in posterior()',x(i),x_info(i,2),x_info(i,3)
   end if
end do
if (violation) then
   posterior = machine_max()
   return
end if 

!update fit_param with current position x
do i = 1, n
   fit_param(x_pos(i,1),x_pos(i,2)) = x(i)
end do

!form E = neg log posterior = neg log likelihood + neg log prior
lhd = likelihood(model_spec, fit_param, vis_data, triple_data)
pri = prior(model_param, model_prior, x, x_pos)
posterior = lhd + pri

end function posterior

!==============================================================================

function likelihood(model_spec, fit_param, vis_data, triple_data)

!computes the negative log likelihood of data given model

!subroutine arguments
double precision :: likelihood
character(len=128), dimension(:,:), allocatable :: model_spec
double precision, dimension(:,:), allocatable :: fit_param
double precision, dimension(:,:), allocatable :: vis_data, triple_data

!local variables
integer :: i
double precision :: lambda, u, v, data_vis, data_vis_err, model_vis
double precision :: u1, v1, u2, v2, data_amp, data_amp_err, data_phase
double precision :: data_phase_err, model_amp, model_phase
double complex :: vis, vis1, vis2, vis3

likelihood = 0D0

!sum over the visibility data points
do i = 1, size(vis_data,1)
   
   !extract points
   lambda = vis_data(i,1)
   u = vis_data(i,2)
   v = vis_data(i,3)
   data_vis = vis_data(i,4)
   data_vis_err = vis_data(i,5)

   if (data_vis_err>0D0) then
      
      !compute model-predicted visibility amplitude squared
      vis = cmplx_vis(model_spec, fit_param, lambda, u, v)
      model_vis = (modulus(vis)**2D0)
      
      !compute contribution to likelihood
      likelihood = likelihood + &
           (((data_vis-model_vis)**2D0)/(2D0*(data_vis_err**2D0)))

   end if

end do

!sum over triple product amplitude and phase data points
do i = 1, size(triple_data,1)

   !extract points
   lambda = triple_data(i,1)
   u1 = triple_data(i,2)
   v1 = triple_data(i,3)
   u2 = triple_data(i,4)
   v2 = triple_data(i,5)
   data_amp = triple_data(i,6)
   data_amp_err = triple_data(i,7)
   data_phase = triple_data(i,8)
   data_phase_err = triple_data(i,9)

   if ((data_amp_err>0D0).or.(data_phase_err>0D0)) then
      !compute model-predicted triple product
      vis1 = cmplx_vis(model_spec, fit_param, lambda, u1, v1)
      vis2 = cmplx_vis(model_spec, fit_param, lambda, u2, v2)
      vis3 = cmplx_vis(model_spec, fit_param, lambda, -(u1+u2), -(v1+v2))
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
           (((modx(data_phase,model_phase))**2D0)/(2D0*(data_phase_err**2D0)))
   end if
   end do

end function likelihood

!==============================================================================

function prior(model_param, model_prior, x, x_pos)

!function arguments
double precision :: prior
double precision, dimension(:,:), allocatable :: model_param, model_prior
double precision, dimension(:), allocatable :: x
integer, dimension(:,:), allocatable :: x_pos

!local variables
integer :: i
double precision :: value, value0, err

prior = 0D0

!loop over the variable parameters
do i = 1, size(x_pos,1)

   value0 = model_param(x_pos(i,1),x_pos(i,2)) !parameter original value
   err = model_prior(x_pos(i,1),x_pos(i,2)) !parameter prior width
   value = x(i) !parameter current value as it is varied

   !calculate contribution to the prior
   prior = prior + (((value-value0)**2D0)/(2D0*(err**2D0)))

end do

end function prior

!==============================================================================

end module Fit
