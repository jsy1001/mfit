program calc

use maths

implicit none

!variables
double precision :: nu, x, y, res, resx, resy, mu, sigma, G05DDF
double precision, dimension(:,:), allocatable :: A
double precision, dimension(2) :: xx
integer :: yyyy, mm, dd, h, m, s, hah, ham
double precision :: sum, sqrdsum, mean, var, stdev, unbvar, unbstdev
double precision :: sqrdmean, fourmean, sixmean, foursum, sixsum
double complex, dimension(:), allocatable :: coutput
double precision, dimension(:), allocatable :: routput
integer :: i, j, n, p, k
character(len=1) :: choice

do i=1, 1000

   print *,' '
   print *,'----------------------------------------------------------------'
   print *,'Select:'
   print *,' '
   print *,' a: bessel functionof the 1st kind'
   print *,' b: gamma function'
   print *,' c: laguerre polynomial'
   print *,' d: atan2'
   print *,' e: factorial'
   print *,' f: combinatorial'
   print *,' g: bess1'
   print *,' h: modulo'
   print *,' i: real symmetric indefinite sqaure matrix inversion'
   print *,' '
   print *,' p: n random numbers drawn from gaussian(mu,sigma)'
   print *,' q: n random numbers drawn from uniform(0,1)'
   print *,'----------------------------------------------------------------'
1  read *,choice

   print *,' '
   select case (choice)
   
     case('a')
        print *,'bessel function: input order, arg re part, arg im part'
        read *,nu,x,y
        call bessel(nu, x, y, 1, coutput)
        print *,'complex result'
        print *,coutput(1)
        if (allocated(coutput)) deallocate(coutput)
        
     case('b')
        print *,'gamma function : input argument'
        read *,x
        res = gamma(x)
        print *,'result'
        print *,res

     case('c')
        print *,'laguerre polynomial: input order, argument'
        read *,n,x
        call laguerre(n, x, routput)
        print *,'result'
        print *,routput(n)
        if (allocated(routput)) deallocate(routput)
        
     case('d')
        print *,'atan2: input cmplx number re part, im part'
        read *,x,y
        print *,'result'
        print *,atan2(y,x)

     case('e')
        print *,'factorial n!: input n'
        read *,n
        print *,'result'
        print *,fact(n)

     case('f')
        print *,'combinatorial nCp: input n, p'
        read *,n,p
        print *,'result'
        print *,comb(n,p)

     case('g')
        print *,'bess1: input x'
        read *,x
        print *,'result'
        print *,bess1(x)
        
     case('h')
        print *,'a mod b: input a, b'
        read *,x,y
        print *,'result'
        print *,modulo(x,y)

     case('i')
        print *,'real symm indef sq matrix A inv, order nxn: input n'
        read *,n
        allocate(A(n,n))
        do j = 1, n
           do k = 1, j
              print *,'input element for row',j,'col',k
              read *,x
              A(j,k) = x
              A(k,j) = x
           end do
        end do
        print *,' '
        print *,'input matrix:'
        do j = 1, n
           print *,A(j,:)
        end do
        call inv_rsis(A)
        print *,' '
        print *,'output matrix:'
        do j = 1, n
           print *,A(j,:)
        end do
        deallocate(A)
           
     case('p')
        print *,'n random numbers from N(mu,sigma): input n, mu, sigma'
        read *,n,mu,sigma
        print *,'numbers'
        call rnd_initialise()
        do j = 1, n
           x = rnd_gauss(mu,sigma)
           print *,j,' ',x
        end do

     case('q')
        print *,'n random numbers from U(0,1): input n'
        read *,n
        print *,'numbers'
        call rnd_initialise()
        do j = 1, n
           x = rnd_unit()
           print *,j,' ',x
        end do

     case default
        print *,'ILLEGAL CHOICE'
        goto 1
        
   end select

end do

end program calc
