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

program calc

use maths

implicit none

!variables
double precision :: nu, x, y, res
double precision, dimension(:,:), allocatable :: A
double precision, dimension(:), allocatable :: routput
integer :: i, j, n, p, k
character(len=1) :: choice

do i=1, 1000

   print *,' '
   print *,'----------------------------------------------------------------'
   print *,'Select:'
   print *,' '
   print *,' a: bessel function of the 1st kind'
   print *,' b: gamma function'
   print *,' c: laguerre polynomial'
   print *,' d: atan2'
   print *,' e: factorial'
   print *,' f: combinatorial'
   print *,' g: bess1'
   print *,' h: modulo'
   print *,' i: (square) matrix inversion'
   print *,' j: machine-dependent constants'
   print *,'----------------------------------------------------------------'
1  read *,choice

   print *,' '
   select case (choice)
   
     case('a')
        print *,'bessel function: input order, real argument'
        read *,nu,x
        call bessel(nu-floor(nu), x, 1+floor(nu), routput)
        print *,'result'
        print *,routput(1+floor(nu))
        if (allocated(routput)) deallocate(routput)
        
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
        print *,'(square) matrix A inv, order nxn: input n'
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
        call inv_mat(A)
        print *,' '
        print *,'output matrix:'
        do j = 1, n
           print *,A(j,:)
        end do
        deallocate(A)

     case('j')
        print *, 'machine_max', machine_max()
        print *, 'machine_min', machine_min()
        print *, 'machine_precision', machine_precision()

     case default
        print *,'ILLEGAL CHOICE'
        goto 1
        
   end select

end do

end program calc
