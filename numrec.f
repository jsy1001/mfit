c     routines from Numerical Recipes in FORTRAN, Second Edition

      subroutine locate(xx,n,x,j)
      implicit none
      integer j,n
      real x,xx(n)
      integer jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end

      subroutine four1(data,nn,isign)
      implicit none
      integer isign,nn
      real data(2*nn)
      integer i,istep,j,m,mmax,n
      real tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
         if(j.gt.i)then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
         m=n/2
 1       if ((m.ge.2).and.(j.gt.m)) then
            j=j-m
            m=m/2
            goto 1
         endif
         j=j+m
 11   continue
      mmax=2
 2    if (n.gt.mmax) then
         istep=2*mmax
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.d0*sin(0.5d0*theta)**2
         wpi=sin(theta)
         wr=1.0d0
         wi=0.0d0
         do 13 m=1,mmax,2
            do 12 i=m,n,istep
               j=i+mmax
               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
 12         continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
 13      continue
         mmax=istep
         goto 2
      endif
      return
      end

      subroutine realft(data,n,isign)
      implicit none
      integer isign,n
      real data(n)
c     uses four1
      integer i,i1,i2,i3,i4,n2p3
      real c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
         c2=-0.5
         call four1(data,n/2,+1)
      else
         c2=0.5
         theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
         i1=2*i-1
         i2=i1+1
         i3=n2p3-i2
         i4=i3+1
         wrs=sngl(wr)
         wis=sngl(wi)
         h1r=c1*(data(i1)+data(i3))
         h1i=c1*(data(i2)-data(i4))
         h2r=-c2*(data(i2)+data(i4))
         h2i=c2*(data(i1)-data(i3))
         data(i1)=h1r+wrs*h2r-wis*h2i
         data(i2)=h1i+wrs*h2i+wis*h2r
         data(i3)=h1r-wrs*h2r+wis*h2i
         data(i4)=-h1i+wrs*h2i+wis*h2r
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
 11   continue
      if (isign.eq.1) then
         h1r=data(1)
         data(1)=h1r+data(2)
         data(2)=h1r-data(2)
      else
         h1r=data(1)
         data(1)=c1*(h1r+data(2))
         data(2)=c1*(h1r-data(2))
         call four1(data,n/2,-1)
      endif
      return
      end

      subroutine cosft1(y,n)
      implicit none
      integer n
      real y(n+1)
c     uses realft
      integer j
      real sum,y1,y2
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/n
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      sum=0.5*(y(1)-y(n+1))
      y(1)=0.5*(y(1)+y(n+1))
      do 11 j=1,n/2-1
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         y1=0.5*(y(j+1)+y(n-j+1))
         y2=(y(j+1)-y(n-j+1))
         y(j+1)=y1-wi*y2
         y(n-j+1)=y1+wi*y2
         sum=sum+wr*y2
 11   continue
      call realft(y,n,+1)
      y(n+1)=y(2)
      y(2)=sum
      do 12 j=4,n,2
         sum=sum+y(j)
         y(j)=sum
 12   continue
      return
      end

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if (fb.gt.fa) then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
      
      FUNCTION golden(ax,bx,cx,f,tol,xmin)
      REAL golden,ax,bx,cx,tol,xmin,f,R,C
      EXTERNAL f
      PARAMETER (R=.61803399,C=1.-R)
      REAL f1,f2,x0,x1,x2,x3
      integer nits
      x0=ax
      x3=cx

      nits=0

      if(abs(cx-bx).gt.abs(bx-ax))then
        x1=bx
        x2=bx+C*(cx-bx)
      else
        x2=bx
        x1=bx-C*(bx-ax)
      endif

      f1=f(x1)
      f2=f(x2)

1     if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then

         nits=nits+1

         if(nits .gt.500) then
            write (6,*) '**********************************************'
            write (6,*) 'Poor starting values. Too many iter. necessary'
            write (6,*) 'Unable to minimise'
            write (6,*) 'Please try again with different guesses'
            write (6,*) '**********************************************'
            return
         endif

        if(f2.lt.f1)then
          x0=x1
          x1=x2
          x2=R*x1+C*x3
          f1=f2
          f2=f(x2)
        else
          x3=x2
          x2=x1
          x1=R*x2+C*x0
          f2=f1
          f1=f(x1)
        endif
      goto 1
      endif

      if(f1.lt.f2)then
        golden=f1
        xmin=x1
      else
        golden=f2
        xmin=x2
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
