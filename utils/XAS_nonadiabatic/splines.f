c
c     spline functions from numerical recpies modified for real(8)
c
c     Mathias Ljungberg 2006-09-21

c----   Subroutine spline is the Numerical Recipes cubic spline
c       routine. It computes the second derivatives at each node
c       for the data points x and y, real vectors of length n.
c       yp1 and ypn are the endpoint slopes of the spline. If
c       they are set to 1.0e30 or greater then a natural spline
c       is formed. MAX is the maximal number of points, can be changed
c
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL(8) yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=100000)
      INTEGER i,k
      REAL(8) p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
 11   continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
 12   continue
      return
      END
       




csubroutine spline(x,y,n,yp1,ypn,y2)
c        parameter (nmax=401)
c        dimension x(n),y(n),y2(n),u(nmax)
c        if (yp1.gt..99e30) then
c           y2(1)=0.
c           u(1)=0.
c        else
cc           y2(1)=-0.5
c           u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
c        endif
c        do 11 i=2,n-1
c           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
c           p=sig*y2(i-1)+2.
c           y2(i)=(sig-1.)/p
c           u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
c     &         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
c  11    continue
c        if (ypn.gt..99e30) then
c           qn=0.
c           un=0.
c        else
c           qn=0.5
c           un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
c        endif
c        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
c        do 12 k=n-1,1,-1
c           y2(k)=y2(k)*y2(k+1)+u(k)
c  12    continue
c        return
c        end



c
c----   Subroutine splint is the Numerical Recipes sister routine
c       for the cubic spline routine, spline. It takes the data
c       points, xa and ya vectors of length n, and the second
c       derivatives computed by spline, y2a, and it computes the
c       y value for the specified x value.
c

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL(8) x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL(8) a,b,h
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      return
      END




c        subroutine splint(xa,ya,y2a,n,x,y)
c        dimension xa(n),ya(n),y2a(n)
c        klo=1
c        khi=n
c  1     if (khi-klo.gt.1) then
c           k=(khi+klo)/2
c           if(xa(k).gt.x)then
c              khi=k
c           else
c              klo=k
c           endif
c           goto 1
c        endif
c        h=xa(khi)-xa(klo)
c        if (h.eq.0.) pause 'bad xa input.'
c        a=(xa(khi)-x)/h
c        b=(x-xa(klo))/h
c        y=a*ya(klo)+b*ya(khi)+
c     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
c        return
c        end



c      SUBROUTINE cspline(x,y,n,yp1,ypn,y2)
c      INTEGER n,NMAX
c      REAL(8) x(n)
c      COMPLEX(16) yp1,ypn,y(n),y2(n)
c      PARAMETER (NMAX=50000)
c      INTEGER i,k
c      COMPLEX(16) p,qn,sig,un,u(NMAX)
c      if (yp1.gt..99d30) then
c        y2(1)=0.d0
c        u(1)=0.d0
c      else
c        y2(1)=-0.5d0
c        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
c      endif
c      do 11 i=2,n-1
c        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
c        p=sig*y2(i-1)+2.d0
c        y2(i)=(sig-1.d0)/p
c        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
c     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
c     *u(i-1))/p
c 11        continue
c      if (ypn.gt..99d30) then
c        qn=0.d0
c        un=0.d0
c      else
c        qn=0.5d0
c        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
c      endif
c      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
c      do 12 k=n-1,1,-1
c         y2(k)=y2(k)*y2(k+1)+u(k)
c 12         continue
c      return
c      END






c bicubic spline functions, 100105

      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      REAL(8) x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=5000)
CU    USES spline                                                                                                                                                                                                 
      INTEGER j,k
      REAL(8) y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
 11             continue
        call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
 12             continue
 13                 continue
      return
      END

 
      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      REAL(8) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=5000)
CU    USES spline,splint                                                                                                                                                                                          
      INTEGER j,k
      REAL(8) y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
 11             continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
 12         continue
      call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END

