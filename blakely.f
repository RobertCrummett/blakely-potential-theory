      program main
      print *, 'Hello Blakely!'
      end

      subroutine sphere(xq,yq,zq,a,rho,xp,yp,zp,gx,gy,gz)
c
c  Subroutine SPHERE calculates the three components of gravitational
c  attraction at a single point due to a uniform sphere.
c
c  Input parameters:
c    Observation point is (xp,yp,zp), and center of sphere is at
c    (xq,yq,zq). Radius of sphere is a and density is rho. Density
c    in units of kg/(m**3). All distance parameters in units of km.
c
c  Output parameters:
c    Gravitatational components (gx,gy,gz) in units of mGal.
c
      real km2m
      data gamma/6.67e-11/,si2mg/1.e5/,pi/3.14159265/,km2m/1.e3/
      ierror=0
      rx=xp-xq
      ry=yp-yq
      rz=zp-zq
      r=sqrt(rx**2+ry**2+rz**2)
      if(r.eq.0.)pause 'SPHERE: Bad argument detected.'
      r3=r**3
      tmass=4.*pi*rho*(a**3)/3.
      gx=-gamma*tmass*rx/r3
      gy=-gamma*tmass*ry/r3
      gz=-gamma*tmass*rz/r3
      gx=gx*si2mg*km2m
      gy=gy*si2mg*km2m
      gz=gz*si2mg*km2m
      return
      end

      subroutine cylind(xq,zq,a,rho,xp,zp,gx,gz)
c 
c  Subroutine CYLINDer calculates x and z components of gravitational
c  attraction due to a cylinder lying parallel to the y axis.
c 
c  Input parameters:
c    Point of observation is (xp,zp). Axis of cylinder penetrates
c    x,z plan at (xq,zq). Radius of cylindar is a and density is
c    rho. Density in kg/(m**3). All distance parameters in km.
c 
c  Output parameters:
c    Components of gravitational attraction (gx,gz) in mGal.
c 
      real km2m
      data gamma/6.67e-11/,si2mg/1.e5/,pi/3.14159265/,km2m/1.e3/
      rx=xp-xq
      rz=zp-zq
      r2=rx**2+rz**2
      if(r2.eq.0.)pause 'CYLIND: Bad argument detected.'
      tmass=pi*(a**2)*rho
      gx=-2.*gam*tmass*rx/r2
      gz=-2.*gam*tmass*rz/r2
      gx=gx*si2mg*km2m
      gz=gz*si2mg*km2m
      return
      end

      subroutine dipole(xq,yq,zq,a,mi,md,m,xp,yp,zp,bx,by,bz)
c
c  Subroutine DIPOLE computes the three components of magentic
c  induction caused by a uniformly magnetized sphere. x axis
c  is north, z axis is down.
c  
c  Input parameters:
c    Observation point located at (xp,yp,zp). Sphere centered
c    at (xq,yq,zq). Magnetization of sphere defined by
c    intensity m, inclination mi, and declination md. Units
c    of distance irrelevant but must be consistent. All angles
c    in degrees. Intensity of magnetization in A/m. Requires
c    subroutine DIRCOS.
c
c  Output parameters:
c    The three components of magnetic induction (bx,by,bz) in
c    units of nT.
c
      real mi,md,m,mx,my,mz,moment
      data pi/3.14159265/,t2nt/1.e9/,cm/1.e-7/
      call dircos(mi,md,0.,mx,my,mz)
      rx=xp-xq
      ry=yp-yq
      rz=zp-zq
      r2=rx**2+ry**2+rz**2
      r=sqrt(r2)
      if(r.eq.0.)pause 'DIPOLE: Bad argument detected.'
      r5=r**5
      dot=rx*mx+ry*my+rz*mz
      moment=4.*pi*(a**3)*m/3
      bx=cm*moment*(3.*dot*rx-r2*mx)/r5
      by=cm*moment*(3.*dot*ry-r2*my)/r5
      bz=cm*moment*(3.*dot*rz-r2*mz)/r5
      bx=bx*t2nt
      by=by*t2nt
      bz=bz*t2nt
      return
      end
      
      function schmidt(n,m,theta)
c
c  Returns Schmidt normalized associated Legendre polynomial.
c  Requires function fac. Modified from Press et al. (1986)
c 
c  Input parameters:
c    Argument of polynomial is theta, in degrees. Degree and
c    order of polynomial are n and m, respectively. Parameter n
c    must be greater than zero, and m must be greater than or
c    equal to n.
c 
      data d2rad/.017453293/ 
      x=cos(theta*d2rad)
      if(m.lt.0.or.m.gt.n)pause 'Schmidt: Bad argument detected'
      pmm=1.
      if(m.gt.0)then
         somx2=sqrt((1.-x)*(1.+x))
         fact=1.
         do 10 i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.
   10       continue
         end if
      if(n.eq.m)then
         schmidt=pmm
         else
            pmmp1=x*(2*m+1)*pmm
            if(n.eq.m+1)then
               schmidt=pmmp1
               else
                  do 11 nn=m+2,n
                     pnn=(x*(2*nn-1)*pmmp1-(nn+m-1)*pmm)/(nn-m)
                     pmm=pmmp1
                     pmmp1=pnn
   11                continue
                  schmidt=pnn
                  end if
            end if
      if(m.ne.0)then
         xnorm=sqrt(2*fac(n-m)/fac(n+m))
         schmidt=xnorm*schmidt
         end if
      return
      end 

      function fac(n)
c 
c  Function FAC calculates n!
c 
      if(n.lt.0)pause 'FAC: Bad argument detected'
      if(n.eq.0.or.n.eq.1)then
         fac=1
         else
            fac=n
            fac2=fac
   30       fac2=fac2-1.
            fac=fac*fac2
            if(fac2.gt.2)go to 30
            end if
      return
      end 
      
      subroutine gbox(x0,y0,z0,x1,y1,z1,x2,y2,z2,rho,g)
c 
c  Subroutine GBOX computes the vertical attraction of a
c  rectangular prism. Sides of prism are parallel to x,y,z axes,
c  and z axis is vertical down.
c 
c  Input parameters:
c    Observation point is (x0,y0,z0). The prism externs from x1
c    to x2, from y1 to y2, and from z1 to z2 in the x, y, and z
c    directions, respectively. Density of prism is rho. All
c    distance parameters in units of km; rho in units of
c    kg/(m**3).
c 
c  Output parameters:
c    Vertical attraction of gravity g, in mGal.
c 
      real km2m
      dimension x(2),y(2),z(2),isign(2)
      data isign/-1,1/,gamma/6.670e-11/,twopi/6.2831853/,
     &     si2mg/1.e5/,km2m/1.e3/
      x(1)=x0-x1
      y(1)=y0-y1
      z(1)=z0-z1
      x(2)=x0-x2
      y(2)=y0-y2
      z(2)=z0-z2
      sum=0.
      do 1 i=1,2
         do 1 j=1,2
            do 1 k=1,2
               rijk=sqrt(x(i)**2+y(j)*2+z(k)**2)
               ijk=isign(i)*isign(j)*isign(k)
               arg1=atan2(x(i)*y(j), z(k)*rijk)
               if (arg1.lt.0.) arg1=arg1+twopi
               arg2=rijk+y(j)
               arg3=rijk+x(i)
               if(arg2.le.0.)pause 'GBOX: Bad field point'
               if(arg3.le.0.)pause 'GBOX: Bad field point'
               arg2=alog(arg2)
               arg3=alog(arg3)
               sum=sum+ijk*(z(k)*arg1-x(i)*arg2-y(j)*arg3)
    1          continue
      g=rho*gamma*sum*si2mg*km2m
      return
      end
      
      subroutine gpoly(x0,z0,xcorn,zcorn,ncorn,rho,g)
c
c  Subroutine GPOLY computes the vertical attraction of a two-
c  dimensional body with polygonal cross section. Axes are
c  right-handed system with y axis parallel to long direction
c  of body and z axis vertical down.
c
c  Input parameters:
c    Observation point is (x0,z0). Arrays xcorn and zcorn (each
c    of length ncorn) contain the coordinates of the polygon
c    corners, arranged in clockwise order when viewed with x axis
c    to right. Density of body is rho. All distance parameters
c    in units of km; rho in units of kg/(m**3).
c
c  Output parameters:
c    Vertical attraction of gravity g, in mGal.
c
      real km2m
      dimension xcorn(ncorn),zcorn(ncorn)
      data gamma/6.67e-11/,si2mg/1.e5/,km2m/1.e3/
      sum=0.
      do 1 n=1,ncorn
         if (n.eq.ncorn) then
            n2=1
            else
               n2=n+1
               end if
         x1=xcorn(n)-x0
         z1=zcorn(n)-z0
         x2=xcorn(n2)-x0
         z2=zcorn(n2)-z0
         r1sq=x1**2+z1**2
         r2sq=x2**2+z2**2
         if(r1sq.eq.0.)pause 'GPOLY: Field point on corner'
         if(r2sq.eq.0.)pause 'GPOLY: Field point on corner'
         denom=z2-z1
         if(denom.eq.0.)denom=1.e-6
         alpha=(x2-x1)/denom
         beta=(x1*z2-x2*z1)/denom
         factor=beta/(1.+alpha**2)
         term1=0.5*(alog(r2sq)-alog(r1sq))
         term2=atan2(z2,x2)-atan2(z1,x1)
         sum=sum+factor*(term1-alpha*term2)
    1    continue
      g=2.*rho*gam*sum*si2mg*km2m
      return
      end  

      subroutine mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t)
c
c  Subroutine MBOX computes the total-field anomaly of an infinitely
c  extended rectangular prism. Sides of prism are parallel to x,y,z
c  axes, and z is vertical down. Bottom of prism extends to infinity.
c  Two calls to mbox can provide the anomaly of a prism with finite
c  thickness; e.g.,
c
c    call mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t1)
c    call mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t2)
c    t=t1-t2
c
c  Requires subroutine DIRCOS. Method from Bhattacharyya (1964).
c
c  Input parameters:
c    Observation point is (x0,y0,z0). Prism extends from x1 to
c    x2, y1 to y2, and z1 to infinity in x, y, and z directions,
c    respectively. Magnetization defined by inclination mi,
c    declination md, intensity m. Ambient field defined by
c    inclination fi and declintaion fd. x axis has declination
c    theta. Distance units are irrelevant but must be consistent.
c    Angles are in degrees, with inclinations positive below
c    horizontal and declinations positive east of true north.
c    Magnetization in A/m.
c
c  Output parameters:
c    Total-field anomaly t, in nT.
c
      real alpha(2),beta(2),mi,md,m,ma,mb,mc
      data cm/1.e-7/,t2nt/1.e9/
      call dircos(mi,md,theta,ma,mb,mc)
      call dircos(fi,fd,theta,fa,fb,fc)
      fm1=ma*fb+mb*fa
      fm2=ma*fc+mc*fa
      fm3=mb*fc+mc*fb
      fm4=ma*fa
      fm5=mb*fb
      fm6=mc*fc
      alpha(1)=x1-x0
      alpha(2)=x2-x0
      beta(1)=y1-y0
      beta(2)=y2-y0
      h=z1-z0
      t=0.
      hsq=h**2
      do 1 i=1,2
         alphasq=alpha(i)**2
         do 1 j=1,2
            sign=1.
            if(i.ne.j)sign=-1.
            r0sq=alphasq+beta(j)**2+hsq
            r0=sqrt(r0sq)
            r0h=r0*h
            alphabeta=alpha(i)*beta(j)
            arg1=(r0-alpha(i))/(r0+alpha(i))
            arg2=(r0-beta(j))/(r0+beta(j))
            arg3=alphasq+r0h+hsq
            arg4=r0sq+r0h-alphasq
            tlog=fm3*alog(arg1)/2.+fm2*alog(arg2)/2.
     &            -fm1*alog(r0+h)
            tatan=-fm4*atan2(alphabeta,arg3)
     &            -fm5*atan2(alphabeta,arg4)
     &            +fm6*atan2(alphabeta,r0h)
    1 t=t+sign*(tlog+tatan)
      t=t*m*cm*t2nt
      return
      end
      
      subroutine dircos(incl,decl,azim,a,b,c)
c  
c  Subroutine DIRCOS computes direction cosines from inclination
c  and declination.
c  
c  Input parameters:
c    incl:  inclination in degrees positive below horizontal.
c    decl:  declinatino in degrees positive east of true north.
c    azim:  azimuth of x axis in degrees positive east of north.
c
c  Output parameters:
c    a,b,c:   the three direction cosines.
c
      real incl
      data d2rad/.017453293/
      xincl=incl*d2rad
      xdecl=decl*d2rad
      xazim=azim*d2rad
      a=cos(xincl)*cos(xdecl-xazim)
      b=cos(xincl)*sin(xdecl-xazim)
      c=sin(xincl)
      return
      end
      
C subroutine facmag(mx,my,mz,x0,y0,z0,x,y,z,n,fx,fy,fz)
C !
C !  Subroutine FACMAG computes the magnetic field due to surface
C !  charge on a polynomial face. Repeated calls can build the
C !  field of an arbitrary polyhedron. x axis is directed north,
C !  z axis vertical down. Requires subroutines CROSS, ROT, LINE,
C !  and PLANE. Algorithm from Bott (1963).
C !
C !  Input parameters:
C !    Observation point is (x0,y0,z0). Polygon corners defined
C !    by arrays x, y, and z of length n. Magnetization given by
C !    mx,my,mz. Polygon limited to 10 corners. Distance units
C !    are irrelevant but must be consistent; magnetization in A/m.
C !
C !  Output parameters:
C !    Three components of magnetic field (fx,fy,fz), in nT.
C !
C implicit none
C real(4), parameter :: CM=1.e-7,T2NT=1.e9,EPS=1.e-20
C real(4), intent(in) :: mx,my,mz,x0,y0,z0
C integer, intent(in) :: n
C real(4), intent(out) :: fx,fy,fz
C integer :: i,j
C real(4) :: nx,ny,nz,dot,rl,rn,px,py,pz,w,u1,v,w1,rk,us,v2s,v1s, &
C            a2,a1,f2,f1,rho2,rho1,r2,r1,fu2,fu1,fv2,fv1, &
C            fw2,fw1,fu,fv,fw
C real(4), dimension(10) :: x,y,z,u,v2,v1,s,xk,yk,zk,xl,yl,zl
C 
C fx=0.
C fy=0.
C fz=0.
C x(n+1)=x(1)
C y(n+1)=y(1)
C z(n+1)=z(1)
C do i=1,n
C         xl(i)=x(i+1)-x(i)
C         yl(i)=y(i+1)-y(i)
C         zl(i)=z(i+1)-z(i)
C         rl=sqrt(xl(i)**2+yl(i)**2+zl(i)**2)
C         xl(i)=xl(i)/rl
C         yl(i)=yl(i)/rl
C         zl(i)=zl(i)/rl
C enddo
C 
C call cross(xl(2),yl(2),zl(2),xl(1),yl(1),zl(1),nx,ny,nz,rn)
C nx=nx/rn
C ny=ny/rn
C nz=nz/rn
C dot=mx*nx+my*ny+mz*nz
C if (dot.eq.0.) return
C 
C call plane(x0,y0,z0,x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3),px,py,pz,w)
C do i=1,n
C         call rot(x(i),y(i),z(i),x(i+1),y(i+1),z(i+1),nx,ny,nz,px,py,pz,s(i))
C         if (s(i).eq.0.) continue
C         call line(px,py,pz,x(i),y(i),z(i),x(i+1),y(i+1),z(i+1),u1,v,w1,v1(i),v2(i),u(i))
C         rk=sqrt((u1-px)**2+(v-py)**2+(w1-pz)**2)
C         xk(i)=(u1-px)/rk
C         yk(i)=(v-py)/rk
C         zk(i)=(w1-pz)/rk
C enddo
C 
C do j=1,n
C         if (s(j).eq.0.) continue
C         us=u(j)**2
C         v2s=v2(j)**2
C         v1s=v1(j)**2
C         a2=v2(j)/u(j)
C         a1=v1(j)/u(j)
C         f2=sqrt(1.+a2*a2)
C         f1=sqrt(1.+a1*a1)
C         rho2=sqrt(us+v2s)
C         rho1=sqrt(us+v1s)
C         r2=sqrt(us+v2s+w**2)
C         r1=sqrt(us+v1s+w**2)
C         if (w.ne.0.) then
C                 fu2=(a2/f2)*alog((r2+rho2)/abs(w)) &
C                     -.5*alog((r2+v2(j))/(r2-v2(j)))
C                 fu1=(a1/f1)*alog((r1+rho1)/abs(w)) &
C                     -.5*alog((r1+v1(j))/(r1-v1(j)))
C                 fv2=(1./f2)*alog((r2+rho2)/abs(w))
C                 fv1=(1./f1)*alog((r1+rho1)/abs(w))
C                 fw2=atan2((a2*(r2-abs(w))),(r2+a2*a2*abs(w)))
C                 fw1=atan2((a1*(r1-abs(w))),(r1+a1*a1*abs(w)))
C                 fu=dot*(fu2-fu1)
C                 fv=-dot*(fv2-fv1)
C                 fw=(-w*dot/abs(w))*(fw2-fw1)
C         else
C                 fu2=(a2/f2)*(1.+alog((r2+rho2)/EPS)) &
C                     -.5*alog((r2+v2(j))/(r2-v2(j)))
C                 fu1=(a1/f1)*(1.+alog((r1+rho1)/EPS)) &
C                     -.5*alog((r1+v1(j))/(r1-v1(j)))
C                 fv2=(1./f2)*(1.+alog((r2+rho2)/EPS))
C                 fv1=(1./f1)*(1.+alog((r1+rho1)/EPS))
C                 fu=dot*(fu2-fu1)
C                 fv=-dot*(fv2-fv1)
C                 fw=0.
C         endif
C         fx=fx-s(j)*(fu*xk(j)+fv*xl(j)+fw*nx)
C         fy=fy-s(j)*(fu*yk(j)+fv*yl(j)+fw*ny)
C         fz=fz-s(j)*(fu*zk(j)+fv*zl(j)+fw*nz)
C enddo
C fx=fx*CM*T2NT
C fy=fy*CM*T2NT
C fz=fz*CM*T2NT
C end subroutine facmag
C 
C subroutine plane(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z,r)
C !
C !  Subroutine PLANE computes the intersection (x,y,z) of a plane
C !  and a perpendicular line. The plane is defined by the three points
C !  (x1,y1,z1), (x2,y2,z2), and (x3,y3,z3). The line passes through
C !  (x0,y0,z0). Computation is done by transformation and inverse
C !  transformation of coordinates systems.
C !
C implicit none
C real(4), intent(in) :: x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3
C real(4), intent(out) :: x,y,z,r
C real(4) :: x2n,y2n,z2n,x0n,y0n,z0n,x3n,y3n,z3n,a,t11,t12,t13, &
C            t21,t22,t23,t31,t32,t33,tx0,tz0,cx,cy,cz,c,dx,dy,dz,d
C 
C x2n=x2-x1
C y2n=y2-y1
C z2n=z2-z1
C x0n=x0-x1
C y0n=y0-y1
C z0n=z0-z1
C x3n=x3-x1
C y3n=y3-y1
C z3n=z3-z1
C call cross(x3n,y3n,z3n,x2n,y2n,z2n,cx,cy,cz,c)
C call cross(x2n,y2n,z2n,cx,cy,cz,dx,dy,dz,d)
C a=sqrt(x2n**2+y2n**2+z2n**2)
C t11=x2n/a
C t12=y2n/a
C t13=z2n/a
C t21=cx/c
C t22=cy/c
C t23=cz/c
C t31=dx/d
C t32=dy/d
C t33=dz/d
C tx0=t11*x0n+t12*y0n+t13*z0n
C tz0=t31*x0n+t32*y0n+t33*z0n
C r=t21*x0n+t22*y0n+t23*z0n
C x=t11*tx0+t31*tz0
C y=t12*tx0+t32*tz0
C z=t13*tx0+t33*tz0
C x=x+x1
C y=y+y1
C z=z+z1
C end subroutine plane
C 
C subroutine line(x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z,v1,v2,r)
C !
C !  Subroutine LINE determines the intersection (x,y,z) of two
C !  lines. First line is defined by points (x1,y1,z1) and
C !  (x2,y2,z2). Second line is perpendicular to the first and
C !  passes through point (x0,y0,z0). Distance between (x,y,z)
C !  and (x0,y0,z0) is returned as r. Computation is done by a
C !  transformation of coordinate systems.
C !
C implicit none
C real(4), intent(in) :: x0,y0,z0,x1,y1,z1,x2,y2,z2
C real(4), intent(out) :: x,y,z,v1,v2,r
C real(4) :: tx0,ty0,tz0,tx2,ty2,tz2,a,cx,cy,cz,c,dx,dy,dz,d, &
C            tt11,tt12,tt13,tt21,tt22,tt23,tt31,tt32,tt33,u0
C 
C tx0=x0-x1
C ty0=y0-y1
C tz0=z0-z1
C tx2=x2-x1
C ty2=y2-y1
C tz2=z2-z1
C a=sqrt(tx2**2+ty2**2+tz2**2)
C call cross(tx2,ty2,tz2,tx0,ty0,tz0,cx,cy,cz,c)
C call cross(cx,cy,cz,tx2,ty2,tz2,dx,dy,dz,d)
C tt11=tx2/a
C tt12=ty2/a
C tt13=tz2/a
C tt21=dx/d
C tt22=dy/d
C tt23=dz/d
C tt31=cx/c
C tt32=cy/c
C tt33=cz/c
C u0=tt11*tx0+tt12*ty0+tt13*tz0
C r=tt21*tx0+tt22*ty0+tt23*tz0
C x=tt11*u0+x1
C y=tt12*u0+y1
C z=tt13*u0+z1
C v1=-u0
C v2=a-u0
C end subroutine line
C 
C subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz,r)
C !
C ! Subroutine CROSS computes the vector product of two vectors: i.e.,
C !
C !               (cx,cy,cz) = (ax, ay, az) X (bx, by, bz)
C !
C implicit none
C real(4), intent(in) :: ax,ay,az,bx,by,bz
C real(4), intent(out) :: cx,cy,cz,r
C cx=ay*bz-az*by
C cy=az*bx-ax*bz
C cz=ax*by-ay*bx
C r=sqrt(cx**2+cy**2+cz**2)
C end subroutine cross
C 
C subroutine rot(ax,ay,az,bx,by,bz,nx,ny,nz,px,py,pz,s)
C !
C !  Subroutine ROT finds the sense of rotation of the vector
C !  from (ax,ay,az) to (bx,by,bz) with respect to a second
C !  vector through point (px,py,pz). The second vector has
C !  components given by (nx,ny,nz). Returned parameter s is
C !  1 if anticlockwise, -1 if clockwise, or 0 if collinear.
C !
C implicit none
C real(4), intent(in) :: ax,ay,az,bx,by,bz,nx,ny,nz,px,py,pz
C real(4), intent(out) :: s
C real(4) :: x,y,z,cx,cy,cz,c,u,v,w,d
C x=bx-ax
C y=by-ay
C z=bz-az
C call cross(nx,ny,nz,x,y,z,cx,cy,cz,c)
C u=px-ax
C v=py-ay
C w=pz-az
C d=u*cx+v*cy+w*cz
C if (d.lt.0.) then
C         s=1.
C else if (d.eq.0.) then
C         s=0.
C else
C         s=-1.
C endif
C end subroutine rot
C 
C subroutine ribbon(x0,z0,x1,z1,x2,z2,mx,mz,fx,fz,ier)
C !
C !  Subroutine RIBBON computes the x and z components of magnetic
C !  induction due to a single side of a two-dimensional prism with
C !  polygonal cross section. The prism is assumed infinitely
C !  extended parallel to the y axis; z axis is vertical down.
C !
C !  Input parameters:
C !    Observation point is (x0,z0). Coordinates (x1,z1) and
C !    (x2,z2) are two consecutive corners of the polygon taken in
C !    clockwise order around the polygon as viewed with the a
C !    axis to the right. The x and z components of magnetization
C !    mx and mz. Distance units are irrelevant but must be
C !    consistent; magnetization in A/m.
C !
C !  Output paramters:
C !    Components of magnetic field fx and fz, in nT.
C !    Errors are recorded by ier:
C !      ier=0, no errors;
C !      ier=1, two corners are too close (no calculation);
C !      ier=2, field point too close to corner (calculation
C !             continues).
C !
C implicit none
C real(4), parameter :: PI=3.14159265,SMALL=1.e-18,CM=1.e-7,T2NT=1.e9
C real(4), intent(in) :: x0,z0,x1,z1,x2,z2,mx,mz
C real(4), intent(out) :: fx,fz
C integer, intent(out) :: ier
C real(4) :: sx,sz,s,qs,rx1,rz1,rx2,rz2,r1,r2,theta1, &
C            theta2,angle,flog,factor
C ier=0
C sx=x2-x1
C sz=z2-z1
C s=sqrt(sx**2+sz**2)
C 
C !
C !  -- If two corners are too close, return
C !
C if (s.lt.SMALL) then
C         ier=1
C         return
C endif
C 
C sx=sx/s
C sz=sz/s
C qs=mx*sz-mz*sx
C rx1=x1-x0
C rz1=z1-z0
C rx2=x2-x0
C rz2=z2-z0
C 
C !
C !  -- If field point is too near a corner, signal error
C !
C if (abs(rx1).lt.SMALL.and.abs(rz1).lt.SMALL) ier=2
C if (abs(rx2).lt.SMALL.and.abs(rz2).lt.SMALL) ier=2
C if (ier.eq.2) then
C         rx1=SMALL
C         rz1=SMALL
C         rx2=SMALL
C         rz2=SMALL
C endif
C r1=sqrt(rx1**2+rz1**2)
C r2=sqrt(rx2**2+rz2**2)
C theta1=atan2(rz1,rx1)
C theta2=atan2(rz2,rx2)
C angle=theta1-theta2
C if (angle.gt.PI) angle=angle-2.*PI
C if (angle.lt.-PI) angle=angle+2.*PI
C 
C !
C !  -- If field point is too close to side, signal error
C !
C if (abs(angle).gt.(.995*PI)) ier=2
C flog=alog(r2)-alog(r1)
C factor=-2.*CM*qs*T2NT
C fx=factor*(sx*flog-sz*angle)
C fz=factor*(sz*flog+sx*angle)
C end subroutine ribbon
C 
C subroutine fork(lx,cx,signi)
C !
C !  Subroutine FORK calculates the Fourier transform of a one-
C !  dimesional array. Algorithm from Claerbout (1976).
C !
C !  Input/output parameters:
C !    Complex array cx of length lx is the input array. Upon
C !    return, cx contains the transformed array. Length of
C !    array must be a power of 2. If signi=-1., then the forward
C !    calculation is performed; if signi=1., the inverse transform
C !    is performed.
C !
C integer, intent(in) :: lx
C complex(4), dimension(lx) :: cx
C integer, intent(out) :: signi
C complex(4) :: carg,cw,ctemp
C real(4) :: sc
C integer :: i,istep,ipl,j=1,l,m
C 
C sc=sqrt(1./lx)
C do i=1,lx
C         if (i.gt.j) goto 2
C         ctemp=cx(j)*sc
C         cx(j)=cx(i)*sc
C         cx(i)=ctemp
C 2       m = lx/2
C 3       if (j.le.m) continue
C         j=j-m
C         m=m/2
C         if (m.ge.1) goto 3
C         j=j+m
C enddo
C 
C l=1
C 6 istep=2*l
C do m=1,l
C         carg=(0.,1.)*(3.14159265*signi*(m-1))/l
C         cw=cexp(carg)
C         do i=m,lx,istep
C                 ipl=i+l
C                 ctemp=cw*cx(ipl)
C                 cx(ipl)=cx(i)-ctemp
C                 cx(i)=cx(i)+ctemp
C         enddo
C enddo
C l=istep
C if (l.lt.lx) goto 6
C end subroutine fork
