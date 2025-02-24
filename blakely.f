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
      
      subroutine facmag(mx,my,mz,x0,y0,z0,x,y,z,n,fx,fy,fz)
c
c  Subroutine FACMAG computes the magnetic field due to surface
c  charge on a polynomial face. Repeated calls can build the
c  field of an arbitrary polyhedron. x axis is directed north,
c  z axis vertical down. Requires subroutines CROSS, ROT, LINE,
c  and PLANE. Algorithm from Bott (1963).
c 
c  Input parameters:
c    Observation point is (x0,y0,z0). Polygon corners defined
c    by arrays x, y, and z of length n. Magnetization given by
c    mx,my,mz. Polygon limited to 10 corners. Distance units
c    are irrelevant but must be consistent; magnetization in A/m.
c 
c  Output parameters:
c    Three components of magnetic field (fx,fy,fz), in nT.
c 
      real mx,my,mz,nx,ny,nz
      dimension u(10),v2(10),v1(10),s(10),xk(10),yk(10),
     &       zk(10),xl(10),yl(10),zl(10),x(10),y(10),z(10)
      data cm/1.e-7/,t2nt/1.e9/,epsilon/1.e-20/
      fx=0.
      fy=0.
      fz=0.
      x(n+1)=x(1)
      y(n+1)=y(1)
      z(n+1)=z(1)
      do 1 i=1,n
         xl(i)=x(i+1)-x(i)
         yl(i)=y(i+1)-y(i)
         zl(i)=z(i+1)-z(i)
         rl=sqrt(xl(i)**2+yl(i)**2+zl(i)**2)
         xl(i)=xl(i)/rl
         yl(i)=yl(i)/rl
    1    zl(i)=zl(i)/rl
      call cross(xl(2),yl(2),zl(2),xl(1),yl(1),zl(1),nx,ny,nz,rn)
      nx=nx/rn
      ny=ny/rn
      nz=nz/rn
      dot=mx*nx+my*ny+mz*nz
      if(dot.eq.0.)return
      call plane(x0,y0,z0,x(1),y(1),z(1),x(2),y(2),z(2),x(3),
     &           y(3),z(3),px,py,pz,w)
      do 2 i=1,n
         call rot(x(i),y(i),z(i),x(i+1),y(i+1),z(i+1),nx,ny,nz,
     &            px,py,pz,s(i))
         if(s(i).eq.0.)go to 2
         call line(px,py,pz,x(i),y(i),z(i),x(i+1),y(i+1),z(i+1),
     &             u1,v,w1,v1(i),v2(i),u(i))
         rk=sqrt((u1-px)**2+(v-py)**2+(w1-pz)**2)
         xk(i)=(u1-px)/rk
         yk(i)=(v-py)/rk
    2    zk(i)=(w1-pz)/rk
      do 3 j=1,n
         if(s(j).eq.0.)go to 3
         us=u(j)**2
         v2s=v2(j)**2
         v1s=v1(j)**2
         a2=v2(j)/u(j)
         a1=v1(j)/u(j)
         f2=sqrt(1.+a2*a2)
         f1=sqrt(1.+a1*a1)
         rho2=sqrt(us+v2s)
         rho1=sqrt(us+v1s)
         r2=sqrt(us+v2s+w**2)
         r1=sqrt(us+v1s+w**2)
         if(w.ne.0.)then
            fu2=(a2/f2)*alog((r2+rho2)/abs(w))
     &          -.5*alog((r2+v2(j))/(r2-v2(j)))
            fu1=(a1/f1)*alog((r1+rho1)/abs(w))-
     &          .5*alog((r1+v1(j))/(r1-v1(j)))
            fv2=(1./f2)*alog((r2+rho2)/abs(w))
            fv1=(1./f1)*alog((r1+rho1)/abs(w))
            fw2=atan2((a2*(r2-abs(w))),(r2+a2*a2*abs(w)))
            fw1=atan2((a1*(r1-abs(w))),(r1+a1*a1*abs(w)))
            fu=dot*(fu2-fu1)
            fv=-dot*(fv2-fv1)
            fw=(-w*dot/abs(w))*(fw2-fw1)
            else
               fu2=(a2/f2)*(1.+alog((r2+rho2)/epsilon))-
     &             .5*alog((r2+v2(j))/(r2-v2(j)))
               fu1=(a1/f1)*(1.+alog((r1+rho1)/epsilon))-
     &             .5*alog((r1+v1(j))/(r1-v1(j)))
               fv2=(1./f2)*(1.+alog((r2+rho2)/epsilon))
               fv1=(1./f1)*(1.+alog((r1+rho1)/epsilon))
               fu=dot*(fu2-fu1)
               fv=-dot*(fv2-fv1)
               fw=0.
               end if
      fx=fx-s(j)*(fu*xk(j)+fv*xl(j)+fw*nx)
      fy=fy-s(j)*(fu*yk(j)+fv*yl(j)+fw*ny)
      fz=fz-s(j)*(fu*zk(j)+fv*zl(j)+fw*nz)
    3 continue
      fx=fx*cm*t2nt
      fy=fy*cm*t2nt
      fz=fz*cm*t2nt
      return
      end

      subroutine plane(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z,r)
c
c  Subroutine PLANE computes the intersection (x,y,z) of a plane
c  and a perpendicular line. The plane is defined by the three points
c  (x1,y1,z1), (x2,y2,z2), and (x3,y3,z3). The line passes through
c  (x0,y0,z0). Computation is done by transformation and inverse
c  transformation of coordinates systems.
c 
      x2n=x2-x1
      y2n=y2-y1
      z2n=z2-z1
      x0n=x0-x1
      y0n=y0-y1
      z0n=z0-z1
      x3n=x3-x1
      y3n=y3-y1
      z3n=z3-z1
      call cross(x3n,y3n,z3n,x2n,y2n,z2n,cx,cy,cz,c)
      call cross(x2n,y2n,z2n,cx,cy,cz,dx,dy,dz,d)
      a=sqrt(x2n**2+y2n**2+z2n**2)
      t11=x2n/a
      t12=y2n/a
      t13=z2n/a
      t21=cx/c
      t22=cy/c
      t23=cz/c
      t31=dx/d
      t32=dy/d
      t33=dz/d
      tx0=t11*x0n+t12*y0n+t13*z0n
      tz0=t31*x0n+t32*y0n+t33*z0n
      r=t21*x0n+t22*y0n+t23*z0n
      x=t11*tx0+t31*tz0
      y=t12*tx0+t32*tz0
      z=t13*tx0+t33*tz0
      x=x+x1
      y=y+y1
      z=z+z1
      return
      end

      subroutine line(x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z,v1,v2,r)
c 
c  Subroutine LINE determines the intersection (x,y,z) of two
c  lines. First line is defined by points (x1,y1,z1) and
c  (x2,y2,z2). Second line is perpendicular to the first and
c  passes through point (x0,y0,z0). Distance between (x,y,z)
c  and (x0,y0,z0) is returned as r. Computation is done by a
c  transformation of coordinate systems.
c 
      tx0=x0-x1
      ty0=y0-y1
      tz0=z0-z1
      tx2=x2-x1
      ty2=y2-y1
      tz2=z2-z1
      a=sqrt(tx2**2+ty2**2+tz2**2)
      call cross(tx2,ty2,tz2,tx0,ty0,tz0,cx,cy,cz,c)
      call cross(cx,cy,cz,tx2,ty2,tz2,dx,dy,dz,d)
      tt11=tx2/a
      tt12=ty2/a
      tt13=tz2/a
      tt21=dx/d
      tt22=dy/d
      tt23=dz/d
      tt31=cx/c
      tt32=cy/c
      tt33=cz/c
      u0=tt11*tx0+tt12*ty0+tt13*tz0
      r=tt21*tx0+tt22*ty0+tt23*tz0
      x=tt11*u0+x1
      y=tt12*u0+y1
      z=tt13*u0+z1
      v1=-u0
      v2=a-u0
      return
      end

      subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz,r)
c 
c Subroutine CROSS computes the vector product of two vectors: i.e.,
c 
c               (cx,cy,cz) = (ax,ay,az) X (bx,by,bz)
c 
      cx=ay*bz-az*by
      cy=az*bx-ax*bz
      cz=ax*by-ay*bx
      r=sqrt(cx**2+cy**2+cz**2)
      return
      end

      subroutine rot(ax,ay,az,bx,by,bz,nx,ny,nz,px,py,pz,s)
c 
c  Subroutine ROT finds the sense of rotation of the vector
c  from (ax,ay,az) to (bx,by,bz) with respect to a second
c  vector through point (px,py,pz). The second vector has
c  components given by (nx,ny,nz). Returned parameter s is
c  1 if anticlockwise, -1 if clockwise, or 0 if collinear.
c 
      real nx,ny,nz
      x=bx-ax
      y=by-ay
      z=bz-az
      call cross(nx,ny,nz,x,y,z,cx,cy,cz,c)
      u=px-ax
      v=py-ay
      w=pz-az
      d=u*cx+v*cy+w*cz
      if(d)2,3,4
    2 s=1.
      go to 1
    3 s=0.
      go to 1
    4 s=-1.
    1 continue
      return
      end

      subroutine ribbon(x0,z0,x1,z1,x2,z2,mx,mz,fx,fz,ier)
c 
c  Subroutine RIBBON computes the x and z components of magnetic
c  induction due to a single side of a two-dimensional prism with
c  polygonal cross section. The prism is assumed infinitely
c  extended parallel to the y axis; z axis is vertical down.
c 
c  Input parameters:
c    Observation point is (x0,z0). Coordinates (x1,z1) and
c    (x2,z2) are two consecutive corners of the polygon taken in
c    clockwise order around the polygon as viewed with the a
c    axis to the right. The x and z components of magnetization
c    mx and mz. Distance units are irrelevant but must be
c    consistent; magnetization in A/m.
c 
c  Output paramters:
c    Components of magnetic field fx and fz, in nT.
c    Errors are recorded by ier:
c        ier=0, no errors;
c        ier=1, two corners are too close (no calculation);
c        ier=2, field point too close to corner (calculation
c               continues).
c 
      real mx,mz
      data pi/3.14159265/,small/1.e-18/,cm/1.e-7/,t2nt/1.e9/
      ier=0
      sx=x2-x1
      sz=z2-z1
      s=sqrt(sx**2+sz**2)
c
c  -- If two corners are too close, return
c
      if (s.lt.SMALL)then
         ier=1
         return
         end if
      sx=sx/s
      sz=sz/s
      qs=mx*sz-mz*sx
      rx1=x1-x0
      rz1=z1-z0
      rx2=x2-x0
      rz2=z2-z0
c
c  -- If field point is too near a corner, signal error
c
      if(abs(rx1).lt.small.and.abs(rz1).lt.small)ier=2
      if(abs(rx2).lt.small.and.abs(rz2).lt.small)ier=2
      if(ier.eq.2)then
         rx1=small
         rz1=small
         rx2=small
         rz2=small
         end if
      r1=sqrt(rx1**2+rz1**2)
      r2=sqrt(rx2**2+rz2**2)
      theta1=atan2(rz1,rx1)
      theta2=atan2(rz2,rx2)
      angle=theta1-theta2
      if (angle.gt.pi)angle=angle-2.*pi
      if (angle.lt.-pi)angle=angle+2.*pi
c
c  -- If field point is too close to side, signal error
c
      if (abs(angle).gt.(.995*pi))ier=2
      flog=alog(r2)-alog(r1)
      factor=-2.*cm*qs*t2nt
      fx=factor*(sx*flog-sz*angle)
      fz=factor*(sz*flog+sx*angle)
      return
      end

      subroutine fork(lx,cx,signi)
c
c  Subroutine FORK calculates the Fourier transform of a one-
c  dimesional array. Algorithm from Claerbout (1976).
c 
c  Input/output parameters:
c    Complex array cx of length lx is the input array. Upon
c    return, cx contains the transformed array. Length of
c    array must be a power of 2. If signi=-1., then the forward
c    calculation is performed; if signi=1., the inverse transform
c    is performed.
c 
      complex cx(2050),carg,cexp,cw,ctemp
      j=1
      sc=sqrt(1./lx)
      do 5 i=1,lx
         if(i.gt.j)go to 2
            ctemp=cx(j)*sc
            cx(j)=cx(i)*sc
            cx(i)=ctemp
    2    m=lx/2
    3    if(j.le.m)go to 5
            j=j-m
            m=m/2
            if(m.ge.1)go to 3
    5    j=j+m
      l=1
    6 istep=2*l
      do 8 m=1,l
         carg=(0.,1.)*(3.14159265*signi*(m-1))/l
         cw=cexp(carg)
         do 8 i=m,lx,istep
            ipl=i+l
            ctemp=cw*cx(ipl)
            cx(ipl)=cx(i)-ctemp
    8       cx(i)=cx(i)+ctemp
      l=istep
      if(l.lt.lx)go to 6
      return
      end

      subroutine fourn(data,nn,ndim,isign)
c
c  Replaces DATA by its NDIM-dimensional discrete Fourier transform,
c  is ISIGN is input as 1.  NN is an integer array of length NDIM,
c  containing the lengths of each dimension (number of complex values),
c  which must all be powers of 2.  DATA is a real array of length twice
c  the product of these lengths, in which the data are stored as in a
c  multidimensional complex Fortran array.  If ISIGN is input as -1,
c  DATA is replaced by its inverse transform times the product of the
c  lengths of all dimensions.  From Press et al. (1986, pp. 451-3).
c
      real*8 wr,wi,wpr,wpi,wtemp,theta
      dimension nn(ndim),data(*)
      ntot=1
      do 11 iidim=1,ndim
        ntot=ntot*nn(iidim)
11    continue
      nprev=1
      do 18 iidim=1,ndim
        n=nn(iidim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          go to 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp1
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        go to 2
        endif
        nprev=n*nprev
18    continue
      return
      end

      subroutine glayer(rho,nx,ny,dx,dy,z1,z2,store)
c
c  Subroutine GLAYER calculates the vertical gravitational
c  attraction on a two-dimensional grid caused by a two-
c  dimensional density confined to a horizontal layer.  The
c  following steps are involved:  (1) Fourier transform the
c  density, (2) multiply by the earth filter, and (3) inverse
c  Fourier transform the product.  Density is specified on a
c  rectangular grid with x and y axes directed north and east,
c  respectively.  z axis is down.  Requires subroutines
c  FOURN, KVALUE, and GFILT.
c
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c         (NOTE: Both nx and ny must be powers of two.)
c    rho - a one-dimensional real array containing the
c         two-dimensional density, in kg/(m**3).  Elements
c         should be in order of west to east, then south to
c         north (i.e., element 1 is the southwest corner,
c         element ny is the southeast corner, element
c         (nx-1)*ny+1 is the northwest corner, and element ny*nx
c         is the northeast corner).
c    store - a one-dimensional real array used internally.
c         It should be dimensioned at least 2*nx*ny.
c    dx - sample interval in the x direction, in km.
c    dy - sample interval in the y direction, in km.
c    z1 - depth to top of layer, in km.  Must be > 0.
c    z2 - depth to bottom of layer, in km.  Must be > z1.
c
c  Output parameters:
c    rho - upon output, rho will contain the gravity anomaly,
c         in mGal, with same orientation as before.
c
      complex crho,cmplx
      real kx,ky,km2m
      dimension rho(nx*ny),store(2*nx*ny),nn(2)
      data pi/3.14159265/,si2mg/1.e5/,km2m/1.e3/
      index(i,j,ncol)=(j-1)*ncol+i
      nn(1)=ny
      nn(2)=nx
      ndim=2
      dkx=2.*pi/(nx*dx)
      dky=2.*pi/(ny*dy)
      do 10 j=1,nx
         do 10 i=1,ny
            ij=index(i,j,ny)
            store(2*ij-1)=rho(ij)
   10       store(2*ij)=0.
      call fourn(store,nn,ndim,-1)
      do 20 j=1,nx
         do 20 i=1,ny
            ij=index(i,j,ny)
            call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
            crho=cmplx(store(2*ij-1),store(2*ij))
            crho=crho*gfilt(kx,ky,z1,z2)
            store(2*ij-1)=real(crho)
   20       store(2*ij)=aimag(crho)
      call fourn(store,nn,ndim,+1)
      do 30 j=1,nx
         do 30 i=1,ny
            ij=index(i,j,ny)
   30       rho(ij)=store(2*ij-1)*si2mg*km2m/(nx*ny)
      return
      end

      subroutine mlayer(mag,nx,ny,dx,dy,z1,z2,mi,md,fi,fd,store)
c
c  Subroutine MLAYER calculates the total-field anomaly on a two-
c  dimensional grid due to a horizontal layer with two-
c  dimensional magnetization.  The following steps are involved:
c  (1) Fourier transform the magnetization, (2) multiply by the
c  earth filter, and (3) inverse Fourier transform the product.
c  Magnetization is specified on a rectangular grid with x and y
c  axes directed north and east, respectively.  z axis is down.
c  Distance units irrelevant but must be consistent.  Requires
c  subroutines FOURN, DIRCOS, KVALUE, and MFILT.
c
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c         (NOTE: Both nx and ny must be powers of two.)
c    mag  - a one-dimensional real array containing the
c         two-dimensional magnetization (in A/m).  Elements should
c         be in order of west to east, then south to north (i.e.,
c         element 1 is southwest corner, element ny is 
c         southeast corner, element (nx-1)*ny+1 is northwest
c         corner, and element ny*nx is northeast corner).
c    store - a one-dimensional real array used internally.
c         It should be dimensioned at least 2*nx*ny.
c    dx - sample interval in the x direction, in km.
c    dy - sample interval in the y direction, in km.
c    z1 - depth to top of layer, in km.  Must be > 0.
c    z2 - depth to bottom of layer, in km.  Must be > z1.
c    mi - inclination of magnetization, in degrees positive below
c         horizontal.
c    md - declination of magnetization, in degrees east of north.
c    fi - inclination of regional field.
c    fd - declination of regional field.
c
c  Output parameters:
c    mag - upon output, mag will contain the total-field anomaly,
c         (in nT) with same orientation as before.
c
      complex cmag,mfilt,cmplx
      real mag,mi,md,mx,my,mz,kx,ky
      dimension mag(nx*ny),store(2*nx*ny),nn(2)
      data pi/3.14159265/,t2nt/1.e9/
      index(i,j,ncol)=(j-1)*ncol+i
      nn(1)=ny
      nn(2)=nx
      ndim=2
      call dircos(mi,md,0.,mx,my,mz)
      call dircos(fi,fd,0.,fx,fy,fz)
      dkx=2.*pi/(nx*dx)
      dky=2.*pi/(ny*dy)
      do 10 j=1,nx
         do 10 i=1,ny
            ij=index(i,j,ny)
            store(2*ij-1)=mag(ij)
   10       store(2*ij)=0.
      call fourn(store,nn,ndim,-1)
      do 20 j=1,nx
         do 20 i=1,ny
            ij=index(i,j,ny)
            call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
            cmag=cmplx(store(2*ij-1),store(2*ij))
            cmag=cmag*mfilt(kx,ky,mx,my,mz,fx,fy,fz,z1,z2)
            store(2*ij-1)=real(cmag)
   20       store(2*ij)=aimag(cmag)
      call fourn(store,nn,ndim,+1)
      do 30 j=1,ny
         do 30 i=1,nx
            ij=index(i,j,ny)
   30       mag(ij)=store(2*ij-1)*t2nt/(nx*ny)
      return
      end

      subroutine kvalue(i,j,nx,ny,dkx,dky,kx,ky)
c  Subroutine KVALUE finds the wavenumber coordinates of one
c  element of a rectangular grid from subroutine FOURN.
c
c  Input parameters:
c    i  - index in the kx direction.
c    j  - index in the ky direction.
c    nx - dimension of grid in ky direction (a power of two).
c    ny - dimension of grid in kx direction (a power of two).
c    dkx - sample interval in the kx direction.
c    dky - sample interval in the ky direction.
c
c  Output parameters:
c    kx - the wavenumber coordinate in the kx direction.
c    ky - the wavenumber coordinate in the ky direction.
c
      real kx,ky
      nyqx=nx/2+1
      nyqy=ny/2+1
      if(j.le.nyqx)then
         kx=(j-1)*dkx
         else
            kx=(j-nx-1)*dkx
            end if
      if(i.le.nyqy)then
         ky=(i-1)*dky
         else
            ky=(i-ny-1)*dky
            end if
      return
      end

      function gfilt(kx,ky,z1,z2)
c
c  Function GFILT calculates the value of the gravitational
c  earth filter at a single (kx,ky) location.
c
c  Input parameters:
c    kx - the wavenumber coordinate in the kx direction, in
c         units of 1/km.
c    ky - the wavenumber coordinate in the ky direction, in
c         units of 1/km.
c    z1 - the depth to the top of the layer, in km.
c    z2 - the depth to the bottom of the layer, in km.
c
c  Output parameters:
c    gfilt - the value of the earth filter
c
      real kx,ky,k
      data pi/3.14159265/,gamma/6.67e-11/
      k=sqrt(kx**2+ky**2)
      if(k.eq.0.)then
         gfilt=2.*pi*gamma*(z2-z1)
         else
            gfilt=2.*pi*gamma*(exp(-k*z1)-exp(-k*z2))/k
            end if
      return
      end

      function mfilt(kx,ky,mx,my,mz,fx,fy,fz,z1,z2)
c
c  Function MFILT calculates the complex value of the earth
c  filter at a single (kx,ky) location.
c
c  Input parameters:
c    kx - the wavenumber coordinate in the kx direction.
c    ky - the wavenumber coordinate in the ky direction.
c    mx - the x direction cosine of the magnetization vector.
c    my - the y direction cosine of the magnetization vector.
c    mz - the z direction cosine of the magnetization vector.
c    fx - the x direction cosine of the regional field vector.
c    fy - the y direction cosine of the regional field vector.
c    fz - the z direction cosine of the regional field vector.
c    z1 - the depth to the top of the layer.
c    z2 - the depth to the bottom of the layer.
c
c  Output parameters:
c    mfilt - the complex value of the earth filter
c
      complex mfilt,thetam,thetaf,cmplx
      real kx,ky,k,mx,my,mz
      data pi/3.14159265/,cm/1.e-7/
      k=sqrt(kx**2+ky**2)
      if(k.eq.0.)then
         mfilt=0.
         else
            thetam=cmplx(mz,(kx*mx+ky*my)/k)
            thetaf=cmplx(mz,(kx*fx+ky*fy)/k)
            mfilt=2.*pi*cm*thetam*thetaf*(exp(-k*z1)-exp(-k*z2))
            end if
      return
      end

      subroutine mtopo(mag,ztop,nx,ny,dx,dy,mi,md,fi,fd,nstop,err
     &                 store,cstore)
c 
c  Subroutine MTOPO calculates the total-field anomaly on a two-
c  dimensional grid due to an infinite half-space with uneven
c  top surface and two-dimensional magnetization.  Method
c  according to Parker (1972, pp. 447-55).
c  Magnetization is specified on a rectangular grid with x and y
c  axes directed north and east, respectively.  z axis is down
c  and anomaly is calculated at z=0; topographic grid should be
c  arranged accordingly.  Units of distance irrelevant but must be
c  consistent.  Requires subroutines FOURN, DIRCOS, KVALUE, and FAC.
c 
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c         (NOTE: Both nx and ny must be powers of two.)
c    mag  - a one-dimensional real array containing the
c         two-dimensional magnetization (in A/m).  Elements should
c         be in order of west to east, then south to north (i.e.,
c         element 1 is southwest corner, element ny is 
c         southeast corner, element (nx-1)*ny+1 is northwest
c         corner, and element ny*nx is northeast corner).
c    ztop - a one-dimensional real array containing the
c         topography of the upper surface, in same units as dx
c         and dy.  Grid layout same as mag.  Note:  z axis is
c         positive down and anomaly is calculated as z=0.  Array
c         ztop is modified by subroutine.
c    store - a one-dimensional real array used internally.
c         Should be dimensioned at least 2*nx*ny.
c    cstore - a one-dimensional complex array used internally.
c         Should be dimensioned at least nx*ny.
c    dx - sample interval in the x direction.
c    dy - sample interval in the y direction.
c    mi - inclination of magnetization, in degrees positive below
c         horizontal.
c    md - declination of magnetization, in degrees east of north.
c    fi - inclination of regional field.
c    fd - declination of regional field.
c    nstop - maximum number of iterations to try.
c    err - convergence criterion.  Iterations stop with the
c         contribution of the last term to the summation is less
c         than err times the contribution of all previous terms.
c 
c  Output parameters:
c    mag - upon output, mag contains the total-field anomaly
c         (in nT) with same grid orientation as before.
c
      dimension mag(nx*ny),ztop(nx*ny),store(2*nx*ny),cstore(nx*ny),
     &                 nn(2)
      real mag,mi,md,mx,my,mz,kx,ky,k
      complex cmplx,cstore,cstep,thetam,thetaf
      data pi/3.14159265/,t2nt/1.e9/,cm/1.e-7/
      index(i,j,ncol)=(j-1)*ncol+i
      nn(1)=ny
      nn(2)=nx
      ndim=2
      call dircos(mi,md,0.,mx,my,mz)
      call dircos(fi,fd,0.,fx,fy,fz)
      dkx=2.*pi/(nx*dx)
      dky=2.*pi/(ny*dy)
      ztpmax=-1.e20
      ztpmin= 1.e20
      do 1 j=1,nx
         do 1 i=1,ny
            ij=index(i,j,ny)
            ztpmax=amax1(ztpmax,ztop(ij))
    1       ztpmin=amin1(ztpmin,ztop(ij))
      ztpmed=ztpmin+(ztpmax-ztpmin)/2.
      do 2 j=1,nx
         do 2 i=1,ny
            ij=index(i,j,ny)
            ztop(ij)=ztop(ij)-ztpmed
    2       cstore(ij)=0.
      write(*,100)
  100 format(/,' Ratio = contribution of Nth term divided by',/,
     &         '         contribution of 0 through N-1 terms',//,
     &         '   N  | Ratio',/,
     &         ' -----|---------')
      n=-1
    3 n=n+1
      do 4 j=1,nx
         do 4 i=1,ny
            ij=index(i,j,ny)
            store(2*ij-1)=mag(ij)*ztop(ij)**n
    4       store(2*ij)=0.
      call fourn(store,nn,ndim,-1)
      abnew=0.
      abold=0.
      do 5 j=1,nx
         do 5 i=1,ny
            ij=index(i,j,ny)
            call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
            k=sqrt(kx**2+ky**2)
            arg=((-k)**n)*exp(-k*ztpmed)/fac(n)
            cstep=arg*cmplx(store(2*ij-1),store(2*ij))
            abnew=abnew+cabs(cstep)
            abold=abold+cabs(cstore(ij))
    5       cstore(ij)=cstore(ij)+cstep
      if(n.eq.0.)go to 3
      ratio=abnew/abold
      write(*,101)n,ratio
  101 format(1x,i5,g10.3)
      if(ratio.gt.err.and.n.lt.nstop)go to 3
      do 6 j=1,nx
         do 6 i=1,ny
            ij=index(i,j,ny)
            if(ij.eq.1)then
                 store(2*ij-1)=0.
               store(2*ij)=0.
               else
                  call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
                  k=sqrt(kx**2+ky**2)
                  thetam=cmplx(mz,(mx*kx+my*ky)/k)
                  thetaf=cmplx(fz,(fx*kx+fy*ky)/k)
                  cstore(ij)=2.*pi*cm*thetam*thetaf*cstore(ij)
                  store(2*ij-1)=real(cstore(ij))
                  store(2*ij)=aimag(cstore(ij))
                  end if
    6       continue
      call fourn(store,nn,ndim,+1)
      do 7 j=1,nx
         do 7 i=1,ny
            ij=index(i,j,ny)
    7       mag(ij)=store(2*ij-1)*t2nt/(nx*ny)
      return
      end

      subroutine contin(grid,nx,ny,dx,dy,dz,store)
c
c  Subroutine CONTIN upward continues gridded potential-field
c  data using the following steps:  (1) Fourier transform the field,
c  (2) multiply by the continuation filter, and (3) inverse Fourier
c  transform the product.  Field values are specified on a
c  rectangular grid with x and y axes directed north and east,
c  respectively.  z axis is down.  North is arbitrary.  Requires
c  subroutines FOURN and KVALUE.
c
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c         (NOTE:  both nx and ny must be powers of two.)
c    grid - a one-dimensional real array containing the
c         two-dimensional magnetization (in A/m).  Elements should
c         be in order of west to east, then south to north (i.e.,
c         element 1 is southwest corner, element ny is 
c         southeast corner, element (nx-1)*ny+1 is northwest
c         corner, and element ny*nx is northeast corner).
c    store - a one-dimensional real array used internally.
c         It should be dimensioned at least 2*nx*ny in the
c         calling progam
c    dx - sample interval in the x direction.
c    dy - sample interval in the y direction.
c    dz - continuation distance,in same units as dx and dy.  Should
c         be greater than zero for upward continuation.
c
c  Output parameters:
c    grid - upon output, grid contains the upward-continued
c         potential field with same orientation as before.
c
      dimension grid(nx*ny),store(2*nx*ny),nn(2)
      real kx,ky,k
      complex cgrid,cmplx
      data pi/3.14159265/
      index(i,j,ncol)=(j-1)*ncol+i
      nn(1)=ny
      nn(2)=nx
      ndim=2
      dkx=2.*pi/(nx*dx)
      dky=2.*pi/(ny*dy)
      do 10 j=1,nx
         do 10 i=1,ny
            ij=index(i,j,ny)
            store(2*ij-1)=grid(ij)
   10       store(2*ij)=0.
      call fourn(store,nn,ndim,-1)
      do 20 j=1,nx
         do 20 i=1,ny
            ij=index(i,j,ny)
            call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
            k=sqrt(kx**2+ky**2)
            cont=exp(-k*dz)
            cgrid=cmplx(store(2*ij-1),store(2*ij))*cont
            store(2*ij-1)=real(cgrid)
   20       store(2*ij)=aimag(cgrid)
      call fourn(store,nn,ndim,+1)
      do 30 j=1,nx
         do 30 i=1,ny
            ij=index(i,j,ny)
   30       grid(ij)=store(2*ij-1)/(nx*ny)
      return
      end

      subroutine verder(grid,nx,ny,dx,dy,norder,store)
c
c  Subroutine VERDER calculates the vertical derivative of
c  gridded potential-field data using the following steps:
c  (1) Fourier transform the field, (2) multiply by the vertical-
c  derivative filter, and (3) inverse Fourier transform the
c  with x and y axes directed north and east, respectively.
c  z axis is down.  North is arbitrarily.  Requires subroutines
c  FOURN and KVALUE.
c
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c         (NOTE:  both nx and ny must be powers of two.)
c    grid - a one-dimensional real array containing the
c         two-dimensional magnetization (in A/m).  Elements should
c         be in order of west to east, then south to north (i.e.,
c         element 1 is southwest corner, element ny is 
c         southeast corner, element (nx-1)*ny+1 is northwest
c         corner, and element ny*nx is northeast corner).
c    store - a one-dimensional real array used internally.
c         It should be dimensioned at least 2*nx*ny in the
c         calling progam
c    dx - sample interval in the x direction, units irrelevant.
c    dy - sample interval in the y direction, units irrelevant.
c    norder - the order of the vertical derivative.
c
c  Output parameters:
c    grid - upon output, grid contains the vertical derivative of
c         the potential field with same orientation as before.
c
      dimension grid(nx*ny),store(2*nx*ny),nn(2)
      complex cgrid,cmplx
      real kx,ky,k
      data pi/3.14159265/
      index(i,j,ncol)=(j-1)*ncol+i
      nn(1)=ny
      nn(2)=nx
      ndim=2
      dkx=2.*pi/(nx*dx)
      dky=2.*pi/(ny*dy)
      do 10 j=1,nx
         do 10 i=1,ny
            ij=index(i,j,ny)
            store(2*ij-1)=grid(ij)
   10       store(2*ij)=0.
      call fourn(store,nn,ndim,-1)
      do 20 j=1,nx
         do 20 i=1,ny
            ij=index(i,j,ny)
            call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
            k=sqrt(kx**2+ky**2)
            cgrid=cmplx(store(2*ij-1),store(2*ij))
            cgrid=cgrid*k**norder
            store(2*ij-1)=real(cgrid)
   20       store(2*ij)=aimag(cgrid)
      call fourn(store,nn,ndim,+1)
      do 30 j=1,nx
         do 30 i=1,ny
            ij=index(i,j,ny)
   30       grid(ij)=store(2*ij-1)/(nx*ny)
      return
      end

      subroutine newvec(grid,nx,ny,dx,dy,fi1,fd1,mi1,md1,fi2,
     &                  fd2,mi2,md2,store)
c
c  Subroutine NEWVEC transforms a gridded total-field anomaly
c  into a new anomaly with new directions of magnetization and
c  ambient field.  NEWVEC uses the following steps:  (1) Foureir
c  transform the field, (2) multiply by the phase filter, and
c  (3) inverse Fourier transform the product.  Anomaly values
c  are specified on a rectangular grid with x and y axes
c  directed north and east, respectively.  z axis is down.
c  Requires subroutines FOURN, DIRCOS, and KVALUE.
c
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c         (NOTE:  both nx and ny must be powers of two.)
c    grid - a one-dimensional real array containing the
c         two-dimensional magnetization (in A/m).  Elements should
c         be in order of west to east, then south to north (i.e.,
c         element 1 is southwest corner, element ny is 
c         southeast corner, element (nx-1)*ny+1 is northwest
c         corner, and element ny*nx is northeast corner).
c    store - a one-dimensional real array used internally.
c         It should be dimensioned at least 2*nx*ny in the
c         calling progam
c    dx - sample interval in the x direction, units irrelevant.
c    dy - sample interval in the y direction, units irrelevant.
c    mi1 - original inclination of magnetization, in degrees.
c    md1 - original declination of magnetization.
c    fi1 - original inclination of ambient field.
c    fd1 - original declination of ambient field.
c    mi2 - new inclination of magnetization.
c    md2 - new declination of magnetization.
c    fi2 - new inclination of ambient field.
c    fd2 - new declination of ambient field.
c
c  Output parameters:
c    grid - upon output, grid contains the transformed total-
c           field anomaly with same orientation as before.
c
      dimension grid(nx*ny),store(2*nx*ny),nn(2)
      complex cgrid,cmplx,thetam1,thetam2,thetaf1,thetaf2
     &        cphase
      real kx,ky,k,mi1,md1,mi2,md2,mx1,my1,mz1,mx2,my2,mz2
      data pi/3.14159265/
      index(i,j,ncol)=(j-1)*ncol+i
      nn(1)=ny
      nn(2)=nx
      ndim=2
      dkx=2.*pi/(nx*dx)
      dky=2.*pi/(ny*dy)
      call dircos(mi1,md1,0.,mx1,my1,mz1)
      call dircos(fi1,fd1,0.,fx1,fy1,fz1)
      call dircos(mi2,md2,0.,mx2,my2,mz2)
      call dircos(fi2,fd2,0.,fx2,fy2,fz2)
      do 10 j=1,nx
         do 10 i=1,ny
            ij=index(i,j,ny)
            store(2*ij-1)=grid(ij)
   10       store(2*ij)=0.
      call fourn(store,nn,ndim,-1)
      do 20 j=1,nx
         do 20 i=1,ny
            ij=index(i,j,ny)
            if(ij.eq.1)then
               cphase=0.
               else
                  call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
                  k=sqrt(kx**2+ky**2)
                  thetam1=cmplx(mz1,(kx*mx1+ky*my1)/k)
                  thetaf1=cmplx(fz1,(kx*fx1+ky*fy1)/k)
                  thetam2=cmplx(mz2,(kx*mx2+ky*my2)/k)
                  thetaf2=cmplx(fz2,(kx*fx2+ky*fy2)/k)
                  cphase=thetam2*thetaf2/(thetam1*thetaf1)
                  end if
            cgrid=cmplx(store(2*ij-1),store(2*ij))
            cgrid=cgrid*cphase
            store(2*ij-1)=real(cgrid)
   20       store(2*ij)=aimag(cgrid)
      call fourn(store,nn,ndim,+1)
      do 30 j=1,nx
         do 30 i=1,ny
            ij=index(i,j,ny)
   30       grid(ij)=store(2*ij-1)/(nx*ny)
      return
      end

      subroutine pseudo(grid,nx,ny,dx,dy,fi,fd,mi,md,store)
c
c  Subroutine PSEUDO transforms a gridded total-field anomaly
c  into a pseudogravity anomaly using the following steps:
c  (1) Fourier transform the field, (2) multiply by the phase
c  filter, and (3) inverse Fourier transform the product.  Anomaly
c  values are specified on a rectangular grid with x and y axes
c  directed north and east, respectively.  z axis is down.  Ratio
c  of density to magnetization assumed to be 100 kg/(m**3) per
c  1 A/m.  Requires subroutine FOURN, DIRCOS, and KVALUE.
c
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c         (NOTE:  both nx and ny must be powers of two.)
c    grid - a one-dimensional real array containing the
c         two-dimensional total-field anomaly in nT.  Elements
c         should be in order of west to east, then south to north 
c         (i.e., element 1 is southwest corner, element ny is 
c         southeast corner, element (nx-1)*ny+1 is northwest
c         corner, and element ny*nx is northeast corner).
c    store - a one-dimensional real array used internally.
c         Should be dimensioned at least 2*nx*ny.
c    dx - sample interval in the x direction, in km.
c    dy - sample interval in the y direction, in km.
c    mi - inclination of magnetization, in degrees.
c    md - declination of magnetization.
c    fi - inclination of ambient field.
c    fd - declination of ambient field.
c
c  Output parameters:
c    grid - upon output, grid contains the pseudogravity anomaly
c         in mGal with same orientation as before.
c
      dimension grid(nx*ny),store(2*nx*ny),nn(2)
      complex cgrid,cmplx,thetam,thetaf,cpsgr
      real kx,ky,k,mi,md,mx,my,mz,mag,km2m
      data cm/1.e-7/,gamma/6.67e-11/,t2nt/1.e9/,si2mg/1.e5/
      data pi/3.14159265/,rho/100./,mag/1./,km2m/1.e3/
      index(i,j,ncol)=(j-1)*ncol+i
      const=gamma*rho*si2mg*km2m/(cm*mag*t2nt)
      nn(1)=ny
      nn(2)=nx
      ndim=2
      dkx=2.*pi/(nx*dx)
      dky=2.*pi/(ny*dy)
      call dircos(mi,md,0.,mx,my,mz)
      call dircos(fi,fd,0.,fx,fy,fz)
      do 10 j=1,nx
         do 10 i=1,ny
            ij=index(i,j,ny)
            store(2*ij-1)=grid(ij)
   10       store(2*ij)=0.
      call fourn(store,nn,ndim,-1)
      do 20 j=1,nx
         do 20 i=1,ny
            ij=index(i,j,ny)
            if(ij.eq.1)then
               cphase=0.
               else
                  call kvalue(i,j,nx,ny,dkx,dky,kx,ky)
                  k=sqrt(kx**2+ky**2)
                  thetam=cmplx(mz,(mx*kx+my*ky)/k)
                  thetaf=cmplx(fz,(fx*kx+fy*ky)/k)
                  cpsgr=1./(thetam*thetaf*k)
                  end if
            cgrid=cmplx(store(2*ij-1),store(2*ij))
            cgrid=cgrid*cpsgr
            store(2*ij-1)=real(cgrid)
   20       store(2*ij)=aimag(cgrid)
      call fourn(store,nn,ndim,+1)
      do 30 j=1,nx
         do 30 i=1,ny
            ij=index(i,j,ny)
   30       grid(ij)=store(2*ij-1)*const/(nx*ny)
      return
      end

      subroutine hgrad(grid,nx,ny,dx,dy,store)
c
c  Subroutine HGRAD calculates the maximum horizontal
c  gradient of a two-dimensional function using simple finite-
c  difference relations.  Function is specified on a
c  rectangular grid with x and y axes directed north and east,
c  respectively.  North is arbitrary.
c
c  Input parameters:
c    nx - number of elements in the south-to-north direction.
c    ny - number of elements in the west-to-east direction.
c    grid - a one-dimensional real array containing the
c         two-dimensional function.  Elements should
c         be in order of west to east, then south to north (i.e.,
c         element 1 is southwest corner, element ny is 
c         southeast corner, element (nx-1)*ny+1 is northwest
c         corner, and element ny*nx is northeast corner).
c    dx - sample interval in the x direction, in km.
c    dy - sample interval in the y direction, in km.
c    store - a one-dimensional real array used internally.
c         Should be dimensioned at least nx*ny in the calling
c         program
c
c  Output parameters:
c    grid - upon output, grid contains the maximum horizontal
c           gradient with same orientation as before.
c
      dimension grid(nx*ny),store(nx*ny)
      index(i,j,ncol)=(j-1)*ncol+i
      dx2=2.*dx
      dy2=2.*dy
      do 10 j=1,nx
         jm1=j-1
         if(jm1.lt.1)jm1=1
         jp1=j+1
         if(jp1.gt.nx)jp1=nx
         do 10 i=1,ny
            im1=i-1
            if(im1.lt.1)im1=1
            ip1=i+1
            if(ip1.gt.ny)imp=ny
            dfdx=(grid(index(ip1,j,ny))-grid(index(im1,j,ny)))/dx2
            dfdy=(grid(index(i,jp1,ny))-grid(index(i,jm1,ny)))/dy2
            store(index(i,j,ny))=sqrt(dfdx**2+dfdy**2)
   10       continue
      do 20 i=1,nx*ny
   20 grid(i)=store(i)
      return
      end

      subroutine expand(grid,ncol,nrow,grid2,ncol2,nrow2)
c
c  Subroutine EXPAND adds "tapered" rows and columns to a
c  grid.  Input grid(i,j), i=1,2,...,ncol, j=1,2,...,nrow, is
c  modified as follows:
c
c  (1) Each row is expanded in length from ncol to ncol2.
c      New elements of row j linearly change in value from
c      grid(ncol,j) to grid(1,j).
c  (2) Each column is expanded in length from nrow to nrow2.
c      New elements of column 1 linearly change in value
c      from grid(i,nrow) to grid(i,1).
c  (3) All elements at i < or = to ncol and j < or = nrow
c      are left unchanged.
c
c  Input parameters:
c    grid -        one-dimensional real array representing
c                  a two-dimensional grid.
c    ncol,nrow -   dimensions of input grid.
c    ncol2,nrow2 - dimensions of output grid (grid2).
c                  ncol2 must be > ncol; nrow2 must be > nrow
c
c  Output parameters:
c    grid2 -       one-dimensional real array representing
c                  a two-dimension grid, as described earlier.
c
      dimension grid(ncol*nrow),grid2(ncol2*nrow2)
      index(i,j,ncol)=(j-1)*ncol+i
      do 10 j=1,nrow
         do 10 i=1,ncol
            ij=index(i,j,ncol)
            ij2=index(i,j,ncol2)
   10       grid2(ij2)=grid(ij)
      if(ncol2.gt.ncol)then
         do 20 j=1,nrow
            i1=index(1,j,ncol2)
            i2=index(ncol,j,ncol2)
            i3=index(ncol2,j,ncol2)
            step=(grid2(i1)-grid2(i2))/(ncol2-ncol+1)
            do 20 i=ncol+1,ncol2
               ij=index(i,j,ncol2)
   20          grid2(ij)=grid2(i2)+step*(i-ncol)
        else
           pause 'EXPAND:  ncol2.le.ncol'
           end if
      if(nrow2.gt.nrow)then
         do 30 i=1,ncol2
            j1=index(i,1,ncol2)
            j2=index(i,nrow,ncol2)
            j3=index(i,nrow2,ncol2)
            step=(grid2(j1)-grid2(j2))/(nrow2-nrow+1)
            do 30 j=nrow+1,nrow2
               ij=index(i,j,ncol2)
   30          grid2(ij)=grid2(j2)+step*(j-nrow)
         else
            pause 'EXPAND:  nrow2.le.nrow'
            end if
      return
      end
