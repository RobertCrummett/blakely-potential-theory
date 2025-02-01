program main
        write (*,*) 'Hello Blakely!'
end program main

subroutine sphere(xq,yq,zq,a,rho,xp,yp,zp,gx,gy,gz)
!
!  Subroutine SPHERE calculates the three components of gravitational
!  attraction at a single point due to a uniform sphere.
!
!  Input parameters:
!    Observation point is (xp,yp,zp), and center of sphere is at
!    (xq,yq,zq). Radius of sphere is a and density is rho. Density
!    in units of kg/(m**3). All distance parameters in units of km.
!
!  Output parameters:
!    Gravitatational components (gx,gy,gz) in units of mGal.
!
implicit none
real(8), parameter :: GAM=6.67e-11,SI2MG=1.e5,PI=3.14159265,KM2M=1.e3
real(8), intent(in) :: xq,yq,zq,a,rho,xp,yp,zp
real(8), intent(out) :: gx,gy,gz
real(8) :: rx,ry,rz,r,r3,tmass

rx=xp-xq
ry=yp-yq
rz=zp-zq
r=sqrt(rx**2+ry**2+rz**2)
if(r.eq.0.)then
        write(*,*) 'SPHERE: Bad argument detected.'
        stop
endif
r3=r**3
tmass=4.*PI*rho*(a**3)/3.
gx=-GAM*tmass*rx/r3
gy=-GAM*tmass*ry/r3
gz=-GAM*tmass*rz/r3
gx=gx*SI2MG*KM2M
gy=gy*SI2MG*KM2M
gz=gz*SI2MG*KM2M
end subroutine sphere

subroutine cylind(xq,zq,a,rho,xp,zp,gx,gz)
!
!  Subroutine CYLINDer calculates x and z components of gravitational
!  attraction due to a cylinder lying parallel to the y axis.
!
!  Input parameters:
!    Point of observation is (xp,zp). Axis of cylinder penetrates
!    x,z plan at (xq,zq). Radius of cylindar is a and density is
!    rho. Density in kg/(m**3). All distance parameters in km.
!
!  Output parameters:
!    Components of gravitational attraction (gx,gz) in mGal.
!
implicit none
real(8), parameter :: GAM=6.67e-11,SI2MG=1.e5,PI=3.14159265,KM2M=1.e3
real(8), intent(in) :: xq,zq,a,rho,xp,zp
real(8), intent(out) :: gx,gz
real(8) :: rx, rz, r2, tmass

rx=xp-xq
rz=zp-zq
r2=rx**2+rz**2
if(r2.eq.0._8)then
        write(*,*) 'CYLIND: Bad argument detected.'
        stop
endif
tmass=PI*(a**2)*rho
gx=-2.*GAM*tmass*rx/r2
gz=-2.*GAM*tmass*rz/r2
gx=gx*SI2MG*KM2M
gz=gz*SI2MG*KM2M
end subroutine cylind

subroutine dipole(xq,yq,zq,a,mi,md,m,xp,yp,zp,bx,by,bz)
!
!  Subroutine DIPOLE computes the three components of magentic
!  induction caused by a uniformly magnetized sphere. x axis
!  is north, z axis is down.
!  
!  Input parameters:
!    Observation point located at (xp,yp,zp). Sphere centered
!    at (xq,yq,zq). Magnetization of sphere defined by
!    intensity m, inclination mi, and declination md. Units
!    of distance irrelevant but must be consistent. All angles
!    in degrees. Intensity of magnetization in A/m. Requires
!    subroutine DIRCOS.
!
!  Output parameters:
!    The three components of magnetic induction (bx,by,bz) in
!    units of nT.
!
implicit none
real(8), parameter :: PI=3.14159265,T2NT=1.e9,CM=1.e-7
real(8), intent(in) :: xq,yq,zq,a,mi,md,m,xp,yp,zp
real(8), intent(out) :: bx,by,bz
real(8) :: rx,ry,rz,r2,r,r5,dot,moment,mx,my,mz

call dircos(mi,md,0._8,mx,my,mz)
rx=xp-xq
ry=yp-yq
rz=zp-zq
r2=rx**2+ry**2+rz**2
r=sqrt(r2)
if(r.eq.0._8)then
        write(*,*) 'DIPOLE: Bad argument detected.'
        stop
endif
r5=r**5
dot=rx*mx+ry*my+rz*mz
moment=4.*PI*(a**3)*m/3
bx=CM*moment*(3.*dot*rx-r2*mx)/r5
by=CM*moment*(3.*dot*ry-r2*my)/r5
bz=CM*moment*(3.*dot*rz-r2*mz)/r5
bx=bx*T2NT
by=by*T2NT
bz=bz*T2NT
end subroutine dipole

real(8) function schmidt(n,m,theta)
!
!  Returns Schmidt normalized associated Legendre polynomial.
!  Requires function fac. Modified from Press et al. (1986)
!
!  Input parameters:
!    Argument of polynomial is theta, in degrees. Degree and
!    order of polynomial are n and m, respectively. Parameter n
!    must be greater than zero, and m must be greater than or
!    equal to n.
!
implicit none
real(8), parameter :: D2RAD=.017453293
integer, intent(in) :: n,m
real(8), intent(in) :: theta
integer :: i,nn
real(8) :: x,fact,somx2,pmm,pmmp1,pnn,xnorm

interface
        integer function fac(n)
                implicit none
                integer, intent(in) :: n
        end function fac
end interface

x=cos(theta*D2RAD)

if (m.lt.0.or.m.gt.n) then
        write(*,*) 'Schmidt: Bad argument detected'
        stop
endif

pmm=1.

if (m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do i=1,m
                pmm=-pmm*fact*somx2
                fact=fact+2
        enddo
endif

if (n.eq.m) then
        schmidt=pmm
else
        pmmp1=x*(2*m+1)*pmm
        if (n.eq.m+1) then
                schmidt=pmmp1
        else
                do nn=m+2,n
                        pnn=(x*(2*nn-1)*pmmp1-(nn+m-1)*pmm)/(nn-m)
                        pmm=pmmp1
                        pmmp1=pnn
                enddo
                schmidt=pnn
        endif
endif

if (m.ne.0) then
        xnorm=sqrt(2.*fac(n-m)/fac(n+m))
        schmidt=xnorm*schmidt
endif
end function schmidt

integer function fac(n)
!
!  Function FAC calculates n!
!
implicit none
integer, intent(in) :: n
integer :: fac2

if (n.lt.0) then
        write (*,*) 'FAC: Bad argument detected.'
        stop
endif

if (n.eq.0.or.n.eq.1) then
        fac=1
        fac2=fac
30      fac2=fac2-1
        fac=fac*fac2
        if (fac2.gt.2) goto 30
endif
end function fac

subroutine gbox(x0,y0,z0,x1,y1,z1,x2,y2,z2,rho,g)
!
!  Subroutine GBOX computes the vertical attraction of a
!  rectangular prism. Sides of prism are parallel to x,y,z axes,
!  and z axis is vertical down.
!
!  Input parameters:
!    Observation point is (x0,y0,z0). The prism externs from x1
!    to x2, from y1 to y2, and from z1 to z2 in the x, y, and z
!    directions, respectively. Density of prism is rho. All
!    distance parameters in units of km; rho in units of
!    kg/(m**3).
!
!  Output parameters:
!    Vertical attraction of gravity g, in mGal.
!
implicit none
real(8), parameter :: GAM=6.670e-11,TWOPI=6.2831853,SI2MG=1.e5,KM2M=1.e3
integer, dimension(1:2), parameter :: SIGN=(/ -1, 1 /)
real(8), intent(in) :: x0,y0,z0,x1,y1,z1,x2,y2,z2,rho
real(8), intent(out) :: g
real(8) :: sum=0.,rijk,arg1,arg2,arg3
real(8), dimension(1:2) :: x,y,z
integer :: ijk,i,j,k

x(1)=x0-x1
y(1)=y0-y1
z(1)=z0-z1
x(2)=x0-x2
y(2)=y0-y2
z(2)=z0-z2
do i=1,2
        do j=1,2
                do k=1,2
                        rijk=sqrt(x(i)**2+y(j)*2+z(k)**2)
                        ijk=SIGN(i)*SIGN(j)*SIGN(k)
                        arg1=atan2(x(i)*y(j), z(k)*rijk)
                        if (arg1.lt.0.) arg1=arg1+twopi
                        arg2=rijk+y(j)
                        arg3=rijk+x(i)
                        if (arg2.le.0.) then
                                write (*,*) 'GBOX: Bad field point'
                                stop
                        endif
                        if (arg3.le.0.) then
                                write (*,*) 'GBOX: Bad field point'
                                stop
                        endif
                        arg2=dlog(arg2)
                        arg3=dlog(arg3)
                        sum=sum+ijk*(z(k)*arg1 - x(i)*arg2 - y(j)*arg3)
                enddo
        enddo
enddo
g=rho*GAM*sum*SI2MG*KM2M
end subroutine gbox

subroutine gpoly(x0,z0,xcorn,zcorn,ncorn,rho,g)
!
!  Subroutine GPOLY computes the vertical attraction of a two-
!  dimensional body with polygonal cross section. Axes are
!  right-handed system with y axis parallel to long direction
!  of body and z axis vertical down.
!
!  Input parameters:
!    Observation point is (x0,z0). Arrays xcorn and zcorn (each
!    of length ncorn) contain the coordinates of the polygon
!    corners, arranged in clockwise order when viewed with x axis
!    to right. Density of body is rho. All distance parameters
!    in units of km; rho in units of kg/(m**3).
!
!  Output parameters:
!    Vertical attraction of gravity g, in mGal.
!
implicit none
real(8), parameter :: GAM=6.67e-11,SI2MG=1.e5,KM2M=1.e3
real(8), intent(in) :: x0,z0,rho
integer, intent(in) :: ncorn
real(8), dimension(ncorn), intent(in) :: xcorn,zcorn
real(8), intent(out) :: g
integer :: n,n2
real(8) :: sum=0.,x1,z1,x2,z2,r1sq,r2sq,denom,alpha,beta,factor,term1,term2

do n=1,ncorn
        if (n.eq.ncorn) then
                n2=1
        else
                n2=n+1
        endif
        
        x1=xcorn(n)-x0
        z1=zcorn(n)-z0
        x2=xcorn(n2)-x0
        z2=zcorn(n2)-z0
        r1sq=x1**2+z1**2
        r2sq=x2**2+z2**2
        if (r1sq.eq.0.) then
                write (*,*) 'GPOLY: Field point on corner'
                stop
        endif
        if (r2sq.eq.0.) then
                write (*,*) 'GPOLY: Field point on corner'
                stop
        endif
        denom=z2-z1
        if (denom.eq.0.) denom=1.e-6
        alpha=(x2-x1)/denom
        beta=(x1*z2-x2*z1)/denom
        factor=beta/(1.+alpha**2)
        term1=0.5*(dlog(r2sq)-dlog(r1sq))
        term2=atan2(z2,x2)-atan2(z1,x1)
        sum=sum+factor*(term1-alpha*term2)
enddo
g=2.*rho*GAM*sum*SI2MG*KM2M
end subroutine gpoly

subroutine mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t)
!
!  Subroutine MBOX computes the total-field anomaly of an infinitely
!  extended rectangular prism. Sides of prism are parallel to x,y,z
!  axes, and z is vertical down. Bottom of prism extends to infinity.
!  Two calls to mbox can provide the anomaly of a prism with finite
!  thickness; e.g.,
!
!    call mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t1)
!    call mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t2)
!    t=t1-t2
!
!  Requires subroutine DIRCOS. Method from Bhattacharyya (1964).
!
!  Input parameters:
!    Observation point is (x0,y0,z0). Prism extends from x1 to
!    x2, y1 to y2, and z1 to infinity in x, y, and z directions,
!    respectively. Magnetization defined by inclination mi,
!    declination md, intensity m. Ambient field defined by
!    inclination fi and declintaion fd. x axis has declination
!    theta. Distance units are irrelevant but must be consistent.
!    Angles are in degrees, with inclinations positive below
!    horizontal and declinations positive east of true north.
!    Magnetization in A/m.
!
!  Output parameters:
!    Total-field anomaly t, in nT.
!
implicit none
real(8), parameter :: CM=1.e-7,T2NT=1.e9
real(8), intent(in) :: x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta
real(8), intent(out) :: t
integer :: i,j
real(8) :: ma,mb,mc,fa,fb,fc,fm1,fm2,fm3,fm4,fm5,fm6, &
           h,hsq,alphasq,sign,r0sq,r0,r0h,alphabeta, &
           arg1,arg2,arg3,arg4,tlog,tatan
real(8), dimension(2) :: alpha,beta

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
hsq=h**2

do i=1,2
        alphasq=alpha(i)**2
        do j=1,2
                sign=1.
                if (i.ne.j) sign=-1.
                r0sq=alphasq+beta(j)**2+hsq
                r0=sqrt(r0sq)
                r0h=r0*h
                alphabeta=alpha(i)*beta(j)
                arg1=(r0-alpha(i))/(r0+alpha(i))
                arg2=(r0-beta(j))/(r0+beta(j))
                arg3=alphasq+r0h+hsq
                arg4=r0sq+r0h-alphasq
                tlog=fm3*dlog(arg1)/2.+fm2*dlog(arg2)/2. &
                    -fm1*dlog(r0+h)
                tatan=-fm4*atan2(alphabeta,arg3) &
                      -fm5*atan2(alphabeta,arg4) &
                      -fm6*atan2(alphabeta,r0h)
                t=t+sign*(tlog+tatan)
        enddo
enddo
t=t*m*CM*T2NT
end subroutine mbox

subroutine dircos(incl,decl,azim,a,b,c)
!  
!  Subroutine DIRCOS computes direction cosines from inclination
!  and declination.
!  
!  Input parameters:

!    incl:  inclination in degrees positive below horizontal.
!    decl:  declinatino in degrees positive east of true north.
!    azim:  aximuth of x axis in degrees positive east of north.
!
!  Output parameters:
!    a,b,c: the three direction cosines.
!
implicit none
real(8), parameter :: D2RAD=.017453293
real(8), intent(in) :: incl,decl,azim
real(8), intent(out) :: a,b,c
real(8) :: xincl,xdecl,xazim

xincl=incl*D2RAD
xdecl=decl*D2RAD
xazim=azim*D2RAD
a=cos(xincl)*cos(xdecl-xazim)
b=cos(xincl)*sin(xdecl-xazim)
c=sin(xincl)
end subroutine dircos

subroutine facmag(mx,my,mz,x0,y0,z0,x,y,z,n,fx,fy,fz)
!
!  Subroutine FACMAG computes the magnetic field due to surface
!  charge on a polynomial face. Repeated calls can build the
!  field of an arbitrary polyhedron. x axis is directed north,
!  z axis vertical down. Requires subroutines CROSS, ROT, LINE,
!  and PLANE. Algorithm from Bott (1963).
!
!  Input parameters:
!    Observation point is (x0,y0,z0). Polygon corners defined
!    by arrays x, y, and z of length n. Magnetization given by
!    mx,my,mz. Polygon limited to 10 corners. Distance units
!    are irrelevant but must be consistent; magnetization in A/m.
!
!  Output parameters:
!    Three components of magnetic field (fx,fy,fz), in nT.
!
implicit none
real(8), parameter :: CM=1.e-7,T2NT=1.e9,EPS=1.e-20
real(8), intent(in) :: mx,my,mz,x0,y0,z0
integer, intent(in) :: n
real(8), intent(out) :: fx,fy,fz
integer :: i,j
real(8) :: nx,ny,nz,dot,rl,rn,px,py,pz,w,u1,v,w1,rk,us,v2s,v1s, &
           a2,a1,f2,f1,rho2,rho1,r2,r1,fu2,fu1,fv2,fv1, &
           fw2,fw1,fu,fv,fw
real(8), dimension(10) :: x,y,z,u,v2,v1,s,xk,yk,zk,xl,yl,zl

fx=0.
fy=0.
fz=0.
x(n+1)=x(1)
y(n+1)=y(1)
z(n+1)=z(1)
do i=1,n
        xl(i)=x(i+1)-x(i)
        yl(i)=y(i+1)-y(i)
        zl(i)=z(i+1)-z(i)
        rl=sqrt(xl(i)**2+yl(i)**2+zl(i)**2)
        xl(i)=xl(i)/rl
        yl(i)=yl(i)/rl
        zl(i)=zl(i)/rl
enddo

call cross(xl(2),yl(2),zl(2),xl(1),yl(1),zl(1),nx,ny,nz,rn)
nx=nx/rn
ny=ny/rn
nz=nz/rn
dot=mx*nx+my*ny+mz*nz
if (dot.eq.0.) return

call plane(x0,y0,z0,x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3),px,py,pz,w)
do i=1,n
        call rot(x(i),y(i),z(i),x(i+1),y(i+1),z(i+1),nx,ny,nz,px,py,pz,s(i))
        if (s(i).eq.0.) continue
        call line(px,py,pz,x(i),y(i),z(i),x(i+1),y(i+1),z(i+1),u1,v,w1,v1(i),v2(i),u(i))
        rk=sqrt((u1-px)**2+(v-py)**2+(w1-pz)**2)
        xk(i)=(u1-px)/rk
        yk(i)=(v-py)/rk
        zk(i)=(w1-pz)/rk
enddo

do j=1,n
        if (s(j).eq.0.) continue
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
        if (w.ne.0.) then
                fu2=(a2/f2)*dlog((r2+rho2)/abs(w)) &
                    -.5*dlog((r2+v2(j))/(r2-v2(j)))
                fu1=(a1/f1)*dlog((r1+rho1)/abs(w)) &
                    -.5*dlog((r1+v1(j))/(r1-v1(j)))
                fv2=(1./f2)*dlog((r2+rho2)/abs(w))
                fv1=(1./f1)*dlog((r1+rho1)/abs(w))
                fw2=atan2((a2*(r2-abs(w))),(r2+a2*a2*abs(w)))
                fw1=atan2((a1*(r1-abs(w))),(r1+a1*a1*abs(w)))
                fu=dot*(fu2-fu1)
                fv=-dot*(fv2-fv1)
                fw=(-w*dot/abs(w))*(fw2-fw1)
        else
                fu2=(a2/f2)*(1.+dlog((r2+rho2)/EPS)) &
                    -.5*dlog((r2+v2(j))/(r2-v2(j)))
                fu1=(a1/f1)*(1.+dlog((r1+rho1)/EPS)) &
                    -.5*dlog((r1+v1(j))/(r1-v1(j)))
                fv2=(1./f2)*(1.+dlog((r2+rho2)/EPS))
                fv1=(1./f1)*(1.+dlog((r1+rho1)/EPS))
                fu=dot*(fu2-fu1)
                fv=-dot*(fv2-fv1)
                fw=0.
        endif
        fx=fx-s(j)*(fu*xk(j)+fv*xl(j)+fw*nx)
        fy=fy-s(j)*(fu*yk(j)+fv*yl(j)+fw*ny)
        fz=fz-s(j)*(fu*zk(j)+fv*zl(j)+fw*nz)
enddo
fx=fx*CM*T2NT
fy=fy*CM*T2NT
fz=fz*CM*T2NT
end subroutine facmag

subroutine plane(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z,r)
!
!  Subroutine PLANE computes the intersection (x,y,z) of a plane
!  and a perpendicular line. The plane is defined by the three points
!  (x1,y1,z1), (x2,y2,z2), and (x3,y3,z3). The line passes through
!  (x0,y0,z0). Computation is done by transformation and inverse
!  transformation of coordinates systems.
!
implicit none
real(8), intent(in) :: x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3
real(8), intent(out) :: x,y,z,r
real(8) :: x2n,y2n,z2n,x0n,y0n,z0n,x3n,y3n,z3n,a,t11,t12,t13, &
           t21,t22,t23,t31,t32,t33,tx0,tz0,cx,cy,cz,c,dx,dy,dz,d

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
end subroutine plane

subroutine line(x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z,v1,v2,r)
!
!  Subroutine LINE determines the intersection (x,y,z) of two
!  lines. First line is defined by points (x1,y1,z1) and
!  (x2,y2,z2). Second line is perpendicular to the first and
!  passes through point (x0,y0,z0). Distance between (x,y,z)
!  and (x0,y0,z0) is returned as r. Computation is done by a
!  transformation of coordinate systems.
!
implicit none
real(8), intent(in) :: x0,y0,z0,x1,y1,z1,x2,y2,z2
real(8), intent(out) :: x,y,z,v1,v2,r
real(8) :: tx0,ty0,tz0,tx2,ty2,tz2,a,cx,cy,cz,c,dx,dy,dz,d, &
           tt11,tt12,tt13,tt21,tt22,tt23,tt31,tt32,tt33,u0

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
end subroutine line

subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz,r)
!
! Subroutine CROSS computes the vector product of two vectors: i.e.,
!
!               (cx,cy,cz) = (ax, ay, az) X (bx, by, bz)
!
implicit none
real(8), intent(in) :: ax,ay,az,bx,by,bz
real(8), intent(out) :: cx,cy,cz,r
cx=ay*bz-az*by
cy=az*bx-ax*bz
cz=ax*by-ay*bx
r=sqrt(cx**2+cy**2+cz**2)
end subroutine cross

subroutine rot(ax,ay,az,bx,by,bz,nx,ny,nz,px,py,pz,s)
!
!  Subroutine ROT finds the sense of rotation of the vector
!  from (ax,ay,az) to (bx,by,bz) with respect to a second
!  vector through point (px,py,pz). The second vector has
!  components given by (nx,ny,nz). Returned parameter s is
!  1 if anticlockwise, -1 if clockwise, or 0 if collinear.
!
implicit none
real(8), intent(in) :: ax,ay,az,bx,by,bz,nx,ny,nz,px,py,pz
real(8), intent(out) :: s
real(8) :: x,y,z,cx,cy,cz,c,u,v,w,d
x=bx-ax
y=by-ay
z=bz-az
call cross(nx,ny,nz,x,y,z,cx,cy,cz,c)
u=px-ax
v=py-ay
w=pz-az
d=u*cx+v*cy+w*cz
if (d.lt.0.) then
        s=1.
else if (d.eq.0.) then
        s=0.
else
        s=-1.
endif
end subroutine rot

subroutine ribbon(x0,z0,x1,z1,x2,z2,mx,mz,fx,fz,ier)
!
!  Subroutine RIBBON computes the x and z components of magnetic
!  induction due to a single side of a two-dimensional prism with
!  polygonal cross section. The prism is assumed infinitely
!  extended parallel to the y axis; z axis is vertical down.
!
!  Input parameters:
!    Observation point is (x0,z0). Coordinates (x1,z1) and
!    (x2,z2) are two consecutive corners of the polygon taken in
!    clockwise order around the polygon as viewed with the a
!    axis to the right. The x and z components of magnetization
!    mx and mz. Distance units are irrelevant but must be
!    consistent; magnetization in A/m.
!
!  Output paramters:
!    Components of magnetic field fx and fz, in nT.
!    Errors are recorded by ier:
!      ier=0, no errors;
!      ier=1, two corners are too close (no calculation);
!      ier=2, field point too close to corner (calculation
!             continues).
!
implicit none
real(8), parameter :: PI=3.14159265,SMALL=1.e-18,CM=1.e-7,T2NT=1.e9
real(8), intent(in) :: x0,z0,x1,z1,x2,z2,mx,mz
real(8), intent(out) :: fx,fz
integer, intent(out) :: ier
real(8) :: sx,sz,s,qs,rx1,rz1,rx2,rz2,r1,r2,theta1, &
           theta2,angle,flog,factor
ier=0
sx=x2-x1
sz=z2-z1
s=sqrt(sx**2+sz**2)

!
!  -- If two corners are too close, return
!
if (s.lt.SMALL) then
        ier=1
        return
endif

sx=sx/s
sz=sz/s
qs=mx*sz-mz*sx
rx1=x1-x0
rz1=z1-z0
rx2=x2-x0
rz2=z2-z0

!
!  -- If field point is too near a corner, signal error
!
if (abs(rx1).lt.SMALL.and.abs(rz1).lt.SMALL) ier=2
if (abs(rx2).lt.SMALL.and.abs(rz2).lt.SMALL) ier=2
if (ier.eq.2) then
        rx1=SMALL
        rz1=SMALL
        rx2=SMALL
        rz2=SMALL
endif
r1=sqrt(rx1**2+rz1**2)
r2=sqrt(rx2**2+rz2**2)
theta1=atan2(rz1,rx1)
theta2=atan2(rz2,rx2)
angle=theta1-theta2
if (angle.gt.PI) angle=angle-2.*PI
if (angle.lt.-PI) angle=angle+2.*PI

!
!  -- If field point is too close to side, signal error
!
if (abs(angle).gt.(.995*PI)) ier=2
flog=dlog(r2)-dlog(r1)
factor=-2.*CM*qs*T2NT
fx=factor*(sx*flog-sz*angle)
fz=factor*(sz*flog+sx*angle)
end subroutine ribbon

subroutine fork(lx,cx,signi)
!
!  Subroutine FORK calculates the Fourier transform of a one-
!  dimesional array. Algorithm from Claerbout (1976).
!
!  Input/output parameters:
!    Complex array cx of length lx is the input array. Upon
!    return, cx contains the transformed array. Length of
!    array must be a power of 2. If signi=-1., then the forward
!    calculation is performed; if signi=1., the inverse transform
!    is performed.
!
integer, intent(in) :: lx
complex(4), dimension(lx) :: cx
integer, intent(out) :: signi
complex(4) :: carg,cw,ctemp
real(4) :: sc
integer :: i,istep,ipl,j=1,l,m

sc=sqrt(1./lx)
do i=1,lx
        if (i.gt.j) goto 2
        ctemp=cx(j)*sc
        cx(j)=cx(i)*sc
        cx(i)=ctemp
2       m = lx/2
3       if (j.le.m) continue
        j=j-m
        m=m/2
        if (m.ge.1) goto 3
        j=j+m
enddo

l=1
6 istep=2*l
do m=1,l
        carg=(0.,1.)*(3.14159265*signi*(m-1))/l
        cw=cexp(carg)
        do i=m,lx,istep
                ipl=i+l
                ctemp=cw*cx(ipl)
                cx(ipl)=cx(i)-ctemp
                cx(i)=cx(i)+ctemp
        enddo
enddo
l=istep
if (l.lt.lx) goto 6
end subroutine fork
