
c
c
c     =================================================================
      subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     =================================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd component of q.
c     ------------------------------------------------
c
c     # Extend the data from the computational region
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)

      dimension mthbc(2)
      
!       dimension p(1:4000),pd(1:4000),time(1:4000)
      dimension timeIC(500), pressure(500), p(0:2)
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      n = 4000
c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      p0 = 101325.0
      c0 = sqrt(gammagas*p0/rhog) 
      p = p0
      ! Read pressure and time data from file
      open (25,file="a-pICtime.dat",action="read",status="old") 
      do i=1,500
          read(25,*) timeIC(i), pressure(i)
      enddo
      close(25)
      
      ! Look for correct value for pressure in data file
      ddt = timeIC(100) - timeIC(101)
      delt = max(ddt,dt)
      do k=1,500
        time2 = (timeIC(k) - 0.018758) 
        if (abs(time2 - t) <  delt) then
          ! Convert PSI to Pascals
          p(0) = p0 + 6894.75729*pressure(k) !*exp(-0.1*ycell**2)
        end if
      end do
      p(2) = p(0)
      p(1) = p(0)
      ! Assign corresponding pressure values to left boudary ghost cells
      do ibc = 1, mbc  
        q(1,1-ibc) = rhog*(p(ibc)/p0)**(1/gammagas)
        q(2,1-ibc) = (2.0/(gammagas - 1.0))*(-c0 + 
     & sqrt(gammagas*p(ibc)/q(1,1-ibc)))
        q(3,1-ibc) = (p(ibc) + gammagas*pinfgas)/(gammagas - 1.0) +
     & (q(2,1-ibc)**2)/(2.0*q(1,1-ibc))
      end do
      
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 m=1,meqn
         do 115 ibc=1,mbc
               q(m,1-ibc) = q(m,1)
  115       continue
      go to 199

  120 continue
c     # periodic:  
      do 125 m=1,meqn
         do 125 ibc=1,mbc
               q(m,1-ibc) = q(m,mx+1-ibc)
  125       continue
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 ibc=1,mbc
               q(m,1-ibc) = q(m,ibc)
  135       continue
c     # negate the normal velocity:
      do 136 ibc=1,mbc
            q(2,1-ibc) = -q(2,ibc)
  136    continue
      go to 199

  199 continue

c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 m=1,meqn
         do 215 ibc=1,mbc
               q(m,mx+ibc) = q(m,mx)
  215       continue
      go to 299

  220 continue
c     # periodic:  
      do 225 m=1,meqn
         do 225 ibc=1,mbc
               q(m,mx+ibc) = q(m,ibc)
  225       continue
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 ibc=1,mbc
               q(m,mx+ibc) = q(m,mx+1-ibc)
  235       continue
      do 236 ibc=1,mbc
            q(2,mx+ibc) = -q(2,mx+1-ibc)
  236    continue
      go to 299

  299 continue
c
      return
      end

