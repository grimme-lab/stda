! This file is part of stda.
!
! Copyright (C) 2013-2019 Stefan Grimme
!
! stda is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! stda is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with stda.  If not, see <https://www.gnu.org/licenses/>.
!
      subroutine sosor(nroots,xmass,x,y)
      implicit none
      real*4 x(*)
      real*8 y(*),xmass
      real*8 xlam(6),r1(6)
      integer nroots

      integer i,j
      real*8 eau,xau,refval,vorfaktor,reffaktor
      real*8 refindex,rau
      logical da

c     data  ams /  1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0,
c    .10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,
c    .20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,
c    .30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,
c    .40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,
c    .54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,
c    .65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,
c    .79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,
c    .91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,
c    .102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,
c    .118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,
c    .132.9054d0, 137.3300d0, 15*0.000d0, 178.4900d0, 180.9479d0,
c    .183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,
c    .196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,
c    .18*0.000d0,   0.0000d0,  5*0.000d0/

************************************************************************
*                        ORD at 6 nm values          
* conversion from R to alpha:
* P. L. Polavarapu and D. K. Chakraborty
* Chem.Phys. 240 (1999) page 1
************************************************************************

c     xmass=0
c     do i=1,n
c        xmass=xmass+ams(idint(xyz(4,i)))
c     enddo

      xlam(1)=632.8 
      xlam(2)=589.3
      xlam(3)=579.
      xlam(4)=546.
      xlam(5)=436.
      xlam(6)=365.

c refractive index of solvent
      refindex=1.4d0
      reffaktor=(refindex**2+2.0d0)/3.0d0

      do j=1,6

      r1(j)=0.0d0
      do i=1,nroots
            xau =x(i)
c measurement point
            eau =1.d+7/xlam(j)/2.19474625d+5 
c R in au (input in 10-40 cgs) taken from TM
            rau =y(i)/64604.8   
c beta
            r1(j)=r1(j)+(2.*137.036/3.)*rau/(xau**2-eau**2)
      enddo
      vorfaktor=(38652./xmass)*(xlam(2)/xlam(j))**2
      r1(j)=r1(j)*vorfaktor*reffaktor
      enddo

      refval=0
      inquire(file='.ref',exist=da)
      if(da)then
      open(unit=33,file='.ref')
      read(33,*)refval
      close(33)
      endif

      write(*,*) 
      write(*,*) 'SOS specific optical rotation '
      write(*,*) 'including Lorentz factor for common solvent (n=1.4)' 
      write(*,*) 'lambda [eV] alpha[grad*cm^3*g^-1*dm^-1]' 
      do j=1,6
         if(j.eq.2)then
         write(*,142) xlam(j),1.d+7/(8065.54093*xlam(j)),
     .   r1(j),refval
         else
         write(*,143) xlam(j),1.d+7/(8065.54093*xlam(j)),
     .   r1(j)
         endif
      enddo
      write(*,*) 

 142  format(f6.1,f6.2,2f12.2,' ##')
 143  format(f6.1,f6.2,2f12.2)

      end

