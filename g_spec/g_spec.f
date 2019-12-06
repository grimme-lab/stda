************************************************************************
*
************************************************************************

      program main

      implicit real*8 (a-h,o-z)

!      parameter (imax=100000)
!      parameter (jmax=2000)

      integer*8 imax,jmax
      dimension xx(50),xlam(6), r1(10),r2(10)
      real*8, allocatable :: x(:),y(:),e(:),r(:), wei(:)
      integer, allocatable :: iff(:),lab(:)
!      dimension x(imax), y(imax), e(jmax), r(jmax), wei(imax), lab(imax)

      character*120  a8
      logical nol
      data pi /3.14159265358979D0/

************************************************************************
*                       defaults
************************************************************************


      nol=.false.
      itype=2
      iuv=0
      icm=0
      inm=0
      ln=0
      iord=0
      ienan=0
      imax=0
      jmax=0
      fau=1.0d0
      rfaktor=1.0d0
      fktl=1.0d0
      sh=0.0d0
      xmass=1.0d0
      ethr=0.0d0
      forb=0.0d0

      ! preread to get number of states
      do
      read (5,*, end=200)
      imax = imax + 1
      enddo
  200 rewind(5)
    
      tmp=dble(imax)/50.0
      jmax=ceiling(tmp)
      if(jmax.le.2000) then
        jmax=2000
      endif

      allocate(x(imax),y(imax),e(jmax),r(jmax),wei(imax),lab(imax),
     .         iff(imax))

      do i=1,imax
         x(i)=0.0d0
         y(i)=0.0d0
      enddo

      do i=1,jmax
         e(i)=0.0d0
         r(i)=0.0d0
      enddo

************************************************************************
*                       read
************************************************************************


  10  read(5,'(a)',end=100) a8
      if(index(a8,'AU ').ne.0) then
         fau=27.21139570d0
         write(*,*) '** energies in au    **'
      endif
      if(index(a8,'ETHR ').ne.0) then
         write(*,*) '** equal energy thr  **'
         read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         ethr=xx(1)
      endif
      if(index(a8,'FORB ').ne.0) then
         write(*,*) '** forbidden trans set to  **'
         read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         forb=xx(1)
      endif
      if(index(a8,'UV').ne.0) then
         iuv=1
         write(*,*)'** uv-spektrum        **'
      endif
      if(index(a8,'ORD').ne.0) then
         iord=1
         write(*,*)'** cal ord (ord.dat)  **'
      endif
      if(index(a8,'NM').ne.0) then
         inm=1
        write(*,*)'** nm-scale           **'
      endif
      if(index(a8,'CM').ne.0) then
         icm=1
         write(*,*)'** cm-scale enabled   **'
      endif
      if(index(a8,'VELO').ne.0) then
         itype=3
         write(*,*)'** read velocity data **'
      endif
      if(index(a8,'MIX').ne.0) then
         itype=4
         write(*,*)'** read mixed l-v data **'
      endif
      if(index(a8,'LOG').ne.0) then
         ln=1
         write(*,*)'** log scale enabled  **'
      endif
      if(index(a8,'INV').ne.0) then
         ienan=1
         write(*,*)'** enantiomer         **'
      endif
      if(index(a8,'MMASS').ne.0) then
         read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         xmass=xx(1)
         write(*,*)'** molar mass = ',xmass,' **'
      endif
      if(index(a8,'RFAKTOR').ne.0) then
         read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         rfaktor=xx(1)
         write(*,*)'** rfaktor = ', rfaktor,' **'
      endif
      if(index(a8,'LFAKTOR').ne.0) then
         read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         fktl=xx(1)
         write(*,*)'** lfaktor = ', fktl,' **'
         if(fktl.lt.1.0d-5) nol=.true.
      endif
      if(index(a8,'SHIFT').ne.0) then
         read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         sh=xx(1)
         write(*,*) '** shift = ',sh,' **'
      endif
      if(index(a8,'WIDTH').ne.0) then
         read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         wi=xx(1)
         write(*,*) '** width = ',wi,' **'
         do kk=1,imax
            wei(kk)=wi
         enddo
      endif
      if(index(a8,'DATAX').ne.0) then
         k=0
  20     read(5,'(a)',end=100) a8
         if(index(a8,'DATAY').ne.0) goto 22
         call readl(a8,xx,n)
         if(n.gt.0.and.k.le.imax) then
            k=k+1
            kk=k
            if(n.gt.1) kk=idint(xx(1))
            lab(kk)=kk
            x(kk)=xx(2)*fau+sh
         endif
         goto 20
      endif

  22  if(index(a8,'DATAY').ne.0) then
         k=0
  25     read(5,'(a)',end=100) a8
         if(index(a8,'DATAW').ne.0) goto 21
         if(index(a8,'END').ne.0) goto 21
         call readl(a8,xx,n)
         if(n.gt.0.and.k.le.imax) then
            k=k+1
            kk=k
            if(n.gt.1) kk=idint(xx(1))
            y(kk)=xx(itype)
            if(ienan.eq.1) y(kk)=-1.0d0*y(kk)
            if(ln.eq.1.and.y(kk).lt.1.d-6)y(kk)=1.d-4
         endif
         goto 25
      endif

  21  if(index(a8,'DATAW').ne.0) then
  27     read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         if(n.eq.2.and.k.le.imax) then
            kk=idint(xx(1))
            wei(kk)=xx(2)
         endif
         goto 27
      endif

  23  if(index(a8,'DATXY').ne.0) then
         k=0
  29     read(5,'(a)',end=100) a8
         call readl(a8,xx,n)
         if(n.gt.0.and.k.le.imax) then
            k=k+1
            kk=k
            lab(kk)=kk
            x(kk)=xx(2)*fau+sh
            if(iuv.eq.0)y(kk)=xx(5+itype-2)
            if(iuv.eq.1)y(kk)=xx(3+itype-2)
            if(ienan.eq.1) y(kk)=-1.0d0*y(kk)
            if(ln.eq.1.and.y(kk).lt.1.d-6)y(kk)=1.d-4
         endif
         goto 29
      endif

      goto 10
 100  continue

      nroots=k


************************************************************************
*                       sort
************************************************************************

      if(fau.lt.1.10d0) then
         do i=1,nroots
            x(i)=-1*x(i)
         enddo
      endif

      call dsort(nroots,lab,y,wei,x)

      if(fau.lt.1.10d0) then
         do i=1,imax
            x(i)=-1*x(i)
         enddo
      endif

      if (fau.gt.1.10d0) then
         do i=2,nroots
            if(abs(x(i)).gt.1.0d-6) then
              x(i)=x(i)-x(1)
            endif
         enddo
         x(1)=0.0d0
      endif

      xma=-1.d10
      xmi= 1.d10
      do i=1,nroots
         if(x(i).gt.xma) xma=x(i)
         if(x(i).lt.xmi) xmi=x(i)
         if(iuv.eq.1.and.forb.gt.1.d-8.and.y(i).lt.1.d-5)y(i)=forb
      enddo

      dx=xma-xmi
      xma=xma+0.40*dx
      xmi=xmi-0.40*dx
      if(xmi.lt.0)xmi=0.5

      dx=xma-xmi

      write(*,*) '** xmin, xmax = ',xmi, xma

      write(*,*) '** roots      = ', k

************************************************************************
*                        output (xmgr.dat)
************************************************************************

c zero line

      open(unit=8,file='rots.dat')
c      xi= 0.1
c      xa=15.0
c      if(inm.eq.1) then
c         xi=1.d+7/(8065.54093*xi)
c         xa=1.d+7/(8065.54093*xa)
c      endif
c      if(icm.eq.1) then
c         xi=8065.54093*xi
c         xa=8065.54093*xa
c      endif
c
c      if(iuv.eq.0.and.(.not.nol)) then
c         write(8,1000)xi,0.0
c         write(8,1000)xa,0.0
c         write(8,*)'&'
c      endif
c
      yma=-10000000.
      ymi= 10000000.
      ymax=-10000000.

************************************************************************
*                        screen and lines
************************************************************************
      
      av=0
      iav=0
      do i=1,nroots
         iff(i)=1
         if(y(i).gt.yma) yma=y(i)
         if(y(i).lt.ymi) ymi=y(i)
         if(abs(y(i)).gt.ymax) ymax=abs(y(i))
         yy=y(i)
         if(i.le.30)then
            av=av+abs(yy)
            iav=iav+1
         endif
         xnm=1.d+7/(8065.54093*x(i))
         xcm=8065.54093*x(i)
         if(x(i).ne.0.0d0) then
            write(*,'(2i3,3f12.4,f12.6,f12.2)')
     .      i,lab(i),x(i),xnm,xcm,yy,wei(i)
            yy=yy*fktl
            if(iuv.eq.1) yy=yy*1.0d5
            if(ethr.gt.0.0001) goto 542
            if(nol)            goto 542
            if(ln.eq.1) then
                        if(yy.lt.0)yy=1.d-6
                        yy=log10(yy)
            endif
            if(inm.eq.0.and.icm.eq.0)write(8,1000) x(i),yy*rfaktor
            if(inm.eq.1) write(8,1000) xnm,yy*rfaktor
            if(icm.eq.1) write(8,1000) xcm,yy*rfaktor
 542        continue
         endif
      enddo
      write(*,*)'AV ',av/float(iav)

      if(ethr.lt.0.0001) goto 642

      do i=1,nroots
         sum=y(i)
         if(iff(i).gt.0) then
         do j=i+1,nroots
            de=abs(x(i)-x(j))
            if(de.lt.ethr) then
               iff(j)=0
               sum=sum+y(j)
            endif
         enddo
         if(inm.eq.0.and.icm.eq.0)write(8,1000) x(i),sum
         if(inm.eq.1) write(8,1000) 1.d+7/(8065.54093*x(i)),sum 
         if(icm.eq.1) write(8,1000) 8065.54093*x(i),sum
         endif
      enddo

 642  continue

************************************************************************
*                        CD and UV gauss curves
************************************************************************

      ymauv=yma
      ymiuv=ymi
c      if(.not.nol) write(8,*)'&'

      np=jmax

      f1=22.97
      if(iuv.eq.1) f1=9.186d-3
      st=(xma-xmi)/np

      do j=1,np
         ee=xmi+j*st
         e(j)=ee
         do i=1,nroots
            if(abs(x(i)).gt.1.0d-6) then
               ! wei (i.e., WIDTH in input) is the half width at 1/e maximum of the Gaussian
               f3=1.0d0/(f1*1.7724539*wei(i)) 
               width=1.0/wei(i)**2
c
c ! uncomment this if: wei (i.e., WIDTH in input) is the half width at 1/2 maximum (HWHM) of the Gaussian
c              f3=1.0d0/(f1*2.5066272*wei(i)) !
c              width=0.5/wei(i)**2
c
               edif=ee-x(i)
               if(iuv.eq.0) then
                  f2=x(i)*y(i)
               else
                  dip2=6.4607*1.5*y(i)/(3.67493d-2*x(i))
                  f2=x(i)*dip2 
               endif
               r(j)=r(j)+f3*f2*exp(-width*edif**2)
c              write(*,*) 1.0d0/(f1*1.7724539)*6.4607*1.5/3.67493d-2
c              r(j)=r(j)+(1.6196e+4/wei(i))*y(i)*exp(-width*edif**2)
c              r(j)=r(j)+1./wei(i)*exp(-width*edif**2)
            endif
         enddo
       enddo

      yma=-10000000.
      ymi= 10000000.
      do k=1,jmax
         if(r(k).gt.yma) yma=r(k)
         if(r(k).lt.ymi) ymi=r(k)
      enddo

      write(*,*) 'writing rots.dat'
      write(*,*)'** ymin, ymax = ',ymi, yma

      f=yma/ymauv

      write(*,*) 'writing spec.dat'
      open(unit=9,file='spec.dat')
      do i=1,np
         r(i)=r(i)*rfaktor
         if(inm.eq.1) e(i)=1.d+7/(8065.54093*e(i))
         if(icm.eq.1) e(i)=8065.54093*e(i)
         if(ln.eq.1) then
            xxx=log10(r(i))
            if(xxx.lt.0) xxx=0
c            write(8,1000) e(i), xxx
            write(9,1000) e(i), xxx
         else
c            write(8,1000) e(i), r(i)
            write(9,1000) e(i), r(i)
         endif
      enddo
     
      close (9)

      goto 742
c      write(8,*)'&'

************************************************************************
*                        lines with low intensity (< 1% of max)
************************************************************************

      write(8,*)'0.0 0.0'
      if(nol) goto 742

      xxx=100.
      yyy=0.0
      write(8,1000) xxx,yyy
      write(*,*) ymax
      do i=1,nroots
         yy=y(i)
         xnm=1.d+7/(8065.54093*x(i))
         xcm=8065.54093*x(i)
         thr=0.02*ymax
         if(x(i).ne.0.0d0) then
            if(abs(yy).lt.thr) then
               if(ln.eq.1) yy=log10(yy)
               if(inm.eq.0.and.icm.eq.0)write(8,1000) x(i),yy
               if(inm.eq.1) write(8,1000) xnm,yy
               if(icm.eq.1) write(8,1000) 0.001*xcm,yy
            endif
         endif
      enddo

      if(iuv.ne.0) goto 742
      if(iord.eq.0) goto 942
     

 1000 format(f16.6,e14.4)

************************************************************************
*                        ORD spectrum (ord.dat)               
************************************************************************

      do i=1,jmax
         e(i)=0.0d0
         r(i)=0.0d0
      enddo

      open(unit=7,file='ord.dat')
      xmi=0.5
      st=(xma-xmi)/np
      do j=1,np
         ee=xmi+j*st
         e(j)=1.d+7/(8065.54093*ee)
      do i=1,nroots
         yy=y(i)
         if(x(i).ne.0.0d0) then
            xnm =1.d+7/(8065.54093*x(i))
            wnm1=1.d+7/(8065.54093*(x(i)+wei(i)))
            wnm2=1.d+7/(8065.54093*(x(i)-wei(i)))
            wnm=abs(wnm2-wnm1)
            pref=yy*9.145d+1*xnm/wnm
            xxx=(e(j)-xnm)/wnm
            r(j)=r(j)+pref*(xintx2(xxx)-0.5*wnm/(e(j)+xnm))
         endif
      enddo
      enddo

      do j=1,np
         if(icm.eq.1) e(j)=1.d+4/e(j)
         write(7,'(2D14.6)') e(j),r(j)
      enddo


942   continue

************************************************************************
*                        ORD at 6 nm values          
* conversion from R to alpha:
* P. L. Polavarapu and D. K. Chakraborty
* Chem.Phys. 240 (1999) page 1
************************************************************************

      xlam(1)=632.8 
      xlam(2)=589.3
      xlam(3)=579.
      xlam(4)=546.
      xlam(5)=436.
      xlam(6)=365.

      pi=3.14159265358979D0
c refractive index of solvent
      refindex=1.0d0
      reffaktor=(refindex**2+2.0d0)/3.0d0
      vorfaktor=reffaktor*1.333333333*0.1343d-3*pi/xmass

      do j=1,6

         r1(j)=0.0d0
      do i=1,nroots
         if(x(i).ne.0.0d0) then
c energy in au (input in eV)
            xau =x(i)/27.2113957
c measurement point
            eau =1.d+7/xlam(j)/2.19474625d+5 
c angular frequency
            eau=eau*2.*pi
            xau=xau*2.*pi
c energy in cm-1, mesaurement point
            ecm=1.d+7/xlam(j)
c R in au (input in 10-40 cgs)
            rau =y(i)/235.730d0 
c
            r1(j)=r1(j)+vorfaktor*ecm**2*rau/(xau**2-eau**2)

         endif
      enddo
      enddo

      write(*,*) '*****************************'
      write(*,*) ' K R O N I G - K R A M E R S'
      write(*,*) '    TRANSFORMATION (delta)'
      write(*,*) '*****************************'
      write(*,*) '  specific optical rotation '
      write(*,*) '  [alp] in grad*cm^3*g^-1*dm^-1'         
      write(*,*) 
      write(*,*) ' lambda  eV        alp    alp*MolW/100' 
      do j=1,6
         write(*,142) xlam(j),1.d+7/(8065.54093*xlam(j)),
     .r1(j),r1(j)*xmass/100.
      enddo

c           xn  =1.d+7/(8065.54093*x(i))
c           wnm1=1.d+7/(8065.54093*(x(i)+wei(i)))
c           wnm2=1.d+7/(8065.54093*(x(i)-wei(i)))
c           wnm=abs(wnm2-wnm1)
c           pref=yy*9.145d+1*xnm/wnm
c           xxx=(e(j)-xnm)/wnm
c           r1(j)=r1(j)+pref*(xintx2(xxx)-0.5*wnm/(e(j)+xnm))
c           r2(j)=r2(j)+9.145d+1*yy*xnm**2/(e(j)**2-xnm**2) 
c           rdebbohr=yy/92.73
c           r2(j)=r2(j)+8.48e+3*rdebbohr*xnm**2/(e(j)**2-xnm**2) 
c     write(*,*) '*****************************'
c     write(*,*) ' K R O N I G - K R A M E R S'
c     write(*,*) '    TRANSFORMATION (delta)'
c     write(*,*) '*****************************'
c     write(*,*) 
c     write(*,*) ' lambda  eV   phi          alp' 
c     do j=1,5
c        write(*,142) e(j),1.d+7/(8065.54093*e(j)),r2(j),
c    .r2(j)*100/xmass,r1(j)/(r2(j)*100/xmass)
c     enddo

 142  format(f6.1,f6.2,3f12.3)

742   continue

      deallocate(r,x,y,iff,lab,e,wei)

      end


************************************************************************
*
************************************************************************

      subroutine dsort(lab,l,x,w,ew)
      implicit real*8 (a-h,o-z)

      dimension l(*),x(*),w(*),ew(*)


c
c     sortieren
c

      do 140   ii = 2, lab
         i = ii - 1
         k = i
         pp= ew(i)
         do 120   j = ii, lab
            if (ew(j) .gt. pp) go to 120
            k = j
            pp= ew(j)
  120    continue
         if (k .eq. i) go to 140
         ew(k) = ew(i)
         ew(i) = pp

         hilf=x(i)
         x(i)=x(k)
         x(k)=hilf

         hilf=w(i)
         w(i)=w(k)
         w(k)=hilf

         ihilf=l(i)
         l(i)=l(k)
         l(k)=ihilf

  140 continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C exp(-x^2)*int_0^x exp(y^2)dy
C s.grimme, aug 1996
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc


      real*8 function xintx2(xli)
      real*8 st, sum, x
      real*8 xli, xen           
      integer n

      xen=abs(xli)

      if(xen.gt.5.0d0) then
      xintx2=0.12499*(4.0d0/xli) + 0.00392*(4.0d0/xen)**3 
     .                           + 0.00031*(4.0d0/xen)**5
      return
      else
      n=1000

      st=xen/n  

      sum=0.0

      x=0.0

      do i=1,n

         x=x+st

         sum=sum+exp((x-0.5*st)**2)*st

      enddo

      xintx2=sum*exp(-xen**2)

      if(xli.lt.0) xintx2=-1.0*xintx2

      return

      endif
 
      end

      
c     implicit real*8 (a-h,o-z)
c     x=0.0
c     y=0.0
c     i=1
c 10  call xintx2f(y,xint)
c     write(*,'(7x,''a('',i3,'')='',3D16.8)') i,x,xintx2(x),xint
c 10  write(*,'(3D16.8)') x,xintx2(x)
c     x=x+0.100
c     y=y+0.272
c     i=i+1
c     if(x.le.12) goto 10
c     end

C     *****************************************************************         
                                                                                
      SUBROUTINE READL(A1,X,N)                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      CHARACTER*(*) A1                                                      
      DIMENSION X(*)                                                            
      I=0                                                                       
      IS=1                                                                      
  10  I=I+1                                                                     
      X(I)=READAA(A1,IS,IB,IE)                                               
      IF(IB.GT.0 .AND. IE.GT.0) THEN                                            
                                IS=IE                                           
                                GOTO 10                                         
      ENDIF                                                                     
      N=I-1                                                                     
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      FUNCTION READAA(A,ISTART,IEND,IEND2)                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 READAA                                                             
      CHARACTER*(*) A                                                      
      NINE=ICHAR('9')                                                           
      IZERO=ICHAR('0')                                                          
      MINUS=ICHAR('-')                                                          
      IDOT=ICHAR('.')                                                           
      ND=ICHAR('D')                                                             
      NE=ICHAR('E')                                                             
      IBL=ICHAR(' ')                                                            
      IEND=0                                                                    
      IEND2=0                                                                   
      IDIG=0                                                                    
      C1=0                                                                      
      C2=0                                                                      
      ONE=1.D0                                                                  
      X = 1.D0                                                                  
      NL=LEN(A) 
      DO 10 J=ISTART,NL-1                                                       
         N=ICHAR(A(J:J))                                                          
         M=ICHAR(A(J+1:J+1)) 
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO                            
     1 .OR. M.EQ.IDOT)) GOTO 20                                                 
   10 CONTINUE                                                                  
      READAA=0.D0                                                               
      RETURN                                                                    
   20 CONTINUE                                                                  
      IEND=J                                                                    
      DO 30 I=J,NL                                                              
         N=ICHAR(A(I:I))                                                          
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C1=C1*10+N-IZERO                                                    
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                                     
            ONE=-1.D0                                                           
         ELSEIF(N.EQ.IDOT) THEN                                                 
            GOTO 40                                                             
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
      IDIG=0                                                                    
      DO 50 II=I+1,NL                                                           
         N=ICHAR(A(II:II))                                                         
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C2=C2*10+N-IZERO                                                    
            X = X /10                                                           
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                                    
            X=-X                                                                
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   50 CONTINUE                                                                  
C                                                                               
C PUT THE PIECES TOGETHER                                                       
C                                                                               
   60 CONTINUE                                                                  
      READAA= ONE * ( C1 + C2 * X)                                              
      DO 55 J=IEND,NL                                                           
         N=ICHAR(A(J:J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
      RETURN                                                                    
                                                                                
   57 C1=0.0D0                                                                  
      ONE=1.0D0                                                                 
      DO 31 I=J+1,NL                                                            
         N=ICHAR(A(I:I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
   31 CONTINUE                                                                  
   61 CONTINUE                                                                  
   70 READAA=READAA*10**(ONE*C1)                                                
      RETURN                                                                    
      END                                                                       
