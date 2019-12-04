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
c**********************************************************************
c
c rlist     : extract a list of numbers from string SCR
c       LS  : length of SCR
c       IOUT: integer array
c       n   : number of integers found
c
c       maxnum is the maximum of numbers to be read
c
c modified by s. grimme, bonn jan. 1993
c
c**********************************************************************


      subroutine rlist(ls,scr,iout,n)
      implicit real*8 (a-h,o-z)
      character*1 scr(*)
c 100 numbers max
      parameter (maxnum=640)
      parameter (mdummy=-9999)

      dimension iout(*)
      dimension idef(maxnum)

c input > min, < max

      data min /0/
      data max /640/

c preset of iout
      do i=1,maxnum
         idef(i)=mdummy
      enddo

      mm=maxnum
      call rdlist(scr,ls,1,igotit,wert,min,max,mm,idef)

      if (igotit.eq.-1) then
          write(*,*) 'SOMETHINGS WRONG IN <RLIST>'
          stop
      endif

      k=0
      do i=1,mm  
         if(idef(i).ne.mdummy) then
            k=k+1
            iout(k)=i
         endif
      enddo

      n=k

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdebbs(scr,n,kint,kreal,kchar,isucc,iwert,rwert,cwert)
      implicit real*8 (a-h,o-z)
c
c     --- this routine reads one integer/real/character variable
c         as specified by the mask kint/kreal/kchar (two values of
c         which have to be zero) from the input string <scr>
c         it is assumed that integer/real numbers do not need
c         more than 32 characters of space (if wanted, pump that up)
c         if the mask value is positive, the input characters
c         will be cleaned after reading
c
c         NOTE : abs(kchar) = length of the output string <cwert>
c
      integer kint,kreal,kchar
      integer isucc,iwert
      real*8 rwert
      character*(*) cwert,scr
      character*32 readit
      logical prtout
c
      prtout=n.lt.0
c
      if(prtout) write(6,601) scr
  601 format(/,' input string = ',/,a,/)
c
c     --- if read succeeds, isucc contains the length of the string
c
c     --- check input parameters
      if((kint.ne.0.and.kreal.ne.0).or.
     1   (kint.ne.0.and.kchar.ne.0).or.
     2   (kreal.ne.0.and.kchar.ne.0)) stop ' abuse of rdebbs '
c
      iwert=0
      rwert=0.d0
      cwert=' '
c
      isucc=iblank(scr)-1
c
      if(isucc.eq.0) return
c
      if((kint.ne.0.or.kreal.ne.0).and.isucc.gt.32) then
        write(6,901)
  901   format(/,' integer/real number with more than 32 digits ',/)
        goto 800
      elseif(kchar.ne.0.and.abs(kchar).lt.isucc) then
        write(6,902) iabs(kchar)
  902       format(/,' i/o-error : input string is longer than ',i3,
     1             ' characters ',/)
        goto 800
      endif
c
      nbl=32-isucc
      if(kint.ne.0) then
        iact=kint
        if(nbl.gt.0) readit(1:nbl)=' '
        readit(nbl+1:32)=scr(1:isucc)
        read(readit,'(i32)',err=701) iwert
      elseif(kreal.ne.0) then
        iact=kreal
        if(nbl.gt.0) readit(1:nbl)=' '
        readit(nbl+1:32)=scr(1:isucc)
        read(readit,'(g32.0)',err=702) rwert
      elseif(kchar.ne.0) then
        iact=kchar
        cwert=scr(1:isucc)
      endif
c     --- delete the string if iact>0 (iact<0 allows second read)
      if(iact.gt.0) scr(1:isucc)=' '
c
      return
c
  701 write(6,751)
  751 format(/,' i/o-error : input variable is not integer ',/)
      goto 800
  702 write(6,752)
  752 format(/,' i/o-error : input variable is not real ',/)
c
  800 write(6,991) scr
  991 format(/,
     1 ' WARNING : <rdebbs> could not read properly from string ',/,a,/)
      isucc=-1
c
      return
      end
c
      function iblank(zeile)
      implicit real*8 (a-h,o-z)
c
      integer iblank
      character*(*) zeile
c
c     --- delete leading blanks and equal signs ("=")
c
      n=len(zeile)
c
      iblank=index(zeile,' ')
      if(iblank.eq.0) then
        iblank=n+1
        return
      endif
c
      i=0
  100 i=i+1
      if(i.gt.n) goto 200
      if(zeile(i:i).eq.' '.or.zeile(i:i).eq.'=') goto 100
  200 i=i-1
      if(i.eq.0) return
c
      do 400 j=1,n-i
  400   zeile(j:j)=zeile(i+j:i+j)
c
      if(i.lt.n) zeile(n-i+1:)=' '
c
      iblank=index(zeile,' ')
c
      return
      end
c

c ======================================================================
c     --- still more subroutines - gee !!!
c ======================================================================
c
      subroutine rdintg(zeile,iwert,iflag,isucc)
      implicit real*8 (a-h,o-z)
c
      character*(*) zeile
      character cdummy
c
      call rdebbs(zeile,len(zeile),iflag,0,0,isucc,iwert,rdummy,cdummy)
c
      return
      end
c
      subroutine rdreal(zeile,rwert,iflag,isucc)
      implicit real*8 (a-h,o-z)
c
      character*(*) zeile
      character cdummy
c
      call rdebbs(zeile,len(zeile),0,iflag,0,isucc,idummy,rwert,cdummy)
c
      return
      end
c
      subroutine rdchar(zeile,cwert,iflag,isucc)
      implicit real*8 (a-h,o-z)
c
      character*(*) zeile,cwert
c
      call rdebbs(zeile,len(zeile),0,0,iflag,isucc,idummy,rdummy,cwert)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdlist(zeile,lzeile,iprint,igotit,wert,
     1                  nmin,nmax,n,idef)
      implicit real*8 (a-h,o-z)
c
c     subroutine reads a list of positive indices from string <zeile>
c     and tries to read a real number enclosed within parentheses
c     at the end of the list
c      igotit = 1     everything has been read
c               0     index list without real has been read
c              -1     something is fucked
c
c     the end of the index list is given by the first word containing
c     a character not contained in the set <allowd>
c
c     the indices have to fall into the region nmin<=i<=nmax ;
c     a warning will be printed if indices have multiple occurrence
c
c     iprint > 1 : print the input string <line>
c
      dimension idef(n)
      character zeile*(*),line*255,scrzl*255,cdummy
c
c     --- characters which may occur within a true index list -> allowd
c
      character allowd*16
      data allowd / '1234567890-' /
c
      if(lzeile.gt.255) stop ' <rdlist> : line too long '
c
      igotit=-1
c
      locpl=0
      locpr=0
c
c     --- input string <zeile> will be copied to string <scrzl> deleting
c         all blanks; <scrzl> will contain only the list of indices ---
c
      scrzl=' '
      nallow=0
c
      i=0
      j=0
      iblrem=0
  100 i=i+1
      if(i.gt.lzeile) goto 109
      if(zeile(i:i).eq.' ') then
        iblrem=1
        goto 100
      elseif(zeile(i:i).eq.',') then
        iblrem=0
        j=j+1
        scrzl(j:j)=' '
        goto 100
      elseif(zeile(i:i).eq.'-') then
        iblrem=0
        j=j+1
        scrzl(j:j)='-'
        goto 100
      endif
      if(index(allowd,zeile(i:i)).gt.0) then
        if(iblrem.eq.1) then
          iblrem=0
          j=j+1
          scrzl(j:j)=' '
        endif
        j=j+1
        scrzl(j:j)=zeile(i:i)
        goto 100
      endif
      if(zeile(i:i).eq.'(') then
        locpl=i
        if(i.lt.lzeile) locpr=i+index(zeile(i+1:),')')
      elseif(iblrem.eq.0) then
        nallow=1
      endif
  109 iende=j
c
c     --- if the last word is no true integer but a mixture, skip it
c
      if(nallow.eq.1.and.iende.gt.0) then
        j=iende
  110   j=j-1
        if(j.eq.0) goto 129
        if(scrzl(j:j).ne.' ') goto 110
  120   j=j-1
        if(j.eq.0) goto 129
        if(scrzl(j:j).eq.' ') goto 120
  129   iende=j
      endif
c
      if(iende.eq.0) then
        write(6,601)
  601   format(/,' <rdlist> : input index list is empty ',/)
        return
      endif
c
c --- try to read the real number between the parentheses (if any)
c
      igotit=0
      wert=0.d0
c
      if(locpl.gt.0) then
        igotit=-1
        if(locpr.le.locpl+1) then
          write(6,901)
  901     format(/,
     1    ' <rdlist> : expected input of type : ( <real> ) ',/)
        else
          line=zeile(locpl+1:locpr-1)
        endif
c       --- read " integer / integer " or " real "
        islash=index(line,'/')
        if(islash.gt.0) then
          line(islash:islash)=' '
          call rdebbs(line,len(line),1,0,0,isucc,iwert,rdummy,cdummy)
          if(isucc.gt.0) then
            wert=dble(iwert)
            call rdebbs(line,len(line),1,0,0,isucc,iwert,rdummy,cdummy)
            if(isucc.gt.0.and.iwert.ne.0) wert=wert/dble(iwert)
          endif
        else
          call rdebbs(line,len(line),0,1,0,isucc,idummy,rwert,cdummy)
          if(isucc.gt.0) wert=rwert
        endif
        if(isucc.gt.0) igotit=1
      endif
c
c     --- now we are ready to digest the index list
c
      line=scrzl(1:iende)
c
c     --- working area
c
  200 call wisch(line,80)
      if(iprint.gt.1) write(6,611) line
  611 format(' input line = ',/,a)
      if(index(line,' ').eq.1) goto 666
c     --- decide whether index pair or single index has to be read
      call rdebbs(line,len(line),0,0,-len(line),
     .            isucc,idummy,rdummy,scrzl)
      if(iprint.gt.1) write(6,611) scrzl
      ibar=index(scrzl,'-')
      if(ibar.gt.0) then
c       --- read index pair
        kbar=index(line,'-')
        if(kbar.eq.1) then
          write(6,9015)
 9015     format(/,' <rdlist> : illegal use of - ',/)
          igotit=-1
          goto 666
        endif
        line(kbar:kbar)=' '
        call rdebbs(line,len(line),1,0,0,isucc,iwert,rdummy,cdummy)
        if(iwert.lt.nmin.or.iwert.gt.nmax) then
          write(6,902) iwert,nmin,nmax
  902     format(/,' <rdlist> : index = ',i5,
     1             ' lies not within (',i5,',',i5,')',/)
          igotit=-1
          goto 666
        endif
        ist=iwert
        call rdebbs(line,len(line),1,0,0,isucc,iwert,rdummy,cdummy)
        if(isucc.eq.0) then
          write(6,9015)
          igotit=-1
          goto 666
        endif
        if(iwert.lt.nmin.or.iwert.gt.nmax) then
          write(6,902) iwert,nmin,nmax
          igotit=-1
          goto 666
        endif
        if(iwert.lt.ist) then
          write(6,903)
  903     format(/,' <rdlist> : illegal order of indices ',/)
          igotit=-1
          goto 666
        endif
        iend=iwert
        do 250 i=ist,iend
          if(idef(i).eq.1) write(6,904) i
  904     format(/,' <rdlist> : multiple use of index ',i5,/)
          idef(i)=1
  250     continue
      else
c       --- read single index
        call rdebbs(line,len(line),1,0,0,isucc,iwert,rdummy,cdummy)
        if(iwert.lt.nmin.or.iwert.gt.nmax) then
          write(6,902) iwert,nmin,nmax
          igotit=-1
          goto 666
        endif
        if(idef(iwert).eq.1) write(6,904) iwert
        idef(iwert)=1
      endif
      goto 200
c
c     --- end of working area
c
  666 continue
      zeile(1:lzeile)=' '
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wisch(line,n)
      implicit real*8 (a-h,o-z)

c ----------------------------------------------------------------------
c       delete leading blanks in string <line>
c       19.09.90  C.M.K.
c ----------------------------------------------------------------------

      character*(*) line

      m=len(line)

      if(line(1:1).ne.' ') return

      i=m+1
c     --- last non-blank
    1 i=i-1
      if(i.gt.0.and.line(i:i).eq.' ') goto 1
      if(i.eq.0) return
c     --- first non-blank
      j=0
    2 j=j+1
      if(j.lt.i.and.line(j:j).eq.' ') goto 2
      l=j
c     --- shift
      k=1
    3 line(k:k)=line(j:j)
      k=k+1
      j=j+1
      if(j.le.i) goto 3
c     --- fill up with blanks
      if(k.gt.l) l=k
    4 line(l:l)=' '
      l=l+1
      if(l.le.i) goto 4

      return
      end
c ---------------------------------------------------------------------
      subroutine delchr(line,chrset)
      implicit real*8 (a-h,o-z)

c ----------------------------------------------------------------------
c |     replace all characters within <line> occurring in string       |
c |     <chrset> by blanks and nothing else                            |
c ----------------------------------------------------------------------

      character*(*) line,chrset

      n=len(line)
      do 100 i=1,n
        if(index(chrset,line(i:i)).gt.0) line(i:i)=' '
  100   continue

      return
      end

      subroutine delmin(line,chrset)
      implicit real*8 (a-h,o-z)

c ----------------------------------------------------------------------
c |     replace all characters within <line> occurring in string       |
c |     <chrset> by blanks if they are neighboured by letters          |
c ----------------------------------------------------------------------

      character*(*) line,chrset

      character abc*52

      data abc /'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      n=len(line)
      if(n.eq.0) return

      if(n.eq.1) then
        if(index(chrset,line(1:1)).gt.0) line(1:1)=' '
        return
      endif

      do 100 i=1,n
        if(index(chrset,line(i:i)).gt.0) then
          if(i.gt.1.and.index(abc,line(i-1:i-1)).gt.0) then
            line(i:i)=' '
          elseif(i.lt.n.and.(index(abc,line(i+1:i+1)).gt.0.or.
     1                       line(i+1:i+1).eq.' ')) then
            line(i:i)=' '
          endif
        endif
  100   continue

      return
      end

C     *****************************************************************

      SUBROUTINE READL(NL,A1,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 A1(*)
      DIMENSION X(*)
      I=0
      IS=1
  10  I=I+1
      X(I)=XREAD(NL,A1,IS,IB,IE)
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END

C     *****************************************************************

      SUBROUTINE IREADL(NL,A1,IX,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 A1(*)
      DIMENSION IX(*)
      I=0
      IS=1
  10  I=I+1
      IX(I)=IDINT(XREAD(NL,A1,IS,IB,IE))
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END

C     *****************************************************************


      FUNCTION XREAD(NL,A,ISTART,IEND,IEND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XREAD

      CHARACTER*1 A(NL)
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
      DO 10 J=ISTART,NL-1
         N=ICHAR(A(J))
         M=ICHAR(A(J+1))
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO
     1 .OR. M.EQ.IDOT)) GOTO 20
   10 CONTINUE
      XREAD=0.D0
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I))
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
         N=ICHAR(A(II))
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
      XREAD= ONE * ( C1 + C2 * X)
      DO 55 J=IEND,NL
         N=ICHAR(A(J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN

   57 C1=0.0D0
      ONE=1.0D0
      DO 31 I=J+1,NL
         N=ICHAR(A(I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
         IF(N.EQ.MINUS)ONE=-1.0D0
   31 CONTINUE
   61 CONTINUE
   70 XREAD=XREAD*10**(ONE*C1)
      RETURN
      END

