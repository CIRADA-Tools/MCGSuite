
      module BasicRanNumGen
      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran0
c
c     This generates a random number using the Park and Miller
c     method of Press et al.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran0(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer IA,IM,IQ,IR,MASK
c
      real AM
      parameter(IA=16897,IM=2147483647, AM=1./IM)
      parameter(IQ=127773, IR=2836, MASK=123459876)

      integer k
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if(idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran1
c
c     This generates a random number using the Park and Miller
c     method with a Bays-Durham shuffle of Press et al.
c
c     Needs a negative idum to initialize
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran1(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer IA,IM,IQ,IR,NTAB,NDIV
c
      real AM,EPS,RNMX
      parameter(IA=16807,IM=2147483647, AM=1./IM)
      parameter(IQ=127773, IR=2836, NTAB=32,NDIV=1+(IM-1)/NTAB)
      parameter(EPS=1.2e-7,RNMX=1.-EPS)

      integer j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv/NTAB*0/, iy/0/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(idum.le.0 .or. iy.eq.0) then
         idum=max(-idum,1)
         do 11 j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if(idum.lt.0) idum=idum+IM
            if(j .le.NTAB) iv(j)=idum
 11      enddo
         iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if(idum .lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      if(ran1 .lt. 0.) print*, 'bug in ran1'
      return
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran2
c
c    This function generates a random number using L'Ecuyer method
c     from Press et al.
c
c     Needs a negative idum to initialize
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran2(idum)
c     
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer IM1, IM2, IMM1, IA1, IA2, IQ1
      integer IQ2,IR1,IR2,NTAB,NDIV
c
      real AM,EPS,RNMX
c
      parameter(IM1=2147483563,IM2=2147483399,AM=1./IM1)
      parameter(IMM1=IM1-1)
      parameter(IA1=40014,IA2=40692, IQ1=53668, IQ2=52774,IR1=12211)
      parameter(IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB)
      parameter(EPS=1.2e-7)
      parameter(RNMX=1.-EPS)
c
      integer idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/,iy/0/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(idum .le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 12 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if(idum.lt.0) idum=idum+IM1
            if(j .le. NTAB) iv(j)=idum
 12      enddo
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ2)-k*IR1
      if(idum .lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if(idum2 .lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy .lt.1) iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      if(ran2 .le. 0.) print*, 'bug in ran2'
      return
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran3
c
c     Another method of generating random numbers from Press et al.
c
c     Needs a negative idum to initialize
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran3(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer MBIG,MSEED,MZ
      real FAC
c
      parameter(MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
c
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(56)
      SAVE iff,inext,inextp,ma
      DATA iff/0/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(idum.lt.0 .or. iff.eq.0) then
         iff=1
         mj=abs(MSEED-abs(idum))
         mj=mod(mj,MBIG)
         ma(55)=mj
         mk=1
         do 13 i=1,54
            ii=mod(21*i,56)
            ma(ii)=mk
            mk=mj-mk
            if(mk .lt. MZ) mk=mk+MBIG
            mj=ma(ii)
 13      enddo
         do 14 k=1,4
            do 15 i=1, 56
               ma(i)=ma(i)-ma(1+mod(i+30,56))
               if(ma(i) .lt. MZ) ma(i)=ma(i)+MBIG
 15         enddo
 14      enddo
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext .eq. 56) inext=1
      inextp=inextp+1
      if(inextp .eq. 56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj .lt. MZ) mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      end function
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Gasdev
c
c     This function calculates a gaussian deviate with a zero mean
c     and a sigma of 1.
c
c     It can be found in Press et al.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function gasdev(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer iset
      real fac, gset, rsq, v1,v2
      save iset, gset
      data iset/0/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (idum.lt. 0) iset=0
      if(iset .eq. 0) then
 1      v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2.+v2**2.
        if(rsq .ge. 1. .or. rsq .eq. 0.) goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      end function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end module

