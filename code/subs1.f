
      Subroutine Ttmean(isteps)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
      complex*16 zmnt(is,js),dmnt(is,js)
      common/tmean/zmnt,dmnt
      complex*16 zmn1(is,js),dmn1(is,js)
      common/time1/zmn1,dmn1
      rsteps=dfloat(isteps)
      do 1 i=1,is
      do 1 j=1,js
           zmnt(i,j)=zmnt(i,j)+zmn1(i,j)/rsteps
           dmnt(i,j)=dmnt(i,j)+dmn1(i,j)/rsteps
1     continue
      return
      end
c
c
c
      Subroutine Whist
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
      complex*16 zmn1(is,js),dmn1(is,js)
      complex*16 zmn2(is,js),dmn2(is,js)
      common/time1/zmn1,dmn1
      common/time2/zmn2,dmn2
      write(23)zmn1
      write(23)dmn1

      write(23)zmn2
      write(23)dmn2

      return
      end

      Subroutine Restart
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
      complex*16 zmn1(is,js),dmn1(is,js)
      complex*16 zmn2(is,js),dmn2(is,js)
      common/time1/zmn1,dmn1
      common/time2/zmn2,dmn2

      read(24)zmn1
      read(24)dmn1

      read(24)zmn2
      read(24)dmn2

      return
      end
c
c
c
c
c

      Subroutine Gauslat(np2,glts,gwts)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
C  CALCULATE SINE OF GAUSSIAN LATITUDES.  ALSO CALCULATE GAUSSIAN
C  WEIGHTS FOR GAUSSIAN QUADRATURE DONE LATER IN VERIFICATION ROUTINE
C
      implicit double precision (a-h,o-z)
      dimension gwts(np2),glts(np2)
      dimension radcol(jx2),gw(jx2)
C
      eps=1.0d-15
      one=1.0d0
      pi=(4.0d0)*datan(one)
      k2=2*np2
      rk2=k2
      scale=2.0d0/(rk2**2)
      k1=k2-1
      dradz=pi/360.0d0
      rad=0.0d0
C
      do 1000 k=1,np2
      iter=0
      drad=dradz
1     CALL POLY2(K2,RAD,P2)
2     p1=p2
      iter=iter+1
      rad=rad+drad
      CALL POLY2(K2,RAD,P2)
      if (sign(one,p1).eq.sign(one,p2)) go to 2
      if (drad.lt.eps) go to 3
      rad=rad-drad
      drad=drad*0.25d0
      go to 1
3     continue
      radcol(k)=rad
      CALL POLY2(K1,RAD,P1)
      x=dcos(rad)
      w=scale*(one-x*x)/(p1*p1)
      gw(k)=w
      CALL POLY2(K2,RAD,P1)
      glts(np2+1-k)=x
      gwts(np2+1-k)=gw(k)
1000  continue
      return
      end
C**********************************************   POLY2   *************
C
C
      Subroutine Poly2(n,rad,p)
C
C  CALLED BY ROUTINE GAUSLAT
C
      implicit double precision (a-h,o-z)
      x=dcos(rad)
      y1=1.0d0
      y2=x
      do 1 i=2,n
      g=x*y2
      y3=g-y1+g-(g-y1)/dble(i)
      y1=y2
      y2=y3
1     continue
      p=y3
      return
      end

      Subroutine Pmns(mphys)
C
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
      common/aleg/pmn(is,jsp),hmn(is,js),pa(is,jsp),ha(is,js),
     1   glat(jx2),gwgt(jx2),gw,y,cs,cs2,acs2
      double precision pmn,hmn,pa,ha,glat,gwgt,gw,y,cs,cs2,acs2
      common/cconst/ms,ns,mp,np,np2
      dimension a1(is),a2(is)
C
      ns1=ns+1
      pmn(1,1)=dsqrt(0.5d0)
C
      do 100 mm=2,ms
      am=dble(mm)-1.0d0
      a1(mm)=cs*dsqrt(1.0d0+1.0d0/(2.0d0*am))
100   continue
      do 101 mm=2,ms
101   pmn(mm,1)=a1(mm)*pmn(mm-1,1)
      do 200 mm=1,ms
      am=dble(mm)-1.0d0
      a1(mm)=y*dsqrt(2.0d0*am+3.0d0)
200   continue
      do 201 mm=1,ms
201   pmn(mm,2)=a1(mm)*pmn(mm,1)
C
      do 400 nj=3,ns1
      do 300 mm=1,ms
      am=dble(mm)-1.0d0
      an=dble(nj+mm)-2.0d0
      anm=an*an-am*am
      a1(mm)=y*dsqrt((4.0d0*an*an-1.0d0)/anm)
300   a2(mm)=dsqrt((2.0d0*an+1.0d0)*(anm-2.0d0*an+1.0d0)
     1  /((2.0d0*an-3.0d0)*anm))
      do 301 mm=1,ms
301   pmn(mm,nj)=a1(mm)*pmn(mm,nj-1)-a2(mm)*pmn(mm,nj-2)
400   continue
C
      if(mphys.ne.0)go to 999
      hmn(1,1)=0.0d0
C
      do 500 mm=2,ms
      am=dble(mm)-1.0d0
500   a1(mm)=dsqrt(am*am*(2.0d0*am+1.0d0)/
     1(4.0d0*am*am+8.0d0*am+3.0d0))
      do 501 mm=2,ms
501   hmn(mm,1)=a1(mm)*pmn(mm,2)
C
      do 650 nj=2,ns
      do 600 mm=1,ms
      an=dble(nj+mm)-2.0d0
      am=dble(mm)-1.0d0
      anm=an*an-am*am
      a1(mm)=an*dsqrt((anm+2.0d0*an+1.0d0)/
     1(4.0d0*an*an+8.0d0*an+3.0d0))
600   a2(mm)=(an+1.)*dsqrt(anm/(4.0d0*an*an-1.0d0))
      do 601 mm=1,ms
601   hmn(mm,nj)=a1(mm)*pmn(mm,nj+1)-a2(mm)*pmn(mm,nj-1)
650   continue
C
      do 800 nj=1,ns
      rnn1=dble(nj)*(dble(nj)-1.0d0)
      if(nj.eq.1)rnn1=1.0d0
      rni=1.0d0/rnn1
      pa(1,nj)=pmn(1,nj)*rni
      ha(1,nj)=hmn(1,nj)*rni
      do 700 mm=2,ms
      an=dble(nj+mm)-2.0d0
700   a1(mm)=1.0d0/(an*(an+1.0d0))
      do 701 mm=2,ms
      pa(mm,nj)=pmn(mm,nj)*a1(mm)
701   ha(mm,nj)=hmn(mm,nj)*a1(mm)
800   continue
C
      pa(1,1)=0.0d0
      ha(1,1)=0.0d0
C
999   continue
C
      return
      end
C
C
C
      Subroutine Bscst
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 ci
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
      common/tconst/dt,alpha,aee,andree

C
C  SET PHYSICAL CONSTANTS
C
      ae=6.371d6
      r=287.05d0
      cp=1005.d0
      grav=9.8d0
      rkappa=r/cp
      ci=(0.0d0,1.0d0)
C  FOR ROTATION OF WIND IN B.L.
C    USES ROTATION OF 20 DEG
      pi9=4.d0*datan(1.d0)/(9.d0)
      rtcs=dcos(pi9)
      rtsn=dsin(pi9)
C  CONST FOR EVAP FROM OCEANS
      epslr=.622d0*2.51d6/r
      escons=.611d0*dexp(epslr/273.16d0)
      esc1=.622d0*escons
      esc2=-.378d0*escons
C  SET COMPUTATIONAL CONSTANTS
C    RESOLUTION
      ms=is
      ns=js
      mp=ix
      np=jx
      np2=np/2
C
      ps=1.0d0
      aee=ae*ae
      rkap1=1.0d0+rkappa
C
C
C
C
      return
      end
C
C
C
C
C
      Subroutine Legsum(fmn,fmnh,fmsh,pmn,ms,ns,ne,no)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 fmn(is,js),fmnh(is),fmsh(is),fme(is),fmo(is)
      double precision pmn(is,jsp)
C
C  SUBROUTINE DOES INVERSE LEGENDRE TRANSFORM BY SUMMING
C    SPHERICAL HARMONIC COEFF OVER N
C  INPUT% FMN SPHER HARM COEFF (ONLY ONE LAYER)
C         PMN CAN BE PMN OR HMN AT N. HEMIS LAT +Y
C         MS NUMBER OF ZONAL WAVES= MAX(M+1)
C         NS NUMBER OF SPECT DEGREES = MAX(N-M+1)
C         NE BEGINNING INDEX FOR SYMMETRIC SUM
C                =1 FOR PMN,  =2 FOR HMN
C        NO BEGINNING INDEX FOR ANTISYM SUM
C               =2 FOR PMN, =1 FOR HMN
C
C  OUTPUT FMNH, FMSH FOURIER COEFF FOR ALL M
C      AT N. HEM +Y  AND AT S. HEM -Y
C
      do 50 mm=1,ms
      fme(mm)=(0.d0,0.d0)
   50 fmo(mm)=(0.d0,0.d0)
      do 100 nj=ne,ns,2
      do 100 mm=1,ms
  100 fme(mm)=fme(mm)+fmn(mm,nj)*(pmn(mm,nj))
      do 200 nj=no,ns,2
      do 200 mm=1,ms
  200 fmo(mm)=fmo(mm)+fmn(mm,nj)*(pmn(mm,nj))
      do 300 mm=1,ms
      fmnh(mm)=fme(mm)+fmo(mm)
      fmsh(mm)=fme(mm)-fmo(mm)
  300 continue
      fmnh(1)=dcmplx(dble(fmnh(1)),0.0d0)
      fmsh(1)=dcmplx(dble(fmsh(1)),0.0d0)
C
      return
      end
C
C
C
      Subroutine Gql(fmnh,fmsh,fmn,pmn,gw,ms,ns)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 fmnh(is),fmsh(is),fmn(is,js)
      double precision pmn(is,jsp),gw
C
C  COMPUTES AND ACCUMULATES CONTRIBUTION OF FMNH(MM)
C     AND FMSH(MM)  TO GAUSSIAN QUADRATURE
C
C  SYMMETRIC PART
      do 50 nj=1,ns,2
      do 50 mm=1,ms
   50 fmn(mm,nj)=fmn(mm,nj)+(pmn(mm,nj)*gw)
     1*(fmnh(mm)+fmsh(mm))
C  ANTISYMMETRIC PART
      do 60 nj=2,ns,2
      do 60 mm=1,ms
   60 fmn(mm,nj)=fmn(mm,nj)+(pmn(mm,nj)*gw)
     1*(fmnh(mm)-fmsh(mm))
C
      return
      end
C
C
C
      Subroutine Lmn(xmnh,xmsh,ymnh,ymsh,dxymn,pmn,hmn,
     1  gw,acs2,ms,ns)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 xmnh(is),xmsh(is),ymnh(is),ymsh(is),dxymn(is,js),
     1  ci,cim
      double precision pmn(is,jsp),hmn(is,js),gw,acs2
C
C  COMPUTES AND ACCUMULATES CONTRIBUTIONS OF
C    XMNH(MM),XMSH(MM) AND YMNH(MM),YMSH(MM) TO THE
C    GAUSSIAN QUADRATURE IN OPERATOR LMN
C
      ci=(0.d0,1.d0)
      gwac=(gw/acs2)
      do 300 mm=1,ms
      cim=ci*(dble(mm)-1.0d0)
C  SYMMETRIC PART
      do 100 nj=1,ns,2
      dxymn(mm,nj)=dxymn(mm,nj)+gwac*(cim*(pmn(mm,nj))*
     1   (xmnh(mm)+xmsh(mm))
     2   +(hmn(mm,nj))*(ymnh(mm)-ymsh(mm)))
  100 continue
C
C  ANTISYMMETRIC PART
      do 200 nj=2,ns,2
      dxymn(mm,nj)=dxymn(mm,nj)+gwac*(cim*(pmn(mm,nj))*
     1   (xmnh(mm)-xmsh(mm))
     2   +(hmn(mm,nj))*(ymnh(mm)+ymsh(mm)))
  200 continue
  300 continue
C
      return
      end
C
C
C
      Subroutine Umvm
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 zmn2(is,js),dmn2(is,js)
      complex*16 zm(is,2),dm(is,2),ci,cim
      complex*16 tm(is,2),wm(is,2)
      common/time2/zmn2,dmn2
      common/fours/zm,dm,tm,wm
      common/aleg/pmn(is,jsp),hmn(is,js),pa(is,jsp),ha(is,js),
     1  glat(jx2),gwgt(jx2),gw,y,cs,cs2,acs2
      double precision pmn,hmn,pa,ha,glat,gwgt,gw,y,cs,cs2,acs2
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
C
C  SUBROUTINE COMPUTES UM,VM DIRECTLY FROM ZMN,DMN
C  INPUT ZMN2,DMN2,PA,HA,MS,NS FROM COMMON
C
C  OUTPUT UM AND VM AT +/- Y AT LAYER K
C    UM HELD IN ZM, VM HELD IN DM
C
      ci=(0.d0,1.d0)
      zmn2(1,2)=zmn2(1,2)-1.190797402d-4
      CALL LEGSUM(ZMN2(1,1),ZM(1,1),ZM(1,2),HA,MS,NS,2,1)
      CALL LEGSUM(DMN2(1,1),DM(1,1),DM(1,2),PA,MS,NS,1,2)
      do 160 mm=1,ms
      cim=ci*(dble(mm)-1.0d0)
      do 160 ihem=1,2
  160 zm(mm,ihem)=-ae*(zm(mm,ihem)+cim*dm(mm,ihem))
      CALL LEGSUM(ZMN2(1,1),DM(1,1),DM(1,2),PA,MS,NS,1,2)
      CALL LEGSUM(DMN2(1,1),TM(1,1),TM(1,2),HA,MS,NS,2,1)
      do 250 mm=1,ms
      cim=ci*(dble(mm)-1.0d0)
      do 250 ihem=1,2
  250 dm(mm,ihem)=-ae*(dm(mm,ihem)*cim-tm(mm,ihem))
      do 260 ihem=1,2
      zm(1,ihem)=dcmplx(dble(zm(1,ihem)),0.0d0)
  260 dm(1,ihem)=dcmplx(dble(dm(1,ihem)),0.0d0)
      zmn2(1,2)=zmn2(1,2)+1.190797402d-4
      return
      end
C
C
C
      Subroutine Fft1(f,w,ms,mp,num,idir)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      dimension f(ix),w(ix)
      common/fft/ifax(10),trigs(ix,2)
C
C  SUBROUTINE DOES MULTIPLE HALF COMPLEX*16 FFTS
C    SYNTHESIS (SPECTRAL TO PHYSICAL) FOR IDIR=+1
C    ANALYSIS (PHYSICAL TO SPECTRAL) FOR IDIR=-1
C
C  USES TEMPERTON'S  ROUTINES DEVELOPED AT ECMWF
C
C  ON INPUT
C
C     F(MP,NUM)   ARRAY CONTAINING *NUM* FIELDS EACH OF LENGTH *MP*
C                  TO BE TRANSFORMED -- FOR SYNTHESIS THEY MUST BE
C                  ARRANGED AS  REAL,IMAG,REAL,IMAG... FOR L=1,MS AND
C                  SET TO ZERO FOR L=MS+1,MP (DONE BELOW IN LOOP 200)
C
C      WORK(MP,NUM)  IS A WORK ARRAY
C
C      MS            NUMBER OF WAVES (INCLUDING WAVENUMBER 0) MUST BE .LT. MP
C
C      MP            LENGTH OF EACH TRANSFORM (NUMBER OF GRID POINTS)
C
C      NUM           NUMBER OF TRANSFORMS TO BE PERFORMED SIMULTANEOUSLY
C
C      IDIR          DIRECTION OF TRANSFORM (+1 FOR SYNTHESIS, -1 FOR ANALYSIS)
C
C  ON OUTPUT
C
C      F(MP,NUM)      CONTAINS TRANFORMED FIELDS
C
C  ON FIRST CALL TO FFT1, MUST SET UP  FACTORS AND TRIG FUNCTIONS
C
      data itst/0/
      if(itst.eq.1)go to 50
C
C   ON FIRST CALL MUST SET UP FACTORS AND TRIG FUNCTIONS
C
      itst=1
      CALL FAX(IFAX,MP,3)
      CALL FFTRIG(TRIGS,MP,3)
C
   50 if(idir.eq.-1)go to 250
C
C   FOR SYNTHESIS MUST ZERO FILL WAVENUMBERS GT MS
C
      llo=2*ms+1
      lhi=mp
      do 200 i=1,num
      do 100 lo=llo,lhi
  100 f(lo)=0.0d0
      llo=llo+mp
      lhi=lhi+mp
  200 continue
C
C   ANALYSIS REQUIRES NO PREPROCESSING
C
  250 CALL FFT99M(F,W,TRIGS,IFAX,1,MP,MP,NUM,IDIR)
      return
      end
C
C
C
      Subroutine Sptogp
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 zmn2(is,js),dmn2(is,js)
      complex*16 zm(is,2),dm(is,2),tm(is,2),wm(is,2),ci,cim
      common/time2/zmn2,dmn2
      common/fours/zm,dm,tm,wm
      common/grdpt/z(ix,2),d(ix,2),
     1  U(ix,2),V(ix,2)
      common/aleg/pmn(is,jsp),hmn(is,js),pa(is,jsp),ha(is,js),
     1  glat(jx2),gwgt(jx2),gw,y,cs,cs2,acs2
      double precision pmn,hmn,pa,ha,glat,gwgt,gw,y,cs,cs2,acs2
      common/fftw/work(ix,4)
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
C
C
C  SUBROUTINE TRANSFORMS DEPENDENT VARIABLES FROM SPECTRAL
C   (SPHERICAL HARMONIC) SPACE TO GRID SPACE
C   OUTPUT OF SUBROUTINE CONSISTS OF GRID POINT VALUES OF
C   Z,D,U,V,AROUND LAT +/- Y
C  ASSUMES PMN,HMN,PA,HA,GW,Y  AT LA ARE AVAILABLE THROUGH COMMON
C
C  CREATE GRID POINT VALUES OF Z,D  AT +/- Y
C
C
C   LEGENDRE TRANSFORM DONE ONE LAYER AT A TIME
C
      CALL LEGSUM(ZMN2(1,1),Z(1,1),Z(1,2),PMN,MS,NS,1,2)
      CALL LEGSUM(DMN2(1,1),D(1,1),D(1,2),PMN,MS,NS,1,2)

C
C  CREATE FOURIER COEFFICIENTS VALUES OF U,V  AT +/-Y
C  UM HELD IN ZM, VM HELD IN DM
C
      CALL UMVM
      do 200 ihem=1,2
      DO 200 mm=1,ms
      u(2*mm-1,ihem)=dble(zm(mm,ihem))
      u(2*mm,ihem)=dimag(zm(mm,ihem))
      v(2*mm-1,ihem)=dble(dm(mm,ihem))
  200 v(2*mm,ihem)=dimag(dm(mm,ihem))
C
C
C  FOURIER TRANSFORMS (NOTE THAT ORDER OF COMMON BLOCK /GRDPT/
C    IS CRUCIAL FOR PROPER FUNCTIONING OF FFT)
C
C    TRANSFORM Z,D
      CALL FFT1(Z,WORK,MS,MP,4,+1)
C
C   TRANSFORM U,V
      CALL FFT1(U,WORK,MS,MP,4,+1)
C
C  END OF SPHER HARM TO GRID POINT TRANSFORM
C
      return
      end
C
C
C
      Subroutine Gptosp
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 dzdt(is,js),dddt(is,js)
      complex*16 zm(is,2),dm(is,2),tm(is,2),wm(is,2),ci
      common/timed/dzdt,dddt
      common/fours/zm,dm,tm,wm
      common/prod/a(ix,2),b(ix,2)
      common/aleg/pmn(is,jsp),hmn(is,js),pa(is,jsp),ha(is,js),
     1  glat(jx2),gwgt(jx2),gw,y,cs,cs2,acs2
      double precision pmn,hmn,pa,ha,glat,gwgt,gw,y,cs,cs2,acs2
      common/fftw/work(ix,4)
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
C
C  SUBROUTINE COMPUTES CONTRIBUTIONS OF NONLINEAR TERMS
C   AT +/- Y TO THE GAUSSIAN QUADRATURES FOR THE SPECTRAL
C   TIME DERIVATIVES
C
C
C    FOURIER TRANSFORMS DONE FOR ALL LAYERS AND BOTH HEMISPHERES
C      SIMULTANEOUSLY (NOTE THAT ORDER OF COMMON BLOCK /PROD/ IS
C      IMPORTANT)
C
C
C   TRANSFORM A,B
      CALL FFT1(A,WORK,MS,MP,4,-1)
C
C    LEGENDRE TRANSFORMS
C
C
      ms2=2*ms
C
C  FOR VORTICITY EQUATION
      
      CALL LMN(A(1,1),A(1,2),B(1,1),B(1,2),DZDT(1,1),
     1  PMN,HMN,GW,ACS2,MS,NS)
C
      return
      end
C
C
C
      Subroutine Laloop
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 zmn2(is,js),dmn2(is,js)
      complex*16 dzdt(is,js),dddt(is,js)
      complex*16 zm(is,2),dm(is,2),tm(is,2),wm(is,2),ci
      common/time2/zmn2,dmn2
      common/timed/dzdt,dddt
      common/fours/zm,dm,tm,wm
      common/aleg/pmn(is,jsp),hmn(is,js),pa(is,jsp),ha(is,js),
     1  glat(jx2),gwgt(jx2),gw,y,cs,cs2,acs2
      double precision pmn,hmn,pa,ha,glat,gwgt,gw,y,cs,cs2,acs2
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
C  LATITUDE LOOP - COMPUTES NONLINEAR TERMS IN SPECTRAL PROGNOSTIC
C        EQUATIONS
C
C  OUTPUT OF SUBROUTINE __  CONTRIBUTIONS OF NONLINEAR TERMS
C   TO THE SPECTRAL TIME DERIVATIVES OF Z
C   FOR ALL MM,NJ,K
C
      ci=(0.0d0,1.0d0)
C
C  ASSUMES SPHER HARM COEFF AVAILABLE THROUGH COMMON
C
C
C  ZEROING ARRAYS TO HOLD SPECTRAL TIME DERIVATIVES
C   FOR ACCUMULATING GAUSSIAN QUADRATURE
C
      do 50 mm=1,ms
      do 50 nj=1,ns
      dzdt(mm,nj)=(0.0d0,0.0d0)
   50 continue
C
***
C********** BEGIN LATITUDE LOOP
C**
      do 6000 la=1,np2
      y=glat(la)
      gw=gwgt(la)
      cs2=1.0d0-y*y
      cs=dsqrt(cs2)
      acs2=ae*cs2
C
C  COMPUTE PMN, HMN, PA, AND HA
C
      mphys=0
      CALL PMNS(MPHYS)
      CALL SPTOGP
      CALL NLPROD(LA)
      CALL GPTOSP
 6000 continue
C**
C******* END LATITUDE LOOP
C
C  CHANGE APPROPRIATE SIGNS SO THAT SUBROUTINE WILL RETURN
C  DZDT=-LMN(A,B)
C
      do 6100 nj=1,ns
      do 6100 mm=1,ms
      dzdt(mm,nj)=-dzdt(mm,nj)
 6100 continue
C
C
      return
      end
C
C
C
      Subroutine Nlprod(la)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 ci
      common/grdpt/z(ix,2),d(ix,2),
     1  u(ix,2),v(ix,2)
      common/prod/a(ix,2),b(ix,2)
      common/aleg/pmn(is,jsp),hmn(is,js),pa(is,jsp),ha(is,js),
     1  glat(jx2),gwgt(jx2),gw,y,cs,cs2,acs2
      double precision pmn,hmn,pa,ha,glat,gwgt,gw,y,cs,cs2,acs2
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
C
C  SUBROUTINE COMPUTES NONLINEAR PRODUCTS IN
C   GRID SPACE AROUND LAT CIRCLES +/- Y FOR ALL K
C
C  ASSUMES LAYERS ARE NUMBERED FROM TOP DOWN
C
C  COMPUTATIONS DONE IN VERTICAL COLUMNS AT EACH LONGITUDE
C
C
      do 510 ihem=1,2
      do 500 lo=1,mp
C
C
C     ADVECTION IN VORT
C
C  HORIZONTAL ADVECTION
C
C  FOR VORTICITY AND DIVERGENCE EQUATIONS
      a(lo,ihem)=z(lo,ihem)*u(lo,ihem)
      b(lo,ihem)=z(lo,ihem)*v(lo,ihem)
C
  500 continue
  510 continue
C
      return
      end
C
C
C
      Subroutine Implct(nt,filtc)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
C    SEMI-IMPLICIT TIME STEPPING WITH ROBERT'S TIME FILTER
C    TO GET FILTERED PROGNOSTIC VALUES AT TIME STEP NT, THIS
C    ROUTINE MUST BE CALLED (NT+1) TIMES
C
      common /time1/zmn1,dmn1
      common /time2/zmn2,dmn2
      common /timed/dzdt,dddt
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
      common /tconst/dt,alpha,aee,andree
      complex*16 zmn1(is,js),dmn1(is,js),ci
      complex*16 zmn2(is,js),dmn2(is,js)
      complex*16 zmn3(is,js),dmn3(is,js)
      complex*16 dzdt(is,js),dddt(is,js)
      equivalence (zmn3,dzdt),(dmn3,dddt)
      data iflg/0/
C
      if(nt.gt.1) go to 5
      dt2=dt
      filtc=0.0d0
      go to 10
    5 if(iflg.gt.0) go to 55
      iflg=1
      filtc=andree
      dt2=dt*2.0d0
   10 ccc=alpha*dt2
      cccs=ccc*ccc
      bbb=1.0d0-alpha
      bb1=1.0d0/alpha
C
   55 do 75 nj=1,ns
      do 73 mm=1,ms
      zmn3(mm,nj)=zmn1(mm,nj)+dt2*dzdt(mm,nj)
   73 continue
   75 continue
C
      CALL TFILT(FILTC,ZMN1(1,1),ZMN2(1,1),ZMN3(1,1))
      return
      end
C
C
C
      Subroutine Tfilt(filtc,fmn1,fmn2,fmn3)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 fmn1(is,js),fmn2(is,js),fmn3(is,js)
      common/cconst/ms,ns,mp,np,np2
C
      do 100 nj=1,ns
      do 100 mm=1,Ms
  100 fmn1(mm,nj)=fmn2(mm,nj)+filtc*(fmn1(mm,nj)-2.d0*
     1  fmn2(mm,nj)+fmn3(mm,nj))
C
      CALL SWITCH2(FMN2,FMN3,0.D0)
C
      return
      end
C
C
C
      Subroutine Switch2(ctime1,ctime2,forc)
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 ctime1(is,js),ctime2(is,js),chold
      common/cconst/ms,ns,mp,np,np2
C
      if(forc.eq.0.d0)go to 200
C
      do 100 mm=1,ms
      do 100 nj=1,ns
      chold=ctime1(mm,nj)
      ctime1(mm,nj)=ctime2(mm,nj)
      ctime2(mm,nj)=chold
  100 continue
C
      return
C
C  SWITCH IN PREPARATION FOR FORWARD STEP
C
  200 continue
      do 300 mm=1,ms
      do 300 nj=1,ns
      ctime1(mm,nj)=ctime2(mm,nj)
  300 continue
C
      return
      end
C
C
C
      Subroutine Diffsn
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 zmn1(is,js),dmn1(is,js)
      complex*16 dzdt(is,js),dddt(is,js),ci
      common/time1/zmn1,dmn1
      common/timed/dzdt,dddt
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2
C
C  SUBROUTINE ADDS  -K*DEL**4   TO Z TENDENCIES
C
      a4=ae*ae*ae*ae
      dkha4=dkh/a4
      ekha4=ekh/a4
C
C  FOR Z ADD DIFF ONLY TO UPPER HALF OF RHOMBOID, FOR N.GT.MS-1
C
      do 2000 mm=2,ms
      mm2=mm-2
      n1=ns-mm2
      do 2000 nj=n1,ns
      n=nj+mm-2
      rn=dble(n*(n+1))
      rnn2=rn*rn
      eka4rn=ekha4*rnn2
      dzdt(mm,nj)=dzdt(mm,nj)-eka4rn*zmn1(mm,nj)
2000  continue
C
      return
      end
C
C
C
C
C
	Subroutine Zero1(x)
        parameter(is=41,js=41,jsp=js+1)
        parameter(ix=128,jx=102,jx2=jx/2)
        implicit double precision (a-h,o-z)
	complex*16 x(is,js)
        do 1 i=1,is
        do 1 j=1,js
           x(i,j)=dcmplx(0.0d0,0.0d0)
1	continue
        return
        end

	Subroutine Zero2(x)
        parameter(is=41,js=41,jsp=js+1)
        parameter(ix=128,jx=102,jx2=jx/2)
        implicit double precision (a-h,o-z)
	complex*16 x(is,js)
        do 1 i=1,is
        do 1 j=1,js
           x(i,j)=dcmplx(0.0d0,0.0d0)
1	continue
        return
        end

	Subroutine Wrtdat(z,d,nt,irec,nfile)
        parameter(is=41,js=41,jsp=js+1)
        parameter(ix=128,jx=102,jx2=jx/2)
        implicit double precision (a-h,o-z)
        complex*8 out(is,js)
        complex*16 z(is,js),d(is,js)
        iday=nt/144
        write(6,*)'wrtdat called',iday,nfile
        z(1,2)=z(1,2)-1.190797402D-4
        do 2 i=1,is
        do 2 j=1,js
           out(i,j)=cmplx(z(i,j))
2	continue
        write(6,*)'get to here 1'
        z(1,2)=z(1,2)+1.190797402D-4
c       write data
        write(nfile)out
1	continue
        write(6,*)'get to here 2'
        do 12 i=1,is
        do 12 j=1,js
           out(i,j)=cmplx(d(i,j))
12	continue
        write(6,*)'get to here 3'
c       write data
        write(nfile)out
11	continue
        return
        end

      Subroutine Damp
      parameter(is=41,js=41,jsp=js+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 zmn1(is,js),dmn1(is,js)
      complex*16 dzdt(is,js),dddt(is,js),ci
      common/time1/zmn1,dmn1
      common/timed/dzdt,dddt
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
      common/cconst/ms,ns,mp,np,np2

      ray=dcmplx(1.0d0/dble(14.7*24*60*60))
      zmn1(1,2)=zmn1(1,2)-1.190797402d-4
      do 2000 mm=1,ms
      do 2000 nj=1,ns
      dzdt(mm,nj)=dzdt(mm,nj)-ray*zmn1(mm,nj)
2000  continue
      zmn1(1,2)=zmn1(1,2)+1.190797402d-4
      return
      end
C
C
