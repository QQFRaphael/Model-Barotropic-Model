      program main
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      COMPLEX         PHISMN(is,js)
      COMMON /TIME2/  PHISMN

ccccccccccccc    in this program the first dimension of phismn is
ccccccccccccc    m(zonal wave number) where as the second dimension
ccccccccccccc    is m-n. In x and y the second dimension corresponds
ccccccccccccc    to n-m and the last dimension to m. n-m is the number
ccccccccccccc    of zeros from equator to pole

      COMPLEX         ZM(is,2),DM(is,2),TM(is,2),WM(is,2),CI
      COMMON /FOURS/  ZM,DM,TM,WM
      COMMON /ALEG/   PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1                GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/PCONST/AE,R,CP,GRAV,RKAPPA,RKAP1,DKH,EKH,RTCS,RTSN,
     1  GS2RD,GS7E5,EPSLR,ESCONS,ESC1,ESC2,CI
      COMMON/CCONST/MS,NS,MP,NP,NP2,KP,KPM1,KPP1,KPH,KPHL,KPHU
      COMMON/FFTW/WORK(ix,2)
      complex div(is,js),vort(is,js)
      real out(ix,jx),u(ix,jx),v(ix,jx)
      character*50 ifile,ifile2
      DATA  MPHYS/0/
      DATA NSATA/0/,NNEGA/0/,RHMAXA/-1.E+10/,RHMINA/+1.E+10/
      DATA NSATB/0/,NNEGB/0/,RHMAXB/-1.E+10/,RHMINB/+1.E+10/
      write(6,*)'enter grid data set name'
      read(5,*)ifile
      open(1,file=ifile,form=
c    1'unformatted',access='direct',recl=jx*ix)
     1'unformatted',access='direct',recl=jx*ix)
      write(6,*)'enter spectral data set name'
      read(5,*)ifile2
      open(66,file=ifile2,form='unformatted')
      CALL BSCST
      CALL GAUSLAT(NP2,GLAT,GWGT)
c     idays=21
      write(6,*)'enter number of days'
      read(5,*)idays
      write(6,*)'enter number of days to skip'
      read(5,*)iskip
      if(iskip.gt.0)then
      iend=(2)*iskip
      do 100 i=1,iend
      read(66)phismn
100   continue
      endif
      call reset1(out,phismn)
      irec=0
      do 3 i=1,idays
c    
c        read vorticity
c
      read(66)phismn
      call switch(vort,phismn)
10    continue
c
c        read divergence
c
      read(66)phismn
      call switch(div,phismn)
11    continue  

      call umvm2(vort,div,u,v)
     
c
c        write u,v
c

      call layer(u,out)
c     call smooth(out)
      irec=irec+1
      write(1,rec=irec)out
      call reset1(out,phismn)
12    continue

      call layer(v,out)
c     call smooth(out)
      irec=irec+1
      write(1,rec=irec)out
      call reset1(out,phismn)
13    continue

c
14    continue
3     continue
99    format(5e16.8)
      stop
      end

      subroutine umvm2(z,d,u,v)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      complex z(is,js),d(is,js)
      dimension um(ix,2),vm(ix,2)
      COMMON /ALEG/   PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1                GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/PCONST/AE,R,CP,GRAV,RKAPPA,RKAP1,DKH,EKH,RTCS,RTSN,
     1  GS2RD,GS7E5,EPSLR,ESCONS,ESC1,ESC2,CI
      COMMON/CCONST/MS,NS,MP,NP,NP2,KP,KPM1,KPP1,KPH,KPHL,KPHU
      COMPLEX ZM(is,2),DM(is,2),TM(is,2),WM(is,2),CI,CIM
      COMMON/FOURS/ZM,DM,TM,WM
      COMMON/FFTW/WORK(ix,2)
      DIMENSION U(ix,jx),V(ix,jx)
      MPHYS=0
      DO 6000 LA=1,NP2
      Y=GLAT(LA)      
      GW=GWGT(LA)
      CS2=1.-Y*Y
      CS=SQRT(CS2)
      ACS2=AE*CS2
      CALL PMNS(MPHYS)
      CALL UMVM(z,d)
      DO 200 IHEM=1,2
      DO 200 MM=1,MS
      UM(2*MM-1,IHEM)=REAL(ZM(MM,IHEM))
      UM(2*MM,IHEM)=AIMAG(ZM(MM,IHEM))
      VM(2*MM-1,IHEM)=REAL(DM(MM,IHEM))
  200 VM(2*MM,IHEM)=AIMAG(DM(MM,IHEM))
      DO 370  IHEM=1,2
      CALL FFT1(um(1,IHEM),WORK(1,IHEM),MS,MP,1,+1)
      CALL FFT1(vm(1,IHEM),WORK(1,IHEM),MS,MP,1,+1)
  370 CONTINUE
      do 344 i=1,ix
        u(i,np2+la)=um(i,1)/cs
        u(i,np2-la+1)=um(i,2)/cs
        v(i,np2+la)=vm(i,1)/cs
        v(i,np2-la+1)=vm(i,2)/cs
344   continue
6000  CONTINUE
      return
      end

      SUBROUTINE UMVM(zmn2,dmn2)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
      complex zmn2(is,js),dmn2(is,js)
      COMPLEX ZM(is,2),DM(is,2),TM(is,2),WM(is,2),CI,CIM
      COMMON/FOURS/ZM,DM,TM,WM
      COMMON/ALEG/PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1  GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/PCONST/AE,R,CP,GRAV,RKAPPA,RKAP1,DKH,EKH,RTCS,RTSN,
     1  GS2RD,GS7E5,EPSLR,ESCONS,ESC1,ESC2,CI
      COMMON/CCONST/MS,NS,MP,NP,NP2,KP,KPM1,KPP1,KPH,KPHL,KPHU
      CI=(0.,1.)
      CALL LEGSUM(ZMN2(1,1),ZM(1,1),ZM(1,2),HA,MS,NS,2,1)
      CALL LEGSUM(DMN2(1,1),DM(1,1),DM(1,2),PA,MS,NS,1,2)
      DO 160 MM=1,MS
      CIM=CI*(MM-1.)
      DO 160 IHEM=1,2
  160 ZM(MM,IHEM)=-AE*(ZM(MM,IHEM)+CIM*DM(MM,IHEM))
      CALL LEGSUM(ZMN2(1,1),DM(1,1),DM(1,2),PA,MS,NS,1,2)
      CALL LEGSUM(DMN2(1,1),TM(1,1),TM(1,2),HA,MS,NS,2,1)
      DO 250 MM=1,MS
      CIM=CI*(MM-1.)
      DO 250 IHEM=1,2
  250 DM(MM,IHEM)=-AE*(DM(MM,IHEM)*CIM-TM(MM,IHEM))
      DO 260 IHEM=1,2
      ZM(1,IHEM)=CMPLX(REAL(ZM(1,IHEM)),0.)
  260 DM(1,IHEM)=CMPLX(REAL(DM(1,IHEM)),0.)
C
      RETURN
      END

	subroutine layer(u,out)
        parameter(is=41,js=41,isp=is+1)
        parameter(ix=128,jx=102,jx2=jx/2)
        dimension u(ix,jx),out(ix,jx)
        do 1 i=1,ix
        do 1 j=1,jx
           out(i,j)=u(i,j)
1	continue
        return
        end

	subroutine switch(u,out)
        parameter(is=41,js=41,isp=is+1)
        parameter(ix=128,jx=102,jx2=jx/2)
        complex u(is,js),out(is,js)
        do 1 i=1,is
        do 1 j=1,js
           u(i,j)=out(i,j)
1	continue
        return
        end

        subroutine smooth(u)
        parameter(is=41,js=41,isp=is+1)
        parameter(ix=128,jx=102,jx2=jx/2)

c       smoothing of two dimensional array

        real u(ix,jx),x(ix,jx)
        do 2 j=1,jx
        do 3 i=2,47
           x(i,j)=0.5*u(i,j)+0.25*(u(i+1,j)+u(i-1,j))
3       continue
        x(1,j)=0.5*u(1,j)+0.25*(u(2,j)+u(ix,j))
        x(ix,j)=0.5*u(ix,j)+0.25*(u(1,j)+u(47,j))
2       continue
        do 4 i=1,48
        do 5 j=2,39
           u(i,j)=0.5*x(i,j)+0.25*(x(i,j+1)+x(i,j-1))
5       continue
        u(i,1)=0.5*(x(i,1)+x(i,2))
        u(i,jx)=0.5*(x(i,jx)+x(i,39))
4       continue
        return
        end

      subroutine reset1(out,phismn)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      real out(ix,jx)
      complex phismn(is,js)
      do 1 i=1,ix
      do 1 j=1,jx
         out(i,j)=0.0
1     continue
      do 2 i=1,is
      do 2 j=1,js
         phismn(i,j)=0.0
2     continue
      return
      end



       subroutine pstprc(ggs2)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      COMPLEX         PHISMN(is,js)
      COMMON /TIME2/  PHISMN

ccccccccccccc    in this program the first dimension of phismn is
ccccccccccccc    m(zonal wave number) where as the second dimension
ccccccccccccc    n is the total wave number
      COMMON /ALEG/   PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1                GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/PCONST/AE,R,CP,GRAV,RKAPPA,RKAP1,DKH,EKH,RTCS,RTSN,
     1  GS2RD,GS7E5,EPSLR,ESCONS,ESC1,ESC2,CI
      COMMON/CCONST/MS,NS,MP,NP,NP2,KP,KPM1,KPP1,KPH,KPHL,KPHU
      COMMON/FFTW/WORK(ix,2)
      DIMENSION       GGS2(ix,jx),GS2(ix,jx)
      DATA  MPHYS/0/
      DATA NSATA/0/,NNEGA/0/,RHMAXA/-1.E+10/,RHMINA/+1.E+10/
      DATA NSATB/0/,NNEGB/0/,RHMAXB/-1.E+10/,RHMINB/+1.E+10/
      DO 6000 LA=1,NP2
      Y=GLAT(LA)      
      GW=GWGT(LA)
      CS2=1.-Y*Y
      CS=SQRT(CS2)
      ACS2=AE*CS2
      CALL PMNS(MPHYS)
      CALL LEGSUM(PHISMN,gs2(1,1),gs2(1,2),PMN,MS,NS,1,2)
      DO 370  IHEM=1,2
      CALL FFT1(gs2(1,IHEM),WORK(1,IHEM),MS,MP,1,+1)
  370 CONTINUE
      do 344 i=1,ix
        ggs2(i,np2+la)=gs2(i,1)
        ggs2(i,np2-la+1)=gs2(i,2)
 344    continue
C
 6000 CONTINUE
      return
      END
C
      SUBROUTINE BSCST
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      COMMON/PCONST/AE,R,CP,GRAV,RKAPPA,RKAP1,DKH,EKH,RTCS,RTSN,
     1  GS2RD,GS7E5,EPSLR,ESCONS,ESC1,ESC2,CI
      COMMON/CCONST/MS,NS,MP,NP,NP2,KP,KPM1,KPP1,KPH,KPHL,KPHU
C
C
C  SET SPATIAL RESOLUTION AND SPECTRAL TRUNCATION PARAMETERS
C  (NOTE NO. OF OUTPUT LATITUDE CIRCLES (NP) MUST BE EVEN.  IN THE CASE
C  OF A REGULAR OUTPUT GRID, THE EQUATORIAL CIRLCE IS OUTPUT TWICE AND
C  COUNTED AS TWO LATITUDE CIRCLES.  E.G. FOR 2.5 X 2.5 DEG
C  LAT/LON GRID, MP=144 AND NP=74)
C
      MS=is
      NS=js
      MP=ix
      NP=jx
      IF(MOD(NP,2).NE.0) STOP 1
      NP2=NP/2
      ICHK=783
      IF(ICHK.LT.(4*MP+15))  STOP 2
      AE=6.371E6
      R=287.05
      CP=1005.
      GRAV=9.8
      RKAPPA=R/CP
      CI=(0.,1.)
      RETURN
      END
      SUBROUTINE GAUSLAT(NP2,GLTS,GWTS)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
C  CALCULATE SINE OF GAUSSIAN LATITUDES.  ALSO CALCULATE GAUSSIAN
C  WEIGHTS FOR GAUSSIAN QUADRATURE DONE LATER IN VERIFICATION ROUTINE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GWTS(NP2),GLTS(NP2)
      DIMENSION RADCOL(jx2),GW(jx2)
C
      EPS=1.D-12
      ONE=1.0D0
      PI=(4.0D0)*DATAN(ONE)
      K2=2*NP2
      RK2=K2
      SCALE=2.D0/(RK2**2)
      K1=K2-1
      DRADZ=PI/360.0D0
      RAD=0.0D0
C
      DO 1000 K=1,NP2
      ITER=0
      DRAD=DRADZ
1     CALL POLY(K2,RAD,P2)
2     P1=P2
      ITER=ITER+1
      RAD=RAD+DRAD
      CALL POLY(K2,RAD,P2)
      IF (SIGN(ONE,P1).EQ.SIGN(ONE,P2)) GO TO 2
      IF (DRAD.LT.EPS) GO TO 3
      RAD=RAD-DRAD
      DRAD=DRAD*0.25D0
      GO TO 1
3     CONTINUE
      RADCOL(K)=RAD
      CALL POLY(K1,RAD,P1)
      X=DCOS(RAD)
      W=SCALE*(ONE-X*X)/(P1*P1)
      GW(K)=W
      CALL POLY(K2,RAD,P1)
      GLTS(NP2+1-K)=X
      GWTS(NP2+1-K)=GW(K)
1000  CONTINUE
C
C  PRINT RESULTS
C
      CFR=180.0D0/PI
C     PRINT 201
C     PRINT 202,  (J,GLTS(J),GWTS(J),ASIN(GLTS(J))*CFR,J=1,NP2)
C 201 FORMAT(1H0//1H ,5X,"RESULTS FROM SUBROUTINE GAUSLAT"/
C    11H0,8X,"LAT INDX",5X,"SINE GAUS LAT",9X,"GAUS WEIGHT",
C    2 4X,"GAUS LAT(DEG)"//)
C 202 FORMAT(1H ,10X,I3,6X,F17.15,4X,F17.15,4X,F8.4)
C
      RETURN
      END
C**********************************************   POLY   *************
C
C
      SUBROUTINE POLY(N,RAD,P)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  CALLED BY ROUTINE GAUSLAT
C
      X=DCOS(RAD)
      Y1=1.0D0
      Y2=X
      DO 1 I=2,N
      G=X*Y2
      Y3=G-Y1+G-(G-Y1)/DBLE(I)
      Y1=Y2
      Y2=Y3
1     CONTINUE
      P=Y3
      RETURN
      END
C*************************************************   PMNS   ***********
C
C
      SUBROUTINE PMNS(MPHYS)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ALEG/PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1   GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/CCONST/MS,NS,MP,NP,NP2
      DIMENSION A1(is),A2(is)
C
      NS1=NS+1
      PMN(1,1)=DSQRT(0.5D0)
      DO 100 MM=2,MS
      AM=MM-1.0D0
      A1(MM)=CS*DSQRT(1.0D0+1.0D0/(2.0D0*AM))
100   CONTINUE
      DO 101 MM=2,MS
101   PMN(MM,1)=A1(MM)*PMN(MM-1,1)
      DO 200 MM=1,MS
      AM=MM-1.0D0
      A1(MM)=Y*DSQRT(2.0D0*AM+3.0D0)
200   CONTINUE
      DO 201 MM=1,MS
201   PMN(MM,2)=A1(MM)*PMN(MM,1)
      DO 400 NJ=3,NS1
      DO 300 MM=1,MS
      AM=MM-1.0D0
      AN=NJ+MM-2
      ANM=AN*AN-AM*AM
      A1(MM)=Y*DSQRT((4.D0*AN*AN-1.D0)/ANM)
300   A2(MM)=DSQRT((2.D0*AN+1.D0)*(ANM-2.D0*AN+1.D0)
     1  /((2.D0*AN-3.D0)*ANM))
      DO 301 MM=1,MS
301   PMN(MM,NJ)=A1(MM)*PMN(MM,NJ-1)-A2(MM)*PMN(MM,NJ-2)
400   CONTINUE
C
      IF(MPHYS.NE.0)GO TO 999
      HMN(1,1)=0.D0
      DO 500 MM=2,MS
      AM=MM-1.D0
500   A1(MM)=DSQRT(AM*AM*(2.D0*AM+1.D0)/
     1(4.D0*AM*AM+8.D0*AM+3.D0))
      DO 501 MM=2,MS
501   HMN(MM,1)=A1(MM)*PMN(MM,2)
      DO 650 NJ=2,NS
      DO 600 MM=1,MS
      AN=NJ+MM-2.D0
      AM=MM-1.D0
      ANM=AN*AN-AM*AM
      A1(MM)=AN*DSQRT((ANM+2.D0*AN+1.D0)/
     1(4.D0*AN*AN+8.D0*AN+3.D0))
600   A2(MM)=(AN+1.D0)*DSQRT(ANM/(4.D0*AN*AN-1.D0))
      DO 601 MM=1,MS
601   HMN(MM,NJ)=A1(MM)*PMN(MM,NJ+1)-A2(MM)*PMN(MM,NJ-1)
650   CONTINUE
      DO 800 NJ=1,NS
      RNN1=NJ*(NJ-1.D0)
      IF(NJ.EQ.1)RNN1=1.D0
      RNI=1.D0/RNN1
      PA(1,NJ)=PMN(1,NJ)*RNI
      HA(1,NJ)=HMN(1,NJ)*RNI
      DO 700 MM=2,MS
      AN=NJ+MM-2.D0
700   A1(MM)=1.D0/(AN*(AN+1.D0))
      DO 701 MM=2,MS
      PA(MM,NJ)=PMN(MM,NJ)*A1(MM)
701   HA(MM,NJ)=HMN(MM,NJ)*A1(MM)
800   CONTINUE
      PA(1,1)=0.D0
      HA(1,1)=0.D0
C
999   CONTINUE
C
      RETURN
      END
C************************************************   LEGSUM   *********
C
C
      SUBROUTINE LEGSUM(FMN,FMNH,FMSH,PMN,MS,NS,NE,NO)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
      COMPLEX FMN(is,js),FMNH(is),FMSH(is),FME(is),FMO(is)
      DOUBLE PRECISION PMN(is,isp)
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
C  OUTPUT FM(MM,IHEM) FOURIER COEFF FOR ALL M
C      AT N. HEM +Y (IHEM=1) AND AT S. HEM -Y (IHEM=2)
C
      DO 50 MM=1,MS
      FME(MM)=(0.,0.)
   50 FMO(MM)=(0.,0.)
      DO 100 NJ=NE,NS,2
      DO 100 MM=1,MS
  100 FME(MM)=FME(MM)+FMN(MM,NJ)*sngl(PMN(MM,NJ))
      DO 200 NJ=NO,NS,2
      DO 200 MM=1,MS
  200 FMO(MM)=FMO(MM)+FMN(MM,NJ)*sngl(PMN(MM,NJ))
      DO 300 MM=1,MS
      FMNH(MM)=FME(MM)+FMO(MM)
      FMSH(MM)=FME(MM)-FMO(MM)
  300 CONTINUE
      FMNH(1)=CMPLX(REAL(FMNH(1)),0.)
      FMSH(1)=CMPLX(REAL(FMSH(1)),0.)
      RETURN
      END

      SUBROUTINE FFT1(F,W,MS,MP,NUM,IDIR)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
      DIMENSION F(ix),W(ix)
      COMMON/FFT/IFAX(10),TRIGS(ix,2)
C
C  SUBROUTINE DOES MULTIPLE HALF COMPLEX FFTS
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
C      MS            NUMBER OF WAVES (INCLUDING WAVENUMBER 0) MUST BE
C                                                             .LT. MP
C
C      MP            LENGTH OF EACH TRANSFORM (NUMBER OF GRID POINTS)
C
C      NUM           NUMBER OF TRANSFORMS TO BE PERFORMED SIMULTANEOUSLY
C
C      IDIR          DIRECTION OF TRANSFORM (+1 FOR SYNTHESIS, -1 FOR
c                                                            ANALYSIS)
C  ON OUTPUT
C
C      F(MP,NUM)      CONTAINS TRANFORMED FIELDS
C
C  ON FIRST CALL TO FFT1, MUST SET UP  FACTORS AND TRIG FUNCTIONS
C
      DATA ITST/0/
      IF(ITST.EQ.1)GO TO 50
C
C   ON FIRST CALL MUST SET UP FACTORS AND TRIG FUNCTIONS
C
      ITST=1
      CALL FAX(IFAX,MP,3)
      CALL FFTRIG(TRIGS,MP,3)
C
   50 IF(IDIR.EQ.-1)GO TO 250
C
C
C   FOR SYNTHESIS MUST ZERO FILL WAVENUMBERS GT MS
C
      LLO=2*MS+1
      LHI=MP
      DO 200 I=1,NUM
      DO 100 LO=LLO,LHI
  100 F(LO)=0.
      LLO=LLO+MP
      LHI=LHI+MP
  200 CONTINUE
C
C   ANALYSIS REQUIRES NO PREPROCESSING
C
  250 CALL FFT99M(F,W,TRIGS,IFAX,1,MP,MP,NUM,IDIR)
      RETURN
      END
