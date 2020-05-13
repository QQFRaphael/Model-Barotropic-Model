      program main
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16         PHISMN(is,js)
      COMMON /TIME2/  PHISMN

ccccccccccccc    in this program the first dimension of phismn is
ccccccccccccc    m(zonal wave number) where as the second dimension
ccccccccccccc    is m-n. In x and y the second dimension corresponds
ccccccccccccc    to n-m and the last dimension to m. n-m is the number
ccccccccccccc    of zeros from equator to pole

      COMPLEX*16         ZM(is,2),DM(is,2),TM(is,2),WM(is,2),CI
      COMMON /FOURS/  ZM,DM,TM,WM
      COMMON /ALEG/   PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1                GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/PCONST/AE,R,CP,GRAV,RKAPPA,RKAP1,DKH,EKH,RTCS,RTSN,
     1  GS2RD,GS7E5,EPSLR,ESCONS,ESC1,ESC2,CI
      COMMON/CCONST/MS,NS,MP,NP,NP2,KP,KPM1,KPP1,KPH,KPHL,KPHU
      COMMON/FFTW/WORK(ix,2)
      complex*16 tfix(is,js),topog(is,js)
      real*8 out(ix,jx)
      real*4 outt(ix,jx)
      DATA  MPHYS/0/
      DATA NSATA/0/,NNEGA/0/,RHMAXA/-1.D+10/,RHMINA/+1.D+10/
      DATA NSATB/0/,NNEGB/0/,RHMAXB/-1.D+10/,RHMINB/+1.D+10/
      open(2,file='../data/vor.grads.frc',access='direct',form=
c     1'unformatted',recl=ix*jx)
     1'unformatted',recl=ix*jx)
      CALL BSCST
      CALL GAUSLAT(NP2,GLAT,GWGT)
      irec=0
      read(2,rec=1)outt
      do 5 ii=1,ix
      do 5 jj=1,jx
           out(ii,jj)=dble(outt(ii,jj))
5     continue 
      write(6,*)'call preprc'
      call preprc(out)
      write(6,*)'exit preprc'
      do mb=1,is
      do nb=1,js
         topog(mb,nb)=phismn(mb,nb)
      enddo
      enddo
      write(88)topog
      stop
      end
     
      subroutine reset1(out,phismn)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      real*8 out(ix,jx)
      complex*16 phismn(is,js)
      do 1 i=1,ix
      do 1 j=1,jx
         out(i,j)=0.0d0
1     continue
      do 2 i=1,is
      do 2 j=1,js
         phismn(i,j)=0.0d0
2     continue
      return
      end




       subroutine preprc(ggs2)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16         PHISMN(is,js)
      COMMON /TIME2/  PHISMN

ccccccccccccc    in this program the first dimension of phismn is
ccccccccccccc    m(zonal wave number) where as the second dimension
ccccccccccccc    n is the total wave number
      COMPLEX*16         ZM(is,2),DM(is,2),TM(is,2),WM(is,2),CI
      COMMON /FOURS/  ZM,DM,TM,WM
      COMMON /ALEG/   PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1                GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/PCONST/AE,R,CP,GRAV,RKAPPA,RKAP1,DKH,EKH,RTCS,RTSN,
     1  GS2RD,GS7E5,EPSLR,ESCONS,ESC1,ESC2,CI
      COMMON/CCONST/MS,NS,MP,NP,NP2,KP,KPM1,KPP1,KPH,KPHL,KPHU
      COMMON/FFTW/WORK(ix,2)
      REAL*8 GGS2(ix,jx),GS2(ix,2)
      MPHYS=0
      do 200 i=1,ix
         gs2(i,1)=0.0d0
         gs2(i,2)=0.0d0
200   continue
      DO 6000 LA=1,NP2
      Y=GLAT(LA)      
      GW=GWGT(LA)
      CS2=1.0D0-Y*Y
      CS=DSQRT(CS2)
      ACS2=AE*CS2
      CALL PMNS(MPHYS)
      do 344 i=1,mp
        gs2(i,1)=ggs2(i,np2+la)
        gs2(i,2)=ggs2(i,np2-la+1)
344   continue
      CALL FFT1(gs2,WORK,MS,MP,2,-1)
      CALL GQL(gs2(1,1),gs2(1,2),PHISMN,PMN,GW,MS,NS)
C
 6000 CONTINUE
      return
      END
C
      SUBROUTINE BSCST
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      AE=6.371D6
      IF(MOD(NP,2).NE.0) STOP 1
      NP2=NP/2
      ICHK=783
      IF(ICHK.LT.(4*MP+15))  STOP 2
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
      REAL*8 GWTS(NP2),GLTS(NP2)
      REAL*8 RADCOL(jx2),GW(jx2)
C
      EPS=1.D-15
      ONE=1.D0
      PI=4.0D0*DATAN(ONE)
      K2=2*NP2
      RK2=K2
      SCALE=2.0D0/(RK2**2)
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
      CFR=180.D0/PI
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/ALEG/PMN(is,isp),HMN(is,js),PA(is,isp),HA(is,js),
     1   GLAT(jx2),GWGT(jx2),GW,Y,CS,CS2,ACS2
      DOUBLE PRECISION PMN,HMN,PA,HA,GLAT,GWGT,GW,Y,CS,CS2,ACS2
      COMMON/CCONST/MS,NS,MP,NP,NP2
      REAL*8 A1(is),A2(is)
C
      NS1=NS+1
      PMN(1,1)=DSQRT(0.50D0)
      DO 100 MM=2,MS
      AM=DFLOAT(MM)-1.0D0
      A1(MM)=CS*SQRT(1.0D0+1.0D0/(2.0D0*AM))
100   CONTINUE
      DO 101 MM=2,MS
101   PMN(MM,1)=A1(MM)*PMN(MM-1,1)
      DO 200 MM=1,MS
      AM=DFLOAT(MM)-1.0D0
      A1(MM)=Y*DSQRT(2.0D0*AM+3.0D0)
200   CONTINUE
      DO 201 MM=1,MS
201   PMN(MM,2)=A1(MM)*PMN(MM,1)
      DO 400 NJ=3,NS1
      DO 300 MM=1,MS
      AM=DFLOAT(MM)-1.0D0
      AN=NJ+MM-2
      ANM=AN*AN-AM*AM
      A1(MM)=Y*DSQRT((4.0D0*AN*AN-1.0D0)/ANM)
300   A2(MM)=DSQRT((2.0D0*AN+1.0D0)*(ANM-2.0D0*AN+1.0D0)
     1  /((2.0D0*AN-3.0D0)*ANM))
      DO 301 MM=1,MS
301   PMN(MM,NJ)=A1(MM)*PMN(MM,NJ-1)-A2(MM)*PMN(MM,NJ-2)
400   CONTINUE
C
      IF(MPHYS.NE.0)GO TO 999
      HMN(1,1)=0.0D0
      DO 500 MM=2,MS
      AM=DFLOAT(MM)-1.0D0
500   A1(MM)=DSQRT(AM*AM*(2.0D0*AM+1.0D0)/
     1(4.0D0*AM*AM+8.0D0*AM+3.0D0))
      DO 501 MM=2,MS
501   HMN(MM,1)=A1(MM)*PMN(MM,2)
      DO 650 NJ=2,NS
      DO 600 MM=1,MS
      AN=DFLOAT(NJ+MM)-2.0D0
      AM=DFLOAT(MM)-1.0D0
      ANM=AN*AN-AM*AM
      A1(MM)=AN*DSQRT((ANM+2.0D0*AN+1.0D0)/
     1(4.0D0*AN*AN+8.0D0*AN+3.0D0))
600   A2(MM)=(AN+1.0D0)*DSQRT(ANM/(4.0D0*AN*AN-1.0D0))
      DO 601 MM=1,MS
601   HMN(MM,NJ)=A1(MM)*PMN(MM,NJ+1)-A2(MM)*PMN(MM,NJ-1)
650   CONTINUE
      DO 800 NJ=1,NS
      RNN1=NJ*(NJ-1.0D0)
      IF(NJ.EQ.1)RNN1=1.0D0
      RNI=1.0D0/RNN1
      PA(1,NJ)=PMN(1,NJ)*RNI
      HA(1,NJ)=HMN(1,NJ)*RNI
      DO 700 MM=2,MS
      AN=DFLOAT(NJ+MM)-2.0D0
700   A1(MM)=1.0D0/(AN*(AN+1.0D0))
      DO 701 MM=2,MS
      PA(MM,NJ)=PMN(MM,NJ)*A1(MM)
701   HA(MM,NJ)=HMN(MM,NJ)*A1(MM)
800   CONTINUE
      PA(1,1)=0.0D0
      HA(1,1)=0.0D0
C
999   CONTINUE
C
      RETURN
      END
      SUBROUTINE GQL(FMNH,FMSH,FMN,PMN,GW,MS,NS)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 FMNH(is),FMSH(is),FMN(is,js)
      REAL*8 PMN(is,isp),GW
C
C  COMPUTES AND ACCUMULATES CONTRIBUTION OF FMNH(MM)
C     AND FMSH(MM)  TO GAUSSIAN QUADRATURE
C
C  SYMMETRIC PART
      DO 50 NJ=1,NS,2
      DO 50 MM=1,MS
      FMN(MM,NJ)=FMN(MM,NJ)+((PMN(MM,NJ)*GW))
     1*(FMNH(MM)+FMSH(MM))
50    CONTINUE
C  ANTISYMMETRIC PART
      DO 60 NJ=2,NS,2
      DO 60 MM=1,MS
      FMN(MM,NJ)=FMN(MM,NJ)+((PMN(MM,NJ)*GW))
     1*(FMNH(MM)-FMSH(MM))
60    CONTINUE
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 FMN(is,js),FMNH(is),FMSH(is),FME(is),FMO(is)
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
      FME(MM)=(0.0D0,0.0D0)
   50 FMO(MM)=(0.0D0,0.0D0)
      DO 100 NJ=NE,NS,2
      DO 100 MM=1,MS
  100 FME(MM)=FME(MM)+FMN(MM,NJ)*(PMN(MM,NJ))
      DO 200 NJ=NO,NS,2
      DO 200 MM=1,MS
  200 FMO(MM)=FMO(MM)+FMN(MM,NJ)*(PMN(MM,NJ))
      DO 300 MM=1,MS
      FMNH(MM)=FME(MM)+FMO(MM)
      FMSH(MM)=FME(MM)-FMO(MM)
  300 CONTINUE
      FMNH(1)=DCMPLX(DBLE(FMNH(1)),0.0D0)
      FMSH(1)=DCMPLX(DBLE(FMSH(1)),0.0D0)
      RETURN
      END

      SUBROUTINE FFT1(F,W,MS,MP,NUM,IDIR)
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 F(ix),W(ix),TRIGS(ix,2)
      INTEGER IFAX(10)
      COMMON/FFT/IFAX,TRIGS
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
  100 F(LO)=0.0D0
      LLO=LLO+MP
      LHI=LHI+MP
  200 CONTINUE
C
C   ANALYSIS REQUIRES NO PREPROCESSING
C
  250 CALL FFT99M(F,W,TRIGS,IFAX,1,MP,MP,NUM,IDIR)
      RETURN
      END
