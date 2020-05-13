C*********************************************************
C********     KIRTMAN, BEN. GSM        *******************
C********     RHOMBOIDAL TRUNCATION    *******************
C********       ACGM                   *******************
C*********************************************************
C*********************************************************
C
C
C       ********** DOUBLE PRECISION VERSION **************
C
C
      Program Main
      parameter(is=41,js=41,isp=is+1)
      parameter(ix=128,jx=102,jx2=jx/2)
      implicit double precision (a-h,o-z)
C
      complex*16 zmn1(is,js),dmn1(is,js)
      complex*16 zmn2(is,js),dmn2(is,js)
      complex*16 dzdt(is,js),dddt(is,js)
      complex*16 zm(is,2),dm(is,2),ci
      complex*16 tm(is,2),wm(is,2)
      common/time1/zmn1,dmn1
      common/time2/zmn2,dmn2
      common/timed/dzdt,dddt
      common/fours/zm,dm,tm,wm
      common/grdpt/z(ix,2),d(ix,2),
     1  u(ix,2),v(ix,2)
      common/prod/a(ix,2),b(ix,2)
      common/aleg/pmn(is,isp),hmn(is,js),pa(is,isp),ha(is,js),
     1  glat(jx2),gwgt(jx2),gw,y,cs,cs2,acs2
      double precision pmn,hmn,pa,ha,glat,gwgt,gw,y,cs,cs2,acs2
      common/pconst/ae,r,cp,grav,rkappa,rkap1,dkh,ekh,rtcs,rtsn,
     1  gs2rd,gs7e5,epslr,escons,esc1,esc2,ci
C
      complex*16 zmnt(is,js),dmnt(is,js)
      complex*16 phismn(is,js)
      common/tmean/zmnt,dmnt
C
      common/cconst/ms,ns,mp,np,np2
      common/tconst/dt,alpha,aee,andree
      common/fft/ifax(10),trigs(ix,2)
      common/clock/nt
C
C
      common /calendr/ imndy, imn
C
C
C
C     NTIME:  NO. OF TIME STEPS FOR THIS RUN (INTEGER)
C    ANDREE:  TIME FILTER COEF (REAL - E.G. 0.04)
C     ALPHA:  SEMI-IMPLICIT PARAM (REAL - BETWEEN 0.5 AND 1.0)
C       DKH:  HORIZONTAL DIFFUSION COEF FOR DIVERGENCE (REAL)
C       EKH:  HORIZONTAL DIFFUSION COEF FOR OTHER VARBLS (REAL)
C       IMN:  INTEGER MONTH OF THE YEAR
C     IMNDY:  INTEGER DAY OF THE MONTH
C
C
C     UNIT 67 = formt.dat TIME MEAN SPECTRAL COEFFICIENTS
C     UNIT 66 = form.dat  ONCE DAILY SPECTRAL COEFFICIENTS
C     UNIT 23 = outre.dat OUTPUT RESTART FILE
C     UNIT 24 = initial.dat INPUT RESTART FILE
C     UNIT 7  = fort.7 RUN DECK INPUT PARAMETERS
C     UNIT 44 = prescribed divergence
C
C
C     
       read (7,*) ntime
       read (7,*) andree,alpha
       read (7,*) dkh,ekh
       read (7,*) imn, imndy
       write(6,*) ntime
       write(6,*) andree,alpha
       write(6,*) dkh,ekh
       write(6,*) imn,imndy
C
C
C  SET VARIOUS CONSTANTS, SIGMA STRUCTURE, AND BASIC STATE
C
      CALL BSCST
C
C
C  GENERATE GAUSSIAN LATITUDES (SIN) AND WEIGHTS
C
      CALL GAUSLAT(NP2,GLAT,GWGT)
C
C
c
c SET INITIAL CONDITIONS
c
      CALL ZERO1(ZMN1)
      CALL ZERO1(DMN1)
      CALL ZERO1(ZMNT)
      CALL ZERO1(DMNT)
      CALL ZERO1(ZMN2)
      CALL ZERO1(DMN2)
      CALL ZERO1(DZDT)
      CALL ZERO1(DDDT)
      CALL RESTART
C
C     Initialize Divergence forcing
C
      read(44)phismn
      do i=1,is
      do j=1,js
         dmn1(i,j)=phismn(i,j)
         dmn2(i,j)=phismn(i,j)
      enddo
      enddo
c
c
C  BEGIN TIME INTEGRATION
C
C
C        NTIME=NUMBER OF DAYS IN THE FORECAST
C
      iend=(144*ntime)
      isteps=iend
      dt=300.0D0

      do 400 nt=2,iend+1
cwrite(6,*) nt
C
C
C  do 1 time step of the model
C
      nstep=nt-1
      time=dble(nstep)*dt/3600.d0
      timefc=dble(nt)*dt/3600.d0
c
c   IF INTERGER DAY UPDATE DIVERGENCE FORCING 
c
      kk=mod(nt,144)
      if(kk.eq.0)then
        imndy=imndy+1
CCC      update phsimn here
      endif
      do i=1,is
      do j=1,js
         dmn1(i,j)=phismn(i,j)
         dmn2(i,j)=phismn(i,j)
      enddo
      enddo
c
c
C
C  LATITUDE LOOP, COMPUTES NONLINEAR TERMS AND BL FLUXES
C
      CALL LALOOP
C
C
C  SUBGRID SCALE DIFFUSION
C
      CALL DIFFSN
C
C  LINEAR DAMPING
C
      CALL DAMP
C
C  TIME STEP THE MODEL
C
      CALL IMPLCT(NT,FILTC)
C
C  DUMP SNAPSHOT DATA
C  
C
      kk=mod(nt,144)
      write(*,*) kk
      
      if(kk.eq.0)then
        write(*,*) "liuyong"        
        write(6,*)'call wrtdat'
        CALL WRTDAT(ZMN1,DMN1,nt,IREC2,66)
        write(6,*)'write data',nt
      endif
      do i=1,is
      do j=1,js
         dmn1(i,j)=phismn(i,j)
         write(*,*) phismn(1,5)
         dmn2(i,j)=phismn(i,j)
      enddo
      enddo
400   continue
      CALL WHIST
      irec2=1
      stop
      end
