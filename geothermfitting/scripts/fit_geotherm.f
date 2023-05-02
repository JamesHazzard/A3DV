C
C***********************************************************************
C
C       PROGRAM TO FIT NODULE DATA WITH A STEADY STATE CONDUCTIVE MODEL
C       OF THE LITHOSPHERE, AND PLOT THE RESULTING FIT AS TEMPERATURE AS A
C       COMPILE WITH: gfortran outplot_shorttle_brent.f -fdefault-real-8 -fno-align-commons -o outplot_shorttle_brent
C       FUNCTION OF DEPTH
C
      SAVE
      OPEN(UNIT=16,FILE='output.dat')
      OPEN(UNIT=19,FILE='geoth.dat')
      OPEN(UNIT=15,FILE='input.dat')
      CALL FITOBS(NOFJ)
C
c     INITIALIZING PLOT
	  CALL OUTPUT(NOFJ)
C
C      IXFIG=0
C      CALL STARTP(IXFIG)
C      CALL PLTOBS(NOFJ)
C
C     ENDING PLOT
C
C      CALL ENGPLT
      END
C
C***********************************************************************
C
      SUBROUTINE FITOBS(NOFJSV)
      SAVE
C
C     FINDS THE  MECHANICAL BOUNDARY
C     LAYER THICKNESS WHICH BEST FITS P,T ESTIMATES FROM NODULES
C
      EXTERNAL FUNCP
      COMMON/NODDAT/PNOD(10000),TNOD(10000),NOD,FITNOD
      COMMON/INFON/TEM,VIS,TLITH,DENS,ALP,CP,G,COND,D,PFAC,
     C DZ,TP,HSCALE,TCRUST,MMOD,CRUSTK,igrid,min,iwrit
     C ,TUCRUST,TLCRUST,HUCRUST,HLCRUST,hmantle,DEPTHMX
      COMMON/CALC/TTHERM,TEMPL
      COMMON/TITL/TITLE,FILEINP
      COMMON/SWAVEFIL/FILEVS,GSIZE,DEPTHMN,IVS
      INTEGER*4 NOFJ,IVS
	  REAL*4 LAT,LON
	  REAL*8 P,FRET
      CHARACTER*60 TITLE
      CHARACTER*60 FILEINP,FILEVS
      PARAMETER(NP=1)
      DIMENSION XI(NP,NP)
      NAMELIST/FIT/TP,MMOD,CRUSTK,TCRUST,HSCALE,VIS,TEM,
     C TITLE,FILEINP,IFIT,TUCRUST,TLCRUST,HUCRUST,HLCRUST,
     C hmantle,FILEVS,IVS,GSIZE,NOFJ,DEPTHMX,VIS,DEPTHMN,igrid,
     c min,iwrit,LON,LAT
c
c       iwrit=1 (=0 default) write out the geotherm as depths and temperatures
c       to output.dat
c
      iwrit=1
c
c       min (=1 default) reads in P,T data from mineralalogy from fileinp,
c       =0 no P,T data input.  If min=0 ifit is set to 0
c
      min=1
c
c       igrid=1 (default) plots grid lines, =0 does not
c
      igrid=1
C
C       DEPTHMX (DEFAULT=250.) AND DEPTHMN SET DEPTHS TO WHICH THE
C       CALCULATIONS AND PLOTS EXTEND
C
      DEPTHMX=400.
      DEPTHMN=0.
C
C       NOFJ=0 (DEFAULT=1) CALCULATES VS AND TEMP FROM FAUL & JACKSON'S
C       EXPRESSIONS
C
      NOFJ=1
C
C       IF IVS=1 (DEFAULT=0) READS IN SWAVE VELOCITIES AND TEMPERATURES
C       FROM FILEVS AND ADDS THEM TO THE PLOTS
C
      IVS=0
C
C       GSIZE IS GRAIN SIZE IN METRES
C
      GSIZE=1.E-3
C
C       IFIT=1 (DEFAULT) FINDS BEST FITTING TEM
C           =0 USES INPUT VALUE
C
      IFIT=1
C
C TP IS THE POTENTIAL TEMPERATURE OF THE BASE OF THE LITHOSPHERE
C IN DEG C,MEASURED AT THE EARTHS SURFACE.
C  potential temp. is TP, Mean interior is 1300, Deg C, Plumes
C  are up to 1550 Deg C
C
      TP=1333.
C
C     initial mechanical b.l. thickness in km
C
      TEM=100.
C
C       MMOD=0 uses constant conductivity
c
C       MMOD=1 uses Hofmeister's conductivity model
c
C       MMOD=2 uses Jaupart et al.s 1998
c
C       MMOD=3 uses Hofmeister's conductivity model for olivine and
c       Whittington et al.'s for the crust
c
C       MMOD=4 uses Osako et al. 2004 PEPI 311-320 for olivine for cond(P,T)
c       and a constant crustal conductivity of 2.5 w/m K
c
C       MMOD=5 uses Osako et al. 2004 PEPI 311-320 for olivine for cond(P,T)
c       and Whittington et al.'s for the crust
c
C       MMOD=6 uses Korenage & Korenaga 2016 JGR for olivine and oceanic crust
C       with cond(P,T)
c
C       MMOD=7 uses Korenage & Korenaga 2016 JGR for olivine and constant crustal
C       conductivity set by user in input.dat (CRUSTK).
c
C       MMOD=8 uses Korenage & Korenaga 2016 JGR for olivine and temperature-dependent crustal
C       conductivity, with k0 (surface temperature conductivity) set by user in input.dat (CRUSTK).
C
      MMOD=1
C
C       CRUSTK is crustal conductivity used in MMOD=7
C
      CRUSTK=0.0
C
C
C       TCRUST IS THE CRUSTAL THICKNESS IN KM
C
C      TCRUST=40.
C
C       TUCRUST,TLCRUST ARE THE THICKNESSES OF THE UPPER AND LOWER CRUST
C       IN KM. TCRUST IS SET TO TUCRUST+TLCRUST
C
      TUCRUST=22.
      TLCRUST=23.
C
C       HUCRUST, HLCRUST ARE THE HEAT GENERATION RATES IN THE UPPER
C       AND LOWER CRUST IN WATTS/M**3. THESE VALUES HAVE BEEN
C       CHOSEN TO GIVE A SURFACE HEAT FLUX OF ABOUT 46 mW/m**2
C       AND A MANTLE HEAT FLOW OF 12 mW/m**2 WHEN THE CRUSTAL THICKNESS IS 45
C       KM, TO MATCH JAUPART+ JGR 103 15,269-15,286 1998
c       hmantle is the heat generation rate in the mantle part of the
c       mechanical boundary layer
C
      HUCRUST=1.12E-6
      HLCRUST=0.4E-6
      hmantle=0.
C
C       HSCALE IS LENGTH SCALE OVER WHICH THE CRUSTAL RADIOACTIVITY
C       DECREASES IN KM.  IF HSCALE.GT.(4.*TCRUST) THE RADIOACTIVITY
C       IS TAKEN TO BE CONSTANT THROUGHOUT THE CRUST
C
      HSCALE=400.
C
C CP IS SPECIFIC HEAT,
C
      DENS=3.3E+3
      ALP=3.E-5
      G=9.81
      CP=1.187E+3
C
c       VIS THE KINEMATIC VISCOSITY OF THE MANTLE IN M2 S-1
C       Post glacial uplift gives 2.E+17, Craig and McKenzie give 4.E+15
C
C       VIS=4.E+15
   	  VIS=9E+16
C
C PFAC CONVERTS DEPTH IN KM TO PRESSURE IN GPa, P=PFAC*DEPTH
C
C       PFAC=DENS*G*1.E-6
C
C IF CONSISTENT WITH DAN'S 2013 PARAMETERISATION...
	  PFAC=1./30.
C
C       DZ IS DEPTH STEP IN KM
C
      DZ=1.
C
C  D IS DEPTH OF CONVECTING LAYER IN KM
C
      D=7.E+5
      READ(15,FIT)
      if(min.eq.0) ifit=0
      NOFJSV=NOFJ
      TCRUST=TUCRUST+TLCRUST
      WRITE(6,FIT)
      WRITE(16,FIT)
C       print *, CRUSTK
C       stop
      nod=0
      if(ifit.eq.1) then
C
C       READS IN DATA FROM NODULES
C
        OPEN(UNIT=19,FILE=FILEINP)
        I=0
2       I=I+1
        READ(19,105,END=4) PNOD(I),TNOD(I)
105     FORMAT(F12.5,1X,F12.5)
        WRITE(6,105) PNOD(I),TNOD(I)
        WRITE(16,105) PNOD(I),TNOD(I)
        GO TO 2
4       CONTINUE
        CLOSE(4)
        NOD=I-1
      endif
      IF(IFIT.EQ.0) RETURN
C
C       FINDS BEST FITTING TEM IF IFIT=1
C
      P=TEM
      N=NP
      FTOL=0.1
C
      DO I=1,NP
        DO J=1,NP
          XI(I,J)=0.
          IF(I.EQ.J) XI(I,J)=1.
        END DO
      END DO
	  PMIN=20.
      PMAX=400.
	  FRET=BRENT(P,FUNCP,PMIN,PMAX)
	  PRINT *, 'Mechanical Boundary Layer Depth is', P, 'km'
      WRITE(6,180) P
180   FORMAT(' TEM=',F5.1)
      RETURN
      END

C c--------------------------------------------------
C c Brent to iteratively minimise misfit between calculated and predicted lithospheric thickness
C c--------------------------------------------------
      function brent(P,f,AX,CX)

      integer*4 iter
      real*8 a,b,dd,e,etemp,fu,fv,fw,fx,bp,q,r,tol1,tol2,u,v,w,x,xm,BX
      external f

	  CGOLD=.3819660
	  ZEPS=1.0e-10
	  tol=1.e-4
	  itmax=500
	  BX=(AX+CX)/2.

!     Given a function f, and given a bracketing triplet of abscissas AX, bx, cx (such that bx is
!     between AX and cx, and f(bx) is less than both f(AX) and f(CX)), this routine isolates
!     the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
!     the minimum is returned as xmin, and the minimum function value is returned as brent,
!     the returned function value.
!     Parameters: MAXimum allowed number of iterations; golden ratio; and a small number that
!     protects against trying to achieve fractional accuracy for a minimum that happens to be
!     exactly zero.

	  a=min(AX,CX)        !a and b must be in ascending order, though the input
	  b=max(AX,CX)    !abscissas need not be.
	  v=BX          !Initializations...
	  w=v
	  x=v
	  e=0.    !This will be the distance moved on the step before last.
	  fx=f(x)
	  fv=fx
	  fw=fx
      do iter=1,itmax   !Main program loop.
!       write(*,*) 'iteration', iter
        xm=0.5*(a+b)
        tol1=(tol*abs(x))+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-0.5*(b-a))) goto 3    !Test for done here.
        if(abs(e).gt.tol1) then        !Construct a trial parabolic fit.
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          bp=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) bp=-bp
          q=abs(q)
          etemp=e
          e=dd
	  if(abs(bp).ge.abs(0.5*q*etemp).or.bp.le.q*(a-x).or.
     &bp.ge.q*(b-x)) goto 1
	  !The above conditions determine the acceptability of the parabolic fit. Here it is o.k.:
	  dd=bp/q !Take the parabolic step.
	  u=x+dd
	  if(u-a.lt.tol2 .or. b-u.lt.tol2) dd=sign(tol1,xm-x)
	  goto 2     !Skip over the golden section step.
        endif
1       if(x.ge.xm) then    !We arrive here for a golden section step, which we take
          e=a-x     !into the larger of the two segments.
        else
          e=b-x
        endif
        dd=CGOLD*e           !Take the golden section step.
2       if(abs(dd).ge.tol1) then     !Arrive here with dd computed either from parabolic fit, or
          u=x+dd              !else from golden section.
        else
          u=x+sign(tol1,dd)
        endif
        fu=f(u) !This is the one function evaluation per iteration,
        if(fu.le.fx) then    !and now we have to decide what to do with our function
          if(u.ge.x) then     !evaluation. Housekeeping follows:
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif

          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif      !Done with housekeeping. Back for another iteration.
      enddo
      write(*,*) 'brent exceed maximum iterations'
3     p=x    !Arrive here ready to exit with best values.
      brent=fx
      return

      end
C
C***********************************************************************
C
      FUNCTION FUNCP(P)
      SAVE
      COMMON/NODDAT/PNOD(10000),TNOD(10000),NOD,FITNOD
      COMMON/INFON/TEM,VIS,TLITH,DENS,ALP,CP,G,COND,D,PFAC,
     C DZ,TP,HSCALE,TCRUST,MMOD,CRUSTK,igrid,min,iwrit
     C ,TUCRUST,TLCRUST,HUCRUST,HLCRUST,HCRUST,DEPTHMX
      COMMON/CALC/TTHERM,TEMPL
	  COMMON/LITHDAT/DTM,TMOHO,FSURF,FMOHO
	  REAL*4 LAT,LON
	  REAL*8 P
      DIMENSION PREF(1000),TREF(1000),DREF(1000),TAD(1000),
     C TEMPLITH(1000),CONDO(1000),DLITH(1000)
      TEM=P
C
C CALCULATING GEOTHERM
C
      IF(MMOD.EQ.0) CALL TINIT(DREF,TREF,TAD,CONDO,NT)
      IF(MMOD.NE.0) CALL TINIT2(DREF,TREF,TAD,TEMPLITH,
     C   CONDO,DLITH,NT,NLITH)
      DO I=1,NT
        PREF(I)=PFAC*DREF(I)
      END DO
      SUM=0.
      N=0
      DO 20 J=1,NOD
        I=1
C
C     INTERPOLATION TO FIND TEMPERATURES AT PNOD(J)
C
5       I=I+1
        IF(I.GT.NT) GO TO 20
10      IF(PREF(I).LT.PNOD(J)) GO TO 5
        TINT=TREF(I-1)+(TREF(I)-TREF(I-1))*(PNOD(J)-PREF(I-1))/
     C   (PREF(I)-PREF(I-1))
        SUM=SUM+(TNOD(J)-TINT)**2
        N=N+1
20    CONTINUE
      FUNCP=SQRT(SUM/FLOAT(N))
      FITNOD=FUNCP
      WRITE(6,200) TEM,TTHERM,TLITH,FITNOD,TMOHO,TEMPL,FSURF,FMOHO
200   FORMAT('MECH.B.L.=',G15.8,' THERM.B.L.=',G15.8,
     C ' LITH.B.L.=',G15.8,' MISFIT=',F14.8,
     C ' MOHO.TEMP=',G15.8,' LAB.TEMP=',G15.8,' HFLUX=',G15.8,
     C ' MHFLUX=',G15.8)

      RETURN
      END
C
C***********************************************************************
C
C
C***********************************************************************
C
      FUNCTION FUNCV(V)
	  COMMON/PRESSDAT/press
      COMMON/NODDAT/PNOD(10000),TNOD(10000),NOD,FITNOD
      COMMON/INFON/TEM,VIS,TLITH,DENS,ALP,CP,G,COND,D,PFAC,
     C DZ,TP,HSCALE,TCRUST,MMOD,CRUSTK,igrid,min,iwrit
     C ,TUCRUST,TLCRUST,HUCRUST,HLCRUST,hmantle,DEPTHMX
      COMMON/CALC/TTHERM,TEMPL
	  REAL*8 V
      DIMENSION PREF(1000),TREF(1000),DREF(1000),TAD(1000),
     C TEMPLITH(1000),CONDO(1000),DLITH(1000)

	  FUNCV=ABS((130.e+9*(3./2.)*(V**(7./3.)-V**(5./3.))*(1.+(((3./4.)
     C     *(4.8-4.))*(V**(2./3.)-1.))))-press)

      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE OUTPUT(NOFJ)
      SAVE
C
C       CALCULATES AND PLOTS GEOTHERMS FOR MECHANICAL AND THERMAL
C             BOUNDARY LAYERS, AND FOR EQUIVALENT PLATE MODEL
C             WRITTEN BY DAN McKENZIE 12/88, AND USED IN
C         McKENZIE AND BICKLE, J. PETROLOGY, 29, 625-79 (1988)
C
C       MODIFIED FOR JAMES JACKSON TO PLOT WHOLE GEOTHERM AND LIST VALUES
C
      DIMENSION DREF(1000),TREF(1000),TAD(1000),TMEL(1000),PREF(1000)
     C ,XNOD(10000),YNOD(10000),REFT(4),DEPREF(4),TEMPLITH(1000),
     C CONDO(1000),DLITH(1000),VS(1000),DEPTHVS(200),SWAVE(200),
     C STEMP(200),STEMPFJ(200),VSFJ(1000),DVS(1000),TRFNOD(20),
     C xplt(2),yplt(2),tdiam(1000),ddiam(1000),PDIAM(1000)
      COMMON/INFON/TEM,VIS,TLITH,DENS,ALP,CP,G,COND,D,PFAC,
     C DZ,TP,HSCALE,TCRUST,MMOD,CRUSTK,igrid,min,iwrit
     C ,TUCRUST,TLCRUST,HUCRUST,HLCRUST,hmantle,DEPTHMX
      COMMON/CALC/TTHERM,TEMPL
      COMMON/NODDAT/PNOD(10000),TNOD(10000),NOD,FITNOD
      COMMON/LITHDAT/DTM,TMOHO,FSURF,FMOHO
      COMMON/TITL/TITLE,FILEINP
      COMMON/SWAVEFIL/FILEVS,GSIZE,DEPTHMN,IVS
      common/symsiz/rsym
      INTEGER*4 NOFJ,IVS
      CHARACTER*60 TITLE
      CHARACTER*60 FILEINP,FILEVS
      DATA REFT/300.,400.,500.,600./
      IF(IVS.EQ.1) THEN
C
C       READING IN SWAVE VELOCITIES AND DEPTHS, CALCULATING T AND TFJ
C       IF IVS=1
C
        OPEN(UNIT=21,FILE=FILEVS)
        PERIOD=100.
        IV=0
3       CONTINUE
        READ(21,160,END=5) DEP,SW
160     FORMAT(15X,F12.5,10X,F12.5)
        IF(DEP.LT.DEPTHMN) GO TO 3
        IV=IV+1
        DEPTHVS(IV)=DEP
        SWAVE(IV)=SW
        CALL VS2TEMP(PERIOD,GSIZE,DEPTHVS(IV),SWAVE(IV),STEMP(IV))
        WRITE(6,165) DEPTHVS(IV),SWAVE(IV),STEMP(IV)
        WRITE(16,165) DEPTHVS(IV),SWAVE(IV),STEMP(IV)
165     FORMAT(' Depth=',F12.5,' Vs=',F12.5,' T(me)=',F12.5)
        IF(DEPTHVS(IV).LE.DEPTHMX) GO TO 3
C
5       CONTINUE
        CLOSE(21)
        IV=IV-1
      ENDIF
C
C CALCULATING GEOTHERM
C

      OPEN(UNIT=14,FILE='geoth.out')
      IF(MMOD.EQ.0) CALL TINIT(DREF,TREF,TAD,CONDO,NT)
      IF(MMOD.NE.0) CALL TINIT2(DREF,TREF,TAD,TEMPLITH,
     C  CONDO,DLITH,NT,NLITH)
      NE=NT
      NMAX=IFIX(DEPTHMX+0.1)
      DO I=1,NMAX
        PREF(I)=PFAC*DREF(I)
        IF(I.EQ.1) THEN
          TMEL(1)=1373.
        ELSE
          TMEL(I)=TMEL(I-1)+273.
        ENDIF
        CALL TSOLIQ(PREF(I),TMEL(I),DTS,TL,DTL)
        TMEL(I)=TMEL(I)-273.1
        IF(IWRIT.EQ.1) THEN
		  WRITE(14,180) DREF(I),TREF(I),TMEL(I),TAD(I),CONDO(I)
180       FORMAT(F12.5,' ',F12.5,' ',F12.5,' ',
     C    F12.5,' ',F12.5)
        ENDIF
      END DO

      RETURN
      END
C
C*******************************************************************
C
      SUBROUTINE TINIT(DREF,TREF,TAD,CONDO,NT)
      SAVE
C
C     CALCULATES THE INITIAL TEMP USING RICHTER & MCK
C     PARAMETERIZATION JGR 86 1738-1744 1981 OF THE
C     THERMAL B.L., GIVEN BY A(5).
C     NEEDS THE THICKNESS OF THE MECH B.L. TEM IN KM
C     THE VISCOSITY VIS,& THE POT. TEMP OF THE MANTLE
C     BELOW IN DEG C, TP.
C     RETURNS THE INITIAL TEMP, TREF AT DEPTH INTERVALS DZ
C     BOTH IN KM,THE THICKNESS OF THE THERMAL B.L. TTHERM
C     & THE APPARENT PLATE THICKNESS TLITH ALSO IN KM.
C     THE CALCS ARE IN RMKS
C
C TLIT IS THE LITHOSPHERE THICKNESS IN KM,MEASURED FROM
C OCEANIC BATHYMETRY.TTHERM IS THE THICKNESS OF THE THERMAL
C BOUNDARY LAYER.TEM IS THE THICKNESS OF THE MECHANICAL B.L.IN km.
C TTHERM AND TLIT ARE CALCULATED FROM B.L. THEORY
C
C
      DIMENSION DREF(1000),TREF(1000),TAD(1000),TADIF(1000),A(5)
	  DIMENSION CONDO(1000)
      COMMON/INFON/TEM,VIS,TLITH,DENS,ALP,CP,G,COND,D,PFAC,
     C DZ,TP,HSCALE,TCRUST,MMOD,CRUSTK,igrid,min,iwrit
     C ,TUCRUST,TLCRUST,HUCRUST,HLCRUST,hmantle,DEPTHMX
      COMMON/CALC/TTHERM,TEMPL
      TMECH=TEM*1.E+3
      T0=273.
      THI=TP+T0
      DIFF=COND/(DENS*CP)
      A(1)=-1.84
      A(2)=-1.16
      A(3)=-0.28
      A(4)=0.15
      A(5)=0.01
      FAC=EXP((ALP*G*TMECH)/CP)
C
C     MATCHES THE CONDUCTIVE HEAT FLUX THROUGH THE MECHANICAL
C     BOUNDARY LAYER TO THE CONDUCTIVE AND CONVECTIVE HEAT FLUX
C     IN THE THERMAL BOUNDARY LAYER BY ITERATIVE ADJUSTMENTS
C     OF DTM, THE TEMPERATURE DROP ACROSS THE MECHANICAL BOUNDARY
C     LAYER
C
C     DTHT IS THE POT. TEMP DIFF ACROSS THE THERMAL B.L..
C     DTM IS THE REAL TEMP DIFF ACROSS THE MECH B.L..
C     F IS THE HEAT FLUX THROUGH THE MECH B.L..
C     FC IS THE HEAT FLUX IN THE CONVECTIVE REGION
C     THAT IS NOT CARRIED BY THE ADIABATIC TEMP GRADIENT.
C     THI IS THE POT TEMP IN DEG K IN THE MANTLE INTERIOR.
C     RA IS THE RAYLEIGH NO..
C
      DTHT=0.
      DTM=(THI-DTHT)*FAC-T0
10    F=COND*DTM/TMECH
      GRADTS=(ALP*G/CP)*THI*FAC
      FC=F-COND*GRADTS
      RA=(G*ALP*D**4*FC)/(DIFF*COND*VIS)
c
C     RELATION BETWEEN TEMPERATURE DROP AND RAYLEIGH NO. FROM
C     McK AND RICHTER
C
      DTHT=FC*D*1.84/(COND*RA**(0.219))
      THM=THI-DTHT
      DTM2=THM*FAC-T0
      ERR=ABS(DTM2-DTM)
      DTM=DTM2
      IF(ERR.GT.0.1) GO TO 10
C
C     END OF ITERATION
C
C---------------------------------------------------------------------
C                   CALCULATION OF GEOTHERM
C---------------------------------------------------------------------
C
C     MECHANICAL BOUNDARY LAYER
C
      IE=TMECH/(DZ*1.E+3)
      DO I=1,IE
        DREF(I)=(I-1)*DZ*1.E+3
        TAD(I)=THI*EXP((ALP*G*DREF(I))/CP)-T0
        TREF(I)=DTM*DREF(I)/TMECH
      END DO
C
C     THERMAL BOUNDARY LAYER AND ADIABATIC INTERIOR
C
C     ZTW(IDDLE) IS A SCALED DEPTH DEFINED IN McK AND RICHTER
C
      DZTW=1.45*(RA**(0.219))*1.E+3*DZ/D
      IT=8./DZTW+1
C
C     DEFINITION OF THERMAL BOUNDARY LAYER THICKNESS
C
      TTHERM=DZ*(IT-1)*5./8.
C
      IT=20./DZTW+1
      IF(IT.GT.200) IT=200
      DO I=1,IT
        ZTW=DZTW*(I-1)
        DREF(I+IE)=(I-1)*1.E+3*DZ+TMECH
        TAD(I+IE)=THI*EXP((ALP*G*DREF(I+IE))/CP)-T0
        SUM=0
        DO J=1,5
          F=1.
          IF(J.GT.1) F=ZTW**(J-1)
          SUM=SUM+A(J)*F
        END DO
        F=1.
        IF(ZTW.GT.0.) F=EXP(-ZTW)
        SUM=1.84+F*SUM
        TREF(I+IE)=EXP((ALP*G*DREF(I+IE))/CP)*(THM+((FC*
     C   D*(1./RA**(0.219)))/COND)*SUM)-T0
      END DO
      NT=IE+IT

C
C       Conductivity and heat flux arrays
C
      DO I=1,NT
        CALL DERIVS(DREF(I),TREF(I),HS,DYDX,COND,FLUX)
		CONDO(I)=COND
      END DO
C
C-----------------------------------------------------------------------
C
C     CONVERTS DEPTH DREF TO KM, CALCULATES DIFFERENCE BETWEEN
C     GEOTHERM TREF AND ADIABATIC TAD TEMPERATURES, TO FIND
C     EQUIVALENT LITHOSPHERIC THICKNESS TLITH
C
      DO I=1,NT
        DREF(I)=DREF(I)/1.E+3
        TADIF(I)=TAD(I)-TREF(I)
      END DO
C
C     INTEGRATION TO FIND TLITH
C
      SUM=0.
      DO I=2,NT
        SUM=SUM+(DREF(I)-DREF(I-1))*(TADIF(I)+TADIF(I-1))/2
      END DO
      TLITH=2.*SUM/TADIF(1)
C
      WRITE(6,200) RA,TTHERM,TLITH,DTM
200   FORMAT(' R.NO.=',G15.8,' THERM.B.L.=',G15.8,' km'/
     C ' LITH.B.L.=',G15.8,' REAL TEMP AT THE BASE OF MECH.L.=',
     C G15.8)
      TEMPL=THI*EXP((ALP*G*TLITH*1.E+3)/CP)-T0
      F=(COND*DTM*1.E+3)/TMECH
      WRITE(6,205) TEMPL,IT,DZ,F
205   FORMAT(' REAL TEMP AT BASE OF LITHOSP.=',F12.5,' IT=',I3,
     C ' DZ=',F12.5,' HFLUX=',F12.5,' mW m**-2')
      RETURN
      END
C
C*******************************************************************
C
      SUBROUTINE TSOLIQ(P,TS,DTS,TL,DTL)
      SAVE
C
C GIVEN P IN GPA FINDS THE SOLIDUS TEMP TS,
C ITS DERIVATIVE WRT PRESSURE DTS WHEN THE PRESSURE
C IS IN GPA,THE LIQUIDUS TEMP TL,& ITS PRESSURE
C DERIVATIVE WHEN THE PRESSURE IS IN GPA
C ALL CONSTANTS FROM McKENZIE AND BICKLE,
C J. PETROLOGY, 29, 625-79 (1988)
C       TEMPERATURES ARE IN KELVIN
C
C     LIQUIDUS TL
C
      TL=2009.4+4.343*P+180.*ATAN(P/2.2169)
      DTL=4.343+(180./2.2169)/(1.+(P/2.2169)**2)
C
C     SOLIDUS TS
C
C FINDS TS BY ITERATION,GIVEN P,THEN FINDS DTS=DTS/DP
C
      T0=1373.1
      B=1.2E-2
      A=4.968E-4
      C=136.05
C
C     USES PREVIOUS VALUE OF TS AS STARTING POINT FOR ITERATION
C
      T1=TS
      N=0
10    PD=(T1-T0)/C+A*EXP(B*(T1-T0))-P
      N=N+1
      DPDODT=1./C+A*B*EXP(B*(T1-T0))
      DT=-PD/DPDODT
      T1=T1+DT
      IF((ABS(DT).GT.1.E-2).AND.(N.LT.10)) GO TO 10
      IF(N.EQ.10) WRITE(6,100) T1,DT,N
100   FORMAT(' T1=',G10.4,' DT=',G10.4,' N=',I3)
      TS=T1
      DPODT=1./C+A*B*EXP(B*(TS-T0))
      DTS=1./DPODT
      RETURN
      END
C
C*******************************************************************
C
      SUBROUTINE TINIT2(DREF,TREF,TAD,TEMPLITH,CONDO,DLITH,NT,NLITH)
      SAVE
C
C       CHECK DZ IS IN KM
C
C CALCULATES THE INITIAL TEMP USING RICHTER & MCK
C PARAMETERIZATION JGR 86 1738-1744 1981 OF THE
C THERMAL B.L..NEEDS THE THICKNESS OF THE MECH B.L. TEM IN KM
C THE VISCOSITY VIS,& THE POT TEMP OF THE MANTLE
C BELOW IN DEG C TP
C BELOW THE MOHO, HSCALE IS THE LENGTH SCALE IN KM OVER WHICH THE
C CRUST RADIOACTIVE HEATING DECAYS, HEAT(Z)=HEAT(SURFACE)*EXP(-Z/HSCALE)
C RETURNS THE INITIAL TEMP TREF AT DEPTH INTERVALS DZ
C BOTH IN KM,THE THICKNESS OF THE THERMAL B.L. TTHERM
C & THE APPARENT PLATE THICKNESS TLIT ALSO IN KM
C THE CALCS ARE IN RMKS
C
      DIMENSION DREF(1000),TREF(1000),TAD(1000),CONDO(1000),
     C  TEMPLITH(1000),HFLUX(1000),A(5),DLITH(1000)
      COMMON/INFON/TEM,VIS,TLITH,DENS,ALP,CP,G,COND,D,PFAC,
     C DZ,TP,HSCALE,TCRUST,MMOD,CRUSTK,igrid,min,iwrit
     C ,TUCRUST,TLCRUST,HUCRUST,HLCRUST,hmantle,DEPTHMX
      COMMON/LITHDAT/DTM,TMOHO,FSURF,FMOHO
      COMMON/CALC/TTHERM,TEMPL
      TMECH=TEM*1.E+3
      T0=273.
      THI=TP+T0
      D=7.E+5
      A(1)=-1.84
      A(2)=-1.16
      A(3)=-0.28
      A(4)=0.15
      A(5)=0.01
      FAC=EXP((ALP*G*TMECH)/CP)
      DTM=THI*FAC-T0
C
C       INITIALISING SEARCH FOR VALUE OF HS, THE SURFACE HEAT FLUX THAT IS
C       COMPATIBLE WITH THE GIVEN VALUES OF  TMECH, TP AND VIS
C
C       CONDUCTIVITY AT TEMP OF DTM/2.
C
      HS=10.e-3
      X=TMECH/2.
      YT=DTM/2.
C      CALL DERIVS(X,YT,HS,DYT,COND,FLUX)
C
C       CALCULATE STARTING VALUE FOR HS
C
      HS=TUCRUST*1.E+3*HUCRUST+TLCRUST*1.E+3*HLCRUST+HS
c      HS=2.*COND*DTM/TMECH
C
C DTHT IS THE POT. TEMP DIFF ACROSS THE THERMAL B. L.
C DTM IS THE TEMP DIFF ACROSS THE MECH B.L.
C FLUX IS THE HEAT FLUX IN THE MBL
C FC IS THE HEAT FLUX CARRIED BY CONVECTION IN THE CONVECTIVE REGION,
C THAT WHICH IS NOT CARRIED BY CONDUCTION DOWN THE ADIABATIC TEMP GRADIENT
C THI IS THE POT TEMP IN DEG K IN THE MANTLE INTERIOR
C
      H=DZ*1.E+3
      N=TMECH/H
      IF(N*H.GT.0.01*H) N=N+1
10    CONTINUE
      X=0.
      Y=0.
      DREF(1)=X
      TREF(1)=Y
C
C       INTEGRATE TO FIND TEMP AND FLUX AT BOTTOM OF MBL
C
      CALL DERIVS(X,Y,HS,DYDX,COND,FLUX)
      DO I=1,N
        H=DZ*1.E+3
        X=DREF(I)
        DREF(I+1)=I*H
        IF(DREF(I+1).GT.TMECH) THEN
          DREF(I+1)=TMECH
          H=DREF(I+1)-DREF(I)
        ENDIF
        HH=H*0.5
        H6=H/6.
        XH=X+HH
        YT=Y+HH*DYDX
        CALL DERIVS(XH,YT,HS,DYT,COND,FLUX)
        YT=Y+HH*DYT
        CALL DERIVS(XH,YT,HS,DYM,COND,FLUX)
        YT=Y+H*DYM
        DYM=DYT+DYM
        CALL DERIVS(X+H,YT,HS,DYT,COND,FLUX)
        YOUT=Y+H6*(DYDX+DYT+2.*DYM)
        Y=YOUT
        TREF(I+1)=YOUT
        DYDX=DYT
      END DO
C      WRITE(6,195) (DREF(I)*1.E-3,TREF(I),I=1,N+1)
C      WRITE(16,195) (DREF(I)*1.E-3,TREF(I),I=1,N+1)
195   FORMAT(' ',18(1X,F12.5))
C
C       COND IS CONDUCTIVITY AT THE TEMPERATURE OF THE BASE OF THE MBL
C
      DIFF=COND/(DENS*CP)
      GRADTS=(ALP*G/CP)*THI*FAC
      FC=FLUX-COND*GRADTS
      IF(FC.GT.0.) THEN
        RA=(G*ALP*D**4*FC)/(DIFF*COND*VIS)
        DTHT=FC*D*1.84/(COND*RA**(0.219))
      ELSE
        FC=0.
        RA=0.
        DTHT=0.
      ENDIF
      THM=THI-DTHT
C
C       DTMTOP IS TEMPERATURE AT TOP OF TBL
C
      DTMTOP=THM*FAC-T0
C
C       DTMBOT IS TEMPERATURE AT BOTTOM OF MBL
C
C      WRITE(6,196) N
C      WRITE(16,196) N
196   FORMAT(' N=',I6)
      DTMBOT=TREF(N+1)
      TERR=DTMTOP-DTMBOT
      HSI=HS
      HS=HS*(1.+0.1*(TERR/DTMBOT))
c      HS=HS*(1.+0.7*(TERR/DTMBOT))
C
C       ITERATING IF ERR.GT.0.1
C
      ERR=ABS(TERR)
C      WRITE(6,198) HSI,HS,TERR
C      WRITE(16,198) HSI,HS,TERR
198   FORMAT(' HSI=',G15.8,' HS=',G15.8,' TERR=',G15.8)
      IF(ERR.GT.0.1) GO TO 10
C
C       CALCULATING ISENTROPIC TEMP IN MBL AND BOTH TEMPS IN TBL
C
      IE=N+1
      DO I=1,IE
        TAD(I)=THI*EXP((ALP*G*DREF(I))/CP)-T0
      END DO
      DZTW=1.45*(RA**(0.219))*1.E+3*DZ/D
      IT=8./DZTW+1
      TTHERM=DZ*(IT-1)*5./8.
      IT=20./DZTW+1
      NMAX=IFIX(DEPTHMX+0.1)
      IF(IT.GT.200) IT=200
      IF((IT+IE).GT.1000) IT=1000-IE
      IF((IT+IE).LT.NMAX) IT=NMAX-IE
      IF(IT.GE.1) THEN
        DO I=1,IT
          ZTW=DZTW*(I-1)
          DREF(I+IE)=(I-1)*1.E+3*DZ+TMECH
          TAD(I+IE)=THI*EXP((ALP*G*DREF(I+IE))/CP)-T0
          SUM=0
          DO J=1,5
            F=1.
            IF(J.GT.1) F=ZTW**(J-1)
            SUM=SUM+A(J)*F
          END DO
          F=1.
          IF(ZTW.GT.0.) F=EXP(-ZTW)
          SUM=1.84+F*SUM
          TREF(I+IE)=EXP((ALP*G*DREF(I+IE))/CP)*(THM+((FC*
     C     D*(1./RA**(0.219)))/COND)*SUM)-T0
        END DO
      ELSE
        IT=0
      ENDIF
      NT=IE+IT

C
C       Conductivity and heat flux arrays
C
      DO I=1,NT
        CALL DERIVS(DREF(I),TREF(I),HS,DYDX,COND,FLUX)
		CONDO(I)=COND
      END DO
c
c       INTEGRATING CONDUCTIVE GEOTHERM BELOW MBL TO FIND
C       LITHOSPHERIC THICKNESS
C
      DLITH(1)=0.
      TEMPLITH(1)=0.
      HFLUX(1)=HS
      X=DLITH(1)
      Y=TEMPLITH(1)
      CALL DERIVS(X,Y,HS,DYDX,COND,FLUX)
	  CONDO(1)=COND
      DO I=1,NT
        H=DZ*1.E+3
        X=DLITH(I)
        HH=H*0.5
        H6=H/6.
        XH=X+HH
        YT=Y+HH*DYDX
        CALL DERIVS(XH,YT,HS,DYT,COND,FLUX)
        YT=Y+HH*DYT
        CALL DERIVS(XH,YT,HS,DYM,COND,FLUX)
        YT=Y+H*DYM
        DYM=DYT+DYM
        CALL DERIVS(X+H,YT,HS,DYT,COND,FLUX)
        YOUT=Y+H6*(DYDX+DYT+2.*DYM)
        Y=YOUT
        TEMPLITH(I+1)=YOUT
        DLITH(I+1)=DLITH(I)+DZ*1.E+3
        HFLUX(I+1)=FLUX
        IF(TEMPLITH(I+1).GT.TAD(I+1)) GO TO 20
        DYDX=DYT
      END DO
C
C       FIND INTERSECTION OF TEMPLITH AND TAD
C
20      CONTINUE
        NLITH=I+1
        TLITH=DLITH(I)+1./(1.-(TEMPLITH(I+1)-TAD(I+1))/
     C   (TEMPLITH(I)-TAD(I)))
C
C CONVERTS TO KM
C
      DO I=1,NT
        DREF(I)=DREF(I)/1.E+3
        DLITH(I)=DLITH(I)/1.E+3
      END DO
      TLITH=TLITH/1.E+3
C
C       MOHO TEMPERATURE
C
      I=IFIX(TCRUST/DZ)+1
      TMOHO=TREF(I)
      FSURF=HS*1.E+3
      FMOHO=HFLUX(I)*1.E+3
      DTM=DTMBOT
      WRITE(6,200) RA,TTHERM,TLITH,DTM,NT,FSURF,FMOHO
200   FORMAT(' RAY #=',G15.8,' THERM.B.L.=',G15.8,
     C ' km, LITH.B.L.=',G15.8/
     C ' REAL TEMP AT BASE MECH.L.=',G15.8,' NT=',I3/
     C ' SURF. HFLUX=',F12.5,' MOHO HEAT FLUX=',F12.5,' mW m**-2')
      TEMPL=THI*EXP((ALP*G*TLITH*1.E+3)/CP)-T0
      WRITE(6,205) TEMPL,IT,DZ
205   FORMAT(' REAL TEMP AT BASE OF LITH.=',F12.5,' IT=',I3,' DZ=',F5.2)
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE DERIVS(Z,T,HS,DT,CON,FLUX)
      SAVE
C
C       DEPTH Z IN METRES, T TEMPERATURE AT DEPTH Z IN DEG C
C       TCRUST IS THE CRUSTAL THICKNESS IN KM
C
      EXTERNAL FUNCV
	  COMMON/PRESSDAT/press
      COMMON/INFON/TEM,VIS,TLITH,DENS,ALP,CP,G,COND,D,PFAC,
     C DZ,TP,HSCALE,TCRUST,MMOD,CRUSTK,igrid,min,iwrit
     C ,TUCRUST,TLCRUST,HUCRUST,HLCRUST,hmantle,DEPTHMX
C
C       TEMPERATURE DEPENDENT CONDUCTIVITY
C
      ABT=T+273.
      IF(MMOD.EQ.0) THEN
        COND=2.5
        IF(Z.GT.(TCRUST*1.E+3)) COND=3.11
      ELSEIF(MMOD.EQ.1) THEN
        COND=2.5
        IF(Z.GT.(TCRUST*1.E+3)) THEN
          RADCON=1.753E-2-1.0365E-4*ABT+2.2451E-7*ABT**2
     C    -3.4071E-11*ABT**3
c
c          COND=4.7*(298./ABT)**0.4
c     C     *EXP(-5.1333*4.E-5*(ABT-298.))+RADCON
C
C       FROM ANNE HOFMEISTER, RADCON FROM HER EXPRESSION,
C       FIRST TERM IN COND CHANGED TO ALLOW INTEGRATION,
C       THE COEFFICIENT OF 5.3 FITS HER PLOT (BY EYE)
C
        COND=5.3/(1.+0.0015*T)+RADCON

        ENDIF
      ELSEIF(MMOD.EQ.2) THEN
        COND=2.5
        IF(Z.GT.(TCRUST*1.E+3)) THEN
          RADCON=0.368E-9*ABT**3
          COND=1./(0.174+0.000265*ABT)+RADCON
        ENDIF
      ELSEIF(MMOD.EQ.3) THEN
        IF(Z.le.(TCRUST*1.E+3)) THEN
c
c       temperature dependence of cond. for crustal rocks
c       Whittington et al Nature 458 319-321 (2009)
c       misprint corrected
c
          if(abt.le.846.) then
            diffk=(567.3/abt-0.062)*1.e-6
            cpk=199.5+0.0857*abt-5.e+6/(abt*abt)
            cpk=cpk*1.e+3/221.78
          else
            diffk=(0.732-0.000135*abt)*1.e-6
            cpk=229.32+0.0323*abt-47.9e-6/(abt*abt)
            cpk=cpk*1.e+3/221.78
          endif
          cond=diffk*cpk*2.7e+3
c
c       I think the average thermal conductivity of the crust is too high
c       because of the proportion of quartz used.  I therefore reduced it
c       to a value closer to that of schist
c
          cond=0.8*cond
        elseif(Z.gt.(TCRUST*1.E+3)) then
          RADCON=1.753E-2-1.0365E-4*ABT+2.2451E-7*ABT**2
     C    -3.4071E-11*ABT**3
c
c          COND=4.7*(298./ABT)**0.4
c     C     *EXP(-5.1333*4.E-5*(ABT-298.))+RADCON
C
C       FROM ANNE HOFMEISTER, RADCON FROM HER EXPRESSION,
C       FIRST TERM IN COND CHANGED TO ALLOW INTEGRATION,
C       THE COEFFICIENT OF 5.3 FITS HER PLOT (BY EYE)
C
          COND=5.3/(1.+0.0015*T)+RADCON
        endif
      ELSEIF(MMOD.EQ.4) THEN
        COND=2.5
        if(Z.gt.(TCRUST*1.E+3)) then
c
c       expression from Osako et al PEPI 143-144 pp 311-320 2004
c        pressure in GPa
c
          press=z/(30.7*1.e+3)
          cond=(1.175+1264./abt)*exp(0.038*press)
        endif
      ELSEIF(MMOD.EQ.5) THEN
        IF(Z.le.(TCRUST*1.E+3)) THEN
c
c       temperature dependence of cond. for crustal rocks
c       Whittington et al Nature 458 319-321 (2009)
c       misprint corrected
c
          if(abt.le.846.) then
            diffk=(567.3/abt-0.062)*1.e-6
            cpk=199.5+0.0857*abt-5.e+6/(abt*abt)
            cpk=cpk*1.e+3/221.78
          else
            diffk=(0.732-0.000135*abt)*1.e-6
            cpk=229.32+0.0323*abt-47.9e-6/(abt*abt)
            cpk=cpk*1.e+3/221.78
          endif
          cond=diffk*cpk*2.7e+3
c
c       I think the average thermal conductivity of the crust is too high
c       because of the proportion of quartz used.  I therefore reduced it
c       to a value closer to that of schist
c
          cond=0.8*cond
        elseif(Z.gt.(TCRUST*1.E+3)) then
c
c       expression from Osako et al PEPI 143-144 pp 311-320 2004
c        pressure in GPa
c
          press=z/(30.7*1.e+3)
          cond=(1.175+1264./abt)*exp(0.038*press)
        endif
        ELSEIF(MMOD.EQ.6) THEN
          IF(Z.le.(TCRUST*1.E+3)) THEN
c
c       pressure and temperature dependence of cond. for crustal rocks
c       from Grose & Afonso 2013
c
  		      press=(z*1.e+6)/(30.)
		      dV0=(-7.334963115431564676e-23*press**2)
     c      +(7.510867653681621105e-12*press)+1.000184023114681908
		    alphaP0=dV0*exp((6.+1.)*((dV0**(-1.))-1.))
		    rhoP0=2950.*dV0
		    alphaT0=(1.639e-5*(abt-273.))+((1.322e-8/2.)
     c      *((abt**2.)-(273.**2.)))
		    densk=rhoP0*(1.-(alphaP0*alphaT0))
		    cpk=((0.15*((1.6108e3)+(-1.24788e4*(abt**(-1./2.)))
     c      +(-1.728477e9*(abt**(-3.)))))+(0.2*(2.1715e3+(-4.555e-1*abt)
     c      +(1.1332e6*(abt**(-2)))+(-2.22716e4*(abt**(-1./2.)))
     c      +(1.299e-4*(abt**2.))))+(0.65*(1.85757e3+(-3.324e-1*abt)
     c      +(-5.061e6*(abt**(-2)))+(-1.64946e4*(abt**(-1./2.)))
     c      +(1.505e-4*(abt**2.)))))
			diffk=(0.432+(0.44*exp(-(abt-273.)/380.))
     c      +(0.305*exp(-(abt-273.)/145.)))*1e-6
			gs=0.5
			dkP=0.05
			Ar=(1.8*(1.-exp(-((gs**(1.3))/0.15))))
     c      -(1.-exp(-((gs**(0.5))/5.)))
			Br=(11.7*exp(-d/0.159))+(6.*exp(-((d**(3.))/10.)))
			Tar=490.+(1850.*exp(-((gs**0.315)/0.825)))+(875.*exp(-gs/0.18))
			Tbr=2700.+(9000.*exp(-((gs**0.5)/0.205)))
			xa=167.5+(505.*exp(-((gs**0.5)/0.85)))
			xb=465.+(1700.*exp(-((gs**0.94)/0.175)))
			RADCON=Ar*exp(-(((abt-Tar)**(2.))/(2*(xa**(2.)))))+Br
     c      *exp(-(((abt-Tbr)**(2.))/(2*(xb**(2.)))))
			cond=(diffk*cpk*densk*exp(dkP*(press/1.e+9)))+RADCON

          elseif(Z.gt.(TCRUST*1.E+3)) then
c
c       pressure and temperature dependence of cond. for mantle rocks
c       from Grose & Afonso 2013
c
  			press=(z*1.e+6)/(30.)
		      dV0=(-7.334963115431564676e-23*press**2)
     c      +(7.510867653681621105e-12*press)+1.000184023114681908
		    alphaP0=dV0*exp((6.+1.)*((dV0**(-1.))-1.))
		    rhoP0=3330.*dV0
		    alphaT0=(2.832e-5*(abt-273.))+
     c      ((0.758e-8/2.)*((abt**2.)-(273.**2.)))
		    densk=rhoP0*(1.-(alphaP0*alphaT0))
		    cpk=(1580.-(12230.*(abt**(-1./2.)))-(1694.e6*(abt**(-3.))))
			diffk=(0.565+(0.67*exp(-(abt-273.)/590.))
     c      +(1.4*exp(-(abt-273.)/135.)))*1e-6
			gs=0.5
			dkP=0.05
			Ar=(1.8*(1.-exp(-((gs**(1.3))/0.15))))
     c      -(1.-exp(-((gs**(0.5))/5.)))
			Br=(11.7*exp(-d/0.159))+(6.*exp(-((d**(3.))/10.)))
			Tar=490.+(1850.*exp(-((gs**0.315)/0.825)))+(875.*exp(-gs/0.18))
			Tbr=2700.+(9000.*exp(-((gs**0.5)/0.205)))
			xa=167.5+(505.*exp(-((gs**0.5)/0.85)))
			xb=465.+(1700.*exp(-((gs**0.94)/0.175)))
			RADCON=Ar*exp(-(((abt-Tar)**(2.))/(2*(xa**(2.)))))+Br
     c      *exp(-(((abt-Tbr)**(2.))/(2*(xb**(2.)))))
			cond=(diffk*cpk*densk*exp(dkP*(press/1.e+9)))+RADCON
          endif
      ELSEIF(MMOD.EQ.7) THEN
        IF(Z.le.(TCRUST*1.E+3)) THEN
            cond=CRUSTK !Variable crustal value
        ELSE
  			press=(z*1.e+6)/(30.)
            ! Fitted function = -alpha*np.exp(-beta*P[GPa])+gamma  beta = 0.021982, and gamma = 1.3437
			dV0=-0.33633*exp(-0.022523*press*1e-9)+1.3364
		    alphaP0=dV0*exp((6.+1.)*((dV0**(-1.))-1.))
		    rhoP0=3330.*dV0
		    alphaT0=(2.832e-5*(abt-273.))+
     c      ((0.758e-8/2.)*((abt**2.)-(273.**2.)))
		    densk=rhoP0*(1.-(alphaP0*alphaT0))
		    cpk=(1580.-(12230.*(abt**(-1./2.)))-(1694.e6*(abt**(-3.))))
			diffk=(0.565+(0.67*exp(-(abt-273.)/590.))
     c      +(1.4*exp(-(abt-273.)/135.)))*1e-6
			gs=0.5
			dkP=0.05
			Ar=(1.8*(1.-exp(-((gs**(1.3))/0.15))))
     c      -(1.-exp(-((gs**(0.5))/5.)))
			Br=(11.7*exp(-d/0.159))+(6.*exp(-((d**(3.))/10.)))
			Tar=490.+(1850.*exp(-((gs**0.315)/0.825)))+(875.*exp(-gs/0.18))
			Tbr=2700.+(9000.*exp(-((gs**0.5)/0.205)))
			xa=167.5+(505.*exp(-((gs**0.5)/0.85)))
			xb=465.+(1700.*exp(-((gs**0.94)/0.175)))
			RADCON=Ar*exp(-(((abt-Tar)**(2.))/(2*(xa**(2.)))))+Br
     c      *exp(-(((abt-Tbr)**(2.))/(2*(xb**(2.)))))
			cond=(diffk*cpk*densk*exp(dkP*(press/1.e+9)))+RADCON
        ENDIF
      ELSEIF(MMOD.EQ.8) THEN
        IF(Z.le.(TCRUST*1.E+3)) THEN
            hasterok_press=(z*1.e-3)/(30.)
            hasterok_n=6.4-(2.3*log(CRUSTK))
            hasterok_beta=0.1
            hasterok_k1=CRUSTK*((hasterok_n-1.)/hasterok_n)
     c       *(1.+hasterok_beta*hasterok_press)
            hasterok_k2=CRUSTK*(exp(-(T - 25.)/300.))
     c       *(1./hasterok_n)*(1.+hasterok_beta*hasterok_press)
            cond=hasterok_k1+hasterok_k2 !Variable temperature-dependent crustal value

        ELSE
  			press=(z*1.e+6)/(30.)
            ! Fitted function = -alpha*np.exp(-beta*P[GPa])+gamma  beta = 0.021982, and gamma = 1.3437
			! dV0=-0.33633*exp(-0.022523*press*1e-9)+1.3364
                  dV0=(-7.334963115431564676e-23*press**2)
     c      +(7.510867653681621105e-12*press)+1.000184023114681908
		    alphaP0=dV0*exp((6.+1.)*((dV0**(-1.))-1.))
		    rhoP0=3330.*dV0
		    alphaT0=(2.832e-5*(abt-273.))+
     c      ((0.758e-8/2.)*((abt**2.)-(273.**2.)))
		    densk=rhoP0*(1.-(alphaP0*alphaT0))
		    cpk=(1580.-(12230.*(abt**(-1./2.)))-(1694.e6*(abt**(-3.))))
			diffk=(0.565+(0.67*exp(-(abt-273.)/590.))
     c      +(1.4*exp(-(abt-273.)/135.)))*1e-6
			gs=0.5
			dkP=0.05
			Ar=(1.8*(1.-exp(-((gs**(1.3))/0.15))))
     c      -(1.-exp(-((gs**(0.5))/5.)))
			Br=(11.7*exp(-d/0.159))+(6.*exp(-((d**(3.))/10.)))
			Tar=490.+(1850.*exp(-((gs**0.315)/0.825)))+(875.*exp(-gs/0.18))
			Tbr=2700.+(9000.*exp(-((gs**0.5)/0.205)))
			xa=167.5+(505.*exp(-((gs**0.5)/0.85)))
			xb=465.+(1700.*exp(-((gs**0.94)/0.175)))
			RADCON=Ar*exp(-(((abt-Tar)**(2.))/(2*(xa**(2.)))))+Br
     c      *exp(-(((abt-Tbr)**(2.))/(2*(xb**(2.)))))
			cond=(diffk*cpk*densk*exp(dkP*(press/1.e+9)))+RADCON
        ENDIF
      ENDIF
      IF(Z.LE.(TCRUST*1.E+3)) THEN
        IF(Z.LE.(TUCRUST*1.E+3)) THEN
          FLUX=HS-Z*HUCRUST
        ELSE
          FLUX=HS-TUCRUST*1.E+3*HUCRUST-(Z-TUCRUST*1.E+3)*HLCRUST
        ENDIF
      ELSEIF(Z.GT.(TCRUST*1.E+3)) THEN
c
c        FLUX=HS-TCRUST*1.E+3*HCRUST
c
c       added mantle heat generation
c
        FLUX=HS-TUCRUST*1.E+3*HUCRUST-TLCRUST*1.E+3*HLCRUST
     c  -(z-tcrust*1.e+3)*hmantle
      ENDIF
      IF(FLUX.LT.0.) FLUX=0.
      DT=FLUX/COND
      CON=COND
c      WRITE(6,100) Z*1.e-3,FLUX,DT,COND
c      WRITE(16,100) Z*1.e-3,FLUX*1.e+3,DT,COND
100   FORMAT(' depth=',G15.8,' km, flux=',f12.5,', DT=',G15.8,
     c ' cond=',G15.8)
      RETURN
      END
C
C*************************************************************************
C
      SUBROUTINE TEMP2VS(PERIN,DIN,ZIN,TIN,VSO)
      SAVE
C
C       THE CONSTANTS IN MY EXPRESSIONS WERE CHANGED ON 22/12/05
C       GIVEN THE PERIOD , PERIN, IN SECS, AND THE GRAIN SIZE DIN IN METRES,
C       THE DEPTH Z IN KM, AND THE TEMPERATURE TIN IN DEG C
C       RETURNS THE SHEAR WAVE VELOCITY, VSO,
C       CALCULATED USING MY EXPRESSIONS
C
      COMMON/VSINFO/PERIOD,D,Z,VSFJ
      Z=ZIN
      PERIOD=PERIN
      D=DIN
      T=TIN
c
c       finds the S wave velocity at depth Z in km,T in deg C using my
c       empirical parameterisation
c
      CALL SVEL(Z,T,VSO)
      RETURN
C
      ENTRY VS2TEMP(PERIN,DIN,ZIN,VSIN,TO)
c
C       THE CONSTANTS IN MY EXPRESSIONS WERE CHANGED ON 12/4/05
c       finds the temperature at which the velocity is VS using
c       my expression, TO,
c
      Z=ZIN
      PERIOD=PERIN
      D=DIN
      VS=VSIN
      VSFJ=VS
      CALL TEMPV(Z,VS,TO)
      TOL=0.1
      T1=400.
      T2=1550.
      RETURN
      END
C
C************************************************************
C
      SUBROUTINE SVEL(Z,T,VS)
      SAVE
C
C       CALCULATES THE SWAVE VELOCITY, GIVEN THE DEPTH Z IN KM AND
C       THE TEMPERATURE T IN DEG C 22/12/05
C
      DATA CONST,GRAD/4.7201,-2.8E-4/
c
c       FIT FOR T>400, fit to depths of 50 and 75 km only
C
      DATA FREQ,ACT,ACTV,STRETCHY/-1.8E+13,409.E+3,
     C  10.0E-6,9.61E-3/
      DATA GASC/8.3145 /
      FACY=1.+STRETCHY*(Z-50.)/25.
      P=Z/30.
      P=P*1.E+9
      VSCL=CONST+GRAD*T+FREQ*EXP(-(ACT+ACTV*P)/(GASC*(T+273.)))
      VS=VSCL*FACY
      RETURN
C
      ENTRY TEMPV(Z,VS,T)
C
C       CALCULATES THE TEMPERATURE, GIVEN THE DEPTH Z IN KM
C       AND THE VELOCITY VS IN KM/SEC
C
      FACY=1.+STRETCHY*(Z-50.)/25.
      VSCL=VS/FACY
      P=Z/30.
      P=P*1.E+9
      T=1000.
      H=ACT+ACTV*P
C
      ITER=0
10    CONTINUE
      ITER=ITER+1
      TA=T+273.
      EX=FREQ*EXP(-H/(GASC*TA))
      VC=CONST+GRAD*T+EX
      DVCDT=GRAD+(H*EX)/(GASC*TA**2)
      DT=(VSCL-VC)/DVCDT
      T=T+DT
      IF(ABS(DT).GT.0.01) GO TO 10
c      WRITE(6,100) ITER
100   FORMAT(' ITER=',I3)
      RETURN
      END
C
C************************************************************
C
      SUBROUTINE INTERSC(DREF,TREF,NT,DDIAM,TDIAM,IDIAM,XI,YI)
C
C       FINDS WHERE DREF,TREF INTERSECTS DDIAM,TDIAM
C       RETURNS THE INTERSECTION AS XI,YI
C
      dimension dref(1000),tref(1000),ddiam(1000),tdiam(1000)
      i1=1
      i2=1
10    i1=i1+1
15    continue
      if(dref(i1).lt.ddiam(i2)) go to 10
      if(tref(i1).lt.tdiam(i2)) go to 20
      i2=i2+1
      go to 15
20    continue
c
c       intersection is between tref(i1-1), tref(i1) and
c       tdiam(i2-1) and tdiam(i2)
c
      a2=(tref(i1)-tref(i1-1))/(dref(i1)-dref(i1-1))
      a4=(tdiam(i2)-tdiam(i2-1))/(ddiam(i2)-ddiam(i2-1))
      xi=(a2*dref(i1-1)-a4*ddiam(i2-1)+tdiam(i2-1)-tref(i1-1))/
     c   (a2-a4)
      yi=a2*(xi-dref(i1-1))+tref(i1-1)
c
      i=1
      if(i.eq.2) then
        write(6,100) i1,dref(i1-1),dref(i1),
     c    tref(i1-1), tref(i1)
        write(16,100) i1,dref(i1-1),dref(i1),
     c  tref(i1-1), tref(i1)
100     format(' i1=',i3,' x1=',G15.8,' x2=',G15.8,' y1=',G15.8,
     c  ' y2=',g10.4)
        write(6,105) i2,ddiam(i2-1),ddiam(i2),tdiam(i2-1),
     c  tdiam(i2)
        write(16,105) i2,ddiam(i2-1),ddiam(i2),tdiam(i2-1),
     c  tdiam(i2)
105     format(' i2=',i3,' x3=',G15.8,' x4=',G15.8,' y3=',G15.8,
     c  ' y4=',G15.8)
        write(6,110) xi,yi
        write(16,110) xi,yi
110     format(' intersection xi=',G15.8,' yi=',G15.8)
      endif
c
      return
      end
