C      PROGRAM NISTL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This program is written to determine the aerosol dynamics        C
C     according to Discrete-Sectional Model considering two species in C
C     the system.  The number of discrete size is limited to MACDISC   C
C     for both species, and the total number of section is limited to  C
C     MAXSEC.  The actual number of discrete size is noted by NDISC,   C
C     and the actual number of section, NSEC, is determined by the     C
C     size increase factor, SIZEFACC, and the maximum size considered, C
C     XMAXSIZE.  In this system, discrete size particles of different  C
C     species do not interact (no coagulation), but particles in the   C
C     sections are not distinguishable.  This program is written       C
C     based on different number concentration, as noted by ETA.  The   C
C     system is based on the experiment performed in NIST, 1994.  The  C
C     subroutine XDOTEQ called by RK45 is used to determine the ODE by C
C     the Discrete-Sectional method.                                   C
C                                                                      C
C     The reactions are assumed to be                                  C
C                                                                      C
C     Fe(CO)5  ----->  Fe(g)      instantly                            C
C                                                                      C
C                   k1                                                 C
C     Fe(g) + O2  ----->  FeO(g) + O                                   C
C                                                                      C
C               k2                                                     C
C     FeO(g)  ----->  (FeO)2(s)                                        C
C                                                                      C
C                                                                      C
C     Si2O(CH3)6 -----> 2Si(g)     instantly                           C
C                                                                      C
C                  k1                                                  C
C     Si(g) + O2 -----> SiO2(g)                                        C
C                                                                      C
C                k2                                                    C
C     SiO2(g)  ----->  (SiO2)2(s)                                      C
C                                                                      C
C     Condensation occurs in the system and the Kelvin effect is       C
C     considered assuming no vaporization.  Mass balance is considered C
C     by removing all the colliding particles and then adding the      C
C     fraction that is newly formed.  The collision coefficient for    C
C     sections is determined by the integration over the section.      C
C     Written by Chang-Yu Wu in the University of Cincinnati, Aug,     C
C     1994 under the guidance of Dr. Pratim Biswas.                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (MAXDISC = 300, MAXSEC = 300)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XM(MAXDISC+MAXDISC+MAXSEC),DM(MAXDISC+MAXDISC+MAXSEC)
      REAL*8 VG(2,MAXDISC),VG2(2,MAXDISC)
      REAL*8 DLNVG(2,MAXDISC),DDCOEF(MAXDISC,MAXDISC)
      REAL*8 VGETA(2,MAXDISC),XMFAC(2,MAXDISC,MAXDISC)
      REAL*8 VL(MAXSEC),VGS(MAXSEC),VGSSQ(MAXSEC),VGSN(MAXSEC)
      REAL*8 DVLNV(MAXSEC),DVK(MAXSEC),DLNVL(MAXSEC),DLNVL2(MAXSEC)
      REAL*8 DV3(MAXSEC),DVK3(MAXSEC),DV2(MAXSEC),DVK2(MAXSEC)
      REAL*8 CONVV(MAXDISC),CONV3(MAXDISC),VGSQ(2,MAXDISC)
      REAL*8 DP(2,MAXDISC),DPS(MAXSEC),XNPD(2,MAXDISC),XNPDS(2,MAXSEC)
      REAL*8 XNPD0(2,MAXDISC),XNPDS0(2,MAXSEC)
      REAL*8 V1(2),XR1(2),RHO(2),XMW(2),XK1(2),XK2(2),CONCPRE(2)
      REAL*8 PSAT(2),SIGMA(2),CONCSP(2),R1(2),CONV2(MAXDISC),PSAT0(2)
      REAL*8 TIME(20),RHOSEC(MAXSEC),DRHO(MAXSEC),PXM(MAXSEC)
      REAL*8 A(2),B(2),A0(2),A1(2),EA(2),B1(2),TN1(2),TP1(2),A2(2),B2(2)
      REAL*8 XK2TMP(2),XK2RATIO(2),TN2(2),XMINIT,DSQRTSF
      REAL C13,S,S1,S2,S3,S4
      INTEGER M,NDISC(2),NSEC,NDISCT,NDISC1,KFSB(2),ETA,NCAL,NCAL1
      INTEGER MAERO,NEXIDIS(2),NEXISEC,NORDER
      INTEGER NFNUCL,NFCOND,NFCOAG
      INTEGER*4 NN(20),NSP1,NSP2
      INTEGER NCOND(2),NCOND1(2),INDB,NCONSCND
      CHARACTER*12 RESULT1,RESULT2
      COMMON /AERO1/ VL,VG,VGS,VGSN
      COMMON /AERO2/ VGSQ
      COMMON /AERO3/ SIZEFAC,VG2,KFSB,XR1
      COMMON /AERO4/ RHOSEC
      COMMON /CONSTANT/ AVO,PI,XKB,ETA
      COMMON /MONOMER/ V1,RHO
      COMMON /RXN/ XK1,XK2,CO2
      COMMON /SYS1/ TEMP,NCOND1
      COMMON /CALCON1/NDISC,NDISCT,NDISC1,NDISCA,M,NSEC,XNCAL,NCAL1,NCAL
      COMMON /PRES/ XNPD,XNPDS,CONCSP
      COMMON /STAST1/ DVLNV,DVK,DLNSF
      COMMON /STAST2/ DLNVG,DLNVL,DLNVL2,VGSSQ
      COMMON /B1/ CONVV,CONV3,C13,CONV2
      COMMON /B2/ R1,DDCOEF
      COMMON /NAMER/ RESULT2
      COMMON /EV/ VGETA,XMFAC,C2E
      COMMON /TMOD/ XTMOD,DRHO
      COMMON /A/ ALPHA
      COMMON /K2R/ XK2RATIO,MAERO
      COMMON /AEROFLAG/ NFNUCL,NFCOND,NFCOAG
      COMMON /TEST/ NCONSCND
      common /ttt/kkk
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     AVO       Avogadro's number (6.023D23)
C     BDD       Coagulation coefficient of discrete-discrete
C     BDS1      Collision coefficient of discrete with section to form 
C               a new section
C     BDS2      Collision coefficient of discrete with section to remove
C               the section
C     BDS4      Collision coefficient of discrete with section to remove
C               the discrete
C     BDS5      Collision coefficient of discrete with section to form
C               the section
C     BSS1      Collision coefficient of 2 smaller sections to form a 
C               new section
C     BSS2      Collision coefficient of the section with a smaller
C               section to remove the section
C     BSS3      Collision coefficient of 2 same sections to remove the
C               section
C     BSS4      Collision coefficient of the section with a larger 
C               section to remove the section
C     BSS5      Collision coefficient of the section with a smaller 
C               section to form the section
C     BSS6      Collision coefficient of 2 same sections to form
C               the section
C     CONCPRE   Initial precusor concentration (mol/cc)
C     CONCSP    Precursor concentration wrt time (mol/cc)
C     CONVV(I)  1.0/XNMON(I)
C     CONV3(I)  XNMON(I)**(1.0/3.0)
C     COUNT     System time variable (s)
C     CO2       Oxygen concentration (mol/cc)
C     C13       1.0/3.0
C     DIST      Axis distance (cm)
C     DLNSF 	Ln(SIZEFAC)
C     DLNVG	Ln(VG)
C     DM        dM/dt, differential equations
C     DP        Mean particle diameter (cm)
C     DT        Differential time (s)
C     DVK(I)    VL(I+1) - VL(I)
C     DVK2(I)   2.0*DVK(I)
C     DVK3(I)   3.0*DVK(I)
C     DVLNV(I)  VL(I+1)*ln(VL(I+1)) - VL(I)*ln(VL(I))
C     DV2       VL(I+1)**2 - vl(I)**2
C     DV3(I)    VL(I+1)**3 - VL(I)**3
C     ETA       Concentration index: 0-number, 1-volume, 2-colume square
C     FLVEL     Flow velocity wrt time (cm/s)
C     FLVEL0    Initial flow velocity (cm/s)
C     INDB	Index for beta calculation (0-First time full calculation; 
C               1-time modification)
C     KFSB      First section boundary for the first section formation
C               by discrete sizes interaction 
C     M         Total number of ODE's (NDISC + NSEC)
C     MAXDISC   Maximum number of discrete sizes
C     MAXSEC    Maximum number of sections
C     NCAL1     XNCAL - 1
C     NCOND     Index for complete condensation (0-No barrier; 1-Kelvin
C			effect)
C     NCOND1    Index for coagulation of the monomer (0-Monomer 
C			coagulation; 1-Vapor condensation)
C     NCONSCND  Index for Constant bulk vapor concentration for 
C			condensation (0-vapor consumption; 1-Constant 
C			concentration)
C     NDISC     Total number of discrete sizes
C     NDISC1    NDISC + 1
C     NDISC2    NDISC + 2
C     NEXIDIS(I)Number of discrete sizes of existing particles for 
C			Species I
C     NEXISEC   Number of sections of existing particles
C     NFNUCL    Flag for nucleation (0-no,1-by stable monomer 
C			coagulation, 2-given kinetic rate)
C     NFCOND    Flag for condensation (0-no,1-yes)
C     NFCOAG    Flag for coagulation (0-no,1-yes)
C     NN(I)     Calculation steps for stage I
C     NORDER    The number of the order of existing particles
C     NSEC      Total number of sections
C     PI        Circumference ratio (3.1415926)
C     PSAT      Saturation pressure of FeO (atm)
C     RESULT1   Concentration & M2 output filename
C     RESULT2   Aerosol size distribution data filename
C     RHO       Aerosol density (g/cm3)
C     SIGMA     Aerosol surface tension (dyne/cm)
C     SIZEFAC   Size increasing factor based on volume
C     SQRTSF    Square root of the size factor
C     TEMP      Temperature (K)
C     TIME(I)   Time for sytem stage I
C     VG        Discrete volume (cm3)
C     VGS       Mean volume of the section
C     VGSN      Mean volume of the section for n-based model (statistics)
C     VGSQ      VG**2
C     VG2       2.0*VG
C     VL        Lower bound aerosol volume of the section (cm3)
C     V1        Monomer volume (cm3)
C     XKB       Boltzman's Constant, 1.38D-16 erg/K
C     XK1       Aerosol vapor formation rate Constant (cm3/mol/s)
C     XK2       Dimer nucleation rate (cm3/#/s)
C     XM        Aerosol and vapor concentration
C       XM(1)                  Vapor molecule concentration of species 1 
C                              (#/cc)
C       XM(2...NDISC(1))       Discrete size aerosols of species 1 (#/cc)
C       XM(NDISC(1)+1)	       Vapor molecule concentration of species 2 
C                              (#/cc)
C       XM(NDISC(1)+2..NDISCT) Discrete size aerosols of species 1 (#/cc)
C       XM(NDISC1...M)         Sectional size aerosols (#/cc)
C     XMAXSIZE  Upper bound of the largest section (A)
C     XMINIT    Initial FeO concentration for Constant condensation
C     XMW       Molecular weight (g/mol)
C     XNCAL     Total dissection number for the coagulation coefficient
C     XNMON     Mean number of molecules in one aerosol in the section
C     XNMON2    2.0*XNMON
C     XNPD(I)   Surface saturation concentration for section I particle 
C               (#/cc)
C     XR1       Aerosol vapor formation rate (#/mol/s)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SYSOPF.DAT   Output file for all the system parameters and
C                  aerosol species properties.  
C     CONCF.DAT    Output file for vapor concentration and M2 wrt time.
C     OUTPTF.DAT   Output file for particle size distribution and M0, 
C                  M1, dpg, dpv, sigma. 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      S = SECNDS(0.0)
      S1 = SECNDS(S)
C***********************************************************************
C     Program Constant declaration
C***********************************************************************
      AVO = 6.023D23
      PI = 3.1415926
      XKB = 1.38D-16
      C13 = 1.0/3.0
C***********************************************************************
C     Reading system and calculation variable values
C***********************************************************************
      OPEN(1,FILE='SYSOPH.DAT')
      DO 31 II = 1,2
        READ(1,*)NDISC(II)
        WRITE(*,5)NDISC(II)
        READ(1,*)XMW(II)
        WRITE(*,8)XMW(II)
        READ(1,*)RHO(II)
        WRITE(*,9)RHO(II)
        READ(1,*)XK1(II)
        WRITE(*,11)XK1(II)
        READ(1,*)XK2(II)
        WRITE(*,12)XK2(II)
        READ(1,*)CONCPRE(II)
        WRITE(*,14)CONCPRE(II)
        READ(1,*)PSAT0(II)
        WRITE(*,18)PSAT0(II)
        READ(1,*)SIGMA(II)
        WRITE(*,19)SIGMA(II)
        READ(1,*)NEXIDIS(II)
        WRITE(*,23)NEXIDIS(II)
        DO 29 IK = 1,NEXIDIS(II)
          if (II == 1) then
            READ(1,*)NORDER,XM(0+NORDER)
          else
            READ(1,*)NORDER,XM(NDISC(II-1)+NORDER)
          end if
 29     CONTINUE
	READ(1,*)NCOND(II)
	READ(1,*)NCOND1(II)
 31   CONTINUE
      XMINIT = XM(1)
      NDISCT = NDISC(1)+ NDISC(2)
      NDISC1 = NDISCT + 1
      NDISCA = NDISC(1) + 1
      READ(1,*)NEXISEC
      DO 32 IK = 1,NEXISEC
        READ(1,*)NORDER,XM(NDISCT+NORDER)
 32   CONTINUE
      READ(1,*)SIZEFAC
      WRITE(*,6)SIZEFAC
      READ(1,*)XMAXSIZE
      WRITE(*,7)XMAXSIZE
      READ(1,*)TEMP0
      WRITE(*,10)TEMP0 
      READ(1,*)FLVEL0 
      READ(1,*)CO2
      WRITE(*,15)CO2
      READ(1,*)NTSEC
      DO 3 I = 1,NTSEC
	  READ(1,*)TIME(I)
	  WRITE(*,16)TIME(I)
	  READ(1,*)NN(I)
	  WRITE(*,17)NN(I)
  3   CONTINUE
      READ(1,*)NSP1
      READ(1,*)NSP2
      READ(1,*)XNCAL
      WRITE(*,21)XNCAL
      READ(1,*)NCAL
      READ(1,*)ETA
      WRITE(*,20)ETA
      READ(1,*)ALPHA
      WRITE(*,22)ALPHA
      READ(1,*)NFNUCL
      READ(1,*)NFCOND
      READ(1,*)NFCOAG
      READ(1,*)NCONSCND
  5   FORMAT(1X,"The total number of discrete size = ",I2)
  6   FORMAT(1X,"The size increase factor = ",F10.5)
  7   FORMAT(1X,"The upper bound of the largest section (A) = ",F10.5)
  8   FORMAT(1X,"The aerosol species molecular weight (g) = ",F10.5)
  9   FORMAT(1X,"The aerosol species density (g/cm3)",F10.5)
 10   FORMAT(1X,"Temperature (K) = ",F10.5)
 11   FORMAT(1X,"The aerosol vapor formation rate Constant, K1 = "
     &      ,1PE10.4)
 12   FORMAT(1X,"The dimer formation rate Constant, K2 = ",1PE10.4)
 14   FORMAT(1X,"Initial precursor concentration (mol/cc) = ",1PE10.4)
 15   FORMAT(1X,"Initial O2 concentration (mol/cc) = ",1PE10.4)
 16   FORMAT(1X,"Time for the stage (s) = ",1PE10.4)
 17   FORMAT(1X,"Calculation steps for the stage = ",I5)
 18   FORMAT(1X,"Saturation concentration (atm) = ",1PE10.4)
 19   FORMAT(1X,"Surface tension (dyne/cm) = ",F10.5)
 20   FORMAT(1X,"Concentration index = ",I2)
 21   FORMAT(1X,"Number of differentials are = ",F6.1)
 22   FORMAT(1X,"The sticking factor is ", F8.3)
 23   FORMAT(1X,"The # of discrete sizes of existing particles is",I3)
      WRITE(*,*)"CONCENTRATION & M2 OUTPUT FILENAME = ?"
      WRITE(*,*)"------1------label added for debugging"
      READ(1,*)RESULT1
      WRITE(*,*)"AEROSOL SIZE DISTRIBUTION DATA FILE = ?"
      WRITE(*,*)"------3------label added for debugging"
      READ(1,*)RESULT2
      CLOSE(1)
      CALL TEMPHIST(TEMP0,FLVEL0)
      CALL MTSAT(A,B)
      CALL RTDATA(A0,EA)
      CALL NUCLDATA(A1,B1,A2,B2,TN1,TN2,TP1)
C***********************************************************************
C     Set initial values and calculation Constants
C***********************************************************************
      TEMP = TEMP0
      NCAL1 = INT(XNCAL)
      DLNSF = DLOG(SIZEFAC)
      DSQRTSF = DSQRT(SIZEFAC)
      C2E = 2.0D0**(ETA-1)
      V1(1) = XMW(1)/RHO(1)/AVO
      V1(2) = XMW(2)/RHO(2)/AVO
      NSEC = INT(DLOG10(XMAXSIZE**3.0/DBLE(NDISC(1)))/DLOG10(SIZEFAC))
      M = NDISCT + NSEC
      MAERO = M + 2
      WRITE(*,*)"The total number of discrete sizes and sections is ",M
      DO 40 I = 2,M
	  DM(I) = 0.0
  40  CONTINUE
      XM(M+1) = CONCPRE(1)
      XM(MAERO) = CONCPRE(2)
C***********************************************************************
C     Set the discrete size and lower bound of the section, and their
C     related variables
C***********************************************************************
C     ###   Discrete sizes ###    
C**********************************
      DO 61 II = 1,2
       IF (NCOND(II).EQ.0) PSAT0(II) = 0.0D0
       DO 60 I = 1,NDISC(II)
	VG(II,I) = V1(II)*DBLE(I)
	VG2(II,I) = 2.0*VG(II,I)
	VGSQ(II,I) = VG(II,I)*VG(II,I)
	VGETA(II,I) = VG(II,I)**ETA
	DLNVG(II,I) = DLOG(VG(II,I))
	XMFAC(II,I,I) = C2E
	DO 55 J = 1,I-1
	  XMFAC(II,J,I) = (VG(II,J)+VG(II,I))**ETA/
     &                    (VGETA(II,I)+VGETA(II,J))
	  XMFAC(II,I,J) = XMFAC(II,J,I)
  55	CONTINUE
	DP(II,I) = (VG(II,I)*6.0/PI)**C13
	XNPD0(II,I) = PSAT0(II)*DEXP(4.0*SIGMA(II)*V1(II)/
     &               (XKB*TEMP*DP(II,I)))/(82.054*TEMP)*AVO*VGETA(II,1)
      IF (NCOND(II).EQ.0) XNPD0(II,I) = 0.0
        XNPD(II,I) = XNPD0(II,I)
  60   CONTINUE
       R1(II) = 0.5*DP(II,1)
  61  CONTINUE
C***************************************
C     ###   The first section   ###
C***************************************
      VL(1) = VG(1,NDISC(1)) + 0.5*V1(1)
      DLNVL(1) = DLOG(VL(1))
C****************************************
C     ###   All other sections   ###
C****************************************
      DO 70 I = 2,NSEC
        IND = I - 1
	VL(I) = VL(IND)*SIZEFAC
	VGS(IND) = VL(IND)*DSQRTSF
 	VGSN(IND) = (VL(IND)+VL(I))/2.0d0
	VGSSQ(IND) = VGS(IND)*VGS(IND)
	DLNVL(I) = DLOG(VL(I))
	DPS(IND) = (VGS(IND)*6.0/PI)**C13
	XNPDS0(1,IND) = PSAT0(1)*DEXP(4.0*SIGMA(1)*V1(1)/
     &               (XKB*TEMP*DPS(IND)))/(82.054*TEMP)*AVO*VGETA(1,1)
      IF (NCOND(1).EQ.0) XNPDS0(1,IND) = 0.0
        XNPDS(1,IND) = XNPDS0(1,IND)
	XNPDS0(2,IND) = PSAT0(2)*DEXP(4.0*SIGMA(2)*V1(2)/
     &               (XKB*TEMP*DPS(IND)))/(82.054*TEMP)*AVO*VGETA(2,1)
      IF (NCOND(2).EQ.0) XNPDS0(2,IND) = 0.0
        XNPDS(2,IND) = XNPDS0(2,IND)
  70  CONTINUE
      I = NSEC + 1
      VL(I) = VL(NSEC)*SIZEFAC
      VGS(NSEC) = VL(NSEC)*DSQRTSF
      VGSN(NSEC) = (VL(NSEC)+VL(I))/2.0d0
      VGSSQ(NSEC) = VGS(NSEC)*VGS(NSEC)
      DLNVL(I) = DLOG(VL(I))
C***********************************************************************
C     ###   Determine the upper bound of the discrete sizes that can 
C           form the first section   
C           Define the reaction rate Constant 
C***********************************************************************
      DO 75 II = 1,2
        XR1(II) = XK1(II)*AVO*VGETA(II,1)  !*CO2!for debug, by ymm *VGETA(II,1) should be VG
        IF ((VG(II,NDISC(II))*2.0D0).GE.VL(2)) THEN
	  KFSB(II) = AINT(VL(2)/V1(II)) !modified by ymm. initial is KFSB(II) = AINT(VL(2)/V1(I))
        ELSE
	  KFSB(II) = NDISC(II)*2
        END IF
  75  CONTINUE
C***********************************************************************
C     ###   Define the statistical variables for the sections   ###
C***********************************************************************
      DO 68 I = 1,NSEC
	DVK(I) = VL(I+1) - VL(I)
  68  CONTINUE
      IF (ETA.EQ.0) THEN
	DO 71 I = 1,NSEC
	  IND = I + 1
	  DVLNV(I) = VL(IND)*DLNVL(IND) - VL(I)*DLNVL(I)
	  DVK3(I) = 3.0*DVK(I)
	  DV3(I) = VL(IND)**3 - VL(I)**3
  71    CONTINUE
      ELSEIF (ETA.EQ.1) THEN
	DO 72 I = 1,NSEC
	  IND = I + 1
	  DVK2(I) = 2.0*DVK(I)
	  DV2(I) = VL(IND)*VL(IND) - VL(I)*VL(I)
	  DLNVL2(I) = (DLNVL(IND)*DLNVL(IND) - DLNVL(I)*DLNVL(I))/2.0
  72    CONTINUE
      END IF
C***********************************************************************
C     ###   Define the beta calculation variables for the discrete size
C***********************************************************************
      LG = MAX(NDISC(1),NDISC(2))
      DO 74 I = 1,LG
	CONVV(I) = 1.0/DBLE(I)
	CONV3(I) = DBLE(I)**C13
	CONV2(I) = CONV3(I)*CONV3(I)
	DO 79 J = 1,I-1 !modified by ymm, initial version is J = 2,I-1
	  DDCOEF(J,I) = DSQRT(CONVV(I)+CONVV(J))*(CONV3(I)+CONV3(J))**2
  79    CONTINUE
  74  CONTINUE
C***********************************************************************
C     Set up the parameters for the differential solver 
C***********************************************************************
      S2 = SECNDS(S)
      CALL CALBETA(SIZEFAC)
      S3 = SECNDS(S)
      INDB = 1
      COUNT = 0.0
      DIST = 0.0
      OPEN(5,FILE=RESULT1)
      CALL OT(COUNT,XM)
      DO 102 KK = 1,NTSEC
      WRITE(*,*) 'Time section:',KK
      DT = TIME(KK)/DBLE(NN(KK))
      CALL INIT0(DT)
C****************************************************
C     ###   Differential equations solving   ###
C****************************************************
      DO 100 I = 1,NN(KK)
	  COUNT = COUNT + DT
C	  FLVEL = FLVEL0*TEMP/TEMP0
	  DIST = DIST + FLVEL*DT
	  TEMPP = TEMP
	  XK2TMP(1) = XK2(1)
	  XK2TMP(2) = XK2(2)
	  CALL TEMPXT(DIST,TEMP,FLVEL)
	  CALL VAPORSAT(TEMP,PSAT,A,B)
	  CALL RXNRT(TEMP,XK1,A0,EA)
	  CALL NUCLRT(TEMP,XK2,A1,B1,A2,B2,TN1,TN2,TP1)
	  TRATIO = TEMP/TEMPP
          CO2 = CO2/TRATIO
	  XTMOD = SQRT(TEMP/TEMP0)
	  kkk = i
	  DO 86 II = 1,2
	    XK2RATIO(II) = XK2TMP(II)/XK2(II)
	    XR1(II) = XK1(II)*AVO*VGETA(II,1) !*CO2
	    IF (NCOND(II).NE.0) THEN 
	     PRATIO = PSAT(II)/PSAT0(II)
	     DO 108 IJ = 1,NDISC(II)
	      XNPD(II,IJ) = XNPD0(II,IJ)*PRATIO
 108         CONTINUE
             DO 109 IJ = 1,NSEC
              XNPDS(II,IJ) = XNPDS0(II,IJ)*PRATIO
 109         CONTINUE
            END IF
  86      CONTINUE
          DO 120 IJ = 1,MAERO
            XM(IJ) = XM(IJ)/TRATIO !the expand due to temperature increase
 120      CONTINUE
          IF (NCONSCND.EQ.1) XM(1) = XMINIT
	  DO 87 J = 1,NSEC
	     PXM(J) = XM(J+NDISCT)
  87      CONTINUE
	  CALL RK45(XM,DM,COUNT,DT,MAERO)
	  IF (MOD(I,NSP1).EQ.0) CALL OT(COUNT,XM)
	  IF (MOD(I,NSP2).EQ.0) THEN
	    XMM2 = 0.0D0
	    IF (ETA.EQ.0) THEN
	      DO 88 J = 2,NDISC(1)
	        XMM2 = XMM2 + XM(J)*VGSQ(1,J)
  88            CONTINUE
	      DO 89 J = 2,NDISC(2)
		XMM2 = XMM2 + XM(NDISC(1)+J)*VGSQ(2,J)
  89          CONTINUE
	      DO 90 J = 1,NSEC
		XMM2 = XMM2 + XM(J+NDISCT)*DV3(J)/DVK3(J)
  90          CONTINUE
	    ELSEIF (ETA.EQ.1) THEN
	      DO 92 J = 2,NDISC(1)
		XMM2 = XMM2 + XM(J)*VG(1,J)
  92          CONTINUE
	      DO 93 J = 2,NDISC(2)
		XMM2 = XMM2 + XM(NDISC(1)+J)*VG(2,J)
  93          CONTINUE
              DO 94 J = 1,NSEC
		XMM2 = XMM2 + XM(J+NDISCT)*DV2(J)/DVK2(J)
  94          CONTINUE
	    ELSE
	      DO 95 J = 2,NDISC(1)
		XMM2 = XMM2 + XM(J)
  95          CONTINUE
	      DO 96 J = NDISC(1)+2,M
		XMM2 = XMM2 + XM(J)
  96          CONTINUE
	    END IF
	    WRITE(5,105)COUNT,DIST,XM(1)/VGETA(1,1),
     &                  XM(NDISCA)/VGETA(2,1),XMM2,XM(M+1)
	  END IF
	  DO 98 J = 1,NSEC
	    RHOSEC(J) = (PXM(J)*RHOSEC(J)+(XM(J)-PXM(J))*DRHO(J))/XM(J)
  98      CONTINUE
      write(*,*)i,'of',NN(KK)
 100  CONTINUE
 102  CONTINUE
 105  FORMAT(1X,1P6E12.4)
      S4 = SECNDS(S)
      WRITE(5,123)S1,S2,S3,S4
 123  FORMAT(1X,1P4E12.3)
      CLOSE(5)
C
 106  STOP
      END
C
C
C
      SUBROUTINE CALBETA(SIZEFAC)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine is written to determine the coagulation and
C     condensation coefficients in the free molecular regime.  The
C     nucleation rate coefficient (k2) following the reaction form is 
C     also determined.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (MAXDISC = 300, MAXSEC = 300)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VG(2,MAXDISC),VGSQ(2,MAXDISC)
      REAL*8 VL(MAXSEC),DVK(MAXSEC),DVLNV(MAXSEC)
      REAL*8 VGS(MAXSEC),VGSN(MAXSEC)
      REAL*8 VGETA(2,MAXDISC),XMFAC(2,MAXDISC,MAXDISC)
      REAL*8 CONVV(MAXDISC),CONV3(MAXDISC),DDCOEF(MAXDISC,MAXDISC)
      REAL*8 BDD(2,MAXDISC,MAXDISC),CONV2(MAXDISC)
      REAL*8 BDS1(2,MAXDISC,MAXSEC,MAXSEC),BDS25(2,MAXDISC,MAXSEC)
      REAL*8 BDS4(2,MAXDISC,MAXSEC)
      REAL*8 BSS1(MAXSEC,MAXSEC,MAXSEC),BSS25(MAXSEC,MAXSEC)
      REAL*8 BSS36(MAXSEC),BSS4(MAXSEC,MAXSEC)
      REAL*8 R1(2),CONST0(2),V1(2),RHO(2),XK1(2),XK2(2),RHOSEC(MAXSEC)
      REAL C8,C6,C13
      INTEGER ETA,NDISC(2),NDISCT,NDISC1,M,NSEC,NCAL,NCAL1,K,IDD
      INTEGER NSTART(2),NCOND1(2),NFNUCL,NFCOND,NFCOAG
      COMMON /AERO1/ VL,VG,VGS,VGSN
      COMMON /AERO4/ RHOSEC
      COMMON /B/ BDD,BDS1,BDS25,BDS4,BSS1,BSS25,BSS36,BSS4
      COMMON /CONSTANT/ AVO,PI,XKB,ETA
      COMMON /MONOMER/ V1,RHO
      COMMON /SYS1/ TEMP,NCOND1
      COMMON /RXN/ XK1,XK2,CO2
      COMMON /B1/ CONVV,CONV3,C13,CONV2
      COMMON /B2/ R1,DDCOEF
      COMMON /STAST1/ DVLNV,DVK,DLNSF
      COMMON /CALCON1/NDISC,NDISCT,NDISC1,NDISCA,M,NSEC,XNCAL,NCAL1,NCAL
      COMMON /FFF/ I,J,K,IDD,II
      COMMON /EV/ VGETA,XMFAC,C2E
      COMMON /AERO2/ VGSQ
      COMMON /A/ ALPHA
      COMMON /AEROFLAG/ NFNUCL,NFCOND,NFCOAG
      EXTERNAL GLEGQ,FX,F,G,H !,DTWODQ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     BETA(I,J) Coagulation coefficients between particle I and J 
C               (cm3/#/s)
C     RHO       Aerosol density (g/cm3)
C     R1        Monomer radius (cm)
C     TEMP      Temperature (K)
C     V1        Monomer volume (cm3)
C     XKB       Boltzman's Constant, 1.38D-16 erg/K
C     XK2       Dimer nucleation rate (cm3/#/s)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      TMP = 6.0*XKB*TEMP
      DO 6 II = 1,2
        CONST0(II) = DSQRT(TMP*R1(II)/RHO(II))
        IF (NCOND1(II).EQ.0) THEN
          NSTART(II) = 1
        ELSE
          NSTART(II) = 2
        END IF
   6  CONTINUE
      EABS = 1.D-6
      ERE = 0.0001
      IR = 3
      C8 = DSQRT(32.0D0)
      C6 = 1.0/6.0
C*********************************
C     ###   Condensation   ###
C*********************************
      IF (NFCOND.EQ.0) GO TO 5
      write(*,*)'Condensation'
      DO 620 II = 1,2
       IF(NCOND1(II).EQ.0) GOTO 620
       DO 600 I = 2,NDISC(II)
        BT = CONST0(II)*CONV2(I)
	BDD(II,1,I) = BT/VGETA(II,I)
	BDD(II,I,1) = BT/VGETA(II,1)
 600   CONTINUE
       DO 610 I = 1,NSEC
        IND = I + 1
        BDS2 = 0.0
	BTEMP = ALPHA*CONST0(II)/DVK(I)
	DDV = DVK(I)/XNCAL
	XL = VL(I)
	XR = XL + DDV
	DO 603 LL = 1,NCAL1
	 IDD = 6
	 BDS2 = BDS2 + QUAD(GLEGQ,NCAL,FX,XL,XR)
	 IDD = 7
	 BDS4(II,1,I) = BDS4(II,1,I) + QUAD(GLEGQ,NCAL,FX,XL,XR)
	 XL = XR
	 XR = XR +DDV
 603    CONTINUE
        BDS4(II,1,I) = BDS4(II,1,I)*BTEMP
        DV = DVK(I) - V1(II)
        BDS5 = 0.0
        IDD = 9
        DDV1 = DV/XNCAL
        XL = VL(I)
        XR = XL + DDV1
        DO 604 LL = 1,NCAL1
         BDS5 = BDS5 + QUAD(GLEGQ,NCAL,FX,XL,XR)
         XL = XR
         XR = XR + DDV1
 604    CONTINUE
        BDS25(II,1,I) = (BDS2/VGETA(II,1) - BDS5)*BTEMP
	DDV = V1(II)/XNCAL
	XL = VL(IND) - V1(II)
	XR = XL + DDV
	DO 605 LL = 1,NCAL1
          BDS1(II,1,I,IND) = BDS1(II,1,I,IND)+QUAD(GLEGQ,NCAL,FX,XL,XR) !for debug, by ymm, not IND
          XL = XR
          XR = XR + DDV
 605    CONTINUE
        BDS1(II,1,I,IND) = BDS1(II,1,I,IND)*BTEMP
 610   CONTINUE
 620  CONTINUE
C*******************************
C     ###   Coagulation   ###
C*******************************
C====================================
C     Discrete-Discrete interaction
C=====
   5  DO 95 II = 1,2
       IF (NFCOAG.EQ.0) GO TO 614
       write(*,*)'BDD'
       DO 20 I = NSTART(II),NDISC(II)
	DO 10 J = I+1,NDISC(II)
	  BDD(II,I,J) = DDCOEF(I,J)*CONST0(II)/VGETA(II,J)
	  BDD(II,J,I) = DDCOEF(I,J)*CONST0(II)/VGETA(II,I)
  10    CONTINUE
	BDD(II,I,I) = CONST0(II)*C8*DBLE(I)**C6/VGETA(II,I)
  20   CONTINUE
C*******************************
C     ###   Nucleation   ###
C*******************************
 614  IF (NFNUCL.EQ.0) GO TO 21
      IF (NFNUCL.EQ.2) GO TO 615
      BT = CONST0(II)*DSQRT(CONVV(1)+CONVV(1))*(CONV3(1)+CONV3(1))**2 !CONV2(1)
      BDD(II,1,1) = BT/VGETA(II,1)
      GO TO 21
 615  BDD(II,1,1) = XK2(II)/AVO/VGETA(II,1)
C=================================================
C    Discrete-Section : formation of new section
C====
  21  IF (NFCOAG.EQ.0) GO TO 95
      write(*,*)'BDS1'
       IDD = 1
       DO 50 I = 2,NSEC
	DO 45 J = 1,I-1
	  IF ((VL(I)-VL(J+1)).GT.VG(II,NDISC(II))) GOTO 45 !the differance is greater than the biggest discrete
	  BTEMP = CONST0(II)/DVK(J)
	  DO 40 K = NSTART(II),NDISC(II)
	    XLB = VG(II,K) + VL(J)
	    XUB = VG(II,K) + VL(J+1)
	    IF ((XLB.GT.VL(I+1)).OR.(XUB.LT.VL(I))) GOTO 40
	    IF (XLB.LE.VL(I)) THEN  ! Here asuming that dvk(j+1)>dvk(j), namely, sizefc>1
	      DV = XUB - VL(I) 
	      DDV = DV/XNCAL
	      XL = VL(I) - VG(II,K)
	      XR = XL + DDV
	      DO 31 LL = 1,NCAL1
                BDS1(II,K,J,I)=BDS1(II,K,J,I)+QUAD(GLEGQ,NCAL,FX,XL,XR)
                XL = XR
                XR = XR + DDV
  31          CONTINUE
            ELSE
	     IF (XUB.LE.VL(I+1)) THEN
	      DDV = DVK(J)/XNCAL
	      XL = VL(J)
	      XR = XL + DDV
	      DO 33 LL = 1,NCAL1
                BDS1(II,K,J,I)=BDS1(II,K,J,I)+QUAD(GLEGQ,NCAL,FX,XL,XR)
                XL = XR
                XR = XR + DDV
  33          CONTINUE
             ELSE
              DV = VL(I+1) - XLB
	      DDV = DV/XNCAL
	      XL = VL(J)
	      XR = XL + DDV
	      DO 35 LL = 1,NCAL1
                BDS1(II,K,J,I)=BDS1(II,K,J,I)+QUAD(GLEGQ,NCAL,FX,XL,XR)
                XL = XR
                XR = XR + DDV
  35          CONTINUE
             END IF
            END IF
	    BDS1(II,K,J,I) = BDS1(II,K,J,I)*BTEMP
  40      CONTINUE
  45    CONTINUE
  50   CONTINUE
C=====================================================================
C    Discrete-Section : removal & formation of the discrete size or 
C    the section
C=====
      write(*,*)'BDS2,BDS4,BDS5'
       DO 100 I = 1,NSEC
	BTEMP = CONST0(II)/DVK(I)
	DDV = DVK(I)/XNCAL
        DO 90 K = NSTART(II),NDISC(II)
	  XL = VL(I)
	  XR = XL + DDV
	  BDS2 = 0.0d0 !for debug, by ymm, should be deleted, otherwise the condensation will be canceled
	  DO 82 LL = 1,NCAL1
	    IDD = 2
	    BDS2 = BDS2 + QUAD(GLEGQ,NCAL,FX,XL,XR)
            IDD = 4
	    BDS4(II,K,I) = BDS4(II,K,I) + QUAD(GLEGQ,NCAL,FX,XL,XR)
	    XL = XR
	    XR = XR + DDV
  82      CONTINUE
          BDS4(II,K,I) = BDS4(II,K,I)*BTEMP
          BDS5 = 0.0D0
	  IF (VG(II,K).LE.DVK(I)) THEN 
	    DV = DVK(I) - VG(II,K)
	    IDD = 1
	    DDV1 = DV/XNCAL
	    XL = VL(I)
	    XR = XL + DDV1
	    DO 84 LL = 1,NCAL1
	      BDS5 = BDS5 + QUAD(GLEGQ,NCAL,FX,XL,XR)
	      XL = XR
	      XR = XR + DDV1
  84        CONTINUE
          END IF
          BDS25(II,K,I) = (BDS2/VGETA(II,K) - BDS5)*BTEMP
  90    CONTINUE
 100   CONTINUE
  95  CONTINUE
C=======================================================================
C     2 Same Sections : removal & formation 
C=====
      IF (NFCOAG.EQ.0) GO TO 621
      write(*,*)'BSS3,BSS6'
      DO 180 J = 1,NSEC
	BTEMP = CONST0(1)/DVK(J)/DVK(J)
        IDD = 13
	CALL DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES3,ERRE)
        IF (SIZEFAC.LE.2.0) THEN
          RES6 = 0.0
        ELSE
          IDD = 16
          CALL DTWODQ(F,VL(J),DVK(J),G,H,EABS,ERE,IR,RES6,ERRE)
        END IF
        BSS36(J) = BTEMP*(RES3-RES6)*0.5 !for debug, by ymm, RES3+RES6? after collision, two become one
 180  CONTINUE  
C=======================================================================
C     Section-Section: Removal and formation of the section by collision
C     with a larger or a smaller section 
C=====
      write(*,*)'BSS2,BSS4,BSS5'
      DO 250 I = 2,NSEC
	DO 245 J = 1,I-1
	  BTEMP = CONST0(1)/DVK(I)/DVK(J)
	  IDD = 12
	  CALL DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES2,ERRE)
          IDD = 14
	  CALL DTWODQ(F,VL(I),VL(I+1),G,H,EABS,ERE,IR,RES,ERRE)
	  BSS4(J,I) = RES*BTEMP
	  IF (DVK(I).LE.VL(J)) THEN
	    RES5 = 0.0
	  ELSE
	    IF (DVK(I).LE.VL(J+1)) THEN
	      IDD = 15
	      CALL DTWODQ(F,VL(J),DVK(I),G,H,EABS,ERE,IR,RES5,ERRE)
	    ELSE
	      IDD = 17
	      CALL DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES5,ERRE)
	    END IF
	  END IF
	  BSS25(I,J) = BTEMP*(RES2 - RES5)
 245    CONTINUE
 250  CONTINUE
C=======================================================================
C     Section-Section: Formation of new section by collisions of smaller
C     sections
C     K: Larger section
C=====
      write(*,*)'BSS1'
      DO 500 I = 2,NSEC
	DO 495 J = 1,I-1
	  DO 490 K = J,I-1
	    IF ((VL(I)-VL(K+1)).GT.VL(K+1)) GOTO 490
	    XVL12 = VL(K)+VL(J)
	    XVU12 = VL(K+1)+VL(J+1)
	    IF ((VL(I+1).LT.XVL12).OR.(VL(I).GT.XVU12)) GOTO 490
	    BTEMP = CONST0(1)/DVK(J)/DVK(K)
	    IF (J.EQ.K) THEN
	      IF (VL(I).LE.XVL12) THEN
	        IF (VL(I+1).LE.(VL(J)+VL(K+1))) THEN
	          IDD = 21
	          CALL DTWODQ(F,VL(K),VL(I+1)-VL(J),G,H,EABS,ERE,IR,
     &                 RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                ELSE
                  IDD = 22
                  CALL DTWODQ(F,VL(K),VL(K+1),G,H,EABS,ERE,IR,RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                  IDD = 23
                  CALL DTWODQ(F,VL(I+1)-VL(K+1),VL(J+1),G,H,EABS,ERE,
     &                 IR,RES,ERRE)
                  BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                END IF
              ELSE
                IF (VL(I).LE.(VL(J)+VL(K+1))) THEN
                  IDD = 24
                  CALL DTWODQ(F,VL(K),VL(K+1),G,H,EABS,ERE,IR,RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                  IDD = 25
                  CALL DTWODQ(F,VL(J),VL(I)-VL(K),G,H,EABS,ERE,IR,
     &                 RES,ERRE)
                  BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                ELSE
                  IDD = 26
                  CALL DTWODQ(F,VL(I)-VL(J+1),VL(K+1),G,H,EABS,ERE,IR,
     &                 RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                END IF
              END IF
              BSS1(I,J,K) = BSS1(I,J,K)/2.0
            ELSE
              IF (VL(I).LE.XVL12) THEN
                IF (VL(I+1).LE.(VL(J+1)+VL(K))) THEN
                  IDD = 31 
                  CALL DTWODQ(F,VL(J),VL(I+1)-VL(K),G,H,EABS,ERE,IR,
     &                 RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                ELSEIF (VL(I+1).LE.(VL(J)+VL(K+1))) THEN
                  IDD = 32
                  CALL DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                ELSE
                  IDD = 33
                  CALL DTWODQ(F,VL(J),VL(I+1)-VL(K+1),G,H,EABS,ERE,IR,
     &                 RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                  IDD = 34
                  CALL DTWODQ(F,VL(I+1)-VL(K+1),VL(J+1),G,H,EABS,ERE,
     &                 IR,RES,ERRE)
                  BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                END IF
              ELSE
                IF (VL(I).LE.(VL(J+1)+VL(K))) THEN
                  IDD = 35
                  CALL DTWODQ(F,VL(I)-VL(K),VL(J+1),G,H,EABS,ERE,IR
     &                 ,RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                  IDD = 36
                  CALL DTWODQ(F,VL(J),VL(I)-VL(K),G,H,EABS,ERE,IR,
     &                 RES,ERRE)
                  BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                ELSEIF (VL(I).LE.(VL(J)+VL(K+1))) THEN
                  IDD = 37
                  CALL DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                ELSE
                  IDD = 38
                  CALL DTWODQ(F,VL(I)-VL(K+1),VL(J+1),G,H,EABS,ERE,IR,
     &                 RES,ERRE)
                  BSS1(I,J,K) = BTEMP*RES
                END IF
              END IF
            END IF
 490      CONTINUE
 495    CONTINUE
 500  CONTINUE
C
 621  continue
      RETURN
      END
C
C
C
      SUBROUTINE OT(COUNT,XM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This program is written to write the results to the output file  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (MAXDISC = 300, MAXSEC = 300)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,NDISC(2),NDISCT,NDISC1,ETA,NSEC,NCAL1,NCAL
      REAL*8 XM(MAXDISC+MAXDISC+MAXSEC)
      REAL*8 VG(2,MAXDISC),DLNVG(2,MAXDISC),VGSQ(2,MAXDISC)
      REAL*8 VL(MAXSEC),VGS(MAXSEC),VGSSQ(MAXSEC),VGSN(MAXSEC)
      REAL*8 DLNVL(MAXSEC),DLNVL2(MAXSEC)
      REAL*8 VGETA(2,MAXDISC),XMFAC(2,MAXDISC,MAXDISC)
      REAL*8 DVLNV(MAXSEC),DVK(MAXSEC)
      REAL*8 XNPD(2,MAXDISC),XNPDS(2,MAXSEC)
      REAL*8 XN(MAXDISC+MAXDISC+MAXSEC),CONCSP(2)
      REAL*8 XNK(MAXDISC+MAXDISC+MAXSEC)
      CHARACTER*12 RESULT2
      COMMON /CALCON1/NDISC,NDISCT,NDISC1,NDISCA,M,NSEC,XNCAL,NCAL1,NCAL
      COMMON /AERO1/ VL,VG,VGS,VGSN
      COMMON /PRES/ XNPD,XNPDS,CONCSP
      COMMON /CONSTANT/ AVO,PI,XKB,ETA
      COMMON /STAST1/ DVLNV,DVK,DLNSF
      COMMON /STAST2/ DLNVG,DLNVL,DLNVL2,VGSSQ
      COMMON /NAMER/ RESULT2
      COMMON /EV/ VGETA,XMFAC,C2E
      COMMON /AERO2/ VGSQ
      REAL*8 TTVOL2,XNN2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DVK(I)    VL(I+1) - VL(I)
C     DVLNV(I)  VL(I+1)*ln(VL(I+1)) - VL(I)*ln(VL(I))
C     TTVOL     Total volume of species (mole/cc)
C     TVG       Geometric mean volume (cm3)
C     XLNVG     Mean logarithm volume (cm3)
C     XLN2SIG   Mean logarithm standard deviation
C     XM        Aerosol concentration (#/cc)
C     XMM1      Total aerosol volume concentration (cm3/cc)
C     XN        Aerosol size distribution function (#/cc/cm)
C     XNN       Total number concentration (#/cc)     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (ETA.EQ.0) THEN
        XLNVG = 0.0
        XMM1 = 0.0
        XNN = 0.0
        DO 10 I = 2,NDISC(1)
	  XN(I) = XM(I)/VG(1,1)
	  XMM1 = XMM1 + VG(1,I)*XM(I)
	  XNN = XNN + XM(I)
	  XLNVG = XLNVG + DLNVG(1,I)*XM(I)
  10    CONTINUE
        DO 11 I = 2,NDISC(2)
          IND = I+NDISC(1)
	  XN(IND) = XM(IND)/VG(2,1)
	  XMM1 = XMM1 + VG(2,I)*XM(IND)
	  XNN = XNN + XM(IND)
	  XLNVG = XLNVG + DLNVG(2,I)*XM(IND)
  11    CONTINUE
        DO 12 I = 1,NSEC
          IND = I+NDISCT
	  XN(IND) = XM(IND)/DVK(I)
	  XMM1 = XMM1 + VGSN(I)*XM(IND)
	  XNN = XNN + XM(IND)
	  XLNVG = XLNVG + XM(IND)*(DVLNV(I)/DVK(I)-1.0)
  12    CONTINUE
C
        XLNVG = XLNVG/XNN         
        TVG = DEXP(XLNVG)
        VM = XMM1/XNN
        TTVOL = XMM1 + XM(1)*VG(1,1) + XM(1+NDISC(1))*VG(2,1) + 
     &        (CONCSP(1)*VG(1,1)+CONCSP(2)*VG(2,1))*AVO
        TTVOL2 = XMM1 + XM(1)*VG(1,1) + XM(1+NDISC(1))*VG(2,1)
        XNN2 = XNN + XM(1) + XM(1+NDISC(1))
        XN(1) = XM(1)/VG(1,1)
C
        XLN2SIG = 0.0
        DO 80 I = 2,NDISC(1)
	  XLN2SIG = XLN2SIG + (DLNVG(1,I)-XLNVG)**2*XM(I)
  80    CONTINUE
        DO 81 I = 2,NDISC(2)
	  XLN2SIG = XLN2SIG + (DLNVG(2,I)-XLNVG)**2*XM(I+NDISC(1))
  81    CONTINUE
        DO 82 I = 1,NSEC
	  XLN2SIG = XLN2SIG + XM(I+NDISCT)*((VL(I+1)*(DLNVL(I+1)-XLNVG)
     &          **2 - VL(I)*(DLNVL(I)-XLNVG)**2 - 2.0*DVLNV(I))/DVK(I)+
     &          2.0D0*(1.0D0+XLNVG))
  82    CONTINUE
        XLN2SIG = XLN2SIG/XNN/9.0
        TXLN2SIG = DEXP(DSQRT(XLN2SIG))
C
      ELSEIF (ETA.EQ.1) THEN
        XNN = 0.0
        XLNVG = 0.0
        XMM1 = 0.0
        DO 210 I = 2,NDISC(1)
          XNK(I) = XM(I)/VG(1,I)
	  XN(I) = XNK(I)/VG(1,1)
	  XNN = XNN + XNK(I)
	  XMM1 = XMM1 + XM(I)
	  XLNVG = XLNVG + DLNVG(1,I)*XNK(I)
 210    CONTINUE
        DO 211 I = 2,NDISC(2)
          IND = I+NDISC(1)
          XNK(IND) = XM(IND)/VG(2,I)
	  XN(IND) = XNK(IND)/VG(2,1)
	  XNN = XNN + XNK(IND)
	  XMM1 = XMM1 + XM(IND)
	  XLNVG = XLNVG + DLNVG(2,I)*XNK(IND)
 211    CONTINUE
        DO 212 I = 1,NSEC
          IND = I+NDISCT
          XNK(IND) = XM(IND)/DVK(I)
	  XN(IND) = XNK(IND)/VGS(I)
	  XNN = XNN + XNK(IND)*DLNSF
	  XMM1 = XMM1 + XM(IND)
	  XLNVG = XLNVG + XNK(IND)*DLNVL2(I)
 212    CONTINUE
C
        XLNVG = XLNVG/XNN         
        TVG = DEXP(XLNVG)
        VM = XMM1/XNN
        TTVOL = XMM1 + XM(1) + XM(1+NDISC(1)) + 
     &        (CONCSP(1)*VG(1,1) + CONCSP(2)*VG(2,1))*AVO
        TTVOL2 = XMM1 + XM(1) + XM(1+NDISC(1))
        XNN2 = XNN + XM(1)/VG(1,1) + XM(1+NDISC(1))/VG(2,1)
        XNK(1) = XM(1)/VG(1,1)
	  XN(1) = XNK(1)/VG(1,1)
C
        XLN2SIG = 0.0
        DO 280 I = 2,NDISC(1)
	  XLN2SIG = XLN2SIG + (DLNVG(1,I)-XLNVG)**2*XNK(I)
 280    CONTINUE
        DO 281 I = 2,NDISC(2)
	XLN2SIG=XLN2SIG+(DLNVG(2,I)-XLNVG)**2*XNK(I+NDISC(1))
 281    CONTINUE
        DO 282 I = 1,NSEC
	  XLN2SIG = XLN2SIG + XM(I+NDISCT)/DVK(I)*((DLNVL(I+1)-XLNVG)**3
     &            - (DLNVL(I)-XLNVG)**3)/3.0D0
 282    CONTINUE
        XLN2SIG = XLN2SIG/XNN/9.0
        TXLN2SIG = DEXP(DSQRT(XLN2SIG))
C
      ELSE
        XNN = 0.0
        XLNVG = 0.0
        XMM1 = 0.0
        DO 310 I = 2,NDISC(1)
          XNK(I) = XM(I)/VGSQ(1,I)
	  XN(I) = XNK(I)/VG(1,1)
	  XNN = XNN + XNK(I)
	  XMM1 = XMM1 + XM(I)/VG(1,I)
	  XLNVG = XLNVG + DLNVG(1,I)*XNK(I)
 310    CONTINUE
        DO 311 I = 2,NDISC(2)
          IND = I+NDISC(1)
          XNK(IND) = XM(IND)/VGSQ(2,I)
	  XN(IND) = XNK(IND)/VG(2,1)
	  XNN = XNN + XNK(IND)
	  XMM1 = XMM1 + XM(IND)/VG(2,I)
	  XLNVG = XLNVG + DLNVG(2,I)*XNK(IND)
 311    CONTINUE
        DO 312 I = 1,NSEC
          IND = I+NDISCT
          XNK(IND) = XM(IND)/VL(I+1)/VL(I)
	  XN(IND) = XM(IND)/VGSSQ(I)/DVK(I)
	  XNN = XNN + XNK(IND)
          XMM1 = XMM1 + XM(IND)*DLNSF/DVK(I)
	  XLNVG = XLNVG + XM(IND)/DVK(I)*((1.0+DLNVL(I))/VL(I)
     &            -(1.0+DLNVL(I+1))/VL(I+1))
 312    CONTINUE
C
      XLNVG = XLNVG/XNN         
      TVG = DEXP(XLNVG)
      VM = XMM1/XNN
      TTVOL = XMM1 + XM(1)/VG(1,1) + XM(NDISCA)/VG(2,1) + 
     &        (CONCSP(1)*VG(1,1)+CONCSP(2)*VG(2,1))*AVO
      TTVOL2 = XMM1 + XM(1)/VG(1,1) + XM(NDISCA)/VG(2,1)
      XNN2 = XNN + XM(1)/VG(1,1)**2 + XM(1+NDISC(1))/VG(2,1)**2
      XNK(1) = XM(1)/VG(1,1)**2
	XN(1) = XNK(1)/VG(1,1)
C
      XLN2SIG = 0.0
      DO 380 I = 2,NDISC(1)
	  XLN2SIG = XLN2SIG + (DLNVG(1,I)-XLNVG)**2*XNK(I)
 380  CONTINUE
      DO 381 I = 2,NDISC(2)
	  XLN2SIG=XLN2SIG+(DLNVG(2,I)-XLNVG)**2*XNK(I+NDISC(1))
 381  CONTINUE
      DO 382 I = 1,NSEC
	  XLN2SIG = XLN2SIG+XM(I+NDISCT)/DVK(I)*(((DLNVL(I)-XLNVG+1.0)**2
     &          +1.0)/VL(I) - ((DLNVL(I+1)-XLNVG+1.0)**2+1.0)
     &          /VL(I+1))
 382  CONTINUE
      XLN2SIG = XLN2SIG/XNN/9.0
      TXLN2SIG = DEXP(DSQRT(XLN2SIG))
      END IF
C
      OPEN(2,FILE=RESULT2)
      WRITE(2,110)COUNT
      WRITE(2,111)ETA
      WRITE(2,112)(VG(1,I),XN(I),XM(I),I=2,NDISC(1))
      WRITE(2,112)(VG(2,I),XN(I+NDISC(1)),XM(I+NDISC(1)),I=2,NDISC(2))
      WRITE(2,112)(VGS(I-NDISCT),XN(I),XM(I),I=NDISC1,M)
      WRITE(2,*)
      WRITE(2,118)XNN,TVG,TXLN2SIG,XMM1,VM
      WRITE(2,119)TTVOL
      WRITE(2,*)
      WRITE(2,*)
 110  FORMAT(1X,1PE16.7,'  Sec')
 111  FORMAT(1X,'ETA = ',I1)
 112  FORMAT(1X,1P3E16.7E3)
 118  FORMAT(1X,1P5E16.7)
 119  FORMAT(1X,1PE16.7)
      OPEN(3,FILE='selpre.dat')
      WRITE(3,110)COUNT
      WRITE(3,111)ETA
      DO I=1,NDISC(1)
      WRITE(3,112)VG(1,I)*XNN2/TTVOL2,XN(I)*TTVOL2/XNN2**2
      ENDDO
      DO I=NDISC1,M
      WRITE(3,112)VGS(I-NDISCT)*XNN2/TTVOL2,XN(I)*TTVOL2/XNN2**2
      ENDDO
      WRITE(3,*)
      WRITE(3,118)XNN,TVG,TXLN2SIG,XMM1,VM
      WRITE(3,119)TTVOL
      WRITE(3,*)
      WRITE(3,*)
C      
      OPEN(4,FILE='sdf.dat')
      WRITE(4,120)COUNT,XNN,TVG,TXLN2SIG,XMM1,VM,TTVOL
 120  FORMAT(1X,1P7E16.7E3)
C
      RETURN
      END
C
C
C
      SUBROUTINE XDOTEQ(XM,DM,COUNT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine is written to determine the differential equations 
C     required by RK45.  It is divided into 5 parts: First discrete size
C     , other discrete sizes, first section, other sections and the 
C     vapor.  In each part, the interactions between the discrete-
C     discrete, section-discrete, section-section and condensation are
C     considered
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (MAXDISC = 300, MAXSEC = 300)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NDISC(2),NDISCT,NDISC1,M,KFSB(2),ETA,NSEC,NCAL1,NCAL
      INTEGER MAERO,NCONSCND
      REAL*8 XM(MAXDISC+MAXDISC+MAXSEC),DM(MAXDISC+MAXDISC+MAXSEC)
      REAL*8 VG(2,MAXDISC),VG2(2,MAXDISC)
      REAL*8 VGETA(2,MAXDISC),XMFAC(2,MAXDISC,MAXDISC)
      REAL*8 VL(MAXSEC),VGS(MAXSEC),VGSN(MAXSEC)
      REAL*8 CONVV(MAXDISC),CONV3(MAXDISC),CONV2(MAXDISC)
      REAL*8 TMP(8),XNPD(2,MAXDISC),XNPDS(2,MAXSEC)
      REAL*8 TMC(2,MAXDISC)
      REAL*8 TMCS(2,MAXSEC),TMSCX(2,MAXSEC)
      REAL*8 TMCX(2,MAXDISC),TMCXS(2,MAXSEC)
      REAL*8 BDD(2,MAXDISC,MAXDISC)
      REAL*8 BDS1(2,MAXDISC,MAXSEC,MAXSEC),BDS25(2,MAXDISC,MAXSEC)
      REAL*8 BDS4(2,MAXDISC,MAXSEC)
      REAL*8 BSS1(MAXSEC,MAXSEC,MAXSEC),BSS25(MAXSEC,MAXSEC)
      REAL*8 BSS36(MAXSEC),BSS4(MAXSEC,MAXSEC)
      REAL*8 XXMM(MAXDISC+MAXDISC+MAXSEC,MAXDISC+MAXDISC+MAXSEC)
      REAL*8 XNUCL(2),CONCSP(2),RHOSEC(MAXSEC),DRHO(MAXSEC)
      REAL*8 XK2RATIO(2)
      REAL C13
      REAL*8 XR1(2),XK1(2),XK2(2)
      COMMON /AERO1/ VL,VG,VGS,VGSN
      COMMON /AERO3/ SIZEFAC,VG2,KFSB,XR1
      COMMON /AERO4/ RHOSEC
      COMMON /B/ BDD,BDS1,BDS25,BDS4,BSS1,BSS25,BSS36,BSS4
      COMMON /RXN/ XK1,XK2,CO2
      COMMON /CALCON1/NDISC,NDISCT,NDISC1,NDISCA,M,NSEC,XNCAL,NCAL1,NCAL
      COMMON /CONSTANT/ AVO,PI,XKB,ETA
      COMMON /PRES/ XNPD,XNPDS,CONCSP
      COMMON /EV/ VGETA,XMFAC,C2E
      COMMON /B1/ CONVV,CONV3,C13,CONV2
      COMMON /TMOD/ XTMOD,DRHO
      COMMON /K2R/ XK2RATIO,MAERO
      COMMON /TEST/ NCONSCND
      common /ttt/ kkk
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     BETA(I,J) Coagulation coefficients between particle I and J 
C               (cm3/#/s)
C     FRAC      Fraction
C     KFSB      First section boundary for the first section formation
C               by discrete sizes interaction 
C     M         Total number of discrete size & sections (NDISC + NSEC)
C     NDISC     Total number of discrete sizes
C     NDISC1    NDISC + 1
C     NDISC2    NDISC + 2
C     SIZEFAC   Size increasing factor based on volume
C     TMC(I)    Removal rate of discrete I by condensation on discrete I 
C               (#/cc/s) 
C     TMCS(I)   Removal & formation rate of section I by condensation on 
C               section I(#/cc/s) 
C     TMCX(I)   Total monomer condensation rate on discrete I (#/cc/s)
C     TMCXS(I)  Total monomer condensation rate on section I (#/cc/s)
C     TMP       Temporary variables
C     TMSCS(I)  Formation rate of section I+1 by condensation on section
C               I (#/cc/s)
C     TTMC      Total condensation  (#/cc)
C     VG        Mean volume of the section (cm3)
C     VL        Lower bound aerosol volume of the section (cm3)
C     XK2       Dimer nucleation rate (cm3/#/s)
C     XM        Aerosol and vapor concentration
C       XM(1)                   Vapor molecule concentration (#/cc)
C       XM(2...NDISC)           Discrete size aerosols (#/cc)
C       XM(NDISC1...M)       Sectional size aerosols (#/cc)
C     XNMON     Mean number of molecules in one aerosol in the section
C     XNMON2    2.0*XNMON
C     XNPD(I)   Surface saturation concentration for section I particle 
C               (#/cc)
C     XNUCL     Dimer nucleation rate (#/cc/s)
C     XR1       Aerosol vapor formation rate (#/cc/s)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C***********************************************************************
C     Determine the concentration multiplication
C***********************************************************************
      DO 10 I = 1,NDISC(1)
        XXMM(I,I) = XM(I)*XM(I)*BDD(1,I,I)*XTMOD
        DO 8 J = I+1,NDISC(1)
          XXMM(I,J) = XM(I)*XM(J)*XTMOD
          XXMM(J,I) = XXMM(I,J)*BDD(1,J,I)
          XXMM(I,J) = XXMM(I,J)*BDD(1,I,J)
   8    CONTINUE
        DO 9 J = NDISC1,M
          XXMM(I,J) = XM(I)*XM(J)*XTMOD
   9    CONTINUE
  10  CONTINUE
      DO 15 I = 1,NDISC(2)
        IND1 = I+NDISC(1)
        XXMM(IND1,IND1) = XM(IND1)*XM(IND1)*BDD(2,I,I)*XTMOD
        DO 12 J = I+1,NDISC(2)
          IND2 = J+NDISC(1)
          XXMM(IND1,IND2) = XM(IND1)*XM(IND2)*XTMOD
          XXMM(IND2,IND1) = XXMM(IND1,IND2)*BDD(2,J,I)
          XXMM(IND1,IND2) = XXMM(IND1,IND2)*BDD(2,I,J)
  12    CONTINUE
        DO 14 J = NDISC1,M
          XXMM(IND1,J) = XM(IND1)*XM(J)*XTMOD
  14    CONTINUE
  15  CONTINUE
      DO 20 I = NDISC1,M
        XXMM(I,I) = XM(I)*XM(I)*XTMOD
        DO 18 J = I+1,M
          XXMM(I,J) = XM(I)*XM(J)*XTMOD
          XXMM(J,I) = XXMM(I,J)
  18    CONTINUE
  20  CONTINUE
C***********************************************************************
C     First aerosol size
C     XNUCL: formation
C     TMP(1): removal
C     TMC(2): condensation (removal)
C***********************************************************************
      TMP(1) = 0.0
      XNUCL(1) = BDD(1,1,1)*XM(1)*XM(1)*XTMOD*XK2RATIO(1)
      DO 50 J = 2,NDISC(1)
	  TMP(1) = TMP(1) + XXMM(2,J)
  50  CONTINUE
      DO 55 J = 1,NSEC
	  TMP(1) = TMP(1) + BDS4(1,2,J)*XXMM(2,J+NDISCT)
  55  CONTINUE
      IF (XM(1).GE.XNPD(1,2)) THEN
        TMP(2) = (XM(1)-XNPD(1,2))*XM(2)*XTMOD
        TMCX(1,2) = BDD(1,1,2)*TMP(2)
        TMC(1,2) = BDD(1,2,1)*TMP(2)
      END IF
      DM(2) = XNUCL(1) - TMP(1) - TMC(1,2)
C
      TMP(1) = 0.0
      IND2 = NDISC(1)+2
      XNUCL(2) = BDD(2,1,1)*XM(NDISCA)*XM(NDISCA)*XTMOD*XK2RATIO(2)
      DO 60 J = 2,NDISC(2)
	  TMP(1) = TMP(1) + XXMM(IND2,NDISC(1)+J)
  60  CONTINUE
      DO 65 J = 1,NSEC
	  TMP(1) = TMP(1) + BDS4(2,2,J)*XXMM(IND2,J+NDISCT)
  65  CONTINUE
      IF (XM(NDISCA).GE.XNPD(2,2)) THEN
        TMP(2) = (XM(NDISCA)-XNPD(2,2))*XM(IND2)*XTMOD
        TMCX(2,2) = BDD(2,1,2)*TMP(2)
        TMC(2,2) = BDD(2,2,1)*TMP(2)
      END IF
C
      DM(IND2) = XNUCL(2) - TMP(1) - TMC(2,2)
C***********************************************************************
C     ###   Other discrete sizes   ###
C***********************************************************************
      DO 200 I = 3,NDISC(1)
      DO 80 LL =1,2
	  TMP(LL) = 0.0
  80  CONTINUE
C=======================================================================
C     Removal due to the collisions with other discrete size or section
C=====
      DO 120 J = 2,NDISC(1)
	 TMP(1) = TMP(1) + XXMM(I,J)
 120  CONTINUE
      DO 125 J = 1,NSEC
	 TMP(1) = TMP(1) + BDS4(1,I,J)*XXMM(I,J+NDISCT)
 125  CONTINUE
C=======================================================================
C     Formation from the collisions of smaller discrete sizes
C=====
      DO 180 J = 2,I-2
	 TMP(2) = TMP(2) + XXMM(J,I-J)*XMFAC(1,J,I-J)
 180  CONTINUE
C=======================================================================
C     Condensation
C=====
      IF (XM(1).GE.XNPD(1,I)) THEN
        TMP(3) = (XM(1)-XNPD(1,I))*XM(I)*XTMOD 
        TMCX(1,I) = BDD(1,1,I)*TMP(3)
        TMC(1,I) = BDD(1,I,1)*TMP(3)
      END IF
C
      DM(I) = - TMP(1) + TMP(2) - TMC(1,I) + (TMC(1,I-1)+TMCX(1,I-1))
     &        *XMFAC(1,I-1,1)
C
 200  CONTINUE
C
      DO 199 I = 3,NDISC(2)
      DO 182 LL =1,2
	  TMP(LL) = 0.0
 182  CONTINUE
      IND = I + NDISC(1)
C=======================================================================
C     Removal due to the collisions with other discrete size or section
C=====
      DO 183 J = 2,NDISC(2)
	 TMP(1) = TMP(1) + XXMM(IND,NDISC(1)+J)
 183  CONTINUE
      DO 184 J = 1,NSEC
	 TMP(1) = TMP(1) + BDS4(2,I,J)*XXMM(IND,J+NDISCT)
 184  CONTINUE
C=======================================================================
C     Formation from the collisions of smaller discrete sizes
C=====
      DO 186 J = 2,I-2
	 TMP(2) = TMP(2) + XXMM(NDISC(1)+J,IND-J)*XMFAC(2,J,I-J)
 186  CONTINUE
C=======================================================================
C     Condensation
C=====
      IF (XM(NDISCA).GE.XNPD(2,I)) THEN
        TMP(3) = (XM(NDISCA)-XNPD(2,I))*XM(IND)*XTMOD
        TMCX(2,I) = BDD(2,1,I)*TMP(3)
        TMC(2,I) = BDD(2,I,1)*TMP(3)
      END IF
C
      DM(IND) = -TMP(1) + TMP(2) - TMC(2,I) + 
     &        (TMC(2,I-1)+TMCX(2,I-1))*XMFAC(2,I-1,1)
C
 199  CONTINUE
C***********************************************************************
C    ###   First section   ###
C***********************************************************************
      DO 210 I = 1,6
	  TMP(I) = 0.0
 210  CONTINUE
C=======================================================================
C     Formation from collisions of 2 discrete sizes
C=====
      DO 325 I = 2,NDISC(1)
        IF ((VG2(1,I).GE.VL(1)).AND.(VG2(1,I).LE.VL(2))) 
     &     TMP(1) = TMP(1) + XXMM(I,I)*C2E
	DO 320 J = I+1,NDISC(1)
	  VTEMP = VG(1,I) + VG(1,J)
	  IF ((VTEMP.GE.VL(1)).AND.(VTEMP.LE.VL(2)))
     &       TMP(1) = TMP(1) + (XXMM(I,J)+XXMM(J,I))*XMFAC(1,J,I)
 320    CONTINUE
 325  CONTINUE
C
      DO 315 I = 2,NDISC(2)
        IND1 = NDISC(1)+I
        IF ((VG2(2,I).GE.VL(1)).AND.(VG2(2,I).LE.VL(2)))
     &    TMP(6) = TMP(6) + XXMM(IND1,IND1)*C2E
	DO 310 J = I+1,NDISC(2)
	  VTEMP = VG(2,I) + VG(2,J)
	  IF ((VTEMP.GE.VL(1)).AND.(VTEMP.LE.VL(2))) THEN
	    IND2 = NDISC(1)+J
            TMP(6)=TMP(6)+(XXMM(IND1,IND2)+XXMM(IND2,IND1))*XMFAC(2,J,I)
          END IF
 310    CONTINUE
 315  CONTINUE
C=======================================================================
C     Formation and removal from collisions of one discrete size and 
C     the first section
C=====
      DO 330 I = 2,NDISC(1)
	  TMP(2) = TMP(2) + BDS25(1,I,1)*XXMM(I,NDISC1)
 330  CONTINUE
C
      DO 340 I = 2,NDISC(2)
	  TMP(2) = TMP(2) + BDS25(2,I,1)*XXMM(NDISC(1)+I,NDISC1)
 340  CONTINUE
C=======================================================================
C     Formation and removal of the first section by the collision of 2
C     first section particles
C=====
      TMP(3) = BSS36(1)*XXMM(NDISC1,NDISC1)
C=======================================================================
C     Removal from collisions of the first section with larger sections
C=====
      DO 400 I = 2,NSEC
	  TMP(4) = TMP(4) + BSS4(1,I)*XXMM(NDISC1,I+NDISCT)
 400  CONTINUE
C=======================================================================
C     Condensation
C=====
      IF (XM(1).GE.XNPDS(1,1)) THEN
	  TMP(5) = (XM(1)-XNPDS(1,1))*XM(NDISC1)*XTMOD
	  TMCXS(1,1) = BDS4(1,1,1)*TMP(5)
	  TMCS(1,1) = BDS25(1,1,1)*TMP(5)
	  TMSCX(1,1) = BDS1(1,1,1,2)*TMP(5)
      END IF
      IF (XM(NDISCA).GE.XNPDS(2,1)) THEN
	  TMP(5) = (XM(NDISCA)-XNPDS(2,1))*XM(NDISC1)*XTMOD
	  TMCXS(2,1) = BDS4(2,1,1)*TMP(5)
	  TMCS(2,1) = BDS25(2,1,1)*TMP(5)
	  TMSCX(2,1) = BDS1(2,1,1,2)*TMP(5)
      END IF
C
      DM(NDISC1) = TMP(1)-TMP(2)-TMP(3)-TMP(4)+TMP(6)-TMCS(1,1)+
     &(TMC(1,NDISC(1))+TMCX(1,NDISC(1)))*XMFAC(1,NDISC(1),1)
      DM(NDISC1) = DM(NDISC1)-TMCS(2,1)+(TMC(2,NDISC(2))+
     &  TMCX(2,NDISC(2)))*XMFAC(2,NDISC(2),1)
C***********************************************************************
C     ###   Other sections   ###
C***********************************************************************
      DO 900 I = 2,NSEC
      DO 450 LL = 1,8
	  TMP(LL) = 0.0
 450  CONTINUE
      IND = I+NDISCT
C=======================================================================
C     Formation from collision of discrete sizes
C=====
      DO 490 J = 2,NDISC(1)
        IF ((VG2(1,J).GE.VL(I)).AND.(VG2(1,J).LE.VL(I+1))) 
     &    TMP(1) = TMP(1) + XXMM(J,J)*C2E
	DO 485 K = J+1,NDISC(1)
	  VTEMP = VG(1,J) + VG(1,K)
	  IF ((VTEMP.GE.VL(I)).AND.(VTEMP.LE.VL(I+1)))
     &	    TMP(1) = TMP(1) + (XXMM(J,K)+XXMM(K,J))*XMFAC(1,K,J)
  485   CONTINUE
  490 CONTINUE
      DO 500 J = 2,NDISC(2)
        IND1 = J + NDISC(1)
        IF ((VG2(2,J).GE.VL(I)).AND.(VG2(2,J).LE.VL(I+1))) 
     &    TMP(1) = TMP(1) + XXMM(IND1,IND1)*C2E
	DO 495 K = J+1,NDISC(2)
	  VTEMP = VG(2,J) + VG(2,K)
	  IF ((VTEMP.GE.VL(I)).AND.(VTEMP.LE.VL(I+1))) THEN
	    IND2 = NDISC(1)+K
     	    TMP(1)=TMP(1)+(XXMM(IND1,IND2)+XXMM(IND2,IND1))*XMFAC(2,K,J)
          END IF
  495   CONTINUE
  500 CONTINUE
C=======================================================================
C     Formation from collisions of one discrete size and one smaller 
C     section
C=====
      DO 520 J = 1,I-1
        IND1 = J+NDISCT
        IF ((VL(I)-VL(J+1)).GT.VG(1,NDISC(1))) GOTO 512
	DO 510 K = 2,NDISC(1)
	  TMP(2) = TMP(2) + BDS1(1,K,J,I)*XXMM(K,IND1)
 510    CONTINUE
 512    IF ((VL(I)-VL(J+1)).GT.VG(2,NDISC(2))) GOTO 525
	DO 515 K = 2,NDISC(2)
	  TMP(2) = TMP(2) + BDS1(2,K,J,I)*XXMM(K+NDISC(1),IND1)
 515    CONTINUE
C=======================================================================
C     Formation from collisions of smaller sections
C     K:  Smaller section
C=====
 525    IF (VL(I).GE.(2.0D0*VL(J+1))) GOTO 520
	TMP(3) = TMP(3) + BSS1(I,J,J)*XXMM(IND1,IND1)
	DO 540 K = 1,J-1
	  TMP(3) = TMP(3) + BSS1(I,K,J)*XXMM(K+NDISCT,IND1)
 540    CONTINUE
 520  CONTINUE
C=======================================================================
C     Removal and formation of the specific section by collision with 
C     discrete size
C=====
 560  DO 600 J = 2,NDISC(1)
	  TMP(4) = TMP(4) + BDS25(1,J,I)*XXMM(J,IND)
 600  CONTINUE
      DO 610 J = 2,NDISC(2)
          TMP(4) = TMP(4) + BDS25(2,J,I)*XXMM(J+NDISC(1),IND)
 610  CONTINUE
C=======================================================================
C     Removal and formation of the specific section by collision with 
C     a smaller section
C     FRAC2:  formation
C=====
      DO 650 J = 1,I-1
	TMP(5) = TMP(5) + BSS25(I,J)*XXMM(IND,J+NDISCT)
 650  CONTINUE
C=======================================================================
C     Removal and formation of the section due to collision of the 2 
C     same sections
C     FRAC2: formation
C=====
      TMP(6) = TMP(6) + BSS36(I)*XXMM(IND,IND)
C=======================================================================
C     Removal of the section due to collsions of the section with a 
C     larger section
C=====
      DO 680 J = I+1,NSEC
	  TMP(7) = TMP(7) + BSS4(I,J)*XXMM(IND,J+NDISCT)
 680  CONTINUE
C=======================================================================
C     Condensation
C=====
      IF (XM(1).GE.XNPDS(1,I)) THEN
	  TMP(8) = (XM(1)-XNPDS(1,I))*XM(IND)*XTMOD
	  TMCXS(1,I) = BDS4(1,1,I)*TMP(8)
	  TMCS(1,I) = BDS25(1,1,I)*TMP(8)
	  TMSCX(1,I) = BDS1(1,1,I,I+1)*TMP(8)
      END IF
      IF (XM(NDISCA).GE.XNPDS(2,I)) THEN
	  TMP(8) = (XM(NDISCA)-XNPDS(2,I))*XM(IND)*XTMOD
	  TMCXS(2,I) = BDS4(2,1,I)*TMP(8)
	  TMCS(2,I) = BDS25(2,1,I)*TMP(8)
	  TMSCX(2,I) = BDS1(2,1,I,I+1)*TMP(8)
      END IF
C
      DM(IND) = TMP(1)+TMP(2)+TMP(3)-(TMP(4)+TMP(5)+TMP(6)+TMP(7))
     &-TMCS(1,I)+TMSCX(1,I-1)-TMCS(2,I)+TMSCX(2,I-1)
 900  CONTINUE
C*********************************************************** 
C     ###   Removal and generation of the FeO vapor   ###
C***********************************************************
      IF (NCONSCND.EQ.1) THEN 
        DM(1) = 0.0 
        DM(NDISCA) = 0.0
      ELSE
        TTMC = 0.0
        DO 950 I = 2,NDISC(1)
	  TTMC = TTMC + TMCX(1,I)
 950    CONTINUE
        DO 952 I = 1,NSEC
	  TTMC = TTMC + TMCXS(1,I)
 952    CONTINUE
        DM(1) = XR1(1)*XM(M+1) - TTMC - XNUCL(1)/C2E
C
        TTMC = 0.0
        DO 960 I = 2,NDISC(2)
	  TTMC = TTMC + TMCX(2,I)
 960    CONTINUE
        DO 962 I = 1,NSEC
	  TTMC = TTMC + TMCXS(2,I)
 962    CONTINUE
        DM(NDISCA) = XR1(2)*XM(MAERO) - TTMC - XNUCL(2)/C2E
      END IF
C*******************************************
C     ###   Precursor decomposition   ###
C*******************************************
      DM(M+1) = -XK1(1)*XM(M+1) !*CO2
      DM(MAERO) = -XK1(2)*XM(MAERO)*CO2
C
      RETURN
      END
C
C
C
      SUBROUTINE RK45(X,XDOT,T,H,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXDISC = 300, MAXSEC = 300)
      COMMON/AABB/HA(6),B(6,5),HC1,HC3,HC4,HC5,HC6
      DIMENSION XNEW(MAXDISC+MAXDISC+MAXSEC),F(MAXDISC+MAXDISC+MAXSEC,6)
      DIMENSION SAVE(MAXDISC+MAXDISC+MAXSEC)
      DIMENSION X(*),XDOT(*)
      DATA ZERO/0.0D0/
C
C        RUNGE-KUTTA INTEGRATION SCHEME DEVELOPED BY DALE BETTIS OF THE
C        UNIVERSITY OF TEXAS AT AUSTIN.  LOCAL TRUNCATION ERROR IS H**6.
C        OPTIMIZED COEFFICIENTS VERSION.
C
C        INIT0 MUST BE CALLED ONCE BEFORE PROGRAM EXECUTION
C        AND EACH TIME THE STEPSIZE IS CHANGED.
C
C        **********  DO NOT CHANGE ANYTHING IN THIS ROUTINE EVER!!! *********
C
C        ALL INFO COMES IN AND OUT THROUGH THE COMMON BLOCK AND THE ARGUMENT
C        LIST.
C
      CALL XDOTEQ(X,F(1,1),T)
      DO 40 K=2,6
      N=K-1
      DO 20 J=1,M
      SAVE(J)=ZERO
      DO 10 LAM=1,N
   10 SAVE(J)=SAVE(J)+B(K,LAM)*F(J,LAM)
   20 SAVE(J)=H*SAVE(J)
	TNEW=T+HA(K)
	DO 30 N=1,M
   30 XNEW(N)=X(N)+SAVE(N)
      CALL XDOTEQ(XNEW,F(1,K),TNEW)
   40 CONTINUE
	DO 50 J=1,M
   50 X(J)=X(J)+HC1*F(J,1)+HC3*F(J,3)+HC4*F(J,4)+HC5*F(J,5)+
     Q         HC6*F(J,6)
      RETURN
      END
C
	SUBROUTINE INIT0(H)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AABB/HA(6),B(6,5),HC1,HC3,HC4,HC5,HC6
C
C        **********  DO NOT CHANGE ANYTHING IN THIS ROUTINE EVER!!! *********
C
      HA(1)=0.0D0
      HA(2)=H*.2D0
      HA(3)=H*.3D0
      HA(4)=H*.6D0
      HA(5)=H
      HA(6)=H*11.D0/12.D0
      HC1=H*59.D0/594.D0
      HC3=H*2750.D0/6993.D0
      HC4=H*125.D0/513.D0
      HC5=H*(-1.D0/14.D0)
      HC6=H*2592.D0/7733.D0
      B(2,1)=.2D0
      B(3,1)=3.D0/40.D0
      B(3,2)=9.D0/40.D0
      B(4,1)=.3D0
      B(4,2)=-.9D0
      B(4,3)=1.2D0
      B(5,1)=-11.D0/54.D0
      B(5,2)=2.5D0
      B(5,3)=-70.D0/27.D0
      B(5,4)=35.D0/27.D0
      B(6,1)=-473.D0/10368.D0
      B(6,2)=1595.D0/1728.D0
      B(6,3)=-34595.D0/54432.D0
      B(6,4)=38665.D0/62208.D0
      B(6,5)=7733.D0/145152.D0
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION FX(V)
C
      PARAMETER (MAXDISC = 300,MAXSEC = 300,C12=0.6666666667D0)
      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      REAL*8 VG(2,MAXDISC),VL(MAXSEC),VGS(MAXSEC),VGSN(MAXSEC)
      INTEGER I,J,K,ETA,IDD
      REAL*8 CONVV(MAXDISC),CONV3(MAXDISC),CONV2(MAXDISC)
      REAL C13
      COMMON /AERO1/ VL,VG,VGS,VGSN
      COMMON /B1/ CONVV,CONV3,C13,CONV2
      COMMON /FFF/ I,J,K,IDD,II
      COMMON /CONSTANT/ AVO,PI,XKB,ETA
C
      IF (IDD.EQ.1) FX = DSQRT(CONVV(K)+VG(II,1)/V)*(CONV3(K)
     &   +(V/VG(II,1))**C13)**2*(1.0/VG(II,K)+1.0/V)**ETA
      IF (IDD.EQ.2) FX = DSQRT(CONVV(K)+VG(II,1)/V)*(CONV3(K)
     &   +(V/VG(II,1))**C13)**2
      IF (IDD.EQ.4) FX = DSQRT(CONVV(K)+VG(II,1)/V)*(CONV3(K)
     &   +(V/VG(II,1))**C13)**2/V**ETA
      IF (IDD.EQ.6) FX = (V/VG(II,1))**C12
      IF (IDD.EQ.7) FX = (V/VG(II,1))**C12/V**ETA
      IF (IDD.EQ.9) FX = (V/VG(II,1))**C12*(1.0/VG(II,1)+1.0/V)**ETA
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION F(U,V)
C
      PARAMETER (MAXDISC = 300,MAXSEC = 300)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CONVV(MAXDISC),CONV3(MAXDISC),V1(2),RHO(2),CONV2(MAXDISC)
      REAL C13
      INTEGER I,J,K,ETA,IDD
      COMMON /MONOMER/ V1,RHO
      COMMON /B1/ CONVV,CONV3,C13,CONV2
      COMMON /FFF/ I,J,K,IDD,II
      COMMON /CONSTANT/ AVO,PI,XKB,ETA
C
      IF (IDD.GE.15) THEN
         F = DSQRT(V1(1)/U+V1(1)/V)*((U/V1(1))**C13
     &   +(V/V1(1))**C13)**2*(1.0/U+1.0/V)**ETA
      ELSEIF (IDD.EQ.13) THEN
         F = DSQRT(V1(1)/U+V1(1)/V)*((U/V1(1))**C13
     &     +(V/V1(1))**C13)**2*(1.0/U**ETA+1.0/V**ETA)
      ELSE
         F = DSQRT(V1(1)/U+V1(1)/V)*((U/V1(1))**C13
     &     +(V/V1(1))**C13)**2/U**ETA
      END IF
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION G(U)
C
      PARAMETER (MAXDISC = 300,MAXSEC = 300)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VG(2,MAXDISC),VL(MAXSEC),VGS(MAXSEC),VGSN(MAXSEC)
      INTEGER I,J,K,IDD
      COMMON /AERO1/ VL,VG,VGS,VGSN
      COMMON /FFF/ I,J,K,IDD,II
      IF (IDD.EQ.12) G = VL(I)
      IF (IDD.EQ.13) G = VL(J)
      IF (IDD.EQ.14) G = VL(J)
      IF (IDD.EQ.15) G = VL(I)
      IF (IDD.EQ.16) G = VL(J)
      IF (IDD.EQ.17) G = VL(I)
      IF (IDD.EQ.21) G = VL(J)
      IF (IDD.EQ.22) G = VL(J)
      IF (IDD.EQ.23) G = VL(K)
      IF (IDD.EQ.24) G = VL(I) - VL(K)
      IF (IDD.EQ.25) G = VL(I) - U
      IF (IDD.EQ.26) G = VL(I) - U
      IF (IDD.EQ.31) G = VL(K)
      IF (IDD.EQ.32) G = VL(K)
      IF (IDD.EQ.33) G = VL(K)
      IF (IDD.EQ.34) G = VL(K)
      IF (IDD.EQ.35) G = VL(K)
      IF (IDD.GE.36) G = VL(I) - U
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION H(U)
C
      PARAMETER (MAXDISC = 300,MAXSEC = 300)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VG(2,MAXDISC),VL(MAXSEC),VGS(MAXSEC),VGSN(MAXSEC)
      INTEGER I,J,K,IDD
      COMMON /AERO1/ VL,VG,VGS,VGSN
      COMMON /FFF/ I,J,K,IDD,II
C
      IF (IDD.EQ.12) H = VL(I+1)
      IF (IDD.EQ.13) H = VL(J+1)
      IF (IDD.EQ.14) H = VL(J+1)
      IF (IDD.EQ.15) H = VL(I+1) - U
      IF (IDD.EQ.16) H = VL(J+1) - U
      IF (IDD.EQ.17) H = VL(I+1) - U
      IF (IDD.EQ.21) H = VL(I+1) - U
      IF (IDD.EQ.22) H = VL(I+1) - VL(K+1)
      IF (IDD.EQ.23) H = VL(I+1) - U
      IF (IDD.EQ.24) H = VL(J+1)
      IF (IDD.EQ.25) H = VL(K+1)
      IF (IDD.EQ.26) H = VL(J+1)
      IF (IDD.EQ.31) H = VL(I+1) - U
      IF (IDD.EQ.32) H = VL(I+1) - U
      IF (IDD.EQ.33) H = VL(K+1)
      IF (IDD.EQ.34) H = VL(I+1) - U
      IF (IDD.GE.35) H = VL(K+1)
C
      RETURN
      END
C
C
C
C PROGRAM QUAD.FOR	
C
       FUNCTION QUAD(GLEGQ,N,FX,A,B)
C
C
C THE DOUBLE PRECISION FUNCTION SUBPROGRAM APPROXIMATES  THE VALUE OF
C AN INTEGRAL I(FX,A,B) OF A FUNCTION FX(X) OVER THE INTERVAL (A,B) 
C USING A QUADRATURE RULE WHICH IS SPECIFIED BY THE PARAMETERS METHOD
C AND N.
C
C
C CALLING SEQUENCE: APPROX= QUAD(METHOD,N,FX,A,B)
C
C
C  PARAMETRS:
C ON ENTRY:
C        METHOD-DOUBLE PRECISION SUBROUTINE METHOD(N,A,B,P,W)
C            USED TO SPECIFY THE TYPE OF QUADRATURE RULE TO BE
C            USED AND TO OBTAIN THE APPROPRIATE POINTS AND WEIGHTS 
C             =NCQ- NEWTON-COTES RULES,
C              GLFGQ-GAUSS-LEGENDRE RULES,
C              GLAGQ-GAUSS-LAGUERRE RULES.E
C       N -   INTEGER           
C             THE NUMBER OF QUADRATURE   POINTS TO BE USED.   LIMITS   
C             ON N WHEN METHOD= NCQ  : 1 <= N <= 21,
C                               GLEGQ: 2 <=N  <=20,
C                               GLAG : 2<= N  <=10,
C       FX -  DOUBLE  PRECISION FUNCTION FX(X)    
C             USED FOR EVALUATING THE INTEGRAND.
C      A,B- REAL*8
C            THE LEFT AND RIGHT ENDPOINTS OF THE INTERVAL OF
C                  INTEGRATION.
C 
C   ON RETURN:
C           QUAD - REAL*8
C            THE APPROXIMATE VALUE OF  I(FX,A,B).
C          SAMPLE CALLING PROGRAM:
C    C********** THIS PROGRAM APPROXIMATES THE VALUE OF I(FX,A,B) USING A
C    C           COMPOSITE 5-POINT GAUSS-LEGENDRE RULE. THE NUMBER OF         
C    C          SUBDIVISIONS OF (A,B) IS MDIV.
C               IMPLICIT REAL*8 (A-H,O-Z)
C               EXTERNAL    FX,GLEGO
C              READ (5,*) A,B,MDIV
C              H=(B-A)/MDIV
C              XL=A
C              XR=A+H
C              APPROX=0.DO
C              DO 10 IDIV=1,MDIV
C              APPROX=APPROX+QUAD(GLEGQ,5,FX,XL,XR)
C               XL=XR      
C               XR=XR +H
C            10 CONTINUE 
C   
C              WRITE(6,*) 'APPROX. VALUE OF THE INTEGRAL=',APPROX
C              STOP
C               END
C
C            FUNCTION FX(X)
C             IMPLICIT REAL*8 (A-H,O-Z)
C              FX = .... FORMULA TO EVALUATE THE INTEGRAND...
C              RETURN
C              END
C             
C    REFERENCE : ELEMENTARY NUMERICAL ANALYSIS
C                BY CONTE AND DEBOOR
C                PAGE 299-306
C
C
               IMPLICIT REAL*8 (A-H,O-Z)
               DIMENSION P(20),W(20)
               COMMON X,Y,C,VAR,NVAR,NN
               EXTERNAL FX,GLEGQ
C
C
C       **********************************************************
C       *
C       * RETRIEVE POINT AND WEIGHTS, TRANSFORMED INTO THE SPECIFIED     
C       *  INTEPVAL (A,B), FROM THE SUBROUTINE METHOD, AND APPLY THE
C       *  QUADRATURE FORMULA.
C       *
C       *       
C       **********************************************************
C
C
C
C        
                CALL GLEGQ (N, A, B, P, W)
                  SUM = 0.D0
                  DO 10 I=1,N
                   SUM = SUM + W(I) * FX(P(I))
 10               CONTINUE
                  QUAD=SUM
                RETURN
                END
C   
C
C**************************************************************************
C
      SUBROUTINE GLEGQ( N, A, B, P, W)
C
C
C   THIS DOUBLE PRECISION SUBROUTINE IS DESIGNED TO BE USED BY EITHER OF 
C   THE ROUTINES QUAD OR ADQUAD FOR APPROXIMATING THE VALUE OF THE 
C   DEFINITE INTEGRAL I(FX,A,B) OF FX(X) OVER THE INTERVAL (A,B) USING AN
C   N-POINT GAUSS-LEGENDRE FORMULA.
C   
C   THE PURPOSE OF GLEGQ IS TO SUPPLY, TO QUAD OR ADQUAD, THE N QUADRA-
C   TURE POINTS AND CORRESPONDING WEIGHTS FOR THE SPECIFIED N-POINT
C   FORMULA FOR THE GIVEN INTERVAL OF INTEGRATION.  THIS IS DONE BY
C   TAKING THE POINTS AND WEIGHTS FOR THE "NORMALIZED" INTERVAL (-1,1)
C   AND TRANSFORMING THEM TO (A,B).
C
C   THIS ROUTINE CAN ALSO BE USED BY THE ROUTINE QUADM TO GENERATE
C   ROUTINE GAUSS-LEGENDRE RULES FOR APPROXIMATING INTEGRALS IN M DIMEN-
C   SIONS, M=2, 3.
C
C
C   CALLING SEQUENCE: CALL GLEGQ(N,A,B,P,W)
C
C   PARAMETERS:
C   ON ENTRY:
C        N    -INTEGER
C               THE NUMBER OF QUADRATURE POINTS TO BE USED.  2<=N<=20
C        A, B -REAL*8
C               THE LEFT AND RIGHT ENDPOINTS OF THE INTERVAL OF
C               INTEGRATION.
C
C   ON RETURN:
C        P    -REAL*8(N)
C               THE POINTS FOR THE N-POINT GUASS-LEGENDRE QUADRATURE
C               RULE ON THE INTERVAL (A,B).
C        W    -REAL*8(N)
C               THE WEIGHTS FOR THE N-POINT GAUSS-LEGENDRE QUADRATURE
C               RULE ON THE INTERVAL (A,B).
C
C
C   SAMPLE CALLING PROGRAM:
C        IMPLICIT REAL*8 (A-H,O-Z)
C        EXTERNAL FX,GLEGQ
C        READ(5,*) A,B,N
C        APPROX=QUAD(GLEGQ,FX,N,A,B)
C        WRITE(6,*) 'THE',N,'-POINT GAUSS-LEGENDRE QUADRATURE APPROX.='
C       +           ,APPROX
C        STOP
C        END
C
C        FUNCTION FX(X)
C        IMPLICIT REAL*8 (A-H,O-Z)
C        FX=... FORMULA TO EVALUATE THE INTEGRAND...
C        RETURN
C        END
C 
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PT(100), WT(100), WZ(9), P(N), W(N)
      COMMON X,Y,C,VAR,NVAR,NN
C
C  *******************
C  *
C  * THE VECTOR PT CONTAINS THE LOCATION OF THE NEGATIVE POINTS IN 
C  * (-1,1) USED IN THE FORMULA. IF N IS EVEN, THE ONLY OTHER POINTS
C  * ARE POSITIVE MIRROR IMAGES OF THE NEGATIVES. FOR N ODD, THERE IS AN
C  * ADDITIONAL POINT AT THE ORIGIN.
C  *
C  *******************
C
      DATA PT( 1),PT( 2),PT( 3),PT( 4),PT( 5),PT( 6),PT( 7),PT( 8),
     +   PT( 9),PT( 10)/-5.773502691896259D-01,-7.745966692414834D-01,
     +   -8.611363115940526D-01,-3.399810435848563D-01,
     +   -9.061798459386640D-01,-5.384693101056832D-01,
     +   -9.324695142031520D-01,-6.612093864662646D-01,
     +   -2.386191860831969D-01,-9.491079123427586D-01/
      DATA PT(11),PT(12),PT(13),PT(14),PT(15),PT(16),PT(17),PT(18),
     +   PT(19),PT( 20)/-7.415311855993942D-01,-4.058451513773972D-01,
     +   -9.602898564975362D-01,-7.966664774136267D-01,
     +   -5.255324099163290D-01,-1.834346424956498D-01,
     +   -9.681602395076261D-01,-8.360311073266358D-01,
     +   -6.133714327005905D-01,-3.242534234038089D-01/
      DATA PT(21),PT(22),PT(23),PT(24),PT(25),PT(26),PT(27),PT(28),
     +   PT(29),PT( 30)/-9.739065285171717D-01,-8.650633666889845D-01,
     +   -6.794095682990245D-01,-4.333953941292472D-01,
     +   -1.488743389816312D-01,-9.782286581460570D-01,
     +   -8.870625997680954D-01,-7.301520055740493D-01,
     +   -5.190961292068119D-01,-2.695431559523450D-01/
      DATA PT(31),PT(32),PT(33),PT(34),PT(35),PT(36),PT(37),PT(38),
     +   PT(39),PT( 40)/-9.815606342467190D-01,-9.041172563704749D-01,
     +   -7.699026741943046D-01,-5.873179542866175D-01,
     +   -3.678314989981802D-01,-1.252334085114689D-01,
     +   -9.841830547185882D-01,-9.175983992229781D-01,
     +   -8.015780907333099D-01,-6.423493394403403D-01/
      DATA PT(41),PT(42),PT(43),PT(44),PT(45),PT(46),PT(47),PT(48),
     +   PT(49),PT( 50)/-4.484927510364468D-01,-2.304583159551348D-01,
     +   -9.862838086968123D-01,-9.284348836635734D-01,
     +   -8.272013150697650D-01,-6.872929048116856D-01,
     +   -5.152486363581541D-01,-3.191123689278898D-01,
     +   -1.080549487073437D-01,-9.879925180204854D-01/
      DATA PT(51),PT(52),PT(53),PT(54),PT(55),PT(56),PT(57),PT(58),
     +   PT(59),PT( 60)/-9.372733924007058D-01,-8.482065834104270D-01,
     +   -7.244177313601699D-01,-5.709721726085388D-01,
     +   -3.941513470775634D-01,-2.011940939974345D-01,
     +   -9.894009349916499D-01,-9.445750230732326D-01,
     +   -8.656312023878317D-01,-7.554044083550030D-01/
      DATA PT(61),PT(62),PT(63),PT(64),PT(65),PT(66),PT(67),PT(68),
     +   PT(69),PT( 70)/-6.178762444026438D-01,-4.580167776572274D-01,
     +   -2.816035507792589D-01,-9.501250983763744D-02,
     +   -9.905754753144173D-01,-9.506755217687678D-01,
     +   -8.802391537269859D-01,-7.815140038968014D-01,
     +   -6.576711592166909D-01,-5.126905370864771D-01/
      DATA PT(71),PT(72),PT(73),PT(74),PT(75),PT(76),PT(77),PT(78),
     +   PT(79),PT( 80)/-3.512317634538763D-01,-1.784841814958478D-01,
     +   -9.915651684209309D-01,-9.558239495713978D-01,
     +   -8.926024664975557D-01,-8.037049589725230D-01,
     +   -6.916870430603533D-01,-5.597708310739476D-01,                      
     +   -4.11751161462826D-01,-2.518862256915055D-01/
      DATA PT(81),PT(82),PT(83),PT(84),PT(85),PT(86),PT(87),PT(88),
     +   PT(89),PT(90)/-8.477501304173527D-02,-9.924068438435845D-01,
     +   -9.602081521348301D-01,-9.031559036148179D-01,
     +   -8.227146565371427D-01,-7.209661773352294D-01,
     +   -6.005453046616811D-01,-4.645707413759609D-01,
     +   -3.165640999636298D-01,-1.603586456402254D-01/
      DATA PT(91),PT(92),PT(93),PT(94),PT(95),PT(96),PT(97),PT(98),
     +   PT(99),PT(100)/-9.931285991850949D-01,-9.639719272779138D-01,
     +   -9.122344282513259D-01,-8.391169718222189D-01,
     +   -7.463319064601507D-01,-6.360536807265151D-01,
     +   -5.108670019508271D-01,-3.737060887154196D-01,
     +   -2.277858511416451D-01,-7.652652113349732D-02/
C
C  ***************
C  *
C  * WT CONTAINS THE WEIGHTS FOR THE VALUES IN PT; THE SAME WEIGHTS ARE
C  * ALSO USED FOR THE 'MIRROR IMAGE' POINTS.
C  *
C  ***************
C
      DATA WT(1),WT(2),WT(3),WT(4),WT(5),WT(6),WT(7),WT(8),
     +   WT(9),WT(10)/ 1.000000000000000D00,5.555555555555557D-01,
     +   3.478548451374538D-01,6.521451548625462D-01,
     +   2.369268850561891D-01,4.786286704993665D-01,
     +   1.713244923791703D-01,3.607615730481386D-01,
     +   4.679139345726910D-01,1.294849661688697D-01/
      DATA WT(11),WT(12),WT(13),WT(14),WT(15),WT(16),WT(17),WT(18),
     +   WT(19),WT(20)/ 2.79705391489276 7D-01,3.8183005050501189D-01,
     +   1.012285362903762D-01,2.223810344533745D-01,
     +   3.137066458778873D-01,3.626837833783620D-01,
     +   8.127438836157441D-02,1.806481606948574D-01,
     +   2.606106964029355D-01,3.123470770400028D-01/
      DATA WT(21),WT(22),WT(23),WT(24),WT(25),WT(26),WT(27),WT(28),
     +   WT(29),WT(30)/ 6.667134430868812D-02, 1.494513491505806D-01,
     +   2.190863625159820D-01,2.692667193099964D-01,
     +   2.955242247147529D-01,5.566856711617367D-02,
     +   1.255803694649046D-01,1.862902109277342D-01,
     +   2.331937645919905D-01,2.628045445102467D-01/
      DATA WT(31),WT(32),WT(33),WT(34),WT(35),WT(36),WT(37),WT(38),
     +   WT(39),WT(40)/ 4.717533638751183D-02,1.069393259953184D-01,
     +   1.600783285433462D-02,2.031674267230659D-01,
     +   2.334925365383548D-01,2.491470458134028D-01,
     +   4.048400476531588D-02,9.212149983772845D-02,
     +   1.388735102197872D-02,1.781459807619457D-01/
      DATA WT(41),WT(42),WT(43),WT(44),WT(45),WT(46),WT(47),WT(48),
     +   WT(49),WT(50)/ 2.078160475368885D-01,2.262831802628972D-01,
     +   3.511946033175186D-02,8.015808715976020D-02,
     +   1.215185706879032D-02,1.572031671581935D-01,
     +   1.855383974779378D-01,2.051984637212956D-01,
     +   2.152638534631578D-01,3.075324199611727D-02/
      DATA WT(51),WT(52),WT(53),WT(54),WT(55),WT(56),WT(57),WT(58),
     +   WT(59),WT(60)/ 7.036604748810809D-02,1.071592204671719D-01,
     +   1.395706779261543D-01,1.662692058169939D-01,
     +   1.861610000155622D-01,1.984314853271116D-01,
     +   2.715245941175409D-02,6.225352393864789D-02,
     +   9.515851168249277D-02,1.246289712555339D-01/
      DATA WT(61),WT(62),WT(63),WT(64),WT(65),WT(66),WT(67),WT(68),
     +   WT(69),WT(70)/ 1.495959888165767D-01, 1.691565193950025D-01,
     +    1.826034150449236D-01, 1.894506104550685D-01,
     +    2.414830286854793D-02, 5.545952937398720D-02,
     +    8.503614831717917D-02, 1.118838471934040D-01,
     +    1.351363684685255D-01, 1.540457610768103D-01/
      DATA WT(71),WT(72),WT(73),WT(74),WT(75),WT(76),WT(77),WT(78),
     +   WT(79),WT(80)/ 1.680041021564500D-01, 1.765627053669926D-01,
     +    2.161601352648331D-02, 4.971454889496981D-02,
     +    7.642573025488905D-02, 1.009420441062872D-01,
     +    1.225552067114785D-01, 1.406429146706506D-01,
     +    1.546846751262652D-01, 1.642764837458327D-01/
      DATA WT(81),WT(82),WT(83),WT(84),WT(85),WT(86),WT(87),WT(88),
     +   WT(89),WT(90)/ 1.691423829631436D-01, 1.946178822972648D-02,
     +    4.481422676569960D-02, 6.904454273764122D-02,
     +    9.149002162244999D-02, 1.115666455473340D-01,
     +    1.287539625393362D-01, 1.426067021736066D-01,
     +    1.527660420658597D-01, 1.589688433939543D-01/
      DATA WT(91),WT(92),WT(93),WT(94),WT(95),WT(96),WT(97),WT(98),
     +   WT(99),WT(100)/ 1.761400713915212D-02, 4.060142980038694D-02,
     +    6.267204833410905D-02, 8.327674157670474D-02,
     +    1.019301198172404D-01, 1.181945319615184D-01,
     +    1.316886384491766D-01, 1.420961093183820D-01,
     +    1.491729864726037D-01, 1.527533871307258D-01/
C
C  ***************
C  *
C  * WZ CONTAINS ALL THE WEIGHTS FOR THE ORIGIN.
C  *
C  ***************
C
      DATA WZ/ 8.888888888888890D-01, 5.688888888888889D-01,
     +    4.179591836734694D-01, 3.302393550012598D-01,
     +    2.729250867779006D-01, 2.325515532308739D-01,
     +    2.025782419255613D-01, 1.794464703562065D-01,
     +    1.610544498487837D-01/
C
C
C  ***************
C  *
C  * NSTART INDICATES WHERE THE DESIRED POINTS AND WEIGHTS ARE LOCATED
C  * RETRIEVE THE POINTS AND WEIGHES. TRANSFORM THE POINTS TO (A,B).
C  *
C  ***************
C
         NSTART = ((N/2)*((N-1)/2))
         ND2 = N/2
         C = (B - A)/2.0D0
         D = (B + A)/2.0D0
         DT2 = D*2
         DO 15 J=1,ND2
           NJ = N-J+1
           NSTJ  = NSTART + J
           P(J)  = PT(NSTJ)*C +D
           P(NJ) = -P(J) + DT2
           W(J)  = WT(NSTJ)*C
           W(NJ) = W(J)
15        CONTINUE
          IF ((ND2)*2.EQ.N) GO TO 25
            P(ND2 + 1) = D
            W(ND2 + 1) = WZ(ND2)*C
25        CONTINUE
       RETURN
       END
C
C
C
	SUBROUTINE TEMPHIST(TEMP0,FLVEL0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       The program is written to retrieve temperature history data.   
C       The max number of datapoints is 30.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	INTEGER N
	REAL*8 POSI(3000),TMP(3000),VEL(3000),TEMP0,FLVEL0
	COMMON /TMPHIST/ N,POSI,TMP,VEL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       N       Number of datapoints
C       POSI    Position of data points (cm)
C       TMP     Temperature of the position specified
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	OPEN(25,FILE='TMPHIST.DAT',STATUS='OLD')
	READ(25,*)N
C
	DO 100 I = 1,N
	    READ(25,*)POSI(I),TMP(I),VEL(I)
 100  CONTINUE
	TEMP0 = TMP(1)
	FLVEL0 = VEL(1)
	CLOSE(25)
C
	RETURN
	END
C
C
C
	SUBROUTINE TEMPXT(X,T,V)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       This program is written to determine the temperature of the    C
C       specified point.  The temperature between datapoints is linear C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	INTEGER N 
	REAL*8 POSI(3000),TMP(3000),VEL(3000),X,T,V
	COMMON /TMPHIST/ N,POSI,TMP,VEL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       N       Number of datapoints
C       POSI    Position of data points (cm)
C       TMP     Temperature of the position specified
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	I = 0
  10  I = I + 1
	IF ((X.GE.POSI(I)).AND.(X.LT.POSI(I+1))) THEN
	    T = ((POSI(I+1) - X)*TMP(I) + (X - POSI(I))*TMP(I+1))
     &      /(POSI(I+1) - POSI(I))
          V = ((POSI(I+1) - X)*VEL(I) + (X - POSI(I))*VEL(I+1))
     &      /(POSI(I+1) - POSI(I))
	ELSE
	    GOTO 10
	END IF
C
	RETURN
	END
C
C
C
	SUBROUTINE VAPORSAT(TEMP,PSAT,A,B)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	This subroutine is to determine the saturation pressure at the C
C	given temperature.                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	REAL*8 TEMP,PSAT(2),A(2),B(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	TEMP	Temperature (K)
C	PSAT	Saturation pressure (atm)
C       A       Temperature coefficient
C       B       Intersection
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DO 50 I = 1,2
           PSAT(I) = 10**(A(I)/TEMP+B(I))/760.0D0
  50	CONTINUE
C
	RETURN
	END
C
C
C
	SUBROUTINE MTSAT(A,B)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	This subroutine is to read the saturation vapor pressure data. C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	REAL*8 A(2),B(2)
	CHARACTER*12 NAMEMT(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	A	Temperature coefficient
C       B       Intersection
C   	NAMEMT	Name of the metal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	OPEN(51,FILE='MTPSAT.DAT')
	DO 50 I = 1,2
	  READ(51,*) NAMEMT(I)
	  READ(51,*) A(I),B(I)
  50	CONTINUE
        CLOSE(51)
C
	RETURN
	END
C
C
C
	SUBROUTINE RXNRT(TEMP,XK1,A0,EA)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	This subroutine is to determine the reaction rate              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	REAL*8 TEMP,XK1(2),A0(2),EA(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	TEMP	Temperature (K)
C	XK1	Reaction rate coefficient
C       A0      Pre-Exponent term of the reaction rate
C       EA      Activation enegry
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DO 50 I = 1,2
           XK1(I) = A0(I)*EXP(-EA(I)/TEMP)
  50	CONTINUE
C
	RETURN
	END
C
C
C
	SUBROUTINE RTDATA(A0,EA)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	This subroutine is to read the reaction rate data.             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	REAL*8 A0(2),EA(2)
	CHARACTER*12 NAMEMT(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	A0	Pre-Exponent term of the reaction rate
C       EA      Activation enegry
C   	NAMEMT	Name of the metal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	OPEN(52,FILE='MTRXN.DAT')
	DO 50 I = 1,2
	  READ(52,*) NAMEMT(I)
	  READ(52,*) A0(I),EA(I)
  50	CONTINUE
        CLOSE(52)
C
	RETURN
	END
C
C
C
	SUBROUTINE NUCLRT(TEMP,XK2,A1,B1,A2,B2,TN1,TN2,TP1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	This subroutine is to determine the nucleation rate            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	REAL*8 TEMP,XK2(2),A1(2),B1(2),A2(2),B2(2),TN1(2),TN2(2),TP1(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	TEMP	Temperature (K)
C	XK2	Nucleation rate coefficient
C       A1      Temperature coefficient
C       B1      Intersection
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DO 50 I = 1,2
          IF (TP1(I).EQ.1) XK2(I) = 10**(A1(I)/TEMP+B1(I))
          IF (TP1(I).EQ.2) XK2(I) = A1(I)*(TEMP**TN1(I))*EXP(B1(I)/TEMP)
     &                            + A2(I)*(TEMP**TN2(I))*EXP(B2(I)/TEMP)
  50	CONTINUE
C
	RETURN
	END
C
C
C
	SUBROUTINE NUCLDATA(A1,B1,A2,B2,TN1,TN2,TP1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	This subroutine is to read the nucleation rate data.           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	REAL*8 A1(2),B1(2),A2(2),B2(2),TN1(2),TN2(2),TP1(2)
	CHARACTER*12 NAMEMT(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	A1	Temperature coefficient
C       B1      Intersection
C   	NAMEMT	Name of the metal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	OPEN(53,FILE='MTNUCL.DAT')
	DO 50 I = 1,2
	  READ(53,*) NAMEMT(I)
	  READ(53,*) TP1(I)
	  IF (TP1(I).EQ.1) READ(53,*) A1(I),B1(I)
	  IF (TP1(I).EQ.2) READ(53,*) A1(I),B1(I),TN1(I),A2(I),B2(I),
     &                     TN2(I)
  50	CONTINUE
        CLOSE(53)
C
	RETURN
	END

