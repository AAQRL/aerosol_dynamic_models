module globals

    ! global constants and variables shared between the main program (discsecpb.f90) and others

    ! precision specification
    use preci, wp => dp

    implicit none

    !***********************************************************************
    !     Global definitions
    !***********************************************************************
    !     AVO       Avogadro's number (6.023D23)
    !     BDD       Coagulation coefficient of discrete-discrete
    !     BDS1      Collision coefficient of discrete with section to form 
    !               a new section
    !     BDS2      Collision coefficient of discrete with section to remove
    !               the section
    !     BDS4      Collision coefficient of discrete with section to remove
    !               the discrete
    !     BDS5      Collision coefficient of discrete with section to form
    !               the section
    !     BSS1      Collision coefficient of 2 smaller sections to form a 
    !               new section
    !     BSS2      Collision coefficient of the section with a smaller
    !               section to remove the section
    !     BSS3      Collision coefficient of 2 same sections to remove the
    !               section
    !     BSS4      Collision coefficient of the section with a larger 
    !               section to remove the section
    !     BSS5      Collision coefficient of the section with a smaller 
    !               section to form the section
    !     BSS6      Collision coefficient of 2 same sections to form
    !               the section
    !     CONCPRE   Initial precusor concentration (mol/cc)
    !     CONCSP    Precursor concentration wrt time (mol/cc)
    !     CONVV(I)  1.0/XNMON(I)
    !     CONV3(I)  XNMON(I)**(1.0/3.0)
    !     COUNT     System time variable (s)
    !     CO2       Oxygen concentration (mol/cc)
    !     C13       1.0/3.0
    !     DIST      Axis distance (cm)
    !     DLNSF 	Ln(SIZEFAC)
    !     DLNVG	Ln(VG)
    !     DM        dM/dt, differential equations
    !     DP        Mean particle diameter (cm)
    !     DT        Differential time (s)
    !     DVK(I)    VL(I+1) - VL(I)
    !     DVK2(I)   2.0*DVK(I)
    !     DVK3(I)   3.0*DVK(I)
    !     DVLNV(I)  VL(I+1)*ln(VL(I+1)) - VL(I)*ln(VL(I))
    !     DV2       VL(I+1)**2 - vl(I)**2
    !     DV3(I)    VL(I+1)**3 - VL(I)**3
    !     ETA       Concentration index: 0-number, 1-volume, 2-colume square
    !     FLVEL     Flow velocity wrt time (cm/s)
    !     FLVEL0    Initial flow velocity (cm/s)
    !     INDB	Index for beta calculation (0-First time full calculation; 
    !               1-time modification)
    !     KFSB      First section boundary for the first section formation
    !               by discrete sizes interaction 
    !     M         Total number of ODE's (NDISC + NSEC)
    !     MAXDISC   Maximum number of discrete sizes
    !     MAXSEC    Maximum number of sections
    !     NCAL1     XNCAL - 1
    !     NCOND     Index for complete condensation (0-No barrier; 1-Kelvin
    !			effect)
    !     NCOND1    Index for coagulation of the monomer (0-Monomer 
    !			coagulation; 1-Vapor condensation)
    !     NCONSCND  Index for Constant bulk vapor concentration for 
    !			condensation (0-vapor consumption; 1-Constant 
    !			concentration)
    !     NDISC     Total number of discrete sizes
    !     NDISC1    NDISC + 1
    !     NDISC2    NDISC + 2
    !     NEXIDIS(I)Number of discrete sizes of existing particles for 
    !			Species I
    !     NEXISEC   Number of sections of existing particles
    !     NFNUCL    Flag for nucleation (0-no,1-by stable monomer 
    !			coagulation, 2-given kinetic rate)
    !     NFCOND    Flag for condensation (0-no,1-yes)
    !     NFCOAG    Flag for coagulation (0-no,1-yes)
    !     NN(I)     Calculation steps for stage I
    !     NORDER    The number of the order of existing particles
    !     NSEC      Total number of sections
    !     PI        Circumference ratio (3.1415926)
    !     PSAT      Saturation pressure of FeO (atm)
    !     RESULT1   Concentration & M2 output filename
    !     RESULT2   Aerosol size distribution data filename
    !     RHO       Aerosol density (g/cm3)
    !     SIGMA     Aerosol surface tension (dyne/cm)
    !     SIZEFAC   Size increasing factor based on volume
    !     SQRTSF    Square root of the size factor
    !     TEMP      Temperature (K)
    !     TIME(I)   Time for sytem stage I
    !     VG        Discrete volume (cm3)
    !     VGS       Mean volume of the section
    !     VGSN      Mean volume of the section for n-based model (statistics)
    !     VGSQ      VG**2
    !     VG2       2.0*VG
    !     VL        Lower bound aerosol volume of the section (cm3)
    !     V1        Monomer volume (cm3)
    !     XKB       Boltzman's Constant, 1.38D-16 erg/K
    !     XK1       Aerosol vapor formation rate Constant (cm3/mol/s)
    !     XK2       Dimer nucleation rate (cm3/#/s)
    !     XM        Aerosol and vapor concentration
    !       XM(1)                  Vapor molecule concentration of species 1 
    !                              (#/cc)
    !       XM(2...NDISC(1))       Discrete size aerosols of species 1 (#/cc)
    !       XM(NDISC(1)+1)	       Vapor molecule concentration of species 2 
    !                              (#/cc)
    !       XM(NDISC(1)+2..NDISCT) Discrete size aerosols of species 1 (#/cc)
    !       XM(NDISC1...M)         Sectional size aerosols (#/cc)
    !     XMAXSIZE  Upper bound of the largest section (A)
    !     XMINIT    Initial FeO concentration for Constant condensation
    !     XMW       Molecular weight (g/mol)
    !     XNCAL     Total dissection number for the coagulation coefficient
    !     XNMON     Mean number of molecules in one aerosol in the section
    !     XNMON2    2.0*XNMON
    !     XNPD(I)   Surface saturation concentration for section I particle 
    !               (#/cc)
    !     XR1       Aerosol vapor formation rate (#/mol/s)
    !***********************************************************************

    !***********************************************************************
    !     Constants declaration
    !***********************************************************************
    integer, parameter :: MAXDISC = 300, MAXSEC = 300
    real (wp), parameter :: AVO = 6.023D23 ! Avogadro constant
    real (wp), parameter :: PI = 3.14159265358979323846
    real (wp), parameter :: XKB = 1.38D-16
    real (wp), parameter :: C13 = 1.0/3.0
    real (wp), parameter :: C12 = 2.0/3.0

    !***********************************************************************
    !     Variables declaration
    !***********************************************************************
    real (wp), dimension(MAXDISC+MAXDISC+MAXSEC) :: XM, DM
    real (wp), dimension(2,MAXDISC) :: VG, VG2, DLNVG, VGETA, VGSQ, DP, XNPD, XNPD0
    real (wp), dimension(MAXDISC,MAXDISC) :: DDCOEF
    real (wp), dimension(2,MAXDISC,MAXDISC) :: XMFAC
    real (wp), dimension(MAXSEC) :: VL,VGS,VGSSQ,VGSN,DVLNV,DVK,DLNVL, &
    DLNVL2,DV3,DVK3,DV2,DVK2,DPS,RHOSEC,DRHO,PXM
    real (wp), dimension(MAXDISC) :: CONVV,CONV3,CONV2
    real (wp), dimension(2,MAXSEC) :: XNPDS,XNPDS0
    real (wp), dimension(2) :: V1,XR1,RHO,XMW,XK1,XK2,CONCPRE,PSAT, &
    SIGMA,CONCSP,R1,PSAT0,A,B,A0,A1,EA,B1,TN1,TP1,A2,B2,XK2TMP,XK2RATIO,TN2
    real (wp), dimension(20) :: TIME
    real (wp) :: XMINIT,DSQRTSF
    integer, dimension(2) :: NDISC,KFSB,NEXIDIS,NCOND,NCOND1
    integer, dimension(20) :: NN
    integer :: M,NSEC,NDISCT,NDISC1,ETA,NCAL,NCAL1,MAERO,NEXISEC, &
    NORDER,NFNUCL,NFCOND,NFCOAG,NSP1,NSP2,INDB,NCONSCND
    character(len=12) :: RESULT1,RESULT2
    integer :: NDISCA
    real (wp) :: SIZEFAC
    real (wp) :: CO2
    real (wp) :: TEMP
    real (wp) :: XNCAL
    real (wp) :: DLNSF
    real (wp) :: C2E
    real (wp) :: XTMOD
    real (wp) :: ALPHA
    real (wp) :: kkk

end module globals
