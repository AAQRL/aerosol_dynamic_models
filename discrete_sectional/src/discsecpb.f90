program discsecpb
!***********************************************************************
!     This program is written to determine the aerosol dynamics        *
!     according to Discrete-Sectional Model considering two species in *
!     the system.  The number of discrete size is limited to MACDISC   *
!     for both species, and the total number of section is limited to  *
!     MAXSEC.  The actual number of discrete size is noted by NDISC,   *
!     and the actual number of section, NSEC, is determined by the     *
!     size increase factor, SIZEFACC, and the maximum size considered, *
!     XMAXSIZE.  In this system, discrete size particles of different  *
!     species do not interact (no coagulation), but particles in the   *
!     sections are not distinguishable.  This program is written       *
!     based on different number concentration, as noted by ETA.  The   *
!     system is based on the experiment performed in NIST, 1994.  The  *
!     subroutine XDOTEQ called by RK45 is used to determine the ODE by *
!     the Discrete-Sectional method.                                   *
!                                                                      *
!     The reactions are assumed to be                                  *
!                                                                      *
!     Fe(CO)5  ----->  Fe(g)      instantly                            *
!                                                                      *
!                   k1                                                 *
!     Fe(g) + O2  ----->  FeO(g) + O                                   *
!                                                                      *
!               k2                                                     *
!     FeO(g)  ----->  (FeO)2(s)                                        *
!                                                                      *
!                                                                      *
!     Si2O(CH3)6 -----> 2Si(g)     instantly                           *
!                                                                      *
!                  k1                                                  *
!     Si(g) + O2 -----> SiO2(g)                                        *
!                                                                      *
!                k2                                                    *
!     SiO2(g)  ----->  (SiO2)2(s)                                      *
!                                                                      *
!     Condensation occurs in the system and the Kelvin effect is       *
!     considered assuming no vaporization.  Mass balance is considered *
!     by removing all the colliding particles and then adding the      *
!     fraction that is newly formed.  The collision coefficient for    *
!     sections is determined by the integration over the section.      *
!     Written by Chang-Yu Wu in the University of Cincinnati, Aug,     *
!     1994 under the guidance of Dr. Pratim Biswas.                    *
!***********************************************************************
!     SYSOPF.DAT   Output file for all the system parameters and       *
!                  aerosol species properties.                         *
!     CONCF.DAT    Output file for vapor concentration and M2 wrt time.*
!     OUTPTF.DAT   Output file for particle size distribution and M0,  *
!                  M1, dpg, dpv, sigma.                                *
!***********************************************************************

!***********************************************************************
!     Modules and precision declarations
!***********************************************************************

    ! precision specification
    use preci, wp => dp

    ! shared variables
    use globals

    ! external functions and subroutines
    use aux_data_funcs
    use calbeta_mod
    use ot_mod
    use rk45_init0

    implicit none


    real (wp) :: COUNT, DIST, DT, FLVEL, FLVEL0, TEMP0, TEMPP, &
    TRATIO, PRATIO, S1, S2, S3, S4, XMM2, XMAXSIZE
    integer :: IND, ki, LG, NTSEC, stat
    integer :: im, jm, km ! appended with "m" to avoid name conflicts with shared variables in calbeta_common.f90

!***********************************************************************
!     Timing system begins
!***********************************************************************
    call CPU_TIME(S1)

!***********************************************************************
!     Reading system and calculation variable values
!***********************************************************************
    open(unit=1, iostat=stat, file='SYSOPH.DAT', status='old')
    do im = 1, 2
        read(1,*) NDISC(im)
        write(*,5) NDISC(im)
        read(1,*) XMW(im)
        write(*,8) XMW(im)
        read(1,*) RHO(im)
        write(*,9) RHO(im)
        read(1,*) XK1(im)
        write(*,11) XK1(im)
        read(1,*) XK2(im)
        write(*,12) XK2(im)
        read(1,*) CONCPRE(im)
        write(*,14) CONCPRE(im)
        read(1,*) PSAT0(im)
        write(*,18) PSAT0(im)
        read(1,*) SIGMA(im)
        write(*,19) SIGMA(im)
        read(1,*) NEXIDIS(im)
        write(*,23) NEXIDIS(im)
        do jm = 1, NEXIDIS(im)
            if (im == 1) then
                read(1,*) NORDER, XM(0 + NORDER)
            else
                read(1,*) NORDER, XM(NDISC(im-1) + NORDER)
            end if
        end do
        read(1,*) NCOND(im)
        read(1,*) NCOND1(im)
    end do
    XMINIT = XM(1)
    NDISCT = NDISC(1) + NDISC(2)
    NDISC1 = NDISCT + 1
    NDISCA = NDISC(1) + 1
    read(1,*) NEXISEC
    do im = 1, NEXISEC
        read(1,*) NORDER, XM(NDISCT + NORDER)
    end do
    read(1,*) SIZEFAC
    write(*,6) SIZEFAC
    read(1,*) XMAXSIZE
    write(*,7) XMAXSIZE
    read(1,*) TEMP0
    write(*,10) TEMP0 
    read(1,*) FLVEL0 
    read(1,*) CO2
    write(*,15) CO2
    read(1,*) NTSEC
    do im = 1, NTSEC
        read(1,*) TIME(im)
        write(*,16) TIME(im)
        read(1,*) NN(im)
        write(*,17) NN(im)
    end do
    read(1,*) NSP1
    read(1,*) NSP2
    read(1,*) XNCAL
    write(*,21) XNCAL
    read(1,*) NCAL
    read(1,*) ETA
    write(*,20) ETA
    read(1,*) ALPHA
    write(*,22) ALPHA
    read(1,*) NFNUCL
    read(1,*) NFCOND
    read(1,*) NFCOAG
    read(1,*) NCONSCND
    5 format(1x,"The total number of discrete size = ",i4)
    6 format(1x,"The size increase factor = ",f10.5)
    7 format(1x,"The upper bound of the largest section (A) = ",f10.5)
    8 format(1x,"The aerosol species molecular weight (g) = ",f10.5)
    9 format(1x,"The aerosol species density (g/cm3)",f10.5)
    10 format(1x,"Temperature (K) = ",f10.5)
    11 format(1x,"The aerosol vapor formation rate Constant, K1 = " &
    ,1PE10.4)
    12 format(1x,"The dimer formation rate Constant, K2 = ",1PE10.4)
    14 format(1x,"Initial precursor concentration (mol/cc) = ",1PE10.4)
    15 format(1x,"Initial O2 concentration (mol/cc) = ",1PE10.4)
    16 format(1x,"Time for the stage (s) = ",1PE10.4)
    17 format(1x,"Calculation steps for the stage = ",i5)
    18 format(1x,"Saturation concentration (atm) = ",1PE10.4)
    19 format(1x,"Surface tension (dyne/cm) = ",f10.5)
    20 format(1x,"Concentration index = ",i2)
    21 format(1x,"Number of differentials are = ",f6.1)
    22 format(1x,"The sticking factor is ", f8.3)
    23 format(1x,"The # of discrete sizes of existing particles is",i4)
    read(1,*) RESULT1
    ! print also works: print *, "CONCENTRATION & M2 OUTPUT FILENAME = ", RESULT1
    write(*,*) "CONCENTRATION & M2 OUTPUT FILENAME = ", RESULT1
    read(1,*) RESULT2
    write(*,*) "AEROSOL SIZE DISTRIBUTION DATA FILE = ", RESULT2
    close(1)

!***********************************************************************
!     call functions
!***********************************************************************
    call temphist(TEMP0,FLVEL0)
    call mtsat(A,B)
    call rtdata(A0,EA)
    call nucldata(A1,B1,A2,B2,TN1,TN2,TP1)

!***********************************************************************
!     Set initial values and calculation constants
!***********************************************************************
    TEMP = TEMP0
    NCAL1 = int(XNCAL)
    DLNSF = dlog(SIZEFAC)
    DSQRTSF = dsqrt(SIZEFAC)
    C2E = 2.0D0**(ETA-1)
    V1(1) = XMW(1)/RHO(1)/AVO
    V1(2) = XMW(2)/RHO(2)/AVO
    NSEC = int(dlog10(XMAXSIZE**3.0/dble(NDISC(1)))/dlog10(SIZEFAC))
    M = NDISCT + NSEC
    MAERO = M + 2
    write(*,*)"The total number of discrete sizes and sections is ",M
    do im = 2, M
        DM(im) = 0.0
    end do
    XM(M+1) = CONCPRE(1)
    XM(MAERO) = CONCPRE(2)

!***********************************************************************
!     Set the discrete size and lower bound of the section, and their
!     related variables
!***********************************************************************
!     ###   Discrete sizes ###    
!**********************************
    do im = 1,2
        if (NCOND(im) == 0) PSAT0(im) = 0.0D0
        do jm = 1,NDISC(im)
            VG(im,jm) = V1(im)*dble(jm)
            VG2(im,jm) = 2.0*VG(im,jm)
            VGSQ(im,jm) = VG(im,jm)*VG(im,jm)
            VGETA(im,jm) = VG(im,jm)**ETA
            DLNVG(im,jm) = dlog(VG(im,jm))
            XMFAC(im,jm,jm) = C2E
            do km = 1,jm-1
                XMFAC(im,km,jm) = (VG(im,km)+VG(im,jm))**ETA/(VGETA(im,jm)+VGETA(im,km))
                XMFAC(im,jm,km) = XMFAC(im,km,jm)
            end do
            DP(im,jm) = (VG(im,jm)*6.0/PI)**C13
            XNPD0(im,jm) = PSAT0(im)*dexp(4.0*SIGMA(im)*V1(im)/(XKB*TEMP*DP(im,jm)))/(82.054*TEMP)*AVO*VGETA(im,1)
            if (NCOND(im) == 0) XNPD0(im,jm) = 0.0
            XNPD(im,jm) = XNPD0(im,jm)
        end do
        R1(im) = 0.5*DP(im,1)
    end do
    
!***************************************
!     ###   The first section   ###
!***************************************
    VL(1) = VG(1,NDISC(1)) + 0.5*V1(1)
    DLNVL(1) = dlog(VL(1))

!****************************************
!     ###   All other sections   ###
!****************************************
    do im = 2,NSEC
        IND = im - 1
        VL(im) = VL(IND)*SIZEFAC
        VGS(IND) = VL(IND)*DSQRTSF
        VGSN(IND) = (VL(IND)+VL(im))/2.0d0
        VGSSQ(IND) = VGS(IND)*VGS(IND)
        DLNVL(im) = dlog(VL(im))
        DPS(IND) = (VGS(IND)*6.0/PI)**C13
        XNPDS0(1,IND) = PSAT0(1)*dexp(4.0*SIGMA(1)*V1(1)/(XKB*TEMP*DPS(IND)))/(82.054*TEMP)*AVO*VGETA(1,1)
        if (NCOND(1) == 0) XNPDS0(1,IND) = 0.0
        XNPDS(1,IND) = XNPDS0(1,IND)
        XNPDS0(2,IND) = PSAT0(2)*dexp(4.0*SIGMA(2)*V1(2)/(XKB*TEMP*DPS(IND)))/(82.054*TEMP)*AVO*VGETA(2,1)
        if (NCOND(2) == 0) XNPDS0(2,IND) = 0.0
        XNPDS(2,IND) = XNPDS0(2,IND)
    end do
    im = NSEC + 1
    VL(im) = VL(NSEC)*SIZEFAC
    VGS(NSEC) = VL(NSEC)*DSQRTSF
    VGSN(NSEC) = (VL(NSEC)+VL(im))/2.0d0
    VGSSQ(NSEC) = VGS(NSEC)*VGS(NSEC)
    DLNVL(im) = dlog(VL(im))

!***********************************************************************
!     ###   Determine the upper bound of the discrete sizes that can 
!           form the first section   
!           Define the reaction rate constant 
!***********************************************************************
    do im = 1, 2
        XR1(im) = XK1(im)*AVO*VGETA(im,1)  !*CO2!for debug, by ymm *VGETA(im,1) should be VG
        if ((VG(im,NDISC(im))*2.0D0) >= VL(2)) then
            KFSB(im) = AINT(VL(2)/V1(im)) !modified by ymm. initial is KFSB(im) = AINT(VL(2)/V1(im))
        else
            KFSB(im) = NDISC(im)*2
        end if
    end do

!***********************************************************************
!     ###   Define the statistical variables for the sections   ###
!***********************************************************************
    do im = 1, NSEC
        DVK(im) = VL(im+1) - VL(im)
    end do
    if (ETA == 0) then
        do im = 1,NSEC
            IND = im + 1
            DVLNV(im) = VL(IND)*DLNVL(IND) - VL(im)*DLNVL(im)
            DVK3(im) = 3.0*DVK(im)
            DV3(im) = VL(IND)**3 - VL(im)**3
        end do
    else if (ETA == 1) then
        do im = 1,NSEC
            IND = im + 1
            DVK2(im) = 2.0*DVK(im)
            DV2(im) = VL(IND)*VL(IND) - VL(im)*VL(im)
            DLNVL2(im) = (DLNVL(IND)*DLNVL(IND) - DLNVL(im)*DLNVL(im))/2.0
        end do
    end if
!***********************************************************************
!     ###   Define the beta calculation variables for the discrete size
!***********************************************************************
    LG = max(NDISC(1),NDISC(2))
    do im = 1,LG
        CONVV(im) = 1.0/dble(im)
        CONV3(im) = dble(im)**C13
        CONV2(im) = CONV3(im)*CONV3(im)
        do jm = 1,im-1 !modified by ymm, initial version is jm = 2,im-1
            DDCOEF(jm,im) = dsqrt(CONVV(im)+CONVV(jm))*(CONV3(im)+CONV3(jm))**2
        end do
    end do
!***********************************************************************
!     Set up the parameters for the differential solver 
!***********************************************************************
    call CPU_TIME(S2)

    call CALBETA(SIZEFAC)

    call CPU_TIME(S3)

    INDB = 1
    COUNT = 0.0
    DIST = 0.0
    open(5, file=RESULT1)
    call OT(COUNT,XM)
    do im = 1,NTSEC
        write(*,*) 'Time section:', im
        DT = TIME(im)/dble(NN(im))
        call INIT0(DT)
    !****************************************************
    !     ###   Differential equations solving   ###
    !****************************************************
        do jm = 1,NN(im)
            COUNT = COUNT + DT
        !	  FLVEL = FLVEL0*TEMP/TEMP0
            DIST = DIST + FLVEL*DT
            TEMPP = TEMP
            XK2TMP(1) = XK2(1)
            XK2TMP(2) = XK2(2)
            call TEMPXT(DIST,TEMP,FLVEL)
            call VAPORSAT(TEMP,PSAT,A,B)
            call RXNRT(TEMP,XK1,A0,EA)
            call NUCLRT(TEMP,XK2,A1,B1,A2,B2,TN1,TN2,TP1)
            TRATIO = TEMP/TEMPP
            CO2 = CO2/TRATIO
            XTMOD = sqrt(TEMP/TEMP0)
            kkk = jm
            do km = 1,2
                XK2RATIO(km) = XK2TMP(km)/XK2(km)
                XR1(km) = XK1(km)*AVO*VGETA(km,1) !*CO2
                if (NCOND(km) /= 0) then 
                    PRATIO = PSAT(km)/PSAT0(km)
                    do ki = 1,NDISC(km)
                        XNPD(km,ki) = XNPD0(km,ki)*PRATIO
                    end do
                    do ki = 1,NSEC
                        XNPDS(km,ki) = XNPDS0(km,ki)*PRATIO
                    end do
                end if
            end do
            do km = 1,MAERO
                XM(km) = XM(km)/TRATIO !the expand due to temperature increase
            end do
            if (NCONSCND == 1) XM(1) = XMINIT
            do km = 1,NSEC
                PXM(km) = XM(km+NDISCT)
            end do
            call RK45(XM,DM,COUNT,DT,MAERO)
            if (MOD(jm,NSP1) == 0) call OT(COUNT,XM)
            if (MOD(jm,NSP2) == 0) then
                XMM2 = 0.0D0
                if (ETA == 0) then
                    do km = 2,NDISC(1)
                        XMM2 = XMM2 + XM(km)*VGSQ(1,km)
                    end do
                    do km = 2,NDISC(2)
                        XMM2 = XMM2 + XM(NDISC(1)+km)*VGSQ(2,km)
                    end do
                    do km = 1,NSEC
                        XMM2 = XMM2 + XM(km+NDISCT)*DV3(km)/DVK3(km)
                    end do
                else if (ETA == 1) then
                    do km = 2,NDISC(1)
                        XMM2 = XMM2 + XM(km)*VG(1,km)
                    end do
                    do km = 2,NDISC(2)
                        XMM2 = XMM2 + XM(NDISC(1)+km)*VG(2,km)
                    end do
                    do km = 1,NSEC
                        XMM2 = XMM2 + XM(km+NDISCT)*DV2(km)/DVK2(km)
                    end do
                else
                    do km = 2,NDISC(1)
                        XMM2 = XMM2 + XM(km)
                    end do
                    do km = NDISC(1)+2,M
                        XMM2 = XMM2 + XM(km)
                    end do
                end if
                write(5,105)COUNT,DIST,XM(1)/VGETA(1,1),XM(NDISCA)/VGETA(2,1),XMM2,XM(M+1)
            end if
            do km = 1,NSEC
                RHOSEC(km) = (PXM(km)*RHOSEC(km)+(XM(km)-PXM(km))*DRHO(km))/XM(km)
            end do
            write(*,*) jm,'of',NN(im)
        end do
    end do
    105 format(1X,1P6E12.4)

    call CPU_TIME(S4)

    write(5,123) S1,S2,S3,S4
    123 format(1X,1P4E12.3)
    close(5)

    106 stop

end program discsecpb
