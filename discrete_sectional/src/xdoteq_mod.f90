module xdoteq_mod
    
    ! precision specification
    use preci, wp => dp

    ! shared variables
    use globals
    use calbeta_common

    implicit none

contains

    subroutine XDOTEQ(XM,DM,COUNT)
!*******************************************************************
!     This subroutine is written to determine the differential equations 
!     required by RK45.  It is divided into 5 parts: First discrete size
!     , other discrete sizes, first section, other sections and the 
!     vapor.  In each part, the interactions between the discrete-
!     discrete, section-discrete, section-section and condensation are
!     considered
!*******************************************************************
    integer :: IND, IND1, IND2, LL, I, J, K
    real (wp) :: COUNT, TTMC, VTEMP
    real (wp), dimension(MAXDISC+MAXDISC+MAXSEC) :: XM, DM
    real (wp), dimension(8) :: TMP
    real (wp), dimension(2,MAXDISC) :: TMC,TMCX
    real (wp), dimension(2,MAXSEC) :: TMCS,TMSCX,TMCXS
    real (wp), dimension(MAXDISC+MAXDISC+MAXSEC,MAXDISC+MAXDISC+MAXSEC) :: XXMM
    real (wp), dimension(2) :: XNUCL
!*******************************************************************
!     BETA(I,J) Coagulation coefficients between particle I and J 
!               (cm3/#/s)
!     FRAC      Fraction
!     KFSB      First section boundary for the first section formation
!               by discrete sizes interaction 
!     M         Total number of discrete size & sections (NDISC + NSEC)
!     NDISC     Total number of discrete sizes
!     NDISC1    NDISC + 1
!     NDISC2    NDISC + 2
!     SIZEFAC   Size increasing factor based on volume
!     TMC(I)    Removal rate of discrete I by condensation on discrete I 
!               (#/cc/s) 
!     TMCS(I)   Removal & formation rate of section I by condensation on 
!               section I(#/cc/s) 
!     TMCX(I)   Total monomer condensation rate on discrete I (#/cc/s)
!     TMCXS(I)  Total monomer condensation rate on section I (#/cc/s)
!     TMP       Temporary variables
!     TMSCS(I)  Formation rate of section I+1 by condensation on section
!               I (#/cc/s)
!     TTMC      Total condensation  (#/cc)
!     VG        Mean volume of the section (cm3)
!     VL        Lower bound aerosol volume of the section (cm3)
!     XK2       Dimer nucleation rate (cm3/#/s)
!     XM        Aerosol and vapor concentration
!       XM(1)                   Vapor molecule concentration (#/cc)
!       XM(2...NDISC)           Discrete size aerosols (#/cc)
!       XM(NDISC1...M)       Sectional size aerosols (#/cc)
!     XNMON     Mean number of molecules in one aerosol in the section
!     XNMON2    2.0*XNMON
!     XNPD(I)   Surface saturation concentration for section I particle 
!               (#/cc)
!     XNUCL     Dimer nucleation rate (#/cc/s)
!     XR1       Aerosol vapor formation rate (#/cc/s)
!*******************************************************************
!
!***********************************************************************
!     Determine the concentration multiplication
!***********************************************************************
    do I = 1,NDISC(1)
        XXMM(I,I) = XM(I)*XM(I)*BDD(1,I,I)*XTMOD
        do J = I+1,NDISC(1)
            XXMM(I,J) = XM(I)*XM(J)*XTMOD
            XXMM(J,I) = XXMM(I,J)*BDD(1,J,I)
            XXMM(I,J) = XXMM(I,J)*BDD(1,I,J)
        end do
        do J = NDISC1,M
            XXMM(I,J) = XM(I)*XM(J)*XTMOD
        end do
    end do
    do I = 1,NDISC(2)
        IND1 = I+NDISC(1)
        XXMM(IND1,IND1) = XM(IND1)*XM(IND1)*BDD(2,I,I)*XTMOD
        do J = I+1,NDISC(2)
            IND2 = J+NDISC(1)
            XXMM(IND1,IND2) = XM(IND1)*XM(IND2)*XTMOD
            XXMM(IND2,IND1) = XXMM(IND1,IND2)*BDD(2,J,I)
            XXMM(IND1,IND2) = XXMM(IND1,IND2)*BDD(2,I,J)
        end do
        do J = NDISC1,M
            XXMM(IND1,J) = XM(IND1)*XM(J)*XTMOD
        end do
    end do
    do I = NDISC1,M
        XXMM(I,I) = XM(I)*XM(I)*XTMOD
        do J = I+1,M
            XXMM(I,J) = XM(I)*XM(J)*XTMOD
            XXMM(J,I) = XXMM(I,J)
        end do
    end do
!***********************************************************************
!     First aerosol size
!     XNUCL: formation
!     TMP(1): removal
!     TMC(2): condensation (removal)
!***********************************************************************
    TMP(1) = 0.0
    XNUCL(1) = BDD(1,1,1)*XM(1)*XM(1)*XTMOD*XK2RATIO(1)
    do J = 2,NDISC(1)
        TMP(1) = TMP(1) + XXMM(2,J)
    end do
    do J = 1,NSEC
        TMP(1) = TMP(1) + BDS4(1,2,J)*XXMM(2,J+NDISCT)
    end do
    if (XM(1) >= XNPD(1,2)) then
        TMP(2) = (XM(1)-XNPD(1,2))*XM(2)*XTMOD
        TMCX(1,2) = BDD(1,1,2)*TMP(2)
        TMC(1,2) = BDD(1,2,1)*TMP(2)
    end if
    DM(2) = XNUCL(1) - TMP(1) - TMC(1,2)

    TMP(1) = 0.0
    IND2 = NDISC(1)+2
    XNUCL(2) = BDD(2,1,1)*XM(NDISCA)*XM(NDISCA)*XTMOD*XK2RATIO(2)
    do J = 2,NDISC(2)
        TMP(1) = TMP(1) + XXMM(IND2,NDISC(1)+J)
    end do
    do J = 1,NSEC
        TMP(1) = TMP(1) + BDS4(2,2,J)*XXMM(IND2,J+NDISCT)
    end do
    if (XM(NDISCA) >= XNPD(2,2)) then
        TMP(2) = (XM(NDISCA)-XNPD(2,2))*XM(IND2)*XTMOD
        TMCX(2,2) = BDD(2,1,2)*TMP(2)
        TMC(2,2) = BDD(2,2,1)*TMP(2)
    end if

    DM(IND2) = XNUCL(2) - TMP(1) - TMC(2,2)
!***********************************************************************
!     ###   Other discrete sizes   ###
!***********************************************************************
    do I = 3,NDISC(1)
        do LL =1,2
            TMP(LL) = 0.0
        end do
        !=======================================================================
        !     Removal due to the collisions with other discrete size or section
        !=====
        do J = 2,NDISC(1)
            TMP(1) = TMP(1) + XXMM(I,J)
        end do
        do J = 1,NSEC
            TMP(1) = TMP(1) + BDS4(1,I,J)*XXMM(I,J+NDISCT)
        end do
        !=======================================================================
        !     Formation from the collisions of smaller discrete sizes
        !=====
        do J = 2,I-2
            TMP(2) = TMP(2) + XXMM(J,I-J)*XMFAC(1,J,I-J)
        end do
        !=======================================================================
        !     Condensation
        !=====
        if (XM(1) >= XNPD(1,I)) then
            TMP(3) = (XM(1)-XNPD(1,I))*XM(I)*XTMOD 
            TMCX(1,I) = BDD(1,1,I)*TMP(3)
            TMC(1,I) = BDD(1,I,1)*TMP(3)
        end if
    
        DM(I) = - TMP(1) + TMP(2) - TMC(1,I) + (TMC(1,I-1)+TMCX(1,I-1)) &
        *XMFAC(1,I-1,1)
    
    end do

    do I = 3,NDISC(2)
        do LL =1,2
            TMP(LL) = 0.0
        end do
        IND = I + NDISC(1)
        !=======================================================================
        !     Removal due to the collisions with other discrete size or section
        !=====
        do J = 2,NDISC(2)
            TMP(1) = TMP(1) + XXMM(IND,NDISC(1)+J)
        end do
        do J = 1,NSEC
            TMP(1) = TMP(1) + BDS4(2,I,J)*XXMM(IND,J+NDISCT)
        end do
        !=======================================================================
        !     Formation from the collisions of smaller discrete sizes
        !=====
        do J = 2,I-2
            TMP(2) = TMP(2) + XXMM(NDISC(1)+J,IND-J)*XMFAC(2,J,I-J)
        end do
        !=======================================================================
        !     Condensation
        !=====
        if (XM(NDISCA) >= XNPD(2,I)) then
            TMP(3) = (XM(NDISCA)-XNPD(2,I))*XM(IND)*XTMOD
            TMCX(2,I) = BDD(2,1,I)*TMP(3)
            TMC(2,I) = BDD(2,I,1)*TMP(3)
        end if
    
        DM(IND) = -TMP(1) + TMP(2) - TMC(2,I) +  &
        (TMC(2,I-1)+TMCX(2,I-1))*XMFAC(2,I-1,1)
    
    end do
!***********************************************************************
!    ###   First section   ###
!***********************************************************************
    do I = 1,6
        TMP(I) = 0.0
    end do
!=======================================================================
!     Formation from collisions of 2 discrete sizes
!=====
    do I = 2,NDISC(1)
        if ((VG2(1,I) >= VL(1)) .AND. (VG2(1,I) <= VL(2)))  &
        TMP(1) = TMP(1) + XXMM(I,I)*C2E
        do J = I+1,NDISC(1)
            VTEMP = VG(1,I) + VG(1,J)
            if ((VTEMP >= VL(1)) .AND. (VTEMP <= VL(2))) &
            TMP(1) = TMP(1) + (XXMM(I,J)+XXMM(J,I))*XMFAC(1,J,I)
        end do
    end do

    do I = 2,NDISC(2)
        IND1 = NDISC(1)+I
        if ((VG2(2,I) >= VL(1)) .AND. (VG2(2,I) <= VL(2))) &
        TMP(6) = TMP(6) + XXMM(IND1,IND1)*C2E
        do J = I+1,NDISC(2)
            VTEMP = VG(2,I) + VG(2,J)
            if ((VTEMP >= VL(1)) .AND. (VTEMP <= VL(2))) then
                IND2 = NDISC(1)+J
                TMP(6)=TMP(6)+(XXMM(IND1,IND2)+XXMM(IND2,IND1))*XMFAC(2,J,I)
            end if
        end do
    end do
!=======================================================================
!     Formation and removal from collisions of one discrete size and 
!     the first section
!=====
    do I = 2,NDISC(1)
        TMP(2) = TMP(2) + BDS25(1,I,1)*XXMM(I,NDISC1)
    end do

    do I = 2,NDISC(2)
        TMP(2) = TMP(2) + BDS25(2,I,1)*XXMM(NDISC(1)+I,NDISC1)
    end do
!=======================================================================
!     Formation and removal of the first section by the collision of 2
!     first section particles
!=====
    TMP(3) = BSS36(1)*XXMM(NDISC1,NDISC1)
!=======================================================================
!     Removal from collisions of the first section with larger sections
!=====
    do I = 2,NSEC
        TMP(4) = TMP(4) + BSS4(1,I)*XXMM(NDISC1,I+NDISCT)
    end do
!=======================================================================
!     Condensation
!=====
    if (XM(1) >= XNPDS(1,1)) then
        TMP(5) = (XM(1)-XNPDS(1,1))*XM(NDISC1)*XTMOD
        TMCXS(1,1) = BDS4(1,1,1)*TMP(5)
        TMCS(1,1) = BDS25(1,1,1)*TMP(5)
        TMSCX(1,1) = BDS1(1,1,1,2)*TMP(5)
    end if
    if (XM(NDISCA) >= XNPDS(2,1)) then
        TMP(5) = (XM(NDISCA)-XNPDS(2,1))*XM(NDISC1)*XTMOD
        TMCXS(2,1) = BDS4(2,1,1)*TMP(5)
        TMCS(2,1) = BDS25(2,1,1)*TMP(5)
        TMSCX(2,1) = BDS1(2,1,1,2)*TMP(5)
    end if

    DM(NDISC1) = TMP(1)-TMP(2)-TMP(3)-TMP(4)+TMP(6)-TMCS(1,1)+ &
    (TMC(1,NDISC(1))+TMCX(1,NDISC(1)))*XMFAC(1,NDISC(1),1)
    DM(NDISC1) = DM(NDISC1)-TMCS(2,1)+(TMC(2,NDISC(2))+ &
    TMCX(2,NDISC(2)))*XMFAC(2,NDISC(2),1)
!***********************************************************************
!     ###   Other sections   ###
!***********************************************************************
    do I = 2,NSEC
        do LL = 1,8
            TMP(LL) = 0.0
        end do
        IND = I+NDISCT
        !=======================================================================
        !     Formation from collision of discrete sizes
        !=====
        do J = 2,NDISC(1)
            if ((VG2(1,J) >= VL(I)) .AND. (VG2(1,J) <= VL(I+1)))  &
            TMP(1) = TMP(1) + XXMM(J,J)*C2E
            do K = J+1,NDISC(1)
                VTEMP = VG(1,J) + VG(1,K)
                if ((VTEMP >= VL(I)) .AND. (VTEMP <= VL(I+1))) &
                TMP(1) = TMP(1) + (XXMM(J,K)+XXMM(K,J))*XMFAC(1,K,J)
            end do
        end do
        do J = 2,NDISC(2)
            IND1 = J + NDISC(1)
            if ((VG2(2,J) >= VL(I)) .AND. (VG2(2,J) <= VL(I+1)))  &
            TMP(1) = TMP(1) + XXMM(IND1,IND1)*C2E
            do K = J+1,NDISC(2)
                VTEMP = VG(2,J) + VG(2,K)
                if ((VTEMP >= VL(I)) .AND. (VTEMP <= VL(I+1))) then
                    IND2 = NDISC(1)+K
                    TMP(1)=TMP(1)+(XXMM(IND1,IND2)+XXMM(IND2,IND1))*XMFAC(2,K,J)
                end if
            end do
        end do
        !=======================================================================
        !     Formation from collisions of one discrete size and one smaller 
        !     section
        !=====
        do J = 1,I-1
            IND1 = J+NDISCT
            if ((VL(I)-VL(J+1)) > VG(1,NDISC(1))) GOTO 512
            do K = 2,NDISC(1)
                TMP(2) = TMP(2) + BDS1(1,K,J,I)*XXMM(K,IND1)
            end do
            512 if ((VL(I)-VL(J+1)) > VG(2,NDISC(2))) GOTO 525
            do K = 2,NDISC(2)
                TMP(2) = TMP(2) + BDS1(2,K,J,I)*XXMM(K+NDISC(1),IND1)
            end do
            !=======================================================================
            !     Formation from collisions of smaller sections
            !     K:  Smaller section
            !=====
            525 if (VL(I) >= (2.0D0*VL(J+1))) cycle
            TMP(3) = TMP(3) + BSS1(I,J,J)*XXMM(IND1,IND1)
            do K = 1,J-1
                TMP(3) = TMP(3) + BSS1(I,K,J)*XXMM(K+NDISCT,IND1)
            end do
        end do
        !=======================================================================
        !     Removal and formation of the specific section by collision with 
        !     discrete size
        !=====
        560 do J = 2,NDISC(1)
            TMP(4) = TMP(4) + BDS25(1,J,I)*XXMM(J,IND)
        end do
        do J = 2,NDISC(2)
            TMP(4) = TMP(4) + BDS25(2,J,I)*XXMM(J+NDISC(1),IND)
        end do
        !=======================================================================
        !     Removal and formation of the specific section by collision with 
        !     a smaller section
        !     FRAC2:  formation
        !=====
        do J = 1,I-1
            TMP(5) = TMP(5) + BSS25(I,J)*XXMM(IND,J+NDISCT)
        end do
        !=======================================================================
        !     Removal and formation of the section due to collision of the 2 
        !     same sections
        !     FRAC2: formation
        !=====
        TMP(6) = TMP(6) + BSS36(I)*XXMM(IND,IND)
        !=======================================================================
        !     Removal of the section due to collsions of the section with a 
        !     larger section
        !=====
        do J = I+1,NSEC
            TMP(7) = TMP(7) + BSS4(I,J)*XXMM(IND,J+NDISCT)
        end do
        !=======================================================================
        !     Condensation
        !=====
        if (XM(1) >= XNPDS(1,I)) then
            TMP(8) = (XM(1)-XNPDS(1,I))*XM(IND)*XTMOD
            TMCXS(1,I) = BDS4(1,1,I)*TMP(8)
            TMCS(1,I) = BDS25(1,1,I)*TMP(8)
            TMSCX(1,I) = BDS1(1,1,I,I+1)*TMP(8)
        end if
        if (XM(NDISCA) >= XNPDS(2,I)) then
            TMP(8) = (XM(NDISCA)-XNPDS(2,I))*XM(IND)*XTMOD
            TMCXS(2,I) = BDS4(2,1,I)*TMP(8)
            TMCS(2,I) = BDS25(2,1,I)*TMP(8)
            TMSCX(2,I) = BDS1(2,1,I,I+1)*TMP(8)
        end if
    
        DM(IND) = TMP(1)+TMP(2)+TMP(3)-(TMP(4)+TMP(5)+TMP(6)+TMP(7)) &
        -TMCS(1,I)+TMSCX(1,I-1)-TMCS(2,I)+TMSCX(2,I-1)
    end do
!*********************************************************** 
!     ###   Removal and generation of the FeO vapor   ###
!***********************************************************
    if (NCONSCND == 1) then 
        DM(1) = 0.0 
        DM(NDISCA) = 0.0
    else
        TTMC = 0.0
        do I = 2,NDISC(1)
            TTMC = TTMC + TMCX(1,I)
        end do
        do I = 1,NSEC
            TTMC = TTMC + TMCXS(1,I)
        end do
        DM(1) = XR1(1)*XM(M+1) - TTMC - XNUCL(1)/C2E
    
        TTMC = 0.0
        do I = 2,NDISC(2)
            TTMC = TTMC + TMCX(2,I)
        end do
        do I = 1,NSEC
            TTMC = TTMC + TMCXS(2,I)
        end do
        DM(NDISCA) = XR1(2)*XM(MAERO) - TTMC - XNUCL(2)/C2E
    end if
!*******************************************
!     ###   Precursor decomposition   ###
!*******************************************
    DM(M+1) = -XK1(1)*XM(M+1) !*CO2
    DM(MAERO) = -XK1(2)*XM(MAERO)*CO2

    return
    end subroutine XDOTEQ


end module xdoteq_mod
