module calbeta_mod
    
    ! precision specification
    use preci, wp => dp

    ! shared variables
    use globals
    use calbeta_common

    ! external functions and subroutines
    use funcs
    use quad_glegq
    use dtwodq_mod

    implicit none

contains

    subroutine CALBETA(SIZEFAC)
!*******************************************************************
!     This subroutine is written to determine the coagulation and
!     condensation coefficients in the free molecular regime.  The
!     nucleation rate coefficient (k2) following the reaction form is 
!     also determined.
!*******************************************************************

    implicit none

    real (wp) :: SIZEFAC, BDS2, BDS5, BT, BTEMP, DDV, DDV1, DV, EABS, &
    ERE, ERRE, RES, RES2, RES3, RES5, RES6, XL, XLB, XR,  &
    XUB, XVL12, XVU12 
    integer :: IR, LL, IND
    real (wp), dimension(2) :: CONST0
    real (wp) :: C8,C6,TMP
    integer :: NSTART(2)
    !real (wp) :: FX,F,G,H
    !external GLEGQ,FX,F,G,H
!*******************************************************************
!     BETA(I,J) Coagulation coefficients between particle I and J 
!               (cm3/#/s)
!     RHO       Aerosol density (g/cm3)
!     R1        Monomer radius (cm)
!     TEMP      Temperature (K)
!     V1        Monomer volume (cm3)
!     XKB       Boltzman's Constant, 1.38D-16 erg/K
!     XK2       Dimer nucleation rate (cm3/#/s)
!*******************************************************************
    TMP = 6.0*XKB*TEMP
    do II = 1,2
        CONST0(II) = dsqrt(TMP*R1(II)/RHO(II))
        if (NCOND1(II) == 0) then
            NSTART(II) = 1
        else
            NSTART(II) = 2
        end if
    end do
    EABS = 1.D-6
    ERE = 0.0001
    IR = 3
    C8 = dsqrt(32.0D0)
    C6 = 1.0/6.0
!*********************************
!     ###   Condensation   ###
!*********************************
    if (NFCOND == 0) GO TO 5
    write(*,*)'Condensation'
    do II = 1,2
        if (NCOND1(II) == 0) cycle
        do I = 2,NDISC(II)
            BT = CONST0(II)*CONV2(I)
            BDD(II,1,I) = BT/VGETA(II,I)
            BDD(II,I,1) = BT/VGETA(II,1)
        end do
        do I = 1,NSEC
            IND = I + 1
            BDS2 = 0.0
            BTEMP = ALPHA*CONST0(II)/DVK(I)
            DDV = DVK(I)/XNCAL
            XL = VL(I)
            XR = XL + DDV
            do LL = 1,NCAL1
                IDD = 6
                BDS2 = BDS2 + QUAD(GLEGQ,NCAL,FX,XL,XR)
                IDD = 7
                BDS4(II,1,I) = BDS4(II,1,I) + QUAD(GLEGQ,NCAL,FX,XL,XR)
                XL = XR
                XR = XR +DDV
            end do
            BDS4(II,1,I) = BDS4(II,1,I)*BTEMP
            DV = DVK(I) - V1(II)
            BDS5 = 0.0
            IDD = 9
            DDV1 = DV/XNCAL
            XL = VL(I)
            XR = XL + DDV1
            do LL = 1,NCAL1
                BDS5 = BDS5 + QUAD(GLEGQ,NCAL,FX,XL,XR)
                XL = XR
                XR = XR + DDV1
            end do
            BDS25(II,1,I) = (BDS2/VGETA(II,1) - BDS5)*BTEMP
            DDV = V1(II)/XNCAL
            XL = VL(IND) - V1(II)
            XR = XL + DDV
            do LL = 1,NCAL1
                BDS1(II,1,I,IND) = BDS1(II,1,I,IND)+QUAD(GLEGQ,NCAL,FX,XL,XR) !for debug, by ymm, not IND
                XL = XR
                XR = XR + DDV
            end do
            BDS1(II,1,I,IND) = BDS1(II,1,I,IND)*BTEMP
        end do
    end do
!*******************************
!     ###   Coagulation   ###
!*******************************
!====================================
!     Discrete-Discrete interaction
!=====
    5 do II = 1,2
        if (NFCOAG == 0) GO TO 614
        write(*,*)'BDD'
        do I = NSTART(II),NDISC(II)
            do J = I+1,NDISC(II)
                BDD(II,I,J) = DDCOEF(I,J)*CONST0(II)/VGETA(II,J)
                BDD(II,J,I) = DDCOEF(I,J)*CONST0(II)/VGETA(II,I)
            end do
            BDD(II,I,I) = CONST0(II)*C8*dble(I)**C6/VGETA(II,I)
        end do
    !*******************************
    !     ###   Nucleation   ###
    !*******************************
        614 if (NFNUCL == 0) GO TO 21
        if (NFNUCL == 2) GO TO 615
        BT = CONST0(II)*dsqrt(CONVV(1)+CONVV(1))*(CONV3(1)+CONV3(1))**2 !CONV2(1)
        BDD(II,1,1) = BT/VGETA(II,1)
        GO TO 21
        615 BDD(II,1,1) = XK2(II)/AVO/VGETA(II,1)
    !=================================================
    !    Discrete-Section : formation of new section
    !====
        21 if (NFCOAG == 0) GO TO 95
        write(*,*)'BDS1'
        IDD = 1
        do I = 2,NSEC
            do J = 1,I-1
                if ((VL(I)-VL(J+1)) > VG(II,NDISC(II))) cycle !the differance is greater than the biggest discrete
                BTEMP = CONST0(II)/DVK(J)
                do K = NSTART(II),NDISC(II)
                    XLB = VG(II,K) + VL(J)
                    XUB = VG(II,K) + VL(J+1)
                    if ((XLB > VL(I+1)) .OR. (XUB < VL(I))) cycle
                    if (XLB <= VL(I)) then  ! Here asuming that dvk(j+1)>dvk(j), namely, sizefc>1
                        DV = XUB - VL(I) 
                        DDV = DV/XNCAL
                        XL = VL(I) - VG(II,K)
                        XR = XL + DDV
                        do LL = 1,NCAL1
                            BDS1(II,K,J,I)=BDS1(II,K,J,I)+QUAD(GLEGQ,NCAL,FX,XL,XR)
                            XL = XR
                            XR = XR + DDV
                        end do
                    else
                        if (XUB <= VL(I+1)) then
                            DDV = DVK(J)/XNCAL
                            XL = VL(J)
                            XR = XL + DDV
                            do LL = 1,NCAL1
                                BDS1(II,K,J,I)=BDS1(II,K,J,I)+QUAD(GLEGQ,NCAL,FX,XL,XR)
                                XL = XR
                                XR = XR + DDV
                            end do
                        else
                            DV = VL(I+1) - XLB
                            DDV = DV/XNCAL
                            XL = VL(J)
                            XR = XL + DDV
                            do LL = 1,NCAL1
                                BDS1(II,K,J,I)=BDS1(II,K,J,I)+QUAD(GLEGQ,NCAL,FX,XL,XR)
                                XL = XR
                                XR = XR + DDV
                            end do
                        end if
                    end if
                    BDS1(II,K,J,I) = BDS1(II,K,J,I)*BTEMP
                end do
            end do
        end do
    !=====================================================================
    !    Discrete-Section : removal & formation of the discrete size or 
    !    the section
    !=====
        write(*,*)'BDS2,BDS4,BDS5'
        do I = 1,NSEC
            BTEMP = CONST0(II)/DVK(I)
            DDV = DVK(I)/XNCAL
            do K = NSTART(II),NDISC(II)
                XL = VL(I)
                XR = XL + DDV
                BDS2 = 0.0d0 !for debug, by ymm, should be deleted, otherwise the condensation will be canceled
                do LL = 1,NCAL1
                    IDD = 2
                    BDS2 = BDS2 + QUAD(GLEGQ,NCAL,FX,XL,XR)
                    IDD = 4
                    BDS4(II,K,I) = BDS4(II,K,I) + QUAD(GLEGQ,NCAL,FX,XL,XR)
                    XL = XR
                    XR = XR + DDV
                end do
                BDS4(II,K,I) = BDS4(II,K,I)*BTEMP
                BDS5 = 0.0D0
                if (VG(II,K) <= DVK(I)) then 
                    DV = DVK(I) - VG(II,K)
                    IDD = 1
                    DDV1 = DV/XNCAL
                    XL = VL(I)
                    XR = XL + DDV1
                    do LL = 1,NCAL1
                        BDS5 = BDS5 + QUAD(GLEGQ,NCAL,FX,XL,XR)
                        XL = XR
                        XR = XR + DDV1
                    end do
                end if
                BDS25(II,K,I) = (BDS2/VGETA(II,K) - BDS5)*BTEMP
            end do
        end do
    95 end do
!=======================================================================
!     2 Same Sections : removal & formation 
!=====
    if (NFCOAG == 0) GO TO 621
    write(*,*)'BSS3,BSS6'
    do J = 1,NSEC
        BTEMP = CONST0(1)/DVK(J)/DVK(J)
        IDD = 13
        call DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES3,ERRE)
        if (SIZEFAC <= 2.0) then
            RES6 = 0.0
        else
            IDD = 16
            call DTWODQ(F,VL(J),DVK(J),G,H,EABS,ERE,IR,RES6,ERRE)
        end if
        BSS36(J) = BTEMP*(RES3-RES6)*0.5 !for debug, by ymm, RES3+RES6? after collision, two become one
    end do
!=======================================================================
!     Section-Section: Removal and formation of the section by collision
!     with a larger or a smaller section 
!=====
    write(*,*)'BSS2,BSS4,BSS5'
    do I = 2,NSEC
        do J = 1,I-1
            BTEMP = CONST0(1)/DVK(I)/DVK(J)
            IDD = 12
            call DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES2,ERRE)
            IDD = 14
            call DTWODQ(F,VL(I),VL(I+1),G,H,EABS,ERE,IR,RES,ERRE)
            BSS4(J,I) = RES*BTEMP
            if (DVK(I) <= VL(J)) then
                RES5 = 0.0
            else
                if (DVK(I) <= VL(J+1)) then
                    IDD = 15
                    call DTWODQ(F,VL(J),DVK(I),G,H,EABS,ERE,IR,RES5,ERRE)
                else
                    IDD = 17
                    call DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES5,ERRE)
                end if
            end if
            BSS25(I,J) = BTEMP*(RES2 - RES5)
        end do
    end do
!=======================================================================
!     Section-Section: Formation of new section by collisions of smaller
!     sections
!     K: Larger section
!=====
    write(*,*)'BSS1'
    do I = 2,NSEC
        do J = 1,I-1
            do K = J,I-1
                if ((VL(I)-VL(K+1)) > VL(K+1)) cycle
                XVL12 = VL(K)+VL(J)
                XVU12 = VL(K+1)+VL(J+1)
                if ((VL(I+1) < XVL12) .OR. (VL(I) > XVU12)) cycle
                BTEMP = CONST0(1)/DVK(J)/DVK(K)
                if (J == K) then
                    if (VL(I) <= XVL12) then
                        if (VL(I+1) <= (VL(J)+VL(K+1))) then
                            IDD = 21
                            call DTWODQ(F,VL(K),VL(I+1)-VL(J),G,H,EABS,ERE,IR, &
                            RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                        else
                            IDD = 22
                            call DTWODQ(F,VL(K),VL(K+1),G,H,EABS,ERE,IR,RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                            IDD = 23
                            call DTWODQ(F,VL(I+1)-VL(K+1),VL(J+1),G,H,EABS,ERE, &
                            IR,RES,ERRE)
                            BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                        end if
                    else
                        if (VL(I) <= (VL(J)+VL(K+1))) then
                            IDD = 24
                            call DTWODQ(F,VL(K),VL(K+1),G,H,EABS,ERE,IR,RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                            IDD = 25
                            call DTWODQ(F,VL(J),VL(I)-VL(K),G,H,EABS,ERE,IR, &
                            RES,ERRE)
                            BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                        else
                            IDD = 26
                            call DTWODQ(F,VL(I)-VL(J+1),VL(K+1),G,H,EABS,ERE,IR, &
                            RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                        end if
                    end if
                    BSS1(I,J,K) = BSS1(I,J,K)/2.0
                else
                    if (VL(I) <= XVL12) then
                        if (VL(I+1) <= (VL(J+1)+VL(K))) then
                            IDD = 31 
                            call DTWODQ(F,VL(J),VL(I+1)-VL(K),G,H,EABS,ERE,IR, &
                            RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                        elseif (VL(I+1) <= (VL(J)+VL(K+1))) then
                            IDD = 32
                            call DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                        else
                            IDD = 33
                            call DTWODQ(F,VL(J),VL(I+1)-VL(K+1),G,H,EABS,ERE,IR, &
                            RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                            IDD = 34
                            call DTWODQ(F,VL(I+1)-VL(K+1),VL(J+1),G,H,EABS,ERE, &
                            IR,RES,ERRE)
                            BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                        end if
                    else
                        if (VL(I) <= (VL(J+1)+VL(K))) then
                            IDD = 35
                            call DTWODQ(F,VL(I)-VL(K),VL(J+1),G,H,EABS,ERE,IR &
                            ,RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                            IDD = 36
                            call DTWODQ(F,VL(J),VL(I)-VL(K),G,H,EABS,ERE,IR, &
                            RES,ERRE)
                            BSS1(I,J,K) = BSS1(I,J,K) + BTEMP*RES
                        elseif (VL(I) <= (VL(J)+VL(K+1))) then
                            IDD = 37
                            call DTWODQ(F,VL(J),VL(J+1),G,H,EABS,ERE,IR,RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                        else
                            IDD = 38
                            call DTWODQ(F,VL(I)-VL(K+1),VL(J+1),G,H,EABS,ERE,IR, &
                            RES,ERRE)
                            BSS1(I,J,K) = BTEMP*RES
                        end if
                    end if
                end if
            end do
        end do
    end do
!
    621 continue
    return
    end subroutine CALBETA


end module calbeta_mod
