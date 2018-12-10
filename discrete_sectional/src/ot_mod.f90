module ot_mod
    
    ! precision specification
    use preci, wp => dp

    ! shared variables
    use globals

    implicit none

contains

    subroutine OT(COUNT,XM)
!*******************************************************************
!     This program is written to write the results to the output file
!*******************************************************************
    real (wp), dimension(MAXDISC+MAXDISC+MAXSEC) :: XN,XNK,XM
    real (wp) :: TTVOL,TTVOL2,XNN2,COUNT,TVG,TXLN2SIG,VM,XLN2SIG,XLNVG, &
    XMM1, XNN
    integer :: IND, I
!*******************************************************************
!     DVK(I)    VL(I+1) - VL(I)
!     DVLNV(I)  VL(I+1)*ln(VL(I+1)) - VL(I)*ln(VL(I))
!     TTVOL     Total volume of species (mole/cc)
!     TVG       Geometric mean volume (cm3)
!     XLNVG     Mean logarithm volume (cm3)
!     XLN2SIG   Mean logarithm standard deviation
!     XM        Aerosol concentration (#/cc)
!     XMM1      Total aerosol volume concentration (cm3/cc)
!     XN        Aerosol size distribution function (#/cc/cm)
!     XNN       Total number concentration (#/cc)     
!*******************************************************************
    if (ETA == 0) then
        XLNVG = 0.0
        XMM1 = 0.0
        XNN = 0.0
        do I = 2,NDISC(1)
            XN(I) = XM(I)/VG(1,1)
            XMM1 = XMM1 + VG(1,I)*XM(I)
            XNN = XNN + XM(I)
            XLNVG = XLNVG + DLNVG(1,I)*XM(I)
        end do
        do I = 2,NDISC(2)
            IND = I+NDISC(1)
            XN(IND) = XM(IND)/VG(2,1)
            XMM1 = XMM1 + VG(2,I)*XM(IND)
            XNN = XNN + XM(IND)
            XLNVG = XLNVG + DLNVG(2,I)*XM(IND)
        end do
        do I = 1,NSEC
            IND = I+NDISCT
            XN(IND) = XM(IND)/DVK(I)
            XMM1 = XMM1 + VGSN(I)*XM(IND)
            XNN = XNN + XM(IND)
            XLNVG = XLNVG + XM(IND)*(DVLNV(I)/DVK(I)-1.0)
        end do
    
        XLNVG = XLNVG/XNN         
        TVG = dexp(XLNVG)
        VM = XMM1/XNN
        TTVOL = XMM1 + XM(1)*VG(1,1) + XM(1+NDISC(1))*VG(2,1) + (CONCSP(1)*VG(1,1)+CONCSP(2)*VG(2,1))*AVO
        TTVOL2 = XMM1 + XM(1)*VG(1,1) + XM(1+NDISC(1))*VG(2,1)
        XNN2 = XNN + XM(1) + XM(1+NDISC(1))
        XN(1) = XM(1)/VG(1,1)
    
        XLN2SIG = 0.0
        do I = 2,NDISC(1)
            XLN2SIG = XLN2SIG + (DLNVG(1,I)-XLNVG)**2*XM(I)
        end do
        do I = 2,NDISC(2)
            XLN2SIG = XLN2SIG + (DLNVG(2,I)-XLNVG)**2*XM(I+NDISC(1))
        end do
        do I = 1,NSEC
            XLN2SIG = XLN2SIG + XM(I+NDISCT)*((VL(I+1)*(DLNVL(I+1)-XLNVG)**2 &
            - VL(I)*(DLNVL(I)-XLNVG)**2 - 2.0*DVLNV(I))/DVK(I)+2.0D0*(1.0D0+XLNVG))
        end do
        XLN2SIG = XLN2SIG/XNN/9.0
        TXLN2SIG = dexp(dsqrt(XLN2SIG))
    
    elseif (ETA == 1) then
        XNN = 0.0
        XLNVG = 0.0
        XMM1 = 0.0
        do I = 2,NDISC(1)
            XNK(I) = XM(I)/VG(1,I)
            XN(I) = XNK(I)/VG(1,1)
            XNN = XNN + XNK(I)
            XMM1 = XMM1 + XM(I)
            XLNVG = XLNVG + DLNVG(1,I)*XNK(I)
        end do
        do I = 2,NDISC(2)
            IND = I+NDISC(1)
            XNK(IND) = XM(IND)/VG(2,I)
            XN(IND) = XNK(IND)/VG(2,1)
            XNN = XNN + XNK(IND)
            XMM1 = XMM1 + XM(IND)
            XLNVG = XLNVG + DLNVG(2,I)*XNK(IND)
        end do
        do I = 1,NSEC
            IND = I+NDISCT
            XNK(IND) = XM(IND)/DVK(I)
            XN(IND) = XNK(IND)/VGS(I)
            XNN = XNN + XNK(IND)*DLNSF
            XMM1 = XMM1 + XM(IND)
            XLNVG = XLNVG + XNK(IND)*DLNVL2(I)
        end do
    
        XLNVG = XLNVG/XNN         
        TVG = dexp(XLNVG)
        VM = XMM1/XNN
        TTVOL = XMM1 + XM(1) + XM(1+NDISC(1)) + (CONCSP(1)*VG(1,1) + CONCSP(2)*VG(2,1))*AVO
        TTVOL2 = XMM1 + XM(1) + XM(1+NDISC(1))
        XNN2 = XNN + XM(1)/VG(1,1) + XM(1+NDISC(1))/VG(2,1)
        XNK(1) = XM(1)/VG(1,1)
        XN(1) = XNK(1)/VG(1,1)
    
        XLN2SIG = 0.0
        do I = 2,NDISC(1)
            XLN2SIG = XLN2SIG + (DLNVG(1,I)-XLNVG)**2*XNK(I)
        end do
        do I = 2,NDISC(2)
            XLN2SIG=XLN2SIG+(DLNVG(2,I)-XLNVG)**2*XNK(I+NDISC(1))
        end do
        do I = 1,NSEC
            XLN2SIG = XLN2SIG + XM(I+NDISCT)/DVK(I)*((DLNVL(I+1)-XLNVG)**3 - (DLNVL(I)-XLNVG)**3)/3.0D0
        end do
        XLN2SIG = XLN2SIG/XNN/9.0
        TXLN2SIG = dexp(dsqrt(XLN2SIG))
    
    else
        XNN = 0.0
        XLNVG = 0.0
        XMM1 = 0.0
        do I = 2,NDISC(1)
            XNK(I) = XM(I)/VGSQ(1,I)
            XN(I) = XNK(I)/VG(1,1)
            XNN = XNN + XNK(I)
            XMM1 = XMM1 + XM(I)/VG(1,I)
            XLNVG = XLNVG + DLNVG(1,I)*XNK(I)
        end do
        do I = 2,NDISC(2)
            IND = I+NDISC(1)
            XNK(IND) = XM(IND)/VGSQ(2,I)
            XN(IND) = XNK(IND)/VG(2,1)
            XNN = XNN + XNK(IND)
            XMM1 = XMM1 + XM(IND)/VG(2,I)
            XLNVG = XLNVG + DLNVG(2,I)*XNK(IND)
        end do
        do I = 1,NSEC
            IND = I+NDISCT
            XNK(IND) = XM(IND)/VL(I+1)/VL(I)
            XN(IND) = XM(IND)/VGSSQ(I)/DVK(I)
            XNN = XNN + XNK(IND)
            XMM1 = XMM1 + XM(IND)*DLNSF/DVK(I)
            XLNVG = XLNVG + XM(IND)/DVK(I)*((1.0+DLNVL(I))/VL(I)-(1.0+DLNVL(I+1))/VL(I+1))
        end do
   
        XLNVG = XLNVG/XNN         
        TVG = dexp(XLNVG)
        VM = XMM1/XNN
        TTVOL = XMM1 + XM(1)/VG(1,1) + XM(NDISCA)/VG(2,1) + (CONCSP(1)*VG(1,1)+CONCSP(2)*VG(2,1))*AVO
        TTVOL2 = XMM1 + XM(1)/VG(1,1) + XM(NDISCA)/VG(2,1)
        XNN2 = XNN + XM(1)/VG(1,1)**2 + XM(1+NDISC(1))/VG(2,1)**2
        XNK(1) = XM(1)/VG(1,1)**2
        XN(1) = XNK(1)/VG(1,1)
    
        XLN2SIG = 0.0
        do I = 2,NDISC(1)
            XLN2SIG = XLN2SIG + (DLNVG(1,I)-XLNVG)**2*XNK(I)
        end do
        do I = 2,NDISC(2)
            XLN2SIG=XLN2SIG+(DLNVG(2,I)-XLNVG)**2*XNK(I+NDISC(1))
        end do
        do I = 1,NSEC
            XLN2SIG = XLN2SIG+XM(I+NDISCT)/DVK(I)*(((DLNVL(I)-XLNVG+1.0)**2+1.0)/VL(I) - ((DLNVL(I+1)-XLNVG+1.0)**2+1.0)/VL(I+1))
        end do
        XLN2SIG = XLN2SIG/XNN/9.0
        TXLN2SIG = dexp(dsqrt(XLN2SIG))
    end if

    open(2, file=RESULT2)
    write(2,110)COUNT
    write(2,111)ETA
    write(2,112)(VG(1,I),XN(I),XM(I),I=2,NDISC(1))
    write(2,112)(VG(2,I),XN(I+NDISC(1)),XM(I+NDISC(1)),I=2,NDISC(2))
    write(2,112)(VGS(I-NDISCT),XN(I),XM(I),I=NDISC1,M)
    write(2,*)
    write(2,118)XNN,TVG,TXLN2SIG,XMM1,VM
    write(2,119)TTVOL
    write(2,*)
    write(2,*)
    110 format(1X,1PE16.7,'  Sec')
    111 format(1X,'ETA = ',I1)
    112 format(1X,1P3E16.7E3)
    118 format(1X,1P5E16.7)
    119 format(1X,1PE16.7)
    open(3, file='selpre.dat')
    write(3,110)COUNT
    write(3,111)ETA
    do I=1,NDISC(1)
        write(3,112)VG(1,I)*XNN2/TTVOL2,XN(I)*TTVOL2/XNN2**2
    end do
    do I=NDISC1,M
        write(3,112)VGS(I-NDISCT)*XNN2/TTVOL2,XN(I)*TTVOL2/XNN2**2
    end do
    write(3,*)
    write(3,118)XNN,TVG,TXLN2SIG,XMM1,VM
    write(3,119)TTVOL
    write(3,*)
    write(3,*)
      
    open(4, file='sdf.dat')
    write(4,120)COUNT,XNN,TVG,TXLN2SIG,XMM1,VM,TTVOL
    120 format(1X,1P7E16.7E3)

    return
    end subroutine OT


end module ot_mod
