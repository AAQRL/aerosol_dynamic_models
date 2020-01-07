module aux_data_funcs

    ! precision specification
    use preci, wp => dp

    ! shared variables
    use globals

    implicit none


    ! used for sharing between TEMPHIST and TEMPXT
    integer :: N
    real (wp), dimension(3000) :: POSI,TMP,VEL

contains

    subroutine TEMPHIST(TEMP0,FLVEL0)
!*******************************************************************
!       The program is written to retrieve temperature history data.   
!       The max number of datapoints is 30.
!*******************************************************************
!
    real (wp) :: TEMP0,FLVEL0
    integer :: I
!*******************************************************************
!       N       Number of datapoints
!       POSI    Position of data points (cm)
!       TMP     Temperature of the position specified
!*******************************************************************
    open(25, file='TMPHIST.DAT',STATUS='OLD')
    read(25,*)N

    do I = 1,N
        read(25,*)POSI(I),TMP(I),VEL(I)
    end do
    TEMP0 = TMP(1)
    FLVEL0 = VEL(1)
    close(25)

    return
    end subroutine TEMPHIST


    subroutine TEMPXT(X,T,V)
!*******************************************************************
!       This program is written to determine the temperature of the    *
!       specified point.  The temperature between datapoints is linear *
!*******************************************************************
    real (wp) :: X,T,V
    integer :: I
!*******************************************************************
!       N       Number of datapoints
!       POSI    Position of data points (cm)
!       TMP     Temperature of the position specified
!*******************************************************************
    I = 0
    10 I = I + 1
    if ((X >= POSI(I)) .AND. (X < POSI(I+1))) then
        T = ((POSI(I+1) - X)*TMP(I) + (X - POSI(I))*TMP(I+1)) &
        /(POSI(I+1) - POSI(I))
        V = ((POSI(I+1) - X)*VEL(I) + (X - POSI(I))*VEL(I+1)) &
        /(POSI(I+1) - POSI(I))
    else
        GOTO 10
    end if

    return
    end subroutine TEMPXT



    subroutine VAPORSAT(TEMP,PSAT,A,B)
!*******************************************************************
!	This subroutine is to determine the saturation pressure at the *
!	given temperature.                                             *
!*******************************************************************
    real (wp), dimension(2) :: PSAT,A,B
    real (wp) :: TEMP
    integer :: I
!*******************************************************************
!	TEMP	Temperature (K)
!	PSAT	Saturation pressure (atm)
!       A       Temperature coefficient
!       B       Intersection
!*******************************************************************
    do I = 1,2
        PSAT(I) = 10**(A(I)/TEMP+B(I))/760.0D0
    end do

    return
    end subroutine VAPORSAT


    subroutine MTSAT(A,B)
!*******************************************************************
!	This subroutine is to read the saturation vapor pressure data. *
!*******************************************************************
    real (wp), dimension(2) :: A,B
    character(len=12), dimension(2) :: NAMEMT
    integer :: I
!*******************************************************************
!	A	Temperature coefficient
!       B       Intersection
!   	NAMEMT	Name of the metal
!*******************************************************************
    open(51, file='MTPSAT.DAT')
    do I = 1,2
        read(51,*) NAMEMT(I)
        read(51,*) A(I),B(I)
    end do
    close(51)

    return
    end subroutine MTSAT



    subroutine RXNRT(TEMP,XK1,A0,EA)
!*******************************************************************
!	This subroutine is to determine the reaction rate              *
!*******************************************************************
    real (wp), dimension(2) :: XK1,A0,EA
    real (wp) :: TEMP
    integer :: I
!*******************************************************************
!	TEMP	Temperature (K)
!	XK1	Reaction rate coefficient
!       A0      Pre-Exponent term of the reaction rate
!       EA      Activation enegry
!*******************************************************************
    do I = 1,2
        XK1(I) = A0(I)*EXP(-EA(I)/TEMP)
    end do

    return
    end subroutine RXNRT



    subroutine RTDATA(A0,EA)
!*******************************************************************
!	This subroutine is to read the reaction rate data.             *
!*******************************************************************
    real (wp), dimension(2) :: A0,EA
    character(len=12), dimension(2) :: NAMEMT
    integer :: I
!*******************************************************************
!	A0	Pre-Exponent term of the reaction rate
!       EA      Activation enegry
!   	NAMEMT	Name of the metal
!*******************************************************************
    open(52, file='MTRXN.DAT')
    do I = 1,2
        read(52,*) NAMEMT(I)
        read(52,*) A0(I),EA(I)
    end do
    close(52)

    return
    end subroutine RTDATA



    subroutine NUCLRT(TEMP,XK2,A1,B1,A2,B2,TN1,TN2,TP1)
!*******************************************************************
!	This subroutine is to determine the nucleation rate            *
!*******************************************************************
    real (wp), dimension(2) :: XK2,A1,B1,A2,B2,TN1,TN2,TP1
    real (wp) :: TEMP
    integer :: I
!*******************************************************************
!	TEMP	Temperature (K)
!	XK2	Nucleation rate coefficient
!       A1      Temperature coefficient
!       B1      Intersection
!*******************************************************************
    do I = 1,2
        if (TP1(I) == 1) XK2(I) = 10**(A1(I)/TEMP+B1(I))
        if (TP1(I) == 2) XK2(I) = A1(I)*(TEMP**TN1(I))*exp(B1(I)/TEMP) &
        + A2(I)*(TEMP**TN2(I))*exp(B2(I)/TEMP)
    end do

    return
    end subroutine NUCLRT


    subroutine NUCLDATA(A1,B1,A2,B2,TN1,TN2,TP1)
!*******************************************************************
!	This subroutine is to read the nucleation rate data.           *
!*******************************************************************
    real (wp), dimension(2) :: A1,B1,A2,B2,TN1,TN2,TP1
    character(len=12), dimension(2) :: NAMEMT
    integer :: I
!*******************************************************************
!	A1	Temperature coefficient
!       B1      Intersection
!   	NAMEMT	Name of the metal
!*******************************************************************
    open(53, file='MTNUCL.DAT')
    do I = 1,2
        read(53,*) NAMEMT(I)
        read(53,*) TP1(I)
        if (TP1(I) == 1) read(53,*) A1(I),B1(I)
        if (TP1(I) == 2) read(53,*) A1(I),B1(I),TN1(I),A2(I),B2(I), &
        TN2(I)
    end do
    close(53)

    return
    end subroutine NUCLDATA

end module aux_data_funcs
