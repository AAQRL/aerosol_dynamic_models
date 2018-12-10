module rk45_init0
    
    ! precision specification
    use preci, wp => dp

    ! external functions and subroutines
    use xdoteq_mod

    implicit none


    ! used for sharing between RK45 and INIT0
    real (wp), dimension(6) :: HA
    real (wp), dimension(6,5) :: BB
    real (wp) :: HC1,HC3,HC4,HC5,HC6

contains

    subroutine RK45(X,XDOT,T,H,M)
    integer :: M, LAM, N
    real (wp) :: T, H, TNEW
    real (wp), dimension(MAXDISC+MAXDISC+MAXSEC) :: XNEW, SA
    real (wp), dimension(MAXDISC+MAXDISC+MAXSEC,6) :: F
    real (wp), dimension(*) :: X(*),XDOT(*)
!
!        RUNGE-KUTTA INTEGRATION SCHEME DEVELOPED BY DALE BETTIS OF THE
!        UNIVERSITY OF TEXAS AT AUSTIN.  LOCAL TRUNCATION ERROR IS H**6.
!        OPTIMIZED COEFFICIENTS VERSION.
!
!        INIT0 MUST BE CALLED ONCE BEFORE PROGRAM EXECUTION
!        AND EACH TIME THE STEPSIZE IS CHANGED.
!
!        **********  DO NOT CHANGE ANYTHING IN THIS ROUTINE EVER!!! *********
!
!        ALL INFO COMES IN AND OUT THROUGH THE COMMON BLOCK AND THE ARGUMENT
!        LIST.
!
    call XDOTEQ(X,F(1,1),T)
    do K=2,6
        N=K-1
        do J=1,M
            SA(J)=0.0D0
            do LAM=1,N
                SA(J)=SA(J)+BB(K,LAM)*F(J,LAM)
            end do
            SA(J)=H*SA(J)
        end do
        TNEW=T+HA(K)
        do N=1,M
            XNEW(N)=X(N)+SA(N)
        end do
        call XDOTEQ(XNEW,F(1,K),TNEW)
    end do
    do J=1,M
        X(J)=X(J)+HC1*F(J,1)+HC3*F(J,3)+HC4*F(J,4)+HC5*F(J,5)+ &
        HC6*F(J,6)
    end do
    return
    end subroutine RK45


    subroutine INIT0(H)

    real (wp) :: H
!
!        **********  DO NOT CHANGE ANYTHING IN THIS ROUTINE EVER!!! *********
!
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
    BB(2,1)=.2D0
    BB(3,1)=3.D0/40.D0
    BB(3,2)=9.D0/40.D0
    BB(4,1)=.3D0
    BB(4,2)=-.9D0
    BB(4,3)=1.2D0
    BB(5,1)=-11.D0/54.D0
    BB(5,2)=2.5D0
    BB(5,3)=-70.D0/27.D0
    BB(5,4)=35.D0/27.D0
    BB(6,1)=-473.D0/10368.D0
    BB(6,2)=1595.D0/1728.D0
    BB(6,3)=-34595.D0/54432.D0
    BB(6,4)=38665.D0/62208.D0
    BB(6,5)=7733.D0/145152.D0
    return
    end subroutine INIT0


end module rk45_init0
