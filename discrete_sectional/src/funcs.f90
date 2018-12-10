module funcs
    
    ! precision specification
    use preci, wp => dp

    ! shared variables
    use globals
    use calbeta_common

    implicit none

contains

    function FX(V)

    real (wp) :: FX, V

    if (IDD == 1) FX = dsqrt(CONVV(K)+VG(II,1)/V)*(CONV3(K)+(V/VG(II,1))**C13)**2*(1.0/VG(II,K)+1.0/V)**ETA
    if (IDD == 2) FX = dsqrt(CONVV(K)+VG(II,1)/V)*(CONV3(K)+(V/VG(II,1))**C13)**2
    if (IDD == 4) FX = dsqrt(CONVV(K)+VG(II,1)/V)*(CONV3(K)+(V/VG(II,1))**C13)**2/V**ETA
    if (IDD == 6) FX = (V/VG(II,1))**C12
    if (IDD == 7) FX = (V/VG(II,1))**C12/V**ETA
    if (IDD == 9) FX = (V/VG(II,1))**C12*(1.0/VG(II,1)+1.0/V)**ETA

    return
    end function FX


    function F(U,V)

    real (wp) :: F, U, V

    if (IDD >= 15) then
        F = dsqrt(V1(1)/U+V1(1)/V)*((U/V1(1))**C13+(V/V1(1))**C13)**2*(1.0/U+1.0/V)**ETA
    elseif (IDD == 13) then
        F = dsqrt(V1(1)/U+V1(1)/V)*((U/V1(1))**C13+(V/V1(1))**C13)**2*(1.0/U**ETA+1.0/V**ETA)
    else
        F = dsqrt(V1(1)/U+V1(1)/V)*((U/V1(1))**C13+(V/V1(1))**C13)**2/U**ETA
    end if

    return
    end function F


    function G(U)

    real (wp) :: G, U

    if (IDD == 12) G = VL(I)
    if (IDD == 13) G = VL(J)
    if (IDD == 14) G = VL(J)
    if (IDD == 15) G = VL(I)
    if (IDD == 16) G = VL(J)
    if (IDD == 17) G = VL(I)
    if (IDD == 21) G = VL(J)
    if (IDD == 22) G = VL(J)
    if (IDD == 23) G = VL(K)
    if (IDD == 24) G = VL(I) - VL(K)
    if (IDD == 25) G = VL(I) - U
    if (IDD == 26) G = VL(I) - U
    if (IDD == 31) G = VL(K)
    if (IDD == 32) G = VL(K)
    if (IDD == 33) G = VL(K)
    if (IDD == 34) G = VL(K)
    if (IDD == 35) G = VL(K)
    if (IDD >= 36) G = VL(I) - U

    return
    end function G


    function H(U)

    real (wp) :: H, U

    if (IDD == 12) H = VL(I+1)
    if (IDD == 13) H = VL(J+1)
    if (IDD == 14) H = VL(J+1)
    if (IDD == 15) H = VL(I+1) - U
    if (IDD == 16) H = VL(J+1) - U
    if (IDD == 17) H = VL(I+1) - U
    if (IDD == 21) H = VL(I+1) - U
    if (IDD == 22) H = VL(I+1) - VL(K+1)
    if (IDD == 23) H = VL(I+1) - U
    if (IDD == 24) H = VL(J+1)
    if (IDD == 25) H = VL(K+1)
    if (IDD == 26) H = VL(J+1)
    if (IDD == 31) H = VL(I+1) - U
    if (IDD == 32) H = VL(I+1) - U
    if (IDD == 33) H = VL(K+1)
    if (IDD == 34) H = VL(I+1) - U
    if (IDD >= 35) H = VL(K+1)

    return
    end function H

end module funcs
