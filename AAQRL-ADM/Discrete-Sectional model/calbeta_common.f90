module calbeta_common

    ! variables shared between the CALBETA subroutine (calbeta_mod.f90) and others

    ! precision specification
    use preci, wp => dp

    ! shared variables
    use globals

    ! used for sharing between CALBETA and XDOTEQ
    real (wp), dimension(2,MAXDISC,MAXDISC) :: BDD
    real (wp), dimension(2,MAXDISC,MAXSEC,MAXSEC) :: BDS1
    real (wp), dimension(2,MAXDISC,MAXSEC) :: BDS25, BDS4
    real (wp), dimension(MAXSEC,MAXSEC,MAXSEC) :: BSS1
    real (wp), dimension(MAXSEC,MAXSEC) :: BSS25,BSS4
    real (wp), dimension(MAXSEC) :: BSS36

    ! used for sharing between CALBETA and FX, F, G, H in funcs
    integer :: IDD, I, J, K, II

end module calbeta_common
