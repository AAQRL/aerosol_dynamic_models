module preci

    implicit none
  
    ! NAG, Cray, Intel, gfortran: all the same
    integer, parameter :: sp &
            = selected_real_kind( 6,  37)
    integer, parameter :: dp &
            = selected_real_kind(15, 307)
  
    ! NAG: (30, 291); Cray, Intel, gfortran: (33, 4931)
    integer, parameter :: qp &
            = selected_real_kind(30, 291)

end module preci
