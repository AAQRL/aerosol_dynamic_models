program main
    use global_variables
    use input
    use derived_parameters
    use rhs
    use DVODE_F90_M
    use calculations
    implicit none      
    
    double precision, allocatable :: ATOL(:), Y(:)
    double precision :: RTOL, T, TOUT, RSTATS
    integer :: NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, R, kk, STEPS
    DIMENSION :: RSTATS(22), ISTATS(31)
    TYPE(VODE_OPTS) :: OPTIONS
    OPEN(UNIT=6, FILE = 'results_raw.dat')
    
    call read_inputfile()
    call derived_param()
    
    NEQ = (N_modes*2)+1
    
    allocate(Y(NEQ))
    allocate(ATOL(NEQ))
    
    n10 = p10/(kB*Temp);
    Y(NEQ) = n10;     
    do kk = 1, N_modes
        Y(2*kk-1) = N_init(kk);
		Y(2*kk-0) = v_init(kk)*Y(2*kk-1);
    end do

	! Error Tolerances
    do kk = 2, N_modes
        ATOL(2*kk-1) = (1D-20)
        ATOL(2*kk-0) = (1D-30)
    end do
    ATOL(NEQ) = (1D-10)
    RTOL = 1D-6
    
    T = 0
    TOUT = dt
    STEPS = t_final/dt
    ITASK = 1
    ISTATE = 1
    OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR_VECTOR=ATOL)    
    
	WRITE(6,63)T,Y
	time(1) = T
	nlines = 1
    do IOUT = 1,STEPS
        call DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
        CALL GET_STATS(RSTATS,ISTATS)
        WRITE(6,63)T,Y
       TOUT = TOUT + dt
       time(IOUT+1) = T
	   nlines = nlines + 1
    end do

63  FORMAT(100E14.6)
call calc()
end program main
