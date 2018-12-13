! This module definbes the global_variables which are used at different places int he other modules
module global_variables
    implicit none
    
    ! CONSTANTS
    double precision, parameter :: kb = 1.38064852e-23, EE = 1.6021766208e-19, PI = 3.14159265, KE = 8.99e9, NA = 6.02214086e23
    double precision, parameter :: P = 101325.0, Mg = 29.2e-3
    ! KB = boltzmann Constant (SI)
    ! EE = elementary charge on an electron
    ! KE = 1/(4.0*PI*eps0)
    ! NA = Avagadro's Constant (mol-1)
    ! P = Pressure (N/m2)
    ! Mg = Molecular Weight of gas (Kg/mol)
    
    ! MOLECULAR WEIGHT AND DENSITY 
    double precision :: M , rho , T , gam , Ps, A, B , r_rate
    integer :: nlines
    double precision, dimension(:), allocatable :: vz_pos, z_pos, time, Temperature, vapor_conc
    !double precision :: time(5), Temperature(5), vapor_conc(5)
    double precision :: num_zero_0 , dpg_zero_0 , sigma_g_zero_0 , vapor_zero_0 
    ! M = Molecular Weight (Kg/mol)
    ! rho = Density (Kg/m3)     
	! T = Temperature (K)
    ! gam = surface tension (N/m)
	! Ps = saturation pressure (N/m2)
	! r_rate = reaction rate
    
    
    ! INITIAL
    double precision :: num_zero , dpg_zero, sigma_g_zero, S_zero 
	! num_zero = Initial number concentration
    ! dpg_zero = Initial geometric mean diameter
	! sigma_g_zero = Initial geometric standard deviation
	! S_zero = Initial saturation ratio
	

    
    ! OTHERS
    integer :: flag_COAG , flag_COND , flag_NUCL , flag_REAC ! defines flags for different inputs
    
    
    ! SOLVER VARIABLES
    double precision ::  t0 = 0.0 , dt = 0.01
	! S = saturation ratio
    !integer :: istate , itask = 1, neq
    
    ! Non-Dimensional Parameters
    !double precision :: K_coag = 8.0*KB*1200.0/(3.0*1.81e-5)
    !double precision :: N_infi_0, tau, 
    !double precision :: tau_mono
    
    ! For beta calculation storage
    !double precision, dimension(:,:), allocatable :: beta_all
    !double precision :: N_infi_0, K_coag, tau, N_track_min
    ! beta_all is the matrix which contains all the values of the beta 
    ! beta_safe is the module
    ! beta is the function which finds the value of beta
    ! fill_beta is the subroutine which fills in all the values of beta
    
    
    
end module global_variables

    
!program test_global_variables
!    use global_variables    
!    implicit none
!    
!    print*, t0
!    
!end program test_global_variables
    

    
    
