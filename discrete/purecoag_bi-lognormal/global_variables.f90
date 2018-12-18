module global_variables
    implicit none
    
    ! CONSTANTS
    real, parameter :: KB = 1.38064852e-23, EE = 1.6021766208e-19, PI = 3.14159265, KE = 8.99e9, NA = 6.02214086e23
    real, parameter :: P = 101325.0, Mg = 29.2e-3
    ! KB = boltzmann Constant (SI)
    ! EE = elementary charge on an electron
    ! KE = 1/(4.0*PI*eps0)
    ! NA = Avagadro's Constant (mol-1)
    ! P = Pressure (N/m2)
    ! Mg = Molecular Weight of gas (Kg/mol)
    
    ! MOLECULAR WEIGHT AND DENSITY 
     double precision :: mw = 100e-3, rho = 1000.0, Temp = 300 !1200.0
    ! mw = Molecular Weight (Kg/mol)
    ! rho = Density (Kg/m3)     
    ! Temp = temperature (K)    
        
    ! NUMBER OF DISCRETE SIZES, CHARGES
    integer :: n_discrete = 100 ! number of discrete sizes eg. 100
    !integer :: n_charge = 1 ! number of possible charges eg. 3
    !integer, dimension(:), allocatable :: charge  ! possible charges eg. 0 -1 +1
    
    ! INITIAL
    double precision :: dp0, v0 = 1e-27  !Inital size of the discrete sizes, dp0 (m)
    double precision, dimension(:), allocatable :: Nc0 ! Initial number concentration for each 
    ! discrete size (n_discrete, n_charge)
    
    ! FINAL
     double precision, dimension(:), allocatable :: Nc ! Number concentration for each  (#/m^3)
    ! discrete size (n_discrete, n_charge) at a fixed t
    
    double precision, dimension(:), allocatable :: YNc !, YNc_mono ! for solver (#/m^3)
    
    double precision, dimension(:), allocatable :: dp_dis !dp in each discrete section (m)
    
    ! OTHERS
    integer :: flag ! defines coagulation theory to be used
    
    
    ! SOLVER VARIABLES
    double precision :: rtol , t0 = 0.0, tf , dt 
    double precision, dimension(:), allocatable :: atol
    integer :: istate , itask = 1, neq
    
    ! FOR TEMPHIST
    !double precision, dimension(:,:), allocatable :: T_data
    !integer :: time_points
    
    ! Non-Dimensional Parameters
    !double precision :: K_coag = 8.0*KB*1200.0/(3.0*1.81e-5)
    !double precision :: N_infi_0, tau, 
    !double precision :: tau_mono
    
    ! For beta calculation storage
    double precision, dimension(:,:), allocatable :: beta_all
    double precision :: N_infi_0, K_coag, tau, N_track_min
    ! beta_all is the matrix which contains all the values of the beta 
    ! beta_safe is the module
    ! beta is the function which finds the value of beta
    ! fill_beta is the subroutine which fills in all the values of beta
        
end module global_variables

    
!program test_global_variables
!    use global_variables    
!    implicit none
!    
!    print*, n_charge
!    
!end program test_global_variables
    

    
    
