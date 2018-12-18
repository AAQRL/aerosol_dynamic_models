module global_variables
    implicit none
    
    !! CONSTANTS
    double precision, parameter :: KB = 1.38064852e-23, EE = 1.6021766208e-19, PI = 3.14159265, KE = 8.99e9, NA = 6.02214086e23
    ! KB = boltzmann Constant (SI)
    ! EE = elementary charge on an electron
    ! KE = 1/(4.0*PI*eps0)
    ! NA = Avagadro's Constant (mol-1)
    
    !! PARTICLE/MATERIAL
    double precision ::  p10, p1, ps, S, n10, n1, m1, vm, v_gen, dp_critical, v_critical, I_nucl
    ! p10 = Initial Vapor Pressure (Pa)
    ! p1 = Vapor Pressure (Pa)
    ! ps = Saturation Pressure (Pa)
    ! S = p1/ps
    ! n10 = Molecules/m3 in Gas Phase
    ! n1 = Molecules/m3 in Gas Phase
    ! m1 = Mass of 1 Molecule (kg/molecule)
	! v_gen = Vapor Generation Rate (#/m3/s)
    ! dp_critical = Critical Nuclei Size (m)
    ! v_critical = Critical Nuclei Volume (m3)
    ! I_nucl = Nucleation Rate (#/m3/s)
    
    !! SPECIES PROPERTIES
    double precision ::  Mw_p, rho_p, sigma
    ! Mw_p = Molecular Weight of the species (kg/mol)
    ! rho_p = Density of the species (kg/m3)
	! sigma = Surface Tension (N/m)

    !! GAS PROPERTIES
    double precision ::  Mw_g, rho_g_ref, mu_g_ref, S_ref, lm_ref
    double precision :: lm, mu ! Derived Parameters
    ! Mw_g = Molecular Weight of gas (kg/mol) 
    ! rho_g_ref = Density of gas at 298 K (kg/m3)
    ! mu_g_ref = Viscosity of gas at Tref = 298 K (Pa.s)
    ! S_ref = Sutherland's Constant
    ! lm_ref = Mean free path of gas molecules at Tref = 298 K (m)
    ! mu = Viscosity of gas at T = Temp (Pa.s)
    ! lm = Mean free path of gas molecules at T = Temp (m)


    !! SYSTEM PROPERTIES
    double precision :: P, Temp
    ! P = Pressure (Pa)
    ! Temp = Temperature (K)
    
    ! SIMULATION TIME
    double precision :: t_final, dt
    double precision :: time(1000)!dummy
    ! tf = Final time in seconds
    ! dt = Time step
    ! Minimum Dimensionless Number Concentration Tracking to find atol and rtol
    
	! MODES
    integer :: N_modes
	! N_modes = number of modes in the system
	
	! AEROSOL PROCESSES (if ON set to 1)
	integer :: flag_coag, flag_nucl, flag_cond
    
	! INITIAL SIZE DISTRIBUTION
    double precision :: N_init(1000), v_init(1000)
	
	! VARIABLES FOR CALCULATIONS
	double precision, dimension(:,:), allocatable :: N_result, v_result
    double precision, dimension(:), allocatable :: M, dp_g, v_g, sigma_g, t_result, N_total
    integer :: nlines
	
end module global_variables