! This module calculates the derived variables from global variables
    
module derived_variables
    use global_variables
    implicit none
    
contains
    subroutine derived(ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, delta_t, n_d)
    implicit none
    double precision :: ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero
    double precision :: n_d
    double precision :: t_start, t_final, delta_t
    
    vg_zero = (22.0/7)/6*(dpg_zero)**(3.0) 
    ns = Ps / (kb*T)        
    mu = 1.716e-5 * (T/273)**(2.0/3)                                           ! viscosity of the medium
    lam = mu/P * sqrt(PI*kb*T/(2.0*Mg/Na))                                   ! mean free path of air (m)
    
    mm1 = M/Na                                                                ! mass of a monomer (Kg)
    v1 = mm1/rho                                                              ! volume of a monomer (m3)
    r1 = (1.0/2) * (6.0*v1/PI)**(1.0/3)                                      ! radius of a monomer (m)
    s1 = 4.0 * PI * r1 * r1                                                  ! area of a monomer (m2)
    sig= gam * v1**(2.0/3) / (kb*T)                                          ! surface tension group (dimensionless)
    tau= (ns*s1*(kb*T/(2.0*PI*mm1))**(1.0/2))**(-1.0)                         ! Characterstic time for particle growth (s)
    K  = (2.0*kb*T/(3.0*mu))*ns*tau                                          ! Coagulation coefficient (dimensionless)   

    !if (num_zero > 1e9) then
    !   n_d = num_zero
    !else
    n_d = ns
    !end if
    

    t_start = 0.0
    delta_t = dt/tau
    
    
    end subroutine derived
    
end module derived_variables

!program test_derived_variables
!    use derived_variables
!    implicit none
!    double precision :: ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, t_final, delta_t, n_d
!    call derived(ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, t_final, delta_t, n_d)
!    
!    print*, ns
!    
!end program test_derived_variables