! All the functions are defines here
! Moment to calculate the moment of the size distribution
! coag to calculate the coagulation coefficient. Harmonic mean of free molecular regime, and continuum regime is used
! cond to calculate the condensation coefficient. Harmonic mean of free molecular regime, and continuum regime is used
! nucl to calculate the critical particle size, and the rate of nucleation
! reac to calculate the rate of reaction
    
module function_list
    use global_variables
    use derived_variables
    implicit none
    double precision :: ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, t_final, delta_t, n_d
    
    
    
contains

    function moment(Nv,vgv,sdv,kv)
        implicit none
        double precision :: Nv,vgv,sdv, moment
        double precision :: kv
        call derived(ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, delta_t, n_d)
        moment = Nv*(vgv)**(kv)*exp(9.0/2*(kv*kv)*log(sdv)*log(sdv))
    end function moment   
    
    function coag(sd,rg,flag)
		
		double precision :: sd,rg,M0
		integer :: flag
        double precision :: b0, b2, c1, c2, lnsig, b5, com0fm, com0cn, com2fm, com2cn, com0, com2
		double precision, dimension(2) :: coag
		
		if ((flag == 0) .OR. (sd <= 1e-5) .OR. (rg <= 1e-15)) then
			coag(1) = 0.0
			coag(2) = 0.0
  
		
		else
			b0 = 0.633 + 0.092*sd*sd - 0.022*sd*sd*sd
			b2 = 0.39 + 0.5*sd - 0.214*sd*sd + 0.029*sd*sd*sd
			c1 = (6.0*kb*T*rg/rho)**(1.0/2)
			c2 = 2.0*kb*T/(3.0*mu)
			
			lnsig = log(sd)*log(sd)
			b5 = 1.257
			
			com0fm = c1 * b0 * ( exp(25.0/8* lnsig) + 2.0*exp(5.0/8*lnsig) + exp(1.0/8*lnsig) ) 
			com0cn = c2 * ( 1.0 + exp(lnsig) + b5 * (lam/rg) * exp(1.0/2 * lnsig) * (1.0 + exp(2.0*lnsig)) )
			
			com2fm = 2.0 * c1 * b2 * exp(3.0/2*lnsig) * ( exp(25.0/8* lnsig) + 2.0*exp(5.0/8*lnsig) + exp(1.0/8*lnsig) )
			com2cn = 2.0 * c2 * ( 1.0 + exp(lnsig) + b5 * (lam/rg) * exp(-1.0/2 * lnsig) * (1.0 + exp(-2.0*lnsig)))
			
			com0 = (com0fm * com0cn) / (com0fm + com0cn) 
			com2 = (com2fm * com2cn) / (com2fm + com2cn)
			coag(1) = com0
			coag(2) = com2
		end if
    end function coag
    
    
    
    
    
    function cond(M0,sd,vg,flag)
        implicit none
        double precision :: lnsig, A, F1, F2, C1, etac1, delc1, psic1, etac2, delc2, psic2, M0
        double precision :: sd,vg
        double precision, dimension(3) :: cond
        integer :: flag
        
        lnsig = log(sd)*log(sd)
        A  = 1.0/3
        F1 = (36.0*PI)**A * v1 * ns * sqrt(kb*T/2.0/PI/mm1)
        F2 = sqrt(8.0*kb*T/(PI*mm1))
        C1 = (48.0*PI*PI)**A*ns*v1*lam/3.0
        
        if ((flag == 0) .or. (vg <= 0.0)) then
            cond(1) = 0.0
            cond(2) = 0.0
            cond(3) = 0.0
            
        elseif (flag== 1) then
            etac1 = v1**A /tau * vg**(2.0*A) * exp(2.0*lnsig)
            delc1 = etac1
            psic1 = 2.0 * v1**A /tau *vg**(2.0*A) * exp(8.0 * lnsig)
            etac2 = 8.0/3/tau/(2.0*r1) *lam * v1**(2.0*A) * vg**A * exp(lnsig/2)
            delc2 = etac2
            psic2 = 16.0/3/tau/(2.0*r1) * lam * v1 **(2.0*A) * vg**(A) * exp(3.5 * lnsig)
            cond(1) = etac1 * etac2 / (etac1 + etac2) 
            cond(2) = delc1 * delc2 / (delc1 + delc2) 
            cond(3) = psic1 * psic2 / (psic1 + psic2) 
        elseif (flag ==2) then
            cond(1) = 0.0 
            cond(2) = 0.0 
            cond(3) = 0.0
        end if
    end function cond
    
    function nucl(S,flag)
        double precision :: S, sigmad, xks, kstar, coef1, pep, coef2, x11
        double precision, dimension(2) :: nucl
        integer :: flag
        
        if ((S < 1.001) .or. (flag == 0)) then
            nucl(1) = 0.0
            nucl(2) = 0.0

        elseif (flag == 1) then
            sigmad = gam * v1 ** (2.0/3) / (kb * T)
            xks = PI/6 * (4.0 * sigmad)**(3.0)
            kstar = xks / (log(S))**(3.0)
            coef1 = (2.0 / 9.0 / PI) ** (1.0/3)
            pep = kstar * (log (S)/2)
            if (pep < 60000.0) then
                coef2 = exp(-pep)
                x11 = ns / tau * S**(2.0) * coef1 * sigmad **(0.5) * coef2
                nucl(1) = kstar
                nucl(2) = x11

            else
                nucl(1) = 0.0
                nucl(2) = 0.0
            end if
        end if
        
        if (flag ==2) then
            nucl(1) = 1.0
            nucl(2) = 1.54*1e6 * exp(-74.03e3/(8.314*T)) * S * ns
        end if
    end function nucl

    function reac(flag)
        implicit none
        integer :: flag
        double precision :: Reac
        
        if (flag == 0) then
            Reac = 0.0
        else
            Reac = r_rate*(ns/tau)
        end if
    end function reac
        

end module function_list

!program test_function_list
!    use function_list  
!    use derived_variables
!    implicit none
!    double precision :: ss1 = 100.0
!    double precision, dimension(2) :: nuclout
!    
!    double precision :: momout
!    double precision,dimension(2) :: coagout
!    double precision,dimension(3) :: condout
!    double precision :: num_zero_0 = 1e11, dpg_zero_0 = 10.0e-9, sigma_g_zero_0 = 1.5, vg_zero_0
!    
!    double precision :: insd = 1.0, momo = 0.0
!    call derived(ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, delta_t, n_d)
!
!    vg_zero_0 = (22.0/7.0/6.0) * (dpg_zero_0)**3.0
!    momout  = moment(num_zero_0,vg_zero_0,sigma_g_zero_0,momo)
!    nuclout = nucl(ss1,1)
!    coagout = coag(insd,r1,1)
!    condout = cond(insd,v1,1)
!    print*, num_zero_0,vg_zero_0,sigma_g_zero_0
!    print*, momout
!    print*, coagout(1), coagout(2)
!    print*, condout(1), condout(2), condout(3)
!    print*, nuclout(1), nuclout(2)
!    print*, ns
!    
!end program test_function_list