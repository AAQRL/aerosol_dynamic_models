module beta_safe
    use global_variables
    implicit none
    
    
contains
    
    function get_beta_ij(ii, jj)
    
      implicit none
      double precision :: get_beta_ij
      integer :: ii, jj
      get_beta_ij = beta(ii*v0, jj*v0, flag)
    
    end function get_beta_ij
    
    
    subroutine fill_beta()
    
        implicit none
        integer :: i,j
        
        do i = 1,n_discrete
            do j = 1,n_discrete
                beta_all(i,j) = beta(i*v0,j*v0,flag) 
            end do
        end do
    end subroutine fill_beta
             
    subroutine KK()
        implicit none
        
        double precision :: mu
        
        mu = 1.716e-5 * (Temp/273) **(2.0/3)
        K_coag = 8.0*KB*Temp/(3.0*mu)
        tau = 2.0/(K_coag*N_infi_0)
        tf = tf/tau
        dt = dt/tau
        !print*, t
    end subroutine KK
    
    
    function Diff(dp)
        implicit none
        double precision :: dp, Kn, t, Diff, mu, lam
        mu = 1.716e-5 * (Temp/273) **(2.0/3)
        lam = mu/P * sqrt(PI*KB*Temp/(2.0*Mg/NA))
        Kn = 2.0*lam/dp
        Diff = KB*Temp/(3.0*PI*mu*dp) * ((5.0 + 4.0*Kn + &  
		6.0*Kn*Kn + 18.0*Kn*Kn*Kn)/(5.0 - Kn + (8.0 + PI)*Kn*Kn))
    end function Diff
    
    function beta(v1,v2,flag)
        implicit none
        double precision  :: v1, v2
        double precision  :: beta
        double precision :: dp0, dp1, dp2, c1, c2, l1, l2, g1, g2, & 
   				coeff1, coeff2
        integer :: flag
        
        dp0 = (6.0/PI * v0)**(1.0/3)
        dp1 = (6.0/PI * v1)**(1.0/3)
        dp2 = (6.0/PI * v2)**(1.0/3)
        
    
        select case (flag)
            
        case(0)
        ! Use Fuchs Sutugin Beta
            c1 = sqrt(8.0*KB*Temp/(PI*(mw/NA)))
            c2 = sqrt(8.0*KB*Temp/(PI*(mw/NA)))
            
            l1 = 8.0*Diff(dp1)/(PI*c1)
            l2 = 8.0*Diff(dp2)/(PI*c2)
            
            g1 = 1.0/(3.0*dp1*l1) * ((dp1+l1)**3-(dp1*dp1+l1*l1)**  &
			(3.0/2)) - dp1
            g2 = 1.0/(3.0*dp2*l2) * ((dp2+l2)**3-(dp2*dp2+l2*l2)**  &
			(3.0/2)) - dp2
    
    
            coeff1 = 2.0*PI*(dp1 + dp2) * (Diff(dp1) + Diff(dp2))
            coeff2 = (dp1 + dp2)/(dp1 + dp2 + 2.0*sqrt(g1*g1 +   &
			g2*g2)) + &
		(8.0*(Diff(dp1)+Diff(dp2)))/((dp1+dp2)*sqrt(c1*c1+c2*c2))
            
            if ((dp1 < dp0) .or. (dp2 < dp0)) then
                beta = 0.0
            else 
                beta =(coeff1/coeff2)
            end if
            
    
        case(1)
        ! Use FM Regime Beta
            
            beta = (3.0/(4.0*PI))**(1.0/6) * (6.0*KB*Temp/rho)**(1.0/2) * &
		           (1.0/v1 + 1.0/v2)**(1.0/2) * (v1**(1.0/3) + v2**(1.0/3))**(2.0) 
            if ((dp1 < dp0) .or. (dp2 < dp0)) then
                beta = 0.0
            
            end if
            
        case(2)
            
            beta = 1.0
            
        end select
    end function beta

    
    
    
    end module beta_safe