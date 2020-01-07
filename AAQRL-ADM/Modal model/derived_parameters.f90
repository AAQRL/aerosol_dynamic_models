module derived_parameters
    use input
    use global_variables
    implicit none 
    contains

	!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine derived_param()            
        implicit none
        integer :: i
        call read_inputfile()
        ! gas derived parameters
        lm = lm_ref*(temp/298)*(101325/p)*(1+(S_ref/298))/(1+(S_ref/temp)) ! mean free path (m)
        mu = mu_g_ref*((298+S_ref)/(temp+S_ref))*((temp/298)**(1.5)) ! viscosity (pa.s)
    
        m1 = Mw_p/NA ! molecular mass
        vm = m1/rho_p ! molecular volume
    end subroutine derived_param
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!

    double precision function beta(v1,v2)
		implicit none
        double precision, intent(in) :: v1, v2
		double precision :: dc1, dc2, Kn1, Kn2, D1, D2, c1, c2, l1, l2, g1, g2
		call read_inputfile()
        
        dc1 = (6*v1/PI)**(1.0/3)
        dc2 = (6*v2/PI)**(1.0/3)
        Kn1 = 2*lm_ref/dc1
		Kn2 = 2*lm_ref/dc2            
        
        D1 = (KB*Temp/(3*PI*mu_g_ref*dc1))*((5+4*Kn1+6*(Kn1**2)+18*(Kn1**3))/(5-Kn1+(8+PI)*(Kn1**2)))
        D2 = (KB*Temp/(3*PI*mu_g_ref*dc2))*((5+4*Kn2+6*(Kn2**2)+18*(Kn2**3))/(5-Kn2+(8+PI)*(Kn2**2)))
		c1 = ((8*KB*Temp)/(PI*rho_p*v1))**(0.5)
        c2 = ((8*KB*Temp)/(PI*rho_p*v2))**(0.5)
              
		l1 = 8*D1/(PI*c1)
        l2 = 8*D2/(PI*c2)

		g1 = (1/(3*dc1*l1))*(((dc1+l1)**3) +((dc1*dc1+l1*l1)**(3/2))) - dc1
        g2 = (1/(3*dc2*l2))*(((dc2+l2)**3) +((dc2*dc2+l2*l2)**(3/2))) - dc2
        
		beta = 2*PI*(D1+D2)*(dc1+dc2)*((((dc1+dc2)/(dc1+dc2 + 2*((g1**2 + g2**2)**(0.5)))) + (8*(D1+D2)/((((c1**2 + c2**2)**(0.5)))*(dc1+dc2))))**(-1))
    end function beta
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!
	
	function fractn(N_modes,vk) ! COAGULATION- REDISTRIBUTION OF PARTICLE BETWEEN BINS
		implicit none
        integer, intent(in) :: N_modes
		double precision, intent(in) :: vk(N_modes)
		double precision :: fractn(N_modes,N_modes,N_modes)
		double precision :: v_new(N_modes,N_modes), sum(N_modes,N_modes), L, R
		integer :: i, j, k, flag
		! v_new(i,j) gives the the new volume of particle when particles in mode i and j collide
		! sum(i,j) gives the mode whose volume is least smaller than v(i) + v(j)
		! L is a dummy variable
		
		sum = 0
		do i = 1, N_modes
			do j = 1, N_modes
				v_new(i,j) = vk(i) + vk(j)
				do k = 1, N_modes
					IF (vk(k) < v_new(i,j)) THEN
						flag = 1
					ELSE
						flag = 0
					END IF
					sum(i,j) = sum(i,j) + flag
				end do
				! Now, sum(i,j) = mode # whose volume is least smaller than v(i) + v(j)
			end do
        end do
		
        
		fractn = 0
		do i = 1, N_modes
			do j = 1, N_modes		
				L = sum(i,j)
				IF (L == N_modes) THEN
					fractn(i,j,L)=1
				ELSE
					R = vk(L+1)/vk(L)
					fractn(i,j,L)=(R*vk(L)-v_new(i,j))/((R-1)*vk(L))
					fractn(i,j,L+1)= (v_new(i,j)-vk(L))/((R-1)*vk(L))
				END IF
			end do
		end do
    end function fractn
	
	!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!    
end module derived_parameters