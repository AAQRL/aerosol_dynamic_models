module rhs 
	use global_variables
    use derived_parameters
	implicit none
    
contains
	
    subroutine FEX(NEQ,T,Y,YDOT)
        implicit none
        integer, intent (in) :: NEQ ! Number of ODEs = 2*Number of Modes+1 (Number and Volume balance for each mode)
        double precision, intent (in) :: T
        double precision, intent (in) :: Y(NEQ)
        double precision, intent (out) :: YDOT(NEQ)
        double precision :: N((NEQ-1)/2), V_mode((NEQ-1)/2), v((NEQ-1)/2), dNdt_coag((NEQ-1)/2), dVdt_coag((NEQ-1)/2), frac((NEQ-1)/2,(NEQ-1)/2,(NEQ-1)/2), dNdt_nucl((NEQ-1)/2), dVdt_nucl((NEQ-1)/2) ! Number of modes = (NEQ-1)/2            
        double precision :: dNdt_coag_term_1, dNdt_coag_term_2, dNdt_coag_term_3, dVdt_coag_term_2, dVdt_coag_term_3
		double precision :: B1, Dm, B3, cond_loss, B4
        double precision :: G_C((NEQ-1)/2), G_FM((NEQ-1)/2), dVdt_cond((NEQ-1)/2)
        
        integer :: i, j, k
		
        call read_inputfile()
		call derived_param()
        
		do k = 1, N_modes	! N_modes = NEQ/2 = Number of Modes
			! Y1 = N1; Y2 = V1; Y3 = N2; Y4 = V2; ...
			N(k) = Y(2*k-1) ! Total # conc. of mode k
			V_mode(k) = Y(2*k) ! Total volume of mode k	
			if (k == 1) then
				v(k) = v_init(k) ! Volm. of a particle in mode 1
			else
				v(k) = V_mode(k)/N(k) ! Volm. of a particle in mode k
			end if 
        end do
        
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!NUCLEATION TERMS!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
		dNdt_nucl = 0
		dVdt_nucl = 0
        if (flag_nucl == 1) then
			n1 = Y(NEQ)
			p1 = n1*(kB*Temp)
			S = p1/ps
			
			! Nucleation Rate
			if (S>1) then
				I_nucl = 2.0*(p1/((2.0*PI*m1*kB*Temp)**(1.0/2.0)))*(n1*vm**(2.0/3.0))*(((sigma*vm**(2.0/3.0))/(kB*Temp))**(1.0/2.0))*exp(-(16.0*PI*(sigma**3.0)*(vm**2.0))/(3*((kB*Temp)**3)*((log(S))**2.0)))
				dp_critical = 4.0*sigma*vm/(kB*Temp*log(S))
				v_critical = PI*(dp_critical**3.0)/6.0
				!v_critical = 2*vm
				!dp_critical = ((6.0/PI)*v_critical)**(1.0/3)
			else
				I_nucl = 0
				dp_critical = 0
			end if
		
			dNdt_nucl(1) = I_nucl*v_critical/v(1)
			dVdt_nucl(1) = dNdt_nucl(1)*v(1)
		end if
			
			
                
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!CONDENSATION TERMS!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
		dVdt_cond = 0
		cond_loss = 0
		if (flag_cond == 1) then
			B1 = ((36*PI)**(1.0/3))*vm*(ps/(kB*Temp))*((kB*Temp)/(2*PI*m1))**(1.0/2)	
			Dm = ((8.0*kB*Temp/(PI*m1))**(1.0/2))
			B3 = ((48.0*PI*PI)**(1.0/3))*vm*(ps/(kB*Temp))*(lm/3)*Dm
		
			do k = 1, N_modes
				G_C(k) = B3*(v(k)**(1.0/3))*(S-1)
				G_FM(k) = B1*(v(k)**(2.0/3))*(S-1)
				
				dVdt_cond(k) = (1/((1/G_C(k))+(1/G_FM(k))))*N(k)
			end do 
			
			cond_loss = 0
			do k = 1, N_modes
				cond_loss = cond_loss + dVdt_cond(k)/vm
			end do
		end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!COAGULATION TERMS!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!        
        dNdt_coag = 0
        dVdt_coag = 0
		
		if (flag_coag == 1) then
			frac = fractn(N_modes,v) ! When particles in mode i and j collide, frac(i,j,k) gives fraction of resulting particle that lies in mode k		

			do k = 1, N_modes
						
				! dN/dt coagulation terms
				dNdt_coag_term_1 = 0 ! Tracking loss of particles from mode k due to collision
				dNdt_coag_term_2 = 0 ! Tracking particles formed by collision of particles in Mode i with themselves which may possibly lie in Mode k
				dNdt_coag_term_3 = 0 ! Tracking particles formed by collision of particles in Mode i with Mode j (i=!j) which may possibly lie in Mode k
				
				! dV/dt coagulation terms
				dVdt_coag_term_2 = 0
				dVdt_coag_term_3 = 0

				if (k < N_modes) then ! All modes except last one
					do i = 1, N_modes
						dNdt_coag_term_1 = dNdt_coag_term_1 + (-beta(v(k),v(i))*N(i)*N(k)) ! Rate of collision of particles from mode k with other modes ! Note: There is no 1/2 factor in k-k collision as it also results in loss of 2 particles 
						dNdt_coag_term_2 = dNdt_coag_term_2 + (1.0/2.0)*beta(v(i),v(i))*N(i)*N(i)*frac(i,i,k)
					end do
				
					do i = 2, N_modes
						do j = 1, i-1
							dNdt_coag_term_3 = dNdt_coag_term_3 + beta(v(j),v(i))*N(i)*N(j)*frac(i,j,k)
						end do	
					end do
				
					dNdt_coag(k) = dNdt_coag_term_1 + dNdt_coag_term_2 + dNdt_coag_term_3
					dVdt_coag(k) = v(k)*dNdt_coag(k) ! Except for the last mode, coagulation will not result in volume change of a particle in a mode
					
				else ! Last Mode
					! dN/dt coagulation terms
					dNdt_coag_term_1 = dNdt_coag_term_1 + (-0.5*beta(v(k),v(k))*N(k)*N(k)) ! Note: There is no 1/2 factor as there is loss of only 1 particle when 2 particle in last mode collide among themselves (the newly formed particle will be in same mode)
					do i = 1, (N_modes-1)
						dNdt_coag_term_2 = dNdt_coag_term_2 + (1.0/2.0)*beta(v(i),v(i))*N(i)*N(i)*frac(i,i,k) ! Tracking particles formed by collision of particles in Mode i with themselves which may possibly lie in last mode.
					end do
					
					! when particles in mode j (j<k) collide with last mode, the resulting particle will lie in last mode. There is no change in # conc in last mode.
					if (N_modes > 2) then
						do i = 2, (N_modes-1)
							do j = 1, i-1
									dNdt_coag_term_3 = dNdt_coag_term_3 + beta(v(j),v(i))*N(i)*N(j)*frac(i,j,k) ! Tracking particles formed by collision of particles in Mode i with Mode j (i=!j) which may possibly lie in Mode k
							end do
						end do
					end if

					dNdt_coag(k) = dNdt_coag_term_1 + dNdt_coag_term_2 + dNdt_coag_term_3
					
					! dV/dt coagulation terms
					do i = 1, N_modes-1
						dVdt_coag_term_2 = dVdt_coag_term_2 + (0.5)*beta(v(i),v(i))*N(i)*N(i)*frac(i,i,k)*(v(k)) ! fracn. going to mode k is frac(i,i,k)
					end do
					
					do i = 2, N_modes
						do j = 1, i-1
							if (i == N_modes) then
								dVdt_coag_term_3 = dVdt_coag_term_3 + beta(v(j),v(i))*N(i)*N(j)*frac(i,j,k)*(v(j)) 
							else
								dVdt_coag_term_3 = dVdt_coag_term_3 + beta(v(i),v(j))*N(i)*N(j)*frac(i,j,k)*(v(k))
							end if
						end do	
					end do
					
					dVdt_coag(k) = dVdt_coag_term_2 + dVdt_coag_term_3		
				end if
			end do
		end if
  
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
            
        do k = 1, N_modes
            YDOT((2*k-1)) =  dNdt_coag(k) + dNdt_nucl(k)
            YDOT((2*k-0)) = dVdt_coag(k) + dVdt_nucl(k) + dVdt_cond(k) 
        end do 
        YDOT(NEQ) =  0 -cond_loss -dNdt_nucl(1)*v_critical/vm
        return
    end subroutine FEX
end module rhs