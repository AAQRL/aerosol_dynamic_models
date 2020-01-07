module calculations
    use input
    use global_variables
    implicit none 
    contains

    subroutine calc             
        implicit none
        double precision :: sum1
        integer :: i,j,k
        
        allocate(N_result(N_modes,nlines))
        allocate(v_result(N_modes,nlines))
        allocate(t_result(nlines))
        allocate(N_total(nlines))
        allocate(dp_g(nlines))
        allocate(v_g(nlines))
        allocate(sigma_g(nlines))
        allocate(M(2*(N_modes+1)))

        
        OPEN(unit = 15, file = "results_raw.dat")
		    DO i = 1,nlines
			    READ(15,*)M
                t_result(i) = M(1)
                DO j = 1,N_modes
                    N_result(j,i) = M(2*j);
                    v_result(j,i) = M(2*j+1)/M(2*j);
                END DO
                
                !!! Geometric Mean Diameter Calculation
                sum1 = 0
                N_total(i) = 0
                DO k = 1,N_modes
                    sum1 = sum1 + N_result(k,i)*log(v_result(k,i))
                    N_total(i) = N_total(i) + N_result(k,i)  
                END DO                
                
                v_g(i) = exp((1/N_total(i))*(sum1))
                dp_g(i) = (6*v_g(i)/PI)**(1.0/3)
                
                !!! Geometric Std. Deviation Calculation
                sum1 = 0
                DO k = 1,N_modes
                    sum1 = sum1 + N_result(k,i)*((log((6*v_result(k,i)/PI)**(1.0/3))-log(dp_g(i)))**2.0)/(N_total(i)-1)
                END DO 
                sigma_g(i) = exp(sum1**(0.5))
		    END DO
		CLOSE (unit = 15)
        
        
        OPEN(UNIT=7, FILE = 'results.dat')
         WRITE(7,*)"  t (s)       ", "N_total (#/m3)", "  dp_g (m)    ", "  sigma_g"
            DO i = 1,nlines
                WRITE(7,67) time(i), N_total(i), dp_g(i), sigma_g(i)
            END DO
        67  FORMAT(100E14.6)
        

    end subroutine calc
end module calculations