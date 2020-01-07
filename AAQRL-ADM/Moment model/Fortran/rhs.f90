! This module caluclates the right hand side by combining the terms for coagulation, condensation, nucleation,
! and reaction.
module rhs 
	
	use global_variables
    use derived_variables
    use function_list
	implicit none	

contains
    subroutine FEX(NEQ,time,Y,YDOT)
        implicit none
        integer, intent (in) :: NEQ
        double precision, intent (in) :: time
        double precision, intent (in) :: Y(NEQ)
        double precision, intent (out) :: YDOT(NEQ)

        double precision :: M0, M1, M2, S, momo1 = 1.0, momo2 = 2.0, momo2by3 = 2.0/3, momo5by3 = 5.0/3, momo8by3 = 8.0/3
        double precision :: vg, rg, linsig, sd
        
        double precision, dimension(2) :: coa
        double precision, dimension(3) :: con
        double precision, dimension(2) :: nuc
        double precision :: rea
             
        double precision :: dxCOA_N, dxCOA_W
        double precision :: dxCON_V, dxCON_W, dxCON_S, CON_cons
        double precision :: dxNUC_N, dxNUC_V, dxNUC_W, dxNUC_S
        double precision :: dxREA_S, dxREA_N, dxREA_V, dxREA_W  
        
        !double precision :: ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, t_final, delta_t, n_d
        
        M0 = y(1)*n_d
        M1 = y(2)*n_d*v1
        M2 = y(3)*n_d*v1*v1
        S  = y(4)
    

    
        if ((y(1) > 1e-50) .and. (y(2) > 1e-50) .and. (y(3) > 1e-50)) then
            vg = M1*M1 / ( M0 ** (1.5) * M2 ** (0.5))
            if (vg > 0) then
                rg = (3.0*vg/(4.0*PI))**(1.0/3)
                linsig = 1.0/9 * log(M0*M2/(M1*M1))
                if (linsig > 0) then
                    sd = exp(linsig**0.5)
                else 
                    sd = 1.0
                end if
            else
                rg = 0.0
                linsig = 0.0
                sd = 1.0
            end if
        else
            vg = 0.0
            rg = 0.0
            linsig = 0.0
            sd = 1.0
        end if
    
        ! Call all the functions that are required here and store the values in appropriate variables
        ! unpack
        coa     =   coag(sd,rg,flag_COAG)
        con     =   cond(M0,sd,vg,flag_COND)
        nuc     =   nucl(S,flag_NUCL)
        rea     =   reac(flag_REAC)

        !######################################## NUCLEATION ############################################ 
        dxNUC_S =   nuc(2) / ns * nuc(1) * tau
        dxNUC_N =   nuc(2) / n_d * tau
        dxNUC_V =   nuc(2) / n_d * nuc(1) * tau
        dxNUC_W =   nuc(2) / n_d * nuc(1) * nuc(1) * tau
        
    
    
        !######################################## REACTION ############################################ 
        dxREA_S =   rea / (ns / tau)
        dxREA_N =   0.0
        dxREA_V =   0.0
        dxREA_W =   0.0
    
    
        !######################################## CONDENSATION ############################################ 
        if ((y(1) > 1e-50) .and. (y(2) > 1e-50)) then
            if ((flag_COND == 0) .or. (flag_COND == 1)) then
                dxCON_S =    con(2) * y(1) * ( S - 1.0 ) * tau /v1       
                dxCON_V =    con(1) * y(1) * ( S - 1.0 ) * tau /v1 * ns/n_d
                dxCON_W =    con(3) * y(2) * ( S - 1.0 ) * tau /v1 * ns/n_d
            elseif (flag_COND == 2) then
                CON_cons =   1e9 * exp(-15155.15/T) * (22.0/7)**(1.0/3) * (6.0)**(2.0/3) 
                dxCON_S =    CON_cons * S * tau  * moment(M0,vg,sd,momo2by3) / (ns * v1)  
                dxCON_V =    CON_cons * S * tau  * moment(M0,vg,sd,momo2by3) / (n_d * v1)
                dxCON_W =    2.0 * CON_cons * S * tau  * moment(M0,vg,sd,momo5by3) / (n_d * v1 * v1)              

            end if
            
        else
            dxCON_S =    0.0       
            dxCON_V =    0.0
            dxCON_W =    0.0
        end if
    
        !######################################## COAGULATION ################################################## 
        if ((y(1) > 1e-50) .and. (y(2) > 1e-50) .and. (y(3) > 1e-50)) then
            dxCOA_N =    coa(1) * y(1) ** 2.0 * n_d * tau
            dxCOA_W =    coa(2) * y(2) ** 2.0 * n_d * tau
        else
            dxCOA_N =    0.0
            dxCOA_W =    0.0
        end if
    
    

        !############ Gathering all the phenomenon together into f(t,y) #########################################
        ydot(1) = dxNUC_N - dxCOA_N
        ydot(2) = dxNUC_V + dxCON_V
        ydot(3) = dxNUC_W + dxCON_W + dxCOA_W
        ydot(4) = 0.0 !dxREA_S - dxNUC_S - dxCON_S
        
        return
    end subroutine FEX
    
end module rhs

!program test_rhs
!    use rhs
!    implicit none
!
!    integer :: NEQ = 4
!    double precision, dimension(4) :: y_in, y_dot_out
!    double precision :: time_in = 0.0
!    call derived(ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, t_final, delta_t, n_d)
!    
!    y_in(1) = num_zero/n_d
!    y_in(2) = moment(num_zero,vg_zero,sigma_g_zero,1)/(v1*n_d)
!    y_in(3) = moment(num_zero,vg_zero,sigma_g_zero,2)/(v1*v1*n_d)
!    y_in(4) = S_zero
!    call FEX(NEQ, time_in, y_in, y_dot_out)
!    print*, y_dot_out
!    print*, tau
!    
!    
!end program test_rhs