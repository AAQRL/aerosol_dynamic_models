! This module solves the 4 differential equations, and calculates the three moments (moment_ode.dat), 
! and superstaturation as a function of time. These values are then used for the calculation of particle
! size distribution (size_distribution.dat)
    
    
module results

      USE DVODE_F90_M
      USE rhs
      USE user_interface
      
      !IMPLICIT NONE
      !INTEGER ITASK, ISTATE, NEQ , IOUT, tfinal
      !DOUBLE PRECISION ATOL, RTOL,  TT, TOUT, Y
      !DIMENSION Y(4), ATOL(4)
      !double precision :: Nc, vg, temp, dp_t, dpmg, sigmag
      !double precision :: momo1 = 1.0, momo2 = 2.0
      !double precision :: y_ini(4)
      
contains
!    
    subroutine output(Nc_initial, dpmg_initial, sigmag_initial, vapor_initial, Nc_final, dpmg_final, sigmag_final, vapor_final,t_final_output)
      implicit none
      INTEGER ITASK, ISTATE, NEQ , IOUT, tfinal
      DOUBLE PRECISION ATOL, RTOL,  TT, TOUT, Y
      DIMENSION Y(4), ATOL(4)
      double precision :: Nc, vg, temp, dp_t, dpmg, sigmag
      double precision :: momo1 = 1.0, momo2 = 2.0
      double precision :: y_ini(4)
      
      double precision :: Nc_final, dpmg_final, sigmag_final, vapor_final, t_final_output
      double precision :: Nc_initial, dpmg_initial, sigmag_initial, vapor_initial
      TYPE (VODE_OPTS) :: OPTIONS
      



      !OPEN (UNIT=6,FILE='moment_ode.dat')
      !OPEN (UNIT=7,FILE='size_disribution.dat')
      
      num_zero = Nc_initial;
      dpg_zero = dpmg_initial;
      sigma_g_zero = sigmag_initial;
      S_zero = vapor_initial/ (Ps / (kb*T) )
      
            !call read_inputfile()
      call derived(ns, mu, lam, mm1, v1, r1, s1, sig, tau, K, vg_zero, t_start, delta_t, n_d)
      !print *, Nc_initial, dpmg_initial, sigmag_initial, vapor_initial, Nc_final, dpmg_final, sigmag_final, vapor_final,t_final_output
      NEQ = 4
      y(1) = num_zero/n_d
      y(2) = moment(num_zero,vg_zero,sigma_g_zero,momo1)/(v1*n_d)
      y(3) = moment(num_zero,vg_zero,sigma_g_zero,momo2)/(v1*v1*n_d)
      y(4) = S_zero
      y_ini(1) = num_zero/n_d
      y_ini(2) = moment(num_zero,vg_zero,sigma_g_zero,momo1)/(v1*n_d)
      y_ini(3) = moment(num_zero,vg_zero,sigma_g_zero,momo2)/(v1*v1*n_d)
      y_ini(4) = S_zero
      
      !print*, y(1), y(2), moment(num_zero,vg_zero,sigma_g_zero,1), num_zero,vg_zero,sigma_g_zero
      TT = t0
      TOUT = delta_t
      
      tfinal = (t_final_output/tau)/delta_t
      print*, tfinal
      RTOL = 1.0D-5
      ATOL(1) = 1.0D-9
      ATOL(2) = 1.0D-9
      ATOL(3) = 1.0D-9
      ATOL(4) = 1.0D-9
      ITASK = 1
      ISTATE = 1

!      write (7,500) '  t(sec)      ','  N(#/cm3)     ', '  dpg (nm)    ', '  sigmag      '
!500   format (a,a,a,a)
!      write (7,600) t0, num_zero*1.0e-6, dpg_zero*1.0e9, sigma_g_zero
!600   FORMAT (ES11.2E3,4ES11.2E3)  
      
      if (tfinal == 0) then
              Nc_final = Nc_initial * 1e-6
              dpmg_final = dpmg_initial * 1e9
              sigmag_final = sigmag_initial
              vapor_final = vapor_initial
      else

      OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR_VECTOR=ATOL)
      DO IOUT = 1, tfinal
        CALL DVODE_F90(FEX,NEQ,Y,TT,TOUT,ITASK,ISTATE,OPTIONS)
!        WRITE (6,90000) TT*tau, Y(1)*n_d, Y(2)*n_d*v1, Y(3)*n_d*v1*v1, Y(4)
!90000   FORMAT (ES11.2E3,4ES11.2E3) 


 
    
        !print >> f2, '[time(min)]',' [Number Conc ( /cm3)]', '[Geomteric mean diameter(nm)]', '[standard deviation]'
        
        Nc = n_d * Y(1)* 1e-6
        !# if else to avoid division by zero
        if ((Y(1) < 1e-40) .or. (Y(3) < 1e-40) .or. (Y(2) < 1e-40)) then
            vg = 0.0
            temp = 0.0
            dp_t = 0.0
        else
            vg = v1 * Y(2)*Y(2)/(Y(1)**(3.0/2)*Y(3)**(1.0/2))
            !pop = v1 * y_ini(1)*V_t[0]/(N_t[0]**(3.0/2)*V2_t[0]**(1.0/2))
            temp = Y(1) * Y(3) / ( Y(2) * Y(2) )
            dp_t = (vg/v1)**(1.0/3)
        end if
            
        !# check if dp_t is less than the diameter for 1 molecule
        if (dp_t < 1) then
            dp_t = 1.0
        !# check if the dp_t value is nan
        else if (dp_t /= dp_t) then
            dp_t = 0.0
        end if
                
            
        !# convert dp_t to dimensional form
        dpmg = 2.0*r1*dp_t*1e9
            
        !# check the value of temp
        if (temp >  1.0) then
            sigmag = exp(sqrt(((1.0/9)*log(temp))))
                
        else
            sigmag = 1.0
        end if
              
        !WRITE (7,90000) TT*tau, Nc, dpmg, sigmag


           
        IF (ISTATE/=2) then
        ISTATE = 2
        end if
        TOUT = TOUT + delta_t
      end do
      
      
      Nc_final = Nc
      dpmg_final = dpmg
      sigmag_final = sigmag
      vapor_final = Y(4) * ns
      
      print*, Nc,dpmg,sigmag,Y(4)*ns
      
      end if

    end subroutine output
      
end module results
    

    