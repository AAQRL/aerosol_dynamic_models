module results
    
    use global_variables
    use rhs
    use DVODE_F90_M
    
    implicit none
    double precision :: t, tout
    integer :: counter, i
    
contains
   
    subroutine output()
        
        double precision :: Ntot, Surf_tot, Vol_tot, dp_g, segma_g
        type (VODE_OPTS) :: OPTIONS
        t = t0
        tout = dt
        istate = 1
    
        open (unit=611,file='output_numb_distr.dat')
        !open (unit=511,file ='stats.dat')
    
          OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,RELERR=rtol, &
                    ABSERR_VECTOR=atol, USER_SUPPLIED_JACOBIAN=.TRUE.)
      
          !Heading of the file
          write(611,*) "dp (m)   "
          write(611,*)"t (in seconds) = ", " " , "Y (in #/m3) = "
                    
          !Write the first row 
          write(unit = 611,fmt = '(D12.4)', advance = "no") 0.0
          do i = 1, n_discrete
                write(unit = 611,fmt = '(D14.6)', advance = "no") dp_dis(i) !(i*v0*6.0/PI)**(1.0/3)      
          end do
                   
          ! Next line (because of advanced no)
          write(611,*) 
        
          ! Write the time t = t0
          write (unit = 611, fmt = '(D12.4)', advance = "no") t0 
        
          !Write the value of y at t = t0
          do i = 1, n_discrete
          
            if (YNc(i) *N_infi_0> 1e0) then
                write(unit = 611,fmt = '(D14.6)', advance = "no") YNc(i)
            else
                write(unit = 611,fmt = '(D14.6)', advance = "no") 0.0
            end if
          
          end do
          
          ! Next line
          write(611,*)

          !!write(511,*)"t = ", " " , "dg = " , " " , " sigmag = "
      
          do counter = 1, int((tf-t0)/dt)
       
            !if (mod(counter, 1) == 0) write(*,*) 'counter is = ', counter
              
            call DVODE_F90(FEX,neq,YNc,t,tout,itask,istate,OPTIONS,J_FCN=JEX)

            tout = tout + dt
                     
          end do
         
          write (unit = 611, fmt = '(D12.4)', advance = "no") t *tau
            
          do i = 1, n_discrete
            
            if (YNc(i)*N_infi_0> 1e-10) then
            
              write(unit = 611,fmt = '(D14.6)', advance = "no") YNc(i)
            
            else
            
              write(unit = 611,fmt = '(D14.6)', advance = "no") 0.0
            
            end if

          end do
                    
        close(unit = 611)
        !close(unit = 511)
        !print *, tau
        
        Ntot = 0. !total number conc. (#/cm^3)
        do i = 1, n_discrete
        
          Ntot = Ntot + YNc(i)
        
        end do
        Ntot = Ntot*N_infi_0/1e+6 !(#/cm^3)
        write (*,*)'Ntot = (#/cm^3)'
        write (*, '(D12.6)') Ntot
               
        dp_g = gmean (YNc) * 1e+9 !geometry mean diameter (nm)
        write(*,*)'dpg = '
        write(*,*)dp_g
        
        segma_g = gstddev (YNc*N_infi_0) !geometry standard deviation 
        write(*,*)'segma_g = '
        write(*,*)segma_g
               
        Surf_tot = 0. !total surface area concentration (m^2/cm^3)
        do i = 1, n_discrete
        
          Surf_tot = Surf_tot + YNc(i)*N_infi_0/1e+6*PI*dp_dis(i)**2
            
        end do
        write (*,*)'Surf_tot = '
        write (*, '(D12.6)') Surf_tot
        
        Vol_tot = 0. !total volume concentration (m^3/cm^3)
        do i = 1, n_discrete
        
          Vol_tot = Vol_tot + YNc(i)*N_infi_0/1e+6*(PI/6.)*dp_dis(i)**3
            
        end do
        write (*,*)'Vol_tot = '
        write (*, '(D12.6)') Vol_tot
        
        open (612, file='output_aver_para.dat')
          
          write (612,*)'Ntot (#/cm^3) '
          write (612,'(D12.6)') Ntot
          
          write(612,*)'dpg (nm) '
          write(612,*)dp_g
          
          write(612,*)'segma_g '
          write(612,*)segma_g
          
          write (612,*)'Surf_tot (m^2/cm^3) '
          write (612, '(D12.6)') Surf_tot
          
          write (612,*)'Vol_tot (m^2/cm^3) '
          write (612, '(D12.6)') Vol_tot
        
        close(612)
      
    end subroutine output
      
end module results