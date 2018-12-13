! This is the main program
program main
    
    use global_variables
    use user_interface
    use results

    implicit none
    double precision :: Nc_final, dpmg_final, sigmag_final, vapor_final, t_final_output
    double precision :: Nc_initial, dpmg_initial, sigmag_initial, vapor_initial
    !double precision :: num_zero_0 = 1e11, dpg_zero_0 = 10.0e-9, sigma_g_zero_0 = 1.5, vapor_zero_0 = 1.85717257286e+16
    !double precision :: num_zero_0 = 0e0, dpg_zero_0 = 0.0e-9, sigma_g_zero_0 = 1.0, vapor_zero_0 = 1.85717257286e+18
    integer :: i
    
        
        call  read_inputfile()   ! To read the inputs
        call  read_temperature_profile()
        
        
        Nc_initial = num_zero_0
        dpmg_initial = dpg_zero_0
        sigmag_initial = sigma_g_zero_0 
        
        OPEN (UNIT=10,FILE='output.txt')        
90010   FORMAT (ES11.2E3,5ES11.2E3) 
        WRITE (10,90010) t0, Temperature(1), Nc_initial*1e-6, dpmg_initial*1e9, sigmag_initial, vapor_initial
        do i = 1,nlines-1
            t_final_output = time(i)
            T = Temperature(i)+273.0
            vapor_initial = vapor_conc(i)
            call  output(Nc_initial , dpmg_initial , sigmag_initial, vapor_initial, Nc_final, dpmg_final, sigmag_final, vapor_final,t_final_output)
            !call output(Nc_final, dpmg_final, sigmag_final, vapor_final, t_final_output)          ! To perform all the calculations, and write the results
            !print*, t_final_output, Nc_final, dpmg_final, sigmag_final, vapor_final
 
            if (Nc_final < 0.0) then 
                Nc_final =0.0 
            end if

            if (dpmg_final < 0.0) then 
                dpmg_final =0.0 
            end if
            
            if (sigmag_final < 1.0) then 
                sigmag_final =1.0 
            end if
            
            if (vapor_final < 0.0) then 
                vapor_final =0.0 
            end if
            
            Nc_initial = Nc_final*1e6
            dpmg_initial = dpmg_final*1e-9
            sigmag_initial = sigmag_final
            
            
            !print*,time(i), T, Nc_final, dpmg_final, sigmag_final, vapor_final
                        
            WRITE (10,90010) z_pos(i+1), T, Nc_final, dpmg_final, sigmag_final, vapor_final

        end do
        close(unit = 10)
end program main