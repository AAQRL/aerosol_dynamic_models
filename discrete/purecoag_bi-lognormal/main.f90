program main
    
    use user_interface
    use results
    !use monodisperse
    use beta_safe
        
    implicit none
    
        double precision :: t1, t2
            !call KK()

            !PRINT*, TAU
            call CPU_TIME(t1)
    
            call read_inputfile()
            
            !allocate(beta_all(n_discrete,n_discrete))
            !call fill_beta()
            call KK()
            call output()
            
            !call print_monodisperse()
            
            !print* , tau
            !deallocate(charge)
            !deallocate(Nc0)
            !deallocate(Nc)
            deallocate(YNc)
            deallocate(atol)
            !deallocate(T_data)
            !deallocate(YNc_mono)
            !deallocate(beta_all)
            
            call CPU_TIME( t2 )
            
            open (613, file='output_CPU&Memory.dat')
            
              write(613,*) 'Elapsed CPU time (s)'
              write(613, '(D12.6)') t2 - t1
              
            close (613)
            !print*, 'Elapsed CPU time = ', t2 - t1
            


end program main
   