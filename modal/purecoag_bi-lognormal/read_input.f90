module input
    use global_variables
    implicit none

contains
    subroutine read_inputfile()    
        implicit none
		integer :: i
        double precision :: dump
	    open(unit = 1234, file = "input_file.txt", status ="old")
            read(1234,*)dump
            read(1234,*)p10
            read(1234,*)ps
			read(1234,*)v_gen
            read(1234,*)dump
		    read(1234,*)Mw_p
		    read(1234,*)rho_p
			read(1234,*)sigma
            read(1234,*)dump
		    read(1234,*)Mw_g
		    read(1234,*)rho_g_ref
		    read(1234,*)mu_g_ref
		    read(1234,*)S_ref
            read(1234,*)lm_ref
            read(1234,*)dump
		    read(1234,*)P
		    read(1234,*)Temp
            read(1234,*)dump
		    read(1234,*)t_final
		    read(1234,*)dt
            read(1234,*)dump
		    read(1234,*)N_modes
			read(1234,*)dump
		    read(1234,*)flag_coag
			read(1234,*)flag_cond
		    read(1234,*)flag_nucl
	    close(unit = 1234)
		
		open(unit = 2345, file = "input_N.txt", status ="old")
			do i = 1, N_modes
				read(2345,*)N_init(i)
			end do 
		close(unit = 2345)
		
		open(unit = 3456, file = "input_v.txt", status ="old")
			do i = 1, N_modes
				read(3456,*)v_init(i)
			end do 
		close(unit = 3456)

    end subroutine read_inputfile
    end module input