! This module is a user-interface for reading different input variables from input.txt and storing them in 
! their respective variables. No calculation is done here. This just interfaces the user with the main code
module user_interface

	use global_variables
	implicit none

    
contains

	subroutine read_inputfile()

	    implicit none
	    integer :: i, j
        !double precision :: num_zero , dpg_zero, sigma_g_zero, S_zero 
		
	    open(unit = 1234, file = "input1.txt", status ="old")
        
            read(1234,*)num_zero_0 , dpg_zero_0 , sigma_g_zero_0
	        read(1234,*)Ps
	        read(1234,*)M, rho, gam
            read(1234,*)flag_COAG	
            read(1234,*)flag_COND	
            read(1234,*)flag_NUCL	
            read(1234,*)flag_REAC	
            read(1234,*)r_rate
            !read(1234,*)num_zero
            !read(1234,*)dpg_zero
            !read(1234,*)sigma_g_zero
            !read(1234,*)S_zero

            !print*, tau
            !print*,dt

	    close(unit = 1234)
    end subroutine read_inputfile
    
 	subroutine read_temperature_profile()

	    implicit none
	    integer :: i
        !double precision :: num_zero , dpg_zero, sigma_g_zero, S_zero 
        
	    
        
       
        OPEN (11, file = "temphist.txt") 
        DO 
            READ (11,*, END=10) 
            nlines = nlines + 1 
        END DO 
        10 CLOSE (11) 

        !print*, nlines
	    allocate(z_pos(nlines))
        allocate(vz_pos(nlines))
	    allocate(Temperature(nlines))
	    allocate(vapor_conc(nlines))
        allocate(time(nlines))
        
        open(unit = 12, file = "temphist.txt", status ="old")
        do i = 1,nlines
	        read(12,*)z_pos(i), vz_pos(i), Temperature(i) , vapor_conc(i)
        end do    
        
        do i = 2,nlines
            time(i-1) = abs((z_pos(i)-z_pos(i-1))*1e-2/((vz_pos(i)+vz_pos(i-1))/2))
        end do     
        
	    close(unit = 12)
    end subroutine read_temperature_profile
	
end module user_interface
    
    
!program test_user_interface
!    use user_interface 
!    implicit none
!    
!    call read_inputfile()
!    call read_temperature_profile()
!    print*, M,time, Temperature
!    
!end program test_user_interface
