module user_interface

	use rhs	
	implicit none
	
	! file name : general system parameters
	character (len = 20) :: filename1, filename2, filename3

	!inital conditions

contains

	subroutine read_inputfile()

	    implicit none
	    integer :: i, j, numb_dx
        double precision :: dpg1, segmag1, dpg2, segmag2, dp_same, &
                            sum_ync, fi, delta_dp, dp1, dp2
        
        dpg1 = 1e-9
        segmag1 = 1.2
        dpg2 = 3e-9
        segmag2 = 1.2
        dp_same = 1.7321e-9 !the same diameter for bi-lognormal distribution
        
        numb_dx = 100
        
	    open(unit = 1234, file = "input.txt", status ="old")
	
	        !reading info from file
            read(1234,*)flag
	        read(1234,*)n_discrete	
	        read(1234,*)mw, rho, dp0
            v0 = (pi/6)*(dp0)**3
	        
            !allocate(nc0(n_discrete))
	        !allocate(nc(n_discrete))
	        allocate(ync(n_discrete))
            allocate(dp_dis(n_discrete)) !Huang Zhang, 2018
		    allocate(atol(n_discrete))
        	
		    neq = n_discrete
        
            do i = 1, n_discrete	
            
              dp_dis(i) = (i*v0*6.0/PI)**(1.0/3)
              if (dp_dis(i) > dp_same) then 
                
                fi = get_lognormal(dp_dis(i)*1e+9, dpg2*1e+9, segmag2)
                dp1 = dp_dis(i)*1e+9 !nm
                dp2 = ((i+1)*v0*6.0/PI)**(1.0/3)*1e+9 !nm
                ync(i) = get_prob_lognormal(dp1, dp2, numb_dx, dpg2*1e+9, segmag2)
                
              else
              
                fi = get_lognormal(dp_dis(i)*1e+9, dpg1*1e+9, segmag1)
                dp1 = dp_dis(i)*1e+9 !nm
                dp2 = ((i+1)*v0*6.0/PI)**(1.0/3)*1e+9 !nm
                ync(i) = get_prob_lognormal(dp1, dp2, numb_dx, dpg1*1e+9, segmag1)
              
              end if
                            
            end do
            !ync(1:n_discrete) = 0.0
            !ync(1) = 1.0
            sum_ync = sum(ync)
                        
            write(*,*)'dp_dis(end) = (nm)'
            write(*, '(D12.4)') dp_dis(n_discrete)*1e+9
	    
		    read(1234,*)N_infi_0
            read(1234,*) tf
            read(1234,*) dt
            read(1234,*) N_track_min
                
            atol(1:n_discrete) = 0.1*N_track_min/N_infi_0
            !rtol = 10**-(log10(N_infi_0/N_track_min) + 1)
            rtol = 10**(-(log10(N_infi_0/N_track_min) + 1))
            
	    close(unit = 1234)

    end subroutine read_inputfile
	
    function get_lognormal(x, x0, sg0)
      
      implicit none
      
      double precision :: get_lognormal
      double precision, intent (in) :: x, x0, sg0
      get_lognormal = 1/((2*PI)**0.5*x*log(sg0))* &
                      exp(-(log(x)-log(x0))**2/(2*(log(sg0))**2))
      return
    
    end function get_lognormal
   
    function get_prob_lognormal(x_begin, x_end, numb_dx, x0, sg0)
    
      implicit none 
      integer :: i_dx
      integer, intent (in) :: numb_dx
      double precision, intent (in) :: x_begin, x_end, x0, sg0
      double precision :: get_prob_lognormal, x1, x2, dx
      dx = (x_end - x_begin)/numb_dx
      
      get_prob_lognormal = 0.
      x1 = x_begin
      do i_dx = 1, numb_dx
        
        x2 = x1 + dx
        get_prob_lognormal = get_prob_lognormal + (get_lognormal(x2, x0, sg0) + get_lognormal(x1, x0, sg0))*dx/2
        x1 = x2
          
      end do
      
      return    
    
    end function get_prob_lognormal
    
end module user_interface
    
    
!program test_user_interface
!    use user_interface 
!    implicit none
!    
!    call read_inputfile()
!    print*, nc0
!    
!end program test_user_interface

    

	

	

	
