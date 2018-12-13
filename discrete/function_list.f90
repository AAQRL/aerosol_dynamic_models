module function_list
    
    use global_variables
    implicit none
    
contains
    
    !subroutine T_hist()
    !    implicit none
    !    
    !    character (len = 20) :: filename7
    !    integer :: i
    ! 
    !    print *, " enter the Temperature file name :"
	   ! read *, filename7
    !    
    !    open(unit = 7, file = filename7, status ="old")
    !    read(7,*) time_points
    !    allocate(T_data(time_points,2))
    !    
    !    do i = 1 , time_points
    !        read(7,*)T_data(i,1:2) 
    !    end do
    !end subroutine T_hist
    
! VARIABLE TEMPERATURE FUNCTION 
  !  function Temp(t)
  !      implicit none
  !      double precision :: t, Temp 
  !      integer :: len_a, i
  !      Temp = 0.0
  !      
  !      len_a = size(T_data(:,1))
  !      
  !      do i = 1, len_a -1
  !          if ((t > T_data(i,1)) .and. (t <= T_data(i+1,1))) then
  !              Temp = T_data(i,2) + (T_data(i+1,2)-T_data(i,2))/  &
		!(T_data(i+1,1)-T_data(i,1)) * (t -T_data(i,1))
  !          end if
  !      end do
  !      
  !      !if (t <= T_data(1,1)) then
  !      !    Temp = T_data(1,2)
  !      !endif
  !      !
  !      !if (t > T_data(len_a,1)) then
  !      !    Temp = T_data(len_a,2)
  !      !endif
  !      !print*, (t <= T_data(1,1))
  !      if (t <= T_data(1,1)) then
  !          Temp = Temp + T_data(1,2) 
  !      end if
  !      
  !      if (t > T_data(len_a,1)) then
		!    Temp = Temp + T_data(len_a,2)
  !      end if
  !
  !  end function Temp

    
    function gmean(array)
        
      implicit none
      
      double precision, dimension(:) :: array
      integer :: i
      double precision :: gmean, dp, sum_array
      
      gmean = 1.0
      sum_array = sum(array)
    
      do i = 1, size(array)
      
        !dp = (6.0*i*v0/PI)**(1.0/3)
        dp = dp_dis(i)
        !gmean = gmean*dp**(array(i)/sum(array))
        gmean = gmean*dp**(array(i)/sum_array)
        
      end do
    
    end function gmean
    
    function gstddev(array)
      
      implicit none
      
      double precision, dimension(:) :: array
      integer :: i
      double precision :: gstddev, dp, dp_gmean, sum_array
      
      gstddev = 0.0
      dp_gmean = gmean(array)      
      sum_array = sum(array)
    
      do i = 1, size(array)
                
        dp = dp_dis(i)
        !dp = (6.0*i*v0/PI)**(1.0/3)
        gstddev = gstddev + (array(i)*(log(dp) - log(dp_gmean))**2)/(sum_array - 1)
        !gstddev = gstddev + (array(i)*(log(dp) - log(gmean(array)))**2)/(sum(array) - 1)
                
      end do
        
      gstddev = exp(gstddev**(1.0/2))
        
    end function gstddev
    
    
!    subroutine calc_beta(q1,q2,v1,v2,beta,flag)
!        integer, intent(in) :: q1,q2 
!        double precision , intent(in) :: v1, v2
!        double precision, intent(inout) :: beta
!        integer, intent(in) :: flag
!        
!        double precision :: dp0, dp1, dp2, c1, c2, l1, l2, g1, g2, & 
!   				coeff1, coeff2, y,  W, t, beta_pre
!        
!        dp0 = (6.0/PI * v0)**(1.0/3)
!        dp1 = (6.0/PI * v1)**(1.0/3)
!        dp2 = (6.0/PI * v2)**(1.0/3)
!        
!    
!        select case (flag)
!            
!        case(0)
!        ! Use Fuchs Sutugin Beta
!            c1 = sqrt(8.0*KB*Temp/(PI*(mw/NA)))
!            c2 = sqrt(8.0*KB*Temp/(PI*(mw/NA)))
!            
!            l1 = 8.0*Diff(dp1)/(PI*c1)
!            l2 = 8.0*Diff(dp2)/(PI*c2)
!            
!            g1 = 1.0/(3.0*dp1*l1) * ((dp1+l1)**3-(dp1*dp1+l1*l1)**  &
!			(3.0/2)) - dp1
!            g2 = 1.0/(3.0*dp2*l2) * ((dp2+l2)**3-(dp2*dp2+l2*l2)**  &
!			(3.0/2)) - dp2
!    
!    
!            coeff1 = 2.0*PI*(dp1 + dp2) * (Diff(dp1) + Diff(dp2))
!            coeff2 = (dp1 + dp2)/(dp1 + dp2 + 2.0*sqrt(g1*g1 +   &
!			g2*g2)) + &
!		(8.0*(Diff(dp1)+Diff(dp2)))/((dp1+dp2)*sqrt(c1*c1+c2*c2))
!            
!            if ((dp1 < dp0) .or. (dp2 < dp0)) then
!                beta = 0.0
!            else 
!                if ((q1 == 0) .or. (q2 == 0)) then
!                    W = 1.0
!                else 
!                    y = KE * EE* EE * q1 * q2/ (KB * Temp* &
!					(dp1 + dp2)/2)
!                    W = 1.0/y * (exp(y) - 1)
!                end if
!                beta = (1.0/W) * (coeff1/coeff2)
!!                beta = beta/K_coag
!            end if
!            
!    
!        case(1)
!        ! Use FM Regime Beta
!            beta_pre = (3.0/(4.0*PI))**(1.0/6) * (6.0*KB*Temp/rho)&
!			**(1.0/2) * &
!		(1.0/v1 + 1.0/v2)**(1.0/2) * (v1**(1.0/3) + v2**(1.0/3))**(2.0) 
!            if ((dp1 < dp0) .or. (dp2 < dp0)) then
!                beta = 0.0
!            else 
!                if ((q1 == 0) .or. (q2 == 0)) then
!                    W = 1.0
!                else 
!                    y = KE * EE* EE * q1 * q2/ (KB * Temp * &
!					(dp1 + dp2)/2)
!                    W = 1.0/y * (exp(y) - 1)
!                end if
!                beta =  beta_pre * (1.0/W)
!!                beta = beta/K_coag
!            end if
!            
!        case(2)
!            beta = 1.0
!            
!        end select
!    end subroutine calc_beta
!        
    end module function_list
    
    
!program test_function_list
!    use function_list  
!    implicit none
!    double precision :: me, t
!    double precision :: beta 
!    !t = 2.5
!    print*, v0
!    call T_hist()
!    print*, T_data
!    print*, " Enter the time "
!    read*, t
!    print*, Temp(t)
!    !call calc_beta(0,0,v0,v0,beta,0)
!    !print*, beta
!    
!end program test_function_list
!    
    

    