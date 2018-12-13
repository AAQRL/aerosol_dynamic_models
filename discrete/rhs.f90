module rhs 
	
	use global_variables
    use function_list
    use beta_safe
	implicit none

contains

    function gain(ndiscrete)
	
	    integer :: ndiscrete
        double precision :: gain
        integer :: i
        gain = 0.0
                
	    do i = 1,ndiscrete
            
            if (( ndiscrete -i > 0 ) .and. (ndiscrete-i < ndiscrete)) then
		      
              gain = gain + 1.0/2*get_beta_ij(i,ndiscrete-i)*Nc(i)*Nc(ndiscrete-i)
              !gain = gain + 1.0/2*beta_all(i,ndiscrete-i)*Nc(i)*Nc(ndiscrete-i)
               
            end if
        
        end do

    end function gain
        
	function loss(ndiscrete)
	
	    integer :: ndiscrete
        double precision :: loss
        integer :: i
        loss = 0.0

	    do i = 1,n_discrete
            
            loss = loss + get_beta_ij(i,ndiscrete)*Nc(ndiscrete)*Nc(i)
			!loss = loss + beta_all(i,ndiscrete)*Nc(ndiscrete)*Nc(i)

	    end do

    end function loss
        
    subroutine FEX(NEQ,T,Y,YDOT)
        implicit none
        integer, intent (in) :: NEQ
        double precision, intent (in) :: T
        double precision, intent (in) :: Y(NEQ)
        double precision, intent (out) :: YDOT(NEQ)
	    integer :: i, j

		Nc = Y
	  
        do j= 1,n_discrete
		
            YDOT(j) =  (2.0/(K_coag))*(gain(j) - loss(j))
        
        end do

    end subroutine FEX
        
    subroutine JEX(NEQ,T,Y,ML,MU,PD,NRPD)
        implicit none
        integer, intent (in) :: NEQ, ML, MU, NRPD
        double precision, intent (in) :: T
        double precision, intent (in) :: y(NEQ)
        double precision, intent (out) :: pd(NRPD,NEQ)
        !double precision, dimension(:,:), allocatable :: pd
        
        integer :: k,j,counter
        
        !pd_gain = 0.0
        !pd_loss = 0.0
        !do k = 1, NEQ
        !    do j = 1, k-1
        !        pd_gain(k,j) = beta_all(j,k-j)*y(j)*y(k-j)
        !    end do
        !end do
        !
        !do k = 1,NEQ
        !    do j = 1, NEQ
        !    pd_loss (k,j) = beta_all(k,j)*y(k)
        !    if (k==j) then
        !        do counter = 1,NEQ
        !            pd_loss(k,j) = pd_loss(k,j) + beta_all(k,j)*y(j)
        !        end do
        !    end if
        !    end do
        !end do
        !
        !
        !
        !pd = (2.0/K_coag)*(pd_gain-pd_loss)
        
        !allocate(pd(n_discrete,n_discrete))
        
        do k = 1,NEQ
            do j = 1, NEQ
                if (j<k) then
                pd(k,j) = get_beta_ij(j,k-j)*y(j)*y(k-j)
                end if
                pd(k,j) = -get_beta_ij(k,j)*y(k)
            if (k==j) then
                do counter = 1,NEQ
                    pd(k,j) = pd(k,j) - get_beta_ij(k,j)*y(j)
                end do
            end if
            end do
        end do
                       
        pd = (2.0/K_coag)*pd
           
     end subroutine JEX
end module rhs
    
!program test_rhs
!    use user_interface
!    implicit none
!    
!    integer :: NEQa
!    double precision :: Ta
!    double precision, dimension(2) :: Ya, YDOTa
!    
!    NEQa = 2
!    Ta = 1.0
!    Ya(1) = 1.0
!    Ya(2) = 2.0
!
!    call read_inputfile()
!    
!    !print*, II(1,1)
!    !print*, loss(1,1), "  ", loss(2,1)
!    !print*, gain(1,1), "  ",  gain(2,1)
!    print*, YDOTa
!    call FEX(NEQa,Ta,Ya,YDOTa)
!    print*, YDOTa
!    
!end program test_rhs