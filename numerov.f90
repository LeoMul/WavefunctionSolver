module numerov
    !integer, parameter :: dp = selected_real_kind(p = 15, r = 307)
    implicit none

    private 
    public numerov_next_step,pointing_func,numerov_whole_interval,numerov_whole_interval_schrodinger

    abstract interface 
        function pointing_func(x)
            real*8,intent (in)::x
            real*8::pointing_func
        end function pointing_func
    end interface

contains
    function numerov_next_step(g_func_ptr,s_func_ptr,h,xcurrent,ycurrent,yprev)
        implicit none
        real*8 :: numerov_next_step
        real*8 :: factor
        real*8 :: hsquared_over_twelve
        real*8, intent (in) :: h,xcurrent,ycurrent,yprev
        procedure(pointing_func),pointer::g_func_ptr
        procedure(pointing_func),pointer::s_func_ptr
        hsquared_over_twelve = h*h/12.0
        factor = (1+hsquared_over_twelve*g_func_ptr(xcurrent+h))
        numerov_next_step = 2*ycurrent*(1-5*hsquared_over_twelve*g_func_ptr(xcurrent))
        numerov_next_step = numerov_next_step -yprev *(1+hsquared_over_twelve*g_func_ptr(xcurrent-h))
        numerov_next_step = numerov_next_step + hsquared_over_twelve*(s_func_ptr(xcurrent+h)+10*s_func_ptr(xcurrent)+s_func_ptr(xcurrent-h))
        numerov_next_step = numerov_next_step/factor
    end function numerov_next_step

    function numerov_next_step_g_s_pre_calculated(h,ycurrent,yprev,gnext,gcurrent,gprev,snext,scurrent,sprev)
        real*8, intent (in) ::h,ycurrent,yprev,gnext,gcurrent,gprev,snext,scurrent,sprev
        real*8 :: hsquared_over_twelve,factor,numerov_next_step_g_s_pre_calculated
        hsquared_over_twelve = h*h/12.0
        factor = (1.0+hsquared_over_twelve*gnext)
        numerov_next_step_g_s_pre_calculated = 2.0*ycurrent*(1.0-5.0*hsquared_over_twelve*gcurrent)
        numerov_next_step_g_s_pre_calculated = numerov_next_step_g_s_pre_calculated -yprev *(1.0+hsquared_over_twelve*gprev)
        numerov_next_step_g_s_pre_calculated = numerov_next_step_g_s_pre_calculated + hsquared_over_twelve*(snext+10.0*scurrent+sprev)
        numerov_next_step_g_s_pre_calculated = numerov_next_step_g_s_pre_calculated/factor
    end function numerov_next_step_g_s_pre_calculated

    subroutine numerov_whole_interval(n,x_0,x_1,y_0,y_1,g_func_ptr,s_func_ptr,x_array,y_array)

        real*8, intent (in) :: x_0,x_1,y_0,y_1
        integer, intent(in) :: n
        procedure(pointing_func),pointer::g_func_ptr
        procedure(pointing_func),pointer::s_func_ptr
        real*8::g_array(n),s_array(n)
        real*8,intent(inout):: x_array(n), y_array(n)
        real*8::h
        integer::i

        h = x_1-x_0

        do i = 1,n 
            x_array(i) = x_0 + (i-1)*h 
            g_array(i) = g_func_ptr(x_array(i))
            s_array(i) = s_func_ptr(x_array(i))
        end do

        y_array(1) = y_0
        y_array(2) = y_1

        do i = 2,n-1
            y_array(i+1) = numerov_next_step_g_s_pre_calculated(h,y_array(i),y_array(i-1),g_array(i+1),g_array(i),g_array(i-1),s_array(i+1),s_array(i),s_array(i-1))
        end do
    end subroutine numerov_whole_interval

    function numerov_next_step_s_is_zero(g_func_ptr,h,xcurrent,ycurrent,yprev)
        implicit none
        real*8 :: numerov_next_step_s_is_zero
        real*8 :: factor
        real*8 :: hsquared_over_twelve
        real*8, intent (in) :: h,xcurrent,ycurrent,yprev
        procedure(pointing_func),pointer::g_func_ptr
        hsquared_over_twelve = h*h/12.0
        factor = (1+hsquared_over_twelve*g_func_ptr(xcurrent+h))
        numerov_next_step_s_is_zero = 2*ycurrent*(1-5*hsquared_over_twelve*g_func_ptr(xcurrent))
        numerov_next_step_s_is_zero = numerov_next_step_s_is_zero -yprev *(1+hsquared_over_twelve*g_func_ptr(xcurrent-h))
        numerov_next_step_s_is_zero = numerov_next_step_s_is_zero/factor
    end function numerov_next_step_s_is_zero

    function numerov_next_step_g_s_pre_calculated_s_is_zero(hsquared_over_twelve,ycurrent,yprev,gnext,gcurrent,gprev)
        real*8, intent (in) ::hsquared_over_twelve,ycurrent,yprev,gnext,gcurrent,gprev
        real*8 :: factor,numerov_next_step_g_s_pre_calculated_s_is_zero
        !hsquared_over_twelve = h*h/12.0
        factor = (1.0+hsquared_over_twelve*gnext)
        numerov_next_step_g_s_pre_calculated_s_is_zero = 2.0*ycurrent*(1.0-5.0*hsquared_over_twelve*gcurrent)
        numerov_next_step_g_s_pre_calculated_s_is_zero = numerov_next_step_g_s_pre_calculated_s_is_zero -yprev *(1.0+hsquared_over_twelve*gprev)
        numerov_next_step_g_s_pre_calculated_s_is_zero = numerov_next_step_g_s_pre_calculated_s_is_zero/factor
    end function numerov_next_step_g_s_pre_calculated_s_is_zero

    subroutine numerov_whole_interval_s_is_zero(n,x_0,x_1,y_0,y_1,g_func_ptr,x_array,y_array)

        real*8, intent (in) :: x_0,x_1,y_0,y_1
        integer, intent(in) :: n
        procedure(pointing_func),pointer::g_func_ptr
        real*8::g_array(n)
        real*8,intent(inout):: x_array(n), y_array(n)
        real*8::h
        integer::i

        h = x_1-x_0

        do i = 1,n 
            x_array(i) = x_0 + (i-1)*h 
            g_array(i) = g_func_ptr(x_array(i))
        end do

        y_array(1) = y_0
        y_array(2) = y_1

        do i = 2,n-1
            y_array(i+1) = numerov_next_step_g_s_pre_calculated_s_is_zero(h*h/12,y_array(i),y_array(i-1),g_array(i+1),g_array(i),g_array(i-1))
        end do
    end subroutine numerov_whole_interval_s_is_zero
    
    function numerov_whole_interval_schrodinger(x_array,h,y_0,y_1,potential_array,E)
        real*8, intent(in) :: x_array(:),potential_array(:),E,h,y_0,y_1
        real*8:: numerov_whole_interval_schrodinger(size(x_array)),energy_minus_potential_array(size(potential_array)),hsquared_over_twelve
        integer::n,i
        
        n = size(x_array)
        numerov_whole_interval_schrodinger(1) = y_0
        numerov_whole_interval_schrodinger(2) = y_1
        energy_minus_potential_array = E - potential_array
        hsquared_over_twelve = h*h/12
        do i = 2,n-1
            numerov_whole_interval_schrodinger(i+1) = numerov_next_step_g_s_pre_calculated_s_is_zero(hsquared_over_twelve,numerov_whole_interval_schrodinger(i),numerov_whole_interval_schrodinger(i-1),energy_minus_potential_array(i+1),energy_minus_potential_array(i),energy_minus_potential_array(i-1))
        end do

    end function numerov_whole_interval_schrodinger


end module numerov

module functions_storage
    use numerov
    implicit none
    public zero_func,constant_func
contains
    function zero_func(x)
            real*8, intent(in)::x
            real*8::zero_func
            zero_func = 0.0
    end function zero_func

    function constant_func(x)
        real*8, intent(in)::x
        real*8::constant_func
        constant_func = 1.0
    end function constant_func

end module functions_storage



!program test_numerov
!    use numerov
!    use functions_storage
!    implicit none
!    real*8::x0,x1,h,y0,y1
!    integer::n
!    real*8,dimension(:),allocatable :: x_array
!    real*8,dimension(:),allocatable :: y_array
!    character(len=*), parameter :: OUT_FILE = 'numerovdata.dat' ! Output file.
!    character(len=*), parameter :: PLT_FILE = 'plot.plt' ! Gnuplot file.
!    integer::i
!
!    
!    procedure (pointing_func),pointer:: g_ptr => constant_func
!    procedure (pointing_func),pointer:: s_ptr => zero_func
!    
!    
!    n = 100000
!    allocate(x_array(n))
!    allocate(y_array(n))
!    !initial conditions and setting x 
!    y0 = 1.0
!    y1 = 0.999999995
!    x0 = 0.0
!    x1 = 0.0001
!    h = x1-x0
!    call numerov_whole_interval(n,x0,x1,y0,y1,g_ptr,s_ptr,x_array,y_array)
!
!    open (1, file = "numerovdata.dat")
!    do i = 1,size(x_array)
!        write(1,*) x_array(i), y_array(i)
!    end do 
!    close(1)
!
!    call execute_command_line('gnuplot -p ' // PLT_FILE)
!
!    deallocate (x_array)
!    deallocate (y_array)
!
!end program test_numerov
!
