!gfortran schrodinger_solution_methods.f90 -c may need to be ran
module schrodinger_solution_mod
    use numerov
    implicit none 
    private 
    public produce_trial_solution,zero_potential
contains

    function zero_potential(x)
        implicit none
        real*8, intent(in)::x
        real*8::zero_potential
        zero_potential = 0.0
    end function zero_potential

    subroutine produce_trial_solution(x0,x1,x_final,psi_0,psi_1,E,x_array,psi_array,V_ptr,N)
        implicit none
        real*8,intent(in)::x0,x1,x_final,psi_0,psi_1,E
        integer,intent(in)::N
        real*8,intent(inout)::x_array(N),psi_array(N)
        real*8::h

        procedure (pointing_func),pointer:: V_ptr

        h = x1-x0
        !N = (x_final-x0)/h
        call numerov_whole_interval_schrodinger(N,E,x0,x1,psi_0,psi_1,V_ptr,x_array,psi_array)


    end subroutine produce_trial_solution

end module schrodinger_solution_mod


program run_solver
    use schrodinger_solution_mod
    use numerov
    implicit none
    
    real*8::x0,x1,h,y0,y1
    real*8:: x_final,E
    integer::n
    real*8,dimension(:),allocatable :: x_array
    real*8,dimension(:),allocatable :: y_array
    procedure (pointing_func),pointer:: V_ptr => zero_potential

    E = 1.0
    y0 = 1.0
    y1 = 0.999999995
    x0 = 0.0
    x1 = 0.0001
    x_final = 1.0
    h = x1-x0
    n = (x_final-x0)/h
    allocate(x_array(n))
    allocate(y_array(n))

    call produce_trial_solution(x0,x1,x_final,y0,y1,E,x_array,y_array,V_ptr,n)
    print*, y_array(n-1)

    deallocate(x_array)
    deallocate(y_array)

end program run_solver




