!gfortran schrodinger_solution_methods.f90 -c may need to be ran
module schrodinger_solution_mod
    use numerov
    implicit none 
    private 
    public produce_trial_solution,lennard_jones_potential,zero_potential,normalise_wave_function,create_potential_array,energy_minus_potential_array,effective_hydrogen_potential,find_bracketing_pair
contains


    function zero_potential(x)
        implicit none
        real*8, intent(in)::x
        real*8::zero_potential
        zero_potential = 0.0
    end function zero_potential

    function lennard_jones_potential(r) result(V)
        implicit none
        real*8, intent(in) :: r
        real*8 :: V
        V = 1.0/(r**12) - 1.0/(r**6)
    end function lennard_jones_potential

    function coulomb_potential(r) result(V) 
        implicit none
        real*8, intent(in) :: r
        real*8 :: V
        V = -1.0/r
    end function coulomb_potential

    function effective_hydrogen_potential(r) result(V) 
        implicit none
        real*8, intent(in) :: r
        integer::l
        real*8 :: V
        l = 1
        V = -1.0/r + 0.5*l*(l+1)/r**2
    end function effective_hydrogen_potential


    function create_potential_array(V_ptr,x_array)
        procedure (pointing_func),pointer:: V_ptr
        real*8,intent(in)::x_array(:)
        real*8::create_potential_array(size(x_array))
        integer::i
        do i=1,size(x_array)
            create_potential_array(i) = V_ptr(x_array(i))
        end do
    end function create_potential_array

    function energy_minus_potential_array(E,V_array)
        implicit none 
        integer::i
        real*8, intent(in)::V_array(:)
        real*8, intent(in)::E
        real*8::energy_minus_potential_array(size(V_array))
        
        !do i = 1,size(V_array
        !    energy_minus_potential_array(i) = E - V_array
        !end do
        energy_minus_potential_array = E-V_array
    end function energy_minus_potential_array

    subroutine produce_trial_solution(x0,x1,x_final,psi_0,psi_1,E,x_array,psi_array,V_ptr,N)
        implicit none
        real*8,intent(in)::x0,x1,x_final,psi_0,psi_1,E
        integer,intent(in)::N
        real*8,intent(inout)::x_array(N),psi_array(N)
        real*8::h,potential_array(N)

        procedure (pointing_func),pointer:: V_ptr

        h = x1-x0
        potential_array = create_potential_array(V_ptr,x_array)
        !N = (x_final-x0)/h
        psi_array = numerov_whole_interval_schrodinger(x_array,psi_0,psi_1,potential_array,E)


    end subroutine produce_trial_solution

    function integrate_trapezium(h,f_array)
        implicit none
        real*8, intent(in)::h,f_array(:)
        real*8::integrate_trapezium
        integer::N,i
        N = size(f_array)
        integrate_trapezium = 0.0
        do i = 2,N-1
            integrate_trapezium = integrate_trapezium + f_array(i)
        end do
        integrate_trapezium = integrate_trapezium + 0.5*(f_array(1)+f_array(N))
        integrate_trapezium = integrate_trapezium*h

    end function integrate_trapezium


    subroutine normalise_wave_function(psi_array,h)
        implicit none
        real*8, intent(in)::h
        real*8, intent(inout)::psi_array(:)
        real*8, dimension(size(psi_array))::psi_sq
        real*8::norm_factor
        integer::i
        psi_sq = psi_array*psi_array
        norm_factor = integrate_trapezium(h,psi_sq)
        do i = 1,size(psi_array)
            psi_array(i) = psi_array(i)/norm_factor
        end do
    end subroutine normalise_wave_function

    function my_linspace(x_0,x_last,N)
        real*8, intent(in)::x_0,x_last
        integer::N,i
        real*8::h
        real*8::my_linspace(N)

        h = (x_last-x_0)/(N-1)

        do i = 1,N
            my_linspace(i) = x_0 + (i-1)*h
            !print*, "hello"
        end do

    end function my_linspace

    function find_bracketing_pair(Estart,deltaE,psi_left_boundary,psi_right_boundary,N,x_0,x_last,V_ptr)
        !future: make this find bracketing pairs which returns a list of tuples?
        real*8, intent(in)::Estart,deltaE,psi_left_boundary,psi_right_boundary,x_0,x_last
        integer, intent(in)::N
        integer::i,max_E_iter
        real*8::psi_array(N),x_array(N),potential_array(N),h
        real*8::currentE,previousE,currentDelta,previousDelta
        real*8::find_bracketing_pair
        real*8::psi_1,PI_squared
        procedure (pointing_func),pointer:: V_ptr
        PI_squared = 9.86960440108935861883449
        psi_1 = 0.1
        currentE = Estart
        x_array = my_linspace(x_0,x_last,N)

        !do i =1,N
        !    print*, x_array(i)
        !end do
        h = x_array(2)-x_array(1)
        potential_array = create_potential_array(V_ptr,x_array)
        max_E_iter = 1000000

        do i = 1,max_E_iter
            currentE = currentE + deltaE 
            psi_array = numerov_whole_interval_schrodinger(x_array,psi_left_boundary,psi_1,potential_array,currentE)
            
            call normalise_wave_function(psi_array,h)
            print*, 'current E units of pi^2/2: ', 2*currentE/PI_squared,' psi at boundary', psi_array(size(psi_array))
        end do 

        find_bracketing_pair = currentE



    end function find_bracketing_pair

    !subroutine find_allowed_solution(Estart,deltaE,psi_right_boundary,N)
    !
    !end subroutine find_allowed_solution

end module schrodinger_solution_mod






program run_solver
    use schrodinger_solution_mod
    use numerov
    implicit none
    character(len=*), parameter :: OUT_FILE = 'solutiondata.dat' ! Output file.
    character(len=*), parameter :: PLT_FILE = 'plot.plt' ! Gnuplot file.
    real*8::x0,y_left,y_right
    real*8:: x_final,E,deltaE
    integer::N,i,quantum_number
    !real*8,dimension(:),allocatable :: x_array
    
    procedure (pointing_func),pointer:: V_ptr => zero_potential

    N = 10000
    !allocate(x_array(N))
    
    y_left = 0.0
    y_right = 0.0
    x0 = 0.0
    x_final = 1.0
    E = 0.0
    deltaE = 0.0001

    E = find_bracketing_pair(E,deltaE,y_left,y_right,N,x0,x_final,V_ptr)

    !call execute_command_line('gnuplot -p ' // PLT_FILE)

    !deallocate(x_array)

end program run_solver




