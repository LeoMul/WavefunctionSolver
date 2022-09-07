!gfortran schrodinger_solution_methods.f90 -c may need to be ran
module schrodinger_solution_mod
    use numerov
    implicit none 
    private 
    public produce_trial_solution,lennard_jones_potential,zero_potential,normalise_wave_function,create_potential_array,energy_minus_potential_array,effective_hydrogen_potential
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
        psi_array = numerov_whole_interval_schrodinger(x_array,h,psi_0,psi_1,potential_array,E)


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

        !print*, size(psi_array)
        norm_factor = integrate_trapezium(h,psi_sq)
        do i = 1,size(psi_array)
            !print*, psi_sq(i)
            psi_array(i) = psi_array(i)/norm_factor
        end do


    end subroutine normalise_wave_function

    function find_bracketing_pair(Estart,deltaE,psi_right_boundary,N)
        !future: make this find bracketing pairs which returns a list of tuples?
        real*8, intent(in)::Estart,deltaE,psi_right_boundary
        integer, intent(in)::N
        integer::i,max_E_iter
        real*8::psi_array(N)
        real*8::currentE,previousE,currentDelta,previousDelta
        real*8::find_bracketing_pair
        currentE = Estart
        do i = 1,max_E_iter


            currentE = currentE + (i-1)*deltaE 
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
    real*8::x0,x1,h,y0,y1
    real*8:: x_final,E
    integer::n,i,quantum_number
    real*8,dimension(:),allocatable :: x_array
    real*8,dimension(:),allocatable :: y_array
    real*8,dimension(:),allocatable :: normed_y_array
    
    real*8::pi,pi_squared,root2
    procedure (pointing_func),pointer:: V_ptr => effective_hydrogen_potential

    pi = 4.D0*DATAN(1.D0)
    pi_squared = pi**2
    root2 = 2**0.5
    x0 = 0.0
    x_final = 2.5
    x1 = 0.01
    quantum_number = 3
    h = x1-x0
    n = (x_final-x0)/h
    E = 10


    allocate(x_array(n))
    allocate(y_array(n))
    allocate(normed_y_array(n))

    do i = 1,size(x_array)
        x_array(i) = x0 + (i-1)*h
        !print *, x_array(i)
        !y_array(i) = root2 * sin(pi*quantum_number*x_array(i))
        !normed_y_array(i) = y_array(i)
    end do
    y_array = create_potential_array(V_ptr,x_array)
    y_array = energy_minus_potential_array(E,y_array)
    open (1, file = "solutiondata.dat", action="write")
        do i = 1,size(x_array)
            write(1,*) x_array(i), y_array(i)!,normed_y_array(i)
        end do 
    close(1)
!
    call execute_command_line('gnuplot -p ' // PLT_FILE)

    deallocate(x_array)
    deallocate(y_array)
    deallocate(normed_y_array)

end program run_solver




