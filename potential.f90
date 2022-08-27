module potentials_module
    implicit none
    integer, parameter :: dp = selected_real_kind(p = 15, r = 307) 

    private 
    public coulomb_potential,lennard_jones_potential,linspace,dp
contains 
    function coulomb_potential(r) result(V) 
        implicit none
        real, intent(in) :: r
        real :: V
        V = -1.0/r
    end function coulomb_potential

    function lennard_jones_potential(r) result(V)
        implicit none
        real(dp), intent(in) :: r
        real(dp) :: V
        V = 1.0/(r**12) - 1.0/(r**6)
    end function lennard_jones_potential

    subroutine linspace(from,final,array)

        real(dp),intent(in) :: from, final
        real(dp), intent(out) :: array(:)
        real(dp):: range 
        integer :: n, i
        n = size(array)
        range = final -from
        if (n == 0 ) return
        if (n == 1) then 
            array(1) = from
            return
        end if
        do i = 1,n
            array(i) = from + range * (i-1)/(n-1)
        end do
    end subroutine linspace
end module potentials_module

    


program potential
    
    use potentials_module

    implicit none

    real:: r
    real :: V

    real(dp)::distances_array(1000)
    real(dp)::potentials_array(size(distances_array))

    integer:: i

    call linspace(from = 0.9_dp, final = 2.5_dp, array = distances_array)

    do i =  1,size(potentials_array)
        potentials_array(i) = lennard_jones_potential(distances_array(i))
    end do
    
    open (1, file = "lennarddata.dat")
    do i = 1,size(potentials_array)
        write(1,*) distances_array(i), potentials_array(i)
    end do 
    close(1)


    !print *, V
end program potential








subroutine coulomb_potential_sub(r)
    implicit none
    real, intent(in):: r
    real :: V
    V = -1.0/r
    print *, V
end subroutine coulomb_potential_sub