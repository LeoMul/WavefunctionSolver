module bisector
    implicit none
    public root_finder
    private 
contains
    function poly(x)
        real*16,intent(in)::x 
        real*16::poly

        poly = (x-2.0)**3   
    end function

    function root_finder(a,b,tol,N)

        real*16,intent(inout)::a,b
        real*16,intent(in)::tol
        integer,intent(in)::N
        real*16::root_finder
        real*16::FA,FP,p
        integer::i

        i = 1
        do while (i<N)
            !print*,"a,p,b ", a,p,b
            FA = poly(a)
            !FB = poly(b)
            p = a + (b-a)/2.0
            !print*,"FA,FP ", FA,FP

            FP = poly(p)

            if (FP == 0.0 .or. (b-a)/2 < tol) then
                root_finder = p
                print*, "method converged, root found ", root_finder
                return
            else 
                IF (FA*FP > 0) THEN
                    a = p
                    FA = FP
                ELSE
                    b = p
                END IF
            end if

            i = i + 1
        end do


        print*, "Method did not converge in allocated ", N, "steps."
        print*, "Final p ",p
        return
    end function root_finder


end module bisector


program root_bisection
    use bisector
    real*16::a,b,tol,root
    integer::N
    a = 0.0
    b = 3.0
    tol = 0.00000000000001
    N = 100
    root = root_finder(a,b,tol,N)
end program root_bisection