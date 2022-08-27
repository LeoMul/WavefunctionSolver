program fortranlearn
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none


    character(len = 3):: first_name
    character(len = 10):: family_name
    character(14)::full_name
    real(dp):: pi
    integer :: amount 
    
    integer:: intarray(0:9)
    
    intarray(:) = 1
    amount = 10
    pi = 3.141592653589793238462643383279

    first_name = 'Leo'
    family_name = 'Mulholland'
    full_name = first_name//' '//family_name
    print *, full_name
    print *, 'Hydrogen Initialised. Ten is: ', amount
    print *, 'Pi is: ', pi
    print *, intarray


end program fortranlearn