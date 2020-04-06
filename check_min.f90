module mod_check_min

contains


function integrate(pdf,min_index,norm,dx) result(ans)
implicit none
real*8, dimension(:), intent(in)            :: pdf
integer, intent(in)                         :: min_index
real*8, intent(in)                          :: dx, norm
real*8                                      :: ans, pi
integer                                     :: i

pi = acos(-1.0d0)
ans= 0.0

do i = 2, min_index
    ans = ans + dx*((pdf(i)+pdf(i-1))/2.0d0)*(4.0*pi*(dx*(i-0.5d0))**2)*norm
enddo   
end function


subroutine check_min ( mat, i_init, min_index )
    implicit none
    real*8, intent(in)   :: mat(:)
    integer, intent(in)  :: i_init
    integer, intent(out) :: min_index
    real*8               :: mat_smooth(size(mat))
    integer              :: i, j, N, N_smooth, N_check
    real*8               :: u, f_min
    logical              :: minimum

    N_smooth = 2
    N_check  = 10

    N = size(mat)

    ! Smoothing
    do i = 1, N
        if (i <= N_smooth ) then
            u = sum( mat(1:i+N_smooth) )/(N_smooth+i-1+1)
        elseif (i+N_smooth >= N ) then
            u  = sum( mat(i-N_smooth:N) )/(N-i+N_smooth+1)
        else
            u = sum( mat(i-N_smooth:i+N_smooth) )/(2*N_smooth+1)
        endif
        mat_smooth(i) = max(u, 1.0d-14)
    enddo

    !mat_smooth = mat

    !Find minimum
    ! f_min = -huge(1.0)
    ! do i_init = 1, N
    !     if ( mat_smooth(i_init) > 0.2 ) then
    !         f_min = mat_smooth(i_init)
    !         min_index = i_init
    !         exit
    !     endif
    ! enddo

    do i = i_init, N

        minimum = .True.
        do j = max(1, i-N_check), i-1!min(N, i+N_check)
            if (i == j) then
                continue
            endif
            if ( mat_smooth(j) <= mat_smooth(i) ) then
                minimum = .false.
                exit
            endif
        enddo

        if (minimum) then
            do j = i+1, min(N, i+N_check)
                if (i == j) then
                    continue
                endif
                if ( mat_smooth(j) < mat_smooth(i) ) then
                    minimum = .false.
                    exit
                endif
            enddo
        endif

        if (minimum) then
            min_index = i
            exit
        endif

    enddo

    if ( .not. minimum ) then
        print*, "ERROR: Could not find any minima"
        min_index = -1
    endif


end subroutine check_min

end module mod_check_min