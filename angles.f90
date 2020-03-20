module mod_angles

implicit none


contains


subroutine get_neighbor_list( dist_matrix, dist_atoms, N_atoms, &
                              N_neighbor, neighbor_list )
    implicit none
    integer, intent(in)  :: N_atoms
    real*8, intent(in)   :: dist_matrix(:,:)
    integer, intent(in)  :: dist_atoms(:,:)
    integer, intent(out) :: N_neighbor(N_atoms)
    integer, allocatable, intent(out) :: neighbor_list(:,:)
    integer :: aux_list(N_atoms, size(dist_matrix,1))
    integer :: i, iat, size_neigh

    aux_list   = 0
    N_neighbor = 0
    do i = 1, size(dist_matrix,1)
        iat = dist_atoms(i,1)
        N_neighbor(iat) = N_neighbor(iat) + 1
        aux_list(iat,N_neighbor(iat)) = i
    enddo

    size_neigh = maxval(N_neighbor)

    allocate(neighbor_list(N_atoms, size_neigh))

    neighbor_list = aux_list(:,:size_neigh)

end subroutine get_neighbor_list



subroutine compute_angle_distr(dist_matrix, dist_atoms, neighbor_list, N_neighbor, &
                               atomtype, type1, type2, bins, angle_distr)
    implicit none
    real*8, intent(in)  :: dist_matrix(:,:)
    integer, intent(in) :: dist_atoms(:,:), neighbor_list(:,:), N_neighbor(:), &
                           atomtype(:), type1, type2, bins
    real*8, intent(out) :: angle_distr(bins)
    real*8, parameter   :: pi = acos(-1.0d0)
    integer :: i, j, k, ind, N_atoms
    real*8  :: rj(3), rk(3), theta, dtheta

    N_atoms = size(N_neighbor)

    dtheta = pi/bins

    angle_distr = 0.0

    do i = 1, N_atoms
        if ( N_neighbor(i) < 2 ) then
            continue
        endif

        if (atomtype(i) == type1 ) then
            do j = 1, N_neighbor(i)
                ind = neighbor_list(i,j)
                rj = dist_matrix(ind,:)
                if ( atomtype(dist_atoms(ind,2)) == type2 ) then
                    do k = j+1, N_neighbor(i)
                        ind = neighbor_list(i,k)
                        rk = dist_matrix(ind,:)

                        if ( atomtype(dist_atoms(ind,2)) == type2 ) then

                            theta =  ( sum(rj*rk) / (norm2(rj)*norm2(rk)) )
                            if (abs(theta) + 1.0D-8 < 1.0d0) then
                                theta = acos(theta)
                                ind = int(theta/dtheta)
                                angle_distr(ind) = angle_distr(ind) + 1
                            endif

                        endif
                    enddo
                endif
            enddo
        endif
    enddo


end subroutine compute_angle_distr

end module mod_angles