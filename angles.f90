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
    real*8  :: aux_dist(N_atoms, size(dist_matrix,1))
    integer :: i, iat, size_neigh
    real*8  :: d

    aux_list   = 0
    N_neighbor = 0
    do i = 1, size(dist_matrix,1)
        iat = dist_atoms(i,1)
        d = norm2( dist_matrix(i,:) )

        ! Insertion sort
        call insert( aux_dist(iat,:), aux_list(iat,:), N_neighbor(iat), d, i )

    enddo
    
    size_neigh = maxval(N_neighbor)

    allocate(neighbor_list(N_atoms, size_neigh))

    neighbor_list = aux_list(:,:size_neigh)

end subroutine get_neighbor_list


subroutine insert(V, W, N, x, y)
    implicit none
    real*8, intent(inout)  :: V(:)
    integer, intent(inout) :: N, W(:)
    real*8, intent(in)     :: x
    integer, intent(in)    :: y
    integer                :: j

    N = N + 1
    V(N) = x
    W(N) = y
    do j = N, 2 , -1
        if ( V(j-1) < V(j) ) then
            exit
        else
            V(j) = V(j-1)
            V(j-1) = x

            W(j) = W(j-1)
            W(j-1) = y
        endif
    enddo

end subroutine insert



subroutine compute_angle_distr(dist_matrix, dist_atoms, neighbor_list, N_neighbor, &
                               atomtype, type1, type2, bins, rcut, angle_distr)
    implicit none
    real*8, intent(in)  :: dist_matrix(:,:), rcut
    integer, intent(in) :: dist_atoms(:,:), neighbor_list(:,:), N_neighbor(:), &
                           atomtype(:), type1, type2, bins
    real*8, intent(out) :: angle_distr(bins)
    real*8, parameter   :: pi = acos(-1.0d0), c = 0.05d0
    integer :: i, j, k, ind, N_atoms, cont
    real*8  :: rj(3), rk(3), theta, dtheta, x

    N_atoms = size(N_neighbor)

    dtheta = pi/bins

    angle_distr = 0.0

    cont = 0
    do i = 1, N_atoms
        if ( N_neighbor(i) < 2 ) then
            continue
        endif

        if (atomtype(i) == type1 ) then
            do j = 1, N_neighbor(i)
                ind = neighbor_list(i,j)
                rj = dist_matrix(ind,:)
                if ( norm2(rj) <= rcut .and. &
                     atomtype(dist_atoms(ind,2)) == type2 ) then
                    do k = j+1, N_neighbor(i)
                        ind = neighbor_list(i,k)
                        rk = dist_matrix(ind,:)

                        if ( norm2(rk) <= rcut .and. &
                             atomtype(dist_atoms(ind,2)) == type2 ) then

                            theta =  ( sum(rj*rk) / (norm2(rj)*norm2(rk)) )
                            if (abs(theta) + 1.0D-8 < 1.0d0) then
                                theta = acos(theta)
                                ind = int(theta/dtheta)
                                angle_distr(ind) = angle_distr(ind) + 1

                                cont = cont + 1

                                !do ind=1, bins
                                !    x=dtheta*ind
                                !    angle_distr(ind)=angle_distr(ind)+exp(-(x-theta)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*pi))
                                !enddo

                            endif

                        endif
                    enddo
                endif
            enddo
        endif
    enddo

    angle_distr = angle_distr/(dtheta*cont)

end subroutine compute_angle_distr


subroutine get_mean_sigma( angle_distr, mean, sigma )
    implicit none
    real*8, intent(in)  :: angle_distr(:)
    real*8, intent(out) :: mean, sigma
    real*8, parameter   :: pi = acos(-1.0d0)
    real*8              :: x, dx
    integer             :: i, N

    N = size(angle_distr)
    dx = pi/N

    mean = 0.0d0
    do i = 1, N
        x = (i-0.5)*dx
        mean = mean + x * angle_distr(i) * dx
    enddo

    sigma = 0.0d0
    do i = 1, N
        x = (i-0.5)*dx
        sigma = sigma + ( x-mean )**2 * angle_distr(i) * dx
    enddo


end subroutine get_mean_sigma


end module mod_angles