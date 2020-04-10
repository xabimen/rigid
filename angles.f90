module mod_angles

implicit none


contains


subroutine get_neighbor_list( dist_matrix, dist_atoms, N_atoms, N_species, atomtype, &
                              N_neighbor, neighbor_list )
    implicit none
    integer, intent(in)  :: N_atoms, N_species, dist_atoms(:,:), atomtype(:)
    real*8, intent(in)   :: dist_matrix(:,:)
    integer, intent(out) :: N_neighbor(N_atoms, N_species)
    integer, allocatable, intent(out) :: neighbor_list(:,:,:,:)
    integer :: aux_list(N_atoms, N_species, size(dist_matrix,1),2)
    real*8  :: aux_dist(N_atoms, size(dist_matrix,1))
    integer :: i, iat, iesp, size_neigh, ins(2)
    real*8  :: d

    aux_list   = 0
    N_neighbor = 0
    do i = 1, size(dist_matrix,1)
        iat  = dist_atoms(i,1)
        iesp = atomtype( dist_atoms(i,2) )
        d = norm2( dist_matrix(i,:) )

        ! Insertion sort
        ins = (/i,dist_atoms(i,2)/)
        !call insert( aux_dist(iat,:), aux_list(iat,:), N_neighbor(iat), d, i )
        call insert2( aux_dist(iat,:), aux_list(iat,iesp,:,:), N_neighbor(iat,iesp), d, ins )

    enddo
    
    size_neigh = maxval(N_neighbor)

    allocate(neighbor_list(N_atoms, N_species, size_neigh,2))

    neighbor_list = aux_list(:,:,:size_neigh,:)

end subroutine get_neighbor_list


! Solo el numero de vecinos que hay hasta el primer minimo
subroutine get_neighbor_list2( dist_matrix, dist_atoms, N_atoms, N_species, atomtype, &
                               mat_neighbor,N_neighbor, neighbor_list )
    implicit none
    integer, intent(in)  :: N_atoms, N_species, dist_atoms(:,:), atomtype(:), &
                            mat_neighbor(:,:)
    real*8, intent(in)   :: dist_matrix(:,:)
    integer, intent(out) :: N_neighbor(N_atoms, N_species)
    integer, allocatable, intent(out) :: neighbor_list(:,:,:,:)
    integer :: aux_list(N_atoms, N_species, size(dist_matrix,1),2)
    real*8  :: aux_dist(N_atoms, N_species, size(dist_matrix,1))
    integer :: i, iat, iesp, size_neigh, ins(2)
    real*8  :: d

    aux_list   = 0
    N_neighbor = 0
    do i = 1, size(dist_matrix,1)
        iat  = dist_atoms(i,1)
        iesp = atomtype( dist_atoms(i,2) )
        d = norm2( dist_matrix(i,:) )

        ! Insertion sort
        ins = (/i,dist_atoms(i,2)/)
        !call insert( aux_dist(iat,:), aux_list(iat,:), N_neighbor(iat), d, i )
        call insert2( aux_dist(iat,iesp,:), aux_list(iat,iesp,:,:), N_neighbor(iat,iesp), d, ins )

    enddo
    
    size_neigh = maxval(N_neighbor)
    allocate(neighbor_list(N_atoms, N_species, size_neigh,2))

    neighbor_list = aux_list(:,:,:size_neigh,:)

    ! N_neighbor = 0
    ! do iat = 1, N_atoms
    !    do iesp = 1, N_species
    !        do i = 1, mat_neighbor(atomtype(iat),iesp)
    !            neighbor_list(iat,iesp,i,:) = aux_list(iat,iesp,i,:)
    !            N_neighbor(iat,iesp) = N_neighbor(iat,iesp) + 1
    !        enddo
    !    enddo
    ! enddo

end subroutine get_neighbor_list2


subroutine insert(V_dist, W_list, N, x_d, y_n)
    implicit none
    real*8, intent(inout)  :: V_dist(:)
    integer, intent(inout) :: N, W_list(:)
    real*8, intent(in)     :: x_d
    integer, intent(in)    :: y_n
    integer                :: j

    N = N + 1
    V_dist(N) = x_d
    W_list(N) = y_n
    do j = N, 2 , -1
        if ( V_dist(j-1) < V_dist(j) ) then
            exit
        else
            V_dist(j) = V_dist(j-1)
            V_dist(j-1) = x_d

            W_list(j) = W_list(j-1)
            W_list(j-1) = y_n
        endif
    enddo

end subroutine insert


subroutine insert2(V_dist, W_list, N, x_d, y_n)
    implicit none
    real*8, intent(inout)  :: V_dist(:)
    integer, intent(inout) :: N, W_list(:,:)
    real*8, intent(in)     :: x_d
    integer, intent(in)    :: y_n(:)
    integer                :: j

    N = N + 1
    V_dist(N)   = x_d
    W_list(N,:) = y_n
    do j = N, 2 , -1
        if ( V_dist(j-1) < V_dist(j) ) then
            exit
        else
            V_dist(j) = V_dist(j-1)
            V_dist(j-1) = x_d

            W_list(j,:)   = W_list(j-1,:)
            W_list(j-1,:) = y_n(:)
        endif
    enddo

end subroutine insert2


subroutine find_pair_index( tag1, tag2, neighbor_order_list, N_neighbor, ind, ierr )
    implicit none
    integer, intent(in)  :: tag1, tag2, neighbor_order_list(:), N_neighbor
    integer, intent(out) :: ind, ierr
    integer :: i, j, ind1, ind2, found1, found2

    found1 = 0
    found2 = 0
    ind1 = -1
    ind2 = -1
    do i = 1, N_neighbor

        if (found1 == 1 .and. found2 == 1) then
            exit
        endif

        if ( neighbor_order_list(i) == tag1 ) then
            ind1 = i
            found1 = found1 + 1
        endif
        if ( neighbor_order_list(i) == tag2 ) then
            ind2 = i
            found2 = found2 + 1
        endif

    enddo

    if (ind1 < 0 .or. ind2 < 0) then
        print*, "ERROR: ", ind1, ind2
        print*, tag1, tag2
        print*, neighbor_order_list(:)
        ierr = -1
        !STOP
    else
        ierr = 0
    endif

    i = min(ind1,ind2)
    j = max(ind1,ind2)

    ind = i + (j-1)*(j-2)/2
    
end subroutine find_pair_index


subroutine update_angle_distr(neighbor_order_list, atomtype, mat_neighbor, &
                              dist_matrix, neighbor_list, angle_distr )
    implicit none
    integer, intent(in)   :: neighbor_order_list(:,:,:), atomtype(:), &
                             neighbor_list(:,:,:,:), mat_neighbor(:,:)
    real*8, intent(in)    :: dist_matrix(:,:)
    real*8, intent(inout) :: angle_distr(:,:,:,:)
    integer               :: N_atoms, N_species, max_pair, bins, N_pair, N_neigh
    integer               :: iat, iesp, ipair, n1, n2, tag1, tag2, ind, ihist, ierr
    real*8                :: rj(3), rk(3), pi, c, theta, x

    pi = acos(-1.0d0)
    c  = 0.03d0

    N_atoms   = size(angle_distr,1)
    N_species = size(angle_distr,2)
    max_pair  = size(angle_distr,3)
    bins      = size(angle_distr,4)


    ! Loop in atoms
    do iat = 1, N_atoms 

        ! Loop in species
        do iesp = 1, N_species 

            if ( atomtype(iat) == iesp ) cycle

            N_neigh = mat_neighbor(atomtype(iat),iesp)
            N_pair = N_neigh*(N_neigh-1)/2

            ! Loop in all posible pairs
            do n1 = 1, mat_neighbor(atomtype(iat),iesp)
                tag1 = neighbor_list(iat,iesp,n1,2)
                rj   = dist_matrix( neighbor_list(iat,iesp,n1,1), : )

                do n2 = n1+1, mat_neighbor(atomtype(iat),iesp)
                    tag2 = neighbor_list(iat,iesp,n2,2)
                    rk   = dist_matrix( neighbor_list(iat,iesp,n2,1), : )

                    ! Find index in the angle_distr corresponding to the pair (n1,n2)
                    call find_pair_index( tag1, tag2, neighbor_order_list(iat,iesp,:), &
                                          N_neigh, ind, ierr )

                    if (ierr == -1) then
                        print*, iat, iesp 
                        write(*,"(100i5)") neighbor_list(iat,iesp,:N_neigh,2)
                        write(*,"(100i5)") neighbor_order_list(iat,iesp,:N_neigh)
                        STOP
                    endif

                    ! Update histogram
                    theta =  ( sum(rj*rk) / (norm2(rj)*norm2(rk)) )
                    if (abs(theta) + 1.0D-8 < 1.0d0) then
                        theta = acos(theta)
                        do ihist=1, bins
                            x=pi/bins*ihist
                            angle_distr(iat,iesp,ind,ihist) = angle_distr(iat,iesp,ind,ihist) + &
                                                              exp(-(x-theta)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*pi))
                        enddo
                    endif

                enddo

            enddo

        enddo

    enddo


end subroutine update_angle_distr


subroutine get_mean_sigma( N_pair, angle_distr, mean, sigma )
    implicit none
    integer, intent(in)   :: N_pair
    real*8, intent(inout) :: angle_distr(:,:)
    real*8, intent(out)   :: mean(:), sigma(:)
    real*8, parameter     :: pi = acos(-1.0d0)
    real*8                :: x, dx, norm
    integer               :: ipair, i, N   

    N = size(angle_distr,2)
    dx = pi/N

    mean = 0.0d0
    sigma = 0.0d0


    do ipair = 1, N_pair
        norm = sum(angle_distr(ipair,:))*dx

        do i = 1, N
            if (angle_distr(ipair,i) > 1.0d-14) then
                angle_distr(ipair,i) = angle_distr(ipair,i)/norm
            endif
        enddo 


        do i = 1, N
            x = (i-0.5)*dx
            mean(ipair) = mean(ipair) + x * angle_distr(ipair,i) * dx
        enddo

        do i = 1, N
            x = (i-0.5)*dx
            sigma(ipair) = sigma(ipair) + ( x-mean(ipair) )**2 * angle_distr(ipair,i) * dx
        enddo

    enddo


    mean = mean*180/pi
    sigma = sqrt(sigma)*180/pi


end subroutine get_mean_sigma

end module mod_angles
