module mod_angles

implicit none


contains


subroutine compare(ref, vec, out)
    implicit none
    integer, intent(in)   :: ref(:), vec(:)
    logical, intent(out) :: out

    integer :: i, j, N, cont

    N = size(vec)

    cont = 0
    do i = 1, size(vec)
        do j = 1, size(ref)
            if (vec(i) == ref(j)) then
                cont = cont + 1
                exit
            endif
        enddo
    enddo

    if (cont == N) then
        out = .true.
    else
        out = .false.
    endif

end subroutine compare


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
subroutine get_neighbor_list2( ext, dist_matrix, dist_atoms, N_atoms, N_species, atomtype, &
                               mat_neighbor,N_neighbor, neighbor_list )
    implicit none
    integer, intent(in)  :: N_atoms, N_species, dist_atoms(:,:), atomtype(:), &
                            mat_neighbor(:,:), ext
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
    
    size_neigh = maxval(N_neighbor)+ext
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


subroutine find_pair_index( ext, tag1, tag2, neighbor_order_list, N_neighbor, ind, ierr )
    implicit none
    integer, intent(in)  :: ext, tag1, tag2, neighbor_order_list(:), N_neighbor
    integer, intent(out) :: ind, ierr
    integer :: i, j, ind1, ind2, found1, found2

    found1 = 0
    found2 = 0
    ind1 = -1
    ind2 = -1
    do i = 1, N_neighbor + ext

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
        ! print*, "ERROR: ", ind1, ind2
        ! print*, tag1, tag2
        ! print*, neighbor_order_list(:)
        ierr = -1
        !STOP
    else
        ierr = 0
    endif

    i = min(ind1,ind2)
    j = max(ind1,ind2)

    ind = i + (j-1)*(j-2)/2
    
end subroutine find_pair_index


subroutine update_angle_distr(ext, neighbor_order_list, atomtype, mat_neighbor, &
                              dist_matrix, neighbor_list, angle_distr, cont_adf )
    implicit none
    integer, intent(in)   :: ext, neighbor_order_list(:,:,:), atomtype(:), &
                             neighbor_list(:,:,:,:), mat_neighbor(:,:)
    real*8, intent(in)    :: dist_matrix(:,:)
    real*8, intent(inout) :: angle_distr(:,:,:,:)
    integer, intent(inout) :: cont_adf(:,:,:)
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
                    call find_pair_index( ext, tag1, tag2, neighbor_order_list(iat,iesp,:), &
                                          N_neigh, ind, ierr )

                    ! if (ierr == -1) then
                    !     print*, iat, iesp 
                    !     write(*,"(100i5)") neighbor_list(iat,iesp,:N_neigh,2)
                    !     write(*,"(100i5)") neighbor_order_list(iat,iesp,:N_neigh)
                    !     STOP
                    ! endif

                    if (ierr == 0) then
                        ! Update histogram
                        theta =  ( sum(rj*rk) / (norm2(rj)*norm2(rk)) )
                        if (abs(theta) + 1.0D-8 < 1.0d0) then
                            cont_adf(iat,iesp,ind) = cont_adf(iat,iesp,ind) + 1
                            theta = acos(theta)
                            ihist = int(theta/pi*bins)
                            angle_distr(iat,iesp,ind,ihist) = angle_distr(iat,iesp,ind,ihist) + 1.0d0
                            ! do ihist=1, bins
                            !     x=pi/bins*ihist
                            !     angle_distr(iat,iesp,ind,ihist) = angle_distr(iat,iesp,ind,ihist) + &
                            !                                       exp(-(x-theta)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*pi))
                            ! enddo
                        endif
                    endif

                enddo

            enddo

        enddo

    enddo


end subroutine update_angle_distr



subroutine apply_smearing_ang(histogram)
    implicit none
    real*8, intent(inout) :: histogram(:)
    integer :: bins, ihist, jhist, N, i
    real*8  :: copy(size(histogram)), theta, x, pi, c, dx, norm

    pi = acos(-1.0d0)
    c  = 0.05d0

    bins = size(histogram)
    copy = 0.0d0

    do ihist = 1, bins
        if (histogram(ihist) < 1.0D-8) cycle
        theta = (pi/bins)*(ihist+0.5d0)
        do jhist = 1, bins
            x = (pi/bins)*jhist
            copy(jhist) = copy(jhist) + histogram(ihist)*exp(-(x-theta)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*pi))
        enddo
    enddo 

    dx = pi/bins
    norm = sum(copy)*dx
    do i = 1, bins
        if (copy(i) > 1.0d-14) then
            histogram(i) = copy(i)/norm
        else
            histogram(i) = copy(i)
        endif
    enddo 


end subroutine apply_smearing_ang


! angle_distr :: Distribution of the distance between each couple of atoms
!               lenght(natoms,N_species,max_neigh,bins)
! contribution_adf :: Decomposition of the total_pdf into the distributions of the first neighbors
!                     lenght(N_species,N_species,max_neigh,bins)
! total_adf  :: total pdf for each pair of species
!               lenght(N_species,N_species,bins)

subroutine total_adf( ext, angle_distr, mat_neighbor, neighbor_list, atomtype, numions, contribution_adf, tot_adf, cont_adf)
    implicit none
    real*8, intent(inout)  :: angle_distr(:,:,:,:)
    integer, intent(in) :: ext, neighbor_list(:,:,:,:), mat_neighbor(:,:), atomtype(:), numions(:), cont_adf(:,:,:)
    real*8, intent(out) :: tot_adf(:,:,:), contribution_adf(:,:,:,:)
    integer :: N_atoms, N_species, max_neigh, bins, iat, iesp, n1, ihist, jesp, N_neigh
    integer, allocatable :: cont(:,:,:)


    N_atoms   = size(angle_distr,1)
    N_species = size(angle_distr,2)
    max_neigh = size(angle_distr,3)
    bins      = size(angle_distr,4)

    tot_adf = 0.0d0
    contribution_adf = 0.0d0

    allocate(cont(N_species, N_species, max_neigh))
    cont = 0

    ! Loop in atoms
    do iat = 1, N_atoms

        ! Loop in species
        do iesp = 1, N_species 

            if ( atomtype(iat) == iesp ) cycle
            N_neigh = mat_neighbor(atomtype(iat),iesp)*(mat_neighbor(atomtype(iat),iesp)-1)/2! + ext
            do n1 = 1, N_neigh

                !call apply_smearing_ang( angle_distr(iat,iesp,n1,:) )

                cont(atomtype(iat),iesp,n1) = cont(atomtype(iat),iesp,n1) + 1
                do ihist = 1, bins

                    !if ( angle_distr(iat,iesp,n1,ihist) > 1.0d-8 ) then
                        !print*, iat, iesp, n1, sum(angle_distr(iat,iesp,n1,:))*acos(-1.0d0)/bins
                        contribution_adf(atomtype(iat),iesp,n1,ihist) = contribution_adf(atomtype(iat),iesp,n1,ihist) + &
                                                                        angle_distr(iat,iesp,n1,ihist)!/cont_adf(iat,iesp,n1)! / &
                                                                        !N_neigh
                    !endif

                enddo

            enddo

        enddo

    enddo

    do iat = 1, N_species
        do iesp = 1, N_species
            if(iat == iesp) cycle
            N_neigh = mat_neighbor(iat,iesp)*(mat_neighbor(iat,iesp)-1)/2! + ext
            do n1 = 1, N_neigh
                contribution_adf(iat,iesp,n1,:) = contribution_adf(iat,iesp,n1,:)/cont(iat,iesp,n1)
                !print*, iat, iesp, n1, sum(contribution_adf(iat,iesp,n1,:))*acos(-1.0d0)/bins
            enddo
        enddo
    enddo

    !tot_adf = sum(contribution_adf,dim=3)

end subroutine total_adf




subroutine get_mean_sigma_angle( N_pair, angle_distr, mean, sigma )
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

        do i = 1, N
            x = (i-0.5)*dx
            mean(ipair) = mean(ipair) + x * angle_distr(ipair,i) * dx
        enddo

        do i = 1, N
            x = (i-0.5)*dx
            sigma(ipair) = sigma(ipair) + ( x-mean(ipair) )**2 * angle_distr(ipair,i) * dx
        enddo

        !print*, sum(angle_distr(ipair,:))*dx


    enddo


    mean = mean*180/pi
    sigma = sqrt(sigma)*180/pi


end subroutine get_mean_sigma_angle

end module mod_angles
