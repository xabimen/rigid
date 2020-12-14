module pair_dist

contains

subroutine find_neigh_index( ext, tag1, neighbor_order_list, N_neighbor, ind, ierr )
    implicit none
    integer, intent(in)  :: ext, tag1, neighbor_order_list(:), N_neighbor
    integer, intent(out) :: ind, ierr
    integer :: i, j, found1

    found1 = 0
    ind = -1
    do i = 1, N_neighbor + ext

        if (found1 == 1 ) then
            exit
        endif

        if ( neighbor_order_list(i) == tag1 ) then
            ind = i
            found1 = found1 + 1
        endif

    enddo

    if (ind < 0 ) then
        ! print*, "ERROR: ", ind1, ind2
        ! print*, tag1, tag2
        ! print*, neighbor_order_list(:)
        ierr = -1
        !STOP
    else
        ierr = 0
    endif

    
end subroutine find_neigh_index


subroutine update_dist_distr(ext, V, neighbor_order_list, atomtype, numIons, mat_neighbor, &
                              dist_matrix, neighbor_list, Rcut, dist_distr, cont_pdf )
    implicit none
    integer, intent(in)   :: ext, neighbor_order_list(:,:,:), atomtype(:), &
                             neighbor_list(:,:,:,:), mat_neighbor(:,:), numIons(:)
    real*8, intent(in)    :: dist_matrix(:,:), Rcut, V
    real*8, intent(inout) :: dist_distr(:,:,:,:)
    integer, intent(inout) :: cont_pdf(:,:,:)
    integer               :: N_atoms, N_species, max_neigh, bins, N_pair, N_neigh
    integer               :: iat, iesp, ipair, n1, n2, tag1, tag2, ind, ihist, ierr
    real*8                :: rj(3), rk(3), pi, c, r, x, norm

    pi = acos(-1.0d0)
    c  = 0.03d0

    N_atoms   = size(dist_distr,1)
    N_species = size(dist_distr,2)
    max_neigh  = size(dist_distr,3)
    bins      = size(dist_distr,4)


    ! Loop in atoms
    do iat = 1, N_atoms 

        ! Loop in species
        do iesp = 1, N_species 

            if ( atomtype(iat) == iesp ) cycle

            N_neigh = mat_neighbor(atomtype(iat),iesp)

            ! Loop in all neighbors
            do n1 = 1, mat_neighbor(atomtype(iat),iesp) + ext
                tag1 = neighbor_list(iat,iesp,n1,2)
                rj   = dist_matrix( neighbor_list(iat,iesp,n1,1), : )

                call find_neigh_index( ext, tag1, neighbor_order_list(iat,iesp,:), N_neigh, ind, ierr )

                if (ierr == 0) then
                    ! Update histogram
                    cont_pdf(iat,iesp,ind) = cont_pdf(iat,iesp,ind) + 1
                    r = norm2(rj)

                    ihist = int(r/rcut*bins)

                    !norm = rcut/bins/( 4.0d0*pi * (rcut/bins)**3 * ( ihist**2 - ihist + 1.0d0/3 )  ) * &
                    !           1.0!/(numIons(iesp)*numIons(atomtype(iat)))!*numions(iesp)

                    dist_distr(iat,iesp,ind,ihist) = dist_distr(iat,iesp,ind,ihist) + 1.0d0


                    ! do ihist=1, bins
                    !     !if ((ihist*3-(ihist-1)**3)<0) print*, ihist, ihist**3, (ihist-1)**3, ihist*3-(ihist-1)**3
                    !     x=rcut/real(bins,8)*ihist

                    !     !print*,  ihist**3 - (ihist-1)**3, 4.0/3 * 2.25d0*(2*ihist-3)**2+6.75d0 
                    !     !print*, 4*pi*x**2*rcut/bins, 4.0d0*pi * (rcut/bins)**3 * ( ihist**2 - ihist + 1.0d0/3 ) 

                    !     norm = rcut/bins/( 4.0d0*pi * (rcut/bins)**3 * ( ihist**2 - ihist + 1.0d0/3 )  ) * &
                    !            1.0!/(numIons(iesp)*numIons(atomtype(iat)))!*numions(iesp)
                    !     !norm = rcut/bins/(4.0d0/3*pi * (rcut/bins)**3 * ( 2.25d0*(2*ihist-3)**2+6.75d0 ) ) !* V

                    !     ! i^3 - (i-1)^3 = 2.25*(2*x-3) + 6.75
                    !     dist_distr(iat,iesp,ind,ihist) = dist_distr(iat,iesp,ind,ihist) + &
                    !                                       exp(-(x-r)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*pi))*norm
                    ! enddo

                endif


            enddo

        enddo

        !if ( atomtype(iat) /= iesp ) print*, iat,iesp, sum(dist_distr(iat,iesp,:,:))
        

    enddo



end subroutine update_dist_distr


function integratee(pdf,min_index,norm,dx) result(ans)
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


subroutine apply_smearing_dist(histogram, rmax, N_file)
    implicit none
    real*8, intent(inout) :: histogram(:)
    real*8, intent(in)    :: rmax
    integer, intent(in)   :: N_file
    integer :: bins, ihist, jhist
    real*8  :: copy(size(histogram)), r, x, pi, c, norm

    pi = acos(-1.0d0)
    c  = 0.05d0

    bins = size(histogram)
    copy = 0.0d0

    do ihist = 1, bins
        if (histogram(ihist) < 1.0D-8) cycle
        r = (rmax/bins)*(ihist-0.5d0)
        do jhist = 1, bins
            x = (rmax/bins)*(jhist-0.5d0)
            copy(jhist) = copy(jhist) + histogram(ihist)*exp(-(x-r)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*pi))
        enddo
    enddo 

    do ihist = 1, bins
        if (copy(ihist) < 1.0D-8) cycle
        histogram(ihist) = copy(ihist)/N_file
    enddo 

    !histogram = copy/N_file
    

end subroutine apply_smearing_dist


subroutine normalize_gdr_dist(histogram, rmax, V_mean, numions_iesp, numions_atomtype_iat)
    implicit none
    real*8, intent(inout) :: histogram(:)
    real*8, intent(in)    :: rmax, V_mean
    integer, intent(in)   :: numions_iesp, numions_atomtype_iat
    integer :: bins, ihist, jhist
    real*8  :: r, x, pi, c, norm

    pi = acos(-1.0d0)
    c  = 0.05d0

    bins = size(histogram)

    do ihist = 1, bins
        if (histogram(ihist) < 1.0D-8) then
            cycle
        endif
        x = (rmax/bins)*ihist
        norm = rmax/bins/( 4.0d0*pi * (rmax/bins)**3 * ( ihist**2 - ihist + 1.0d0/3 )  ) * &
               V_mean/(numions_iesp*numions_atomtype_iat)

        histogram(ihist) = histogram(ihist)*norm
    enddo
    

end subroutine normalize_gdr_dist



! dist_distr :: Distribution of the distance between each couple of atoms
!               lenght(natoms,N_species,max_neigh,bins)
! contribution_pdf :: Decomposition of the total_pdf into the distributions of the first neighbors
!                     lenght(N_species,N_species,max_neigh,bins)
! total_pdf  :: total pdf for each pair of species
!               lenght(N_species,N_species,bins)
subroutine total_pdf( ext, V_mean, N_file, rmax, dist_distr, mat_neighbor, neighbor_list, &
                      atomtype, numions, contribution_pdf, tot_pdf, cont_pdf)
    implicit none
    real*8, intent(in)    :: rmax, V_mean
    integer, intent(in)   :: ext, N_file, neighbor_list(:,:,:,:), mat_neighbor(:,:), atomtype(:), numions(:), cont_pdf(:,:,:)
    real*8, intent(out)   :: tot_pdf(:,:,:), contribution_pdf(:,:,:,:)
    real*8, intent(inout) :: dist_distr(:,:,:,:)
    integer :: N_atoms, N_species, max_neigh, bins, iat, iesp, n1, ihist, jesp


    N_atoms   = size(dist_distr,1)
    N_species = size(dist_distr,2)
    max_neigh = size(dist_distr,3)
    bins      = size(dist_distr,4)

    tot_pdf = 0.0d0
    contribution_pdf = 0.0d0


    ! Loop in atoms
    do iat = 1, N_atoms

        ! Loop in species
        do iesp = 1, N_species 

            if ( atomtype(iat) == iesp ) cycle

            do n1 = 1, mat_neighbor(atomtype(iat),iesp) + ext

                ! call apply_smearing_dist(dist_distr(iat,iesp,n1,:), rmax, V_mean, cont_pdf(iat,iesp,n1),&
                !                          numions(iesp), numIons(atomtype(iat)))

                !print*, iat, iesp, n1, &
                !      integratee(dist_distr(iat,iesp,n1,:),bins,numions(iesp)/V_mean,rmax/real(bins)) 

                do ihist = 1, bins
                    
                    if (dist_distr(iat,iesp,n1,ihist)<1.0d-8) cycle

                    contribution_pdf(atomtype(iat),iesp,n1,ihist) = contribution_pdf(atomtype(iat),iesp,n1,ihist) + &
                                                                    dist_distr(iat,iesp,n1,ihist)!/cont_pdf(iat,iesp,n1)! / &
                                                                    !(mat_neighbor(atomtype(iat),iesp)+ext)
                    tot_pdf(atomtype(iat),iesp,ihist) = tot_pdf(atomtype(iat),iesp,ihist) + &
                                                        dist_distr(iat,iesp,n1,ihist)!/mat_neighbor(atomtype(iat),iesp)

                enddo

            enddo

        enddo

    enddo

    !  do iat = 1, N_atoms

    !     ! Loop in species
    !     do iesp = 1, N_species 

    !         if ( atomtype(iat) == iesp ) cycle

    !         do n1 = 1, mat_neighbor(atomtype(iat),iesp)


    !             print*, iat, iesp, n1, cont_pdf(iat,iesp,n1), &
    !                  integratee(contribution_pdf(atomtype(iat),iesp,n1,:),bins,numions(iesp)/V_mean,rmax/real(bins)) 

    !         enddo

    !     enddo
    !     print*, ""

    ! enddo

    ! tot_pdf = sum(contribution_pdf,dim=3)
    ! do iesp = 1, N_species
    !     do jesp = 1, N_species
    !         tot_pdf(iesp,jesp,:) = tot_pdf(iesp,jesp,:)!mat_neighbor(atomtype(iat),iesp)
    !     enddo
    ! enddo


end subroutine total_pdf


subroutine get_mean_sigma_dist( N_pair, dist_distr, mean, sigma, rmax, norm, normalize )
    implicit none
    integer, intent(in)   :: N_pair
    real*8, intent(in)    :: rmax, norm
    real*8, intent(inout) :: dist_distr(:,:)
    real*8, intent(out)   :: mean(:), sigma(:)
    logical, intent(in)   :: normalize
    real*8, parameter     :: pi = acos(-1.0d0)
    real*8                :: x, dx, aux
    integer               :: ipair, i, N   

    N = size(dist_distr,2)
    dx = rmax/N

    mean = 0.0d0
    sigma = 0.0d0


    do ipair = 1, N_pair

        if (normalize) then
            aux = 0.0d0
            do i = 1, N
                x = (i-0.5)*dx
                aux = aux+  dist_distr(ipair,i) * dx!*(4.0*pi*(dx*(i-0.5d0))**2)*norm  !norm = numions(iesp)/V_mean
            enddo
        else
            aux = 1.0d0
        endif

        do i = 1, N
            x = (i-0.5)*dx
            mean(ipair) = mean(ipair) +  dist_distr(ipair,i) * x/aux * dx !*(4.0*pi*(dx*(i-0.5d0))**2)*norm  !norm = numions(iesp)/V_mean
        enddo

        do i = 1, N
            x = (i-0.5)*dx
            sigma(ipair) = sigma(ipair) + ( x/aux-mean(ipair) )**2 * dist_distr(ipair,i) * dx !*(4.0*pi*(dx*(i-0.5d0))**2)*norm 
        enddo
        sigma = sqrt(sigma)

    enddo



end subroutine get_mean_sigma_dist



!************************************
!************************************
subroutine compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,numIons,V,bins,rcut,gdr2)
implicit none
real*8, dimension(:,:), intent(in)      :: dist_matrix
integer, dimension(:,:), intent(in)     :: dist_atoms
integer, dimension(:), intent(in)       :: atomType, numIons
integer, intent(in)                     :: bins, type1, type2
real*8, intent(in)                      :: rcut, V
real*8, dimension(bins), intent(out)    :: gdr2
real*8, dimension(bins)                 :: gdr
integer                                 :: i, j, k
integer                                 :: len
real*8                                  :: x, c, r, pi


c=0.05d0
gdr=0.0d0
pi=acos(-1.0d0)
len = size(dist_matrix,1)
k=0
do i = 1, len
    if (atomType(dist_atoms(i,1))==Type1) then
        if (atomType(dist_atoms(i,2))==Type2) then
            r=norm2(dist_matrix(i,:))
            k=k+1
            do j=1, bins
                x=rcut/real(bins,8)*j
                gdr(j)=gdr(j)+exp(-(x-r)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*pi))
            enddo
        endif
    endif
enddo



gdr(:)=gdr(:)/real(numIons(type1),8)
gdr2=0.0
do i = 2, bins
    !if (i==1) then
        !gdr2(i) = gdr(i)*rcut/real(bins,8)/(4.0*3.1416*((rcut/real(bins,8))**3*(i**3-(i-1)**3))/3.0d0)*V/(norm)
    !else
        gdr2(i) = ((gdr(i)+gdr(i-1))/2.0d0)*(rcut/real(bins,8))/ &
        (4.0*pi*((rcut/real(bins,8))**3*(i**3-(i-1)**3))/3.0d0)*1/(numIons(type2))   ! iria V en lugar del 1
    !endif
enddo

end subroutine
!************************************
!************************************
subroutine makeMatrices(cell,coordinates,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)

implicit none

real*8, dimension(:,:), intent(in)              :: cell, coordinates
integer, dimension(:), intent(in)               :: atomtype, numIons
real*8, intent(in)                              :: Rmax
integer, dimension(:), allocatable, intent(out) :: N
real*8, intent(out)                             :: V
real*8, dimension(:,:), allocatable, intent(out):: dist_matrix
integer                                         :: species, i, j, c, lengthX,&
                                                   lengthY, lengthZ,&
                                                   quit_marker_x, quit_marker_y, &
                                                   quit_marker_z, k, quadrants,&
                                                   current_cell, natoms,&
                                                   marker, basic_cell, Nfull, &
                                                   c2, c3
integer, dimension(sum(numIons))                :: atomType1
real*8, dimension(sum(numIons))                 :: distances
real*8, dimension(13,3)                         :: vect
real*8, dimension(13)                           :: abs_vect
integer, dimension(8,3)                         :: condition, signum
real*8                                          :: x, y, z, Rij
integer, dimension(100000)                      :: N2
integer, dimension(100000,2)                    :: dist_atoms2
integer, dimension(:,:), allocatable, intent(out)   :: dist_atoms
real*8, dimension(100000,3)                     :: dist_matrix2

V = abs(det(cell))
species = size(numIons)
natoms = sum(numIons)
N2 = 0
distances = 0
dist_matrix2 = 0
c = 0
c2 = 1

!do i = 1, species
!    do j = 1, numIons(i)
!        c = c + 1
!        atomType1(c) = atomType(i)
!    enddo
!enddo

c = 0

condition(1,:) = (/0, 0, 0/)
condition(2,:) = (/1, 0, 0/)
condition(3,:) = (/0, 1, 0/)
condition(4,:) = (/0, 0, 1/)
condition(5,:) = (/1, 1, 0/)
condition(6,:) = (/0, 1, 1/)
condition(7,:) = (/1, 0, 1/)
condition(8,:) = (/1, 1, 1/)
signum(1,:) = (/1, 1, 1/)
signum(2,:) = (/-1, 1, 1/)
signum(3,:) = (/1, -1, 1/)
signum(4,:) = (/1, 1, -1/)
signum(5,:) = (/-1, -1, 1/)
signum(6,:) = (/1, -1, -1/)
signum(7,:) = (/-1, 1, -1/)
signum(8,:) = (/-1, -1, -1/)
vect(1,:) = cell(1,:);
vect(2,:) = cell(2,:);
vect(3,:) = cell(3,:);
vect(4,:) = vect(1,:)+vect(2,:);
vect(5,:) = vect(1,:)-vect(2,:);
vect(6,:) = vect(1,:)+vect(3,:);
vect(7,:) = vect(1,:)-vect(3,:);
vect(8,:) = vect(3,:)+vect(2,:);
vect(9,:) = vect(3,:)-vect(2,:);
vect(10,:) = vect(1,:)+vect(2,:)+vect(3,:);
vect(11,:) = vect(1,:)+vect(2,:)-vect(3,:);
vect(12,:) = vect(1,:)-vect(2,:)+vect(3,:);
vect(13,:) = -vect(1,:)+vect(2,:)+vect(3,:);

do i = 1, 13
    abs_vect(i) = sqrt(dot_product(vect(i,:),vect(i,:)))
enddo

lengthX = ceiling((Rmax+maxval(abs_vect))/minval(abs_vect)) + 1
lengthY = lengthX
lengthZ = lengthX
do i = 0, lengthX
    quit_marker_x = 1
    do j = 0, lengthY
        quit_marker_y = 1
        do k = 0, lengthZ
            quit_marker_z = 1
            do quadrants = 1, 8
                if (condition(quadrants,1)*iszero(i) + condition(quadrants,2)&
                    *iszero(j) + condition(quadrants,3)*iszero(k) == 0) then
                    do current_cell = 1, natoms
                        distances = 0.0
                        marker = 0
                        do basic_cell = 1, natoms
                            x = coordinates(current_cell,1) + &
                            signum(quadrants,1)*i - &
                            coordinates(basic_cell,1)
                            y = coordinates(current_cell,2) + &
                            signum(quadrants,2)*j - &
                            coordinates(basic_cell,2)
                            z = coordinates(current_cell,3) + &
                            signum(quadrants,3)*k - &
                            coordinates(basic_cell,3)

                            Rij = (x*cell(1,1)+y*cell(2,1)+z*cell(3,1))**2
                            Rij = Rij + (x*cell(1,2)+y*cell(2,2)+z*cell(3,2))**2
                            Rij = Rij + (x*cell(1,3)+y*cell(2,3)+z*cell(3,3))**2
                            !print*, k, quadrants, current_cell, basic_cell, Rij
                            !print*, x, y, z
                            if (Rij < Rmax**2 .and. Rij > 0.0001) then
                                quit_marker_z = 0
                                quit_marker_y = 0
                                quit_marker_x = 0


                                !print*, i, j, k, current_cell, basic_cell, matmul((/x,y,z/),cell), atomType1(current_cell)
                                !if (marker == 0 .and. Nfull >= natoms) then
                                c = c + 1
                                    !print*, "barruan"
                                N2(c) = atomType(current_cell)
                                dist_atoms2(c,1) = current_cell
                                dist_atoms2(c,2) = basic_cell
                                dist_matrix2(c,:) = matmul((/x,y,z/),cell)

                                !endif
                                !Nfull = Nfull + 1 - marker
                                marker = 1
                                distances(basic_cell) = sqrt(Rij)
                            endif
                        enddo
                        !if (i+j+k+current_cell == 1) then
                        !    dist_matrix2(1,:) = distances(:)
                        !    typ_j2(1) = typ_j_coef
                        !    c3 = 1
                        !elseif (marker == 1) then
                        !    c2 = c2 + 1
                        !    dist_matrix2(c2,:) = distances(:)
                        !    c3 = c3 + 1
                        !    typ_j2(c3) = typ_j_coef
                        !endif
                    enddo
                endif
            enddo
            if (quit_marker_z == 1) then
                exit
            endif
        enddo
        if (quit_marker_y == 1) then
            exit
        endif
    enddo
    if (quit_marker_x == 1) then
        exit
    endif
enddo

allocate(N(c),dist_matrix(c,3),dist_atoms(c,2))

do i = 1, c
    N(i) = N2(i)
    dist_matrix(i,:) = dist_matrix2(i,:)
    dist_atoms(i,:) = dist_atoms2(i,:)
enddo


end subroutine
!************************************
!************************************
FUNCTION DET (A) RESULT (ans)
IMPLICIT NONE
REAL*8, DIMENSION(:,:), INTENT(IN)  :: A
REAL*8  :: ans

ans =   A(1,1)*A(2,2)*A(3,3)  &
    - A(1,1)*A(2,3)*A(3,2)  &
    - A(1,2)*A(2,1)*A(3,3)  &
    + A(1,2)*A(2,3)*A(3,1)  &
    + A(1,3)*A(2,1)*A(3,2)  &
    - A(1,3)*A(2,2)*A(3,1)

END FUNCTION
!************************************
!************************************
FUNCTION ISZERO (N) RESULT (ANS)
IMPLICIT NONE
INTEGER, INTENT(IN)     :: N
INTEGER                 :: ANS

IF (N == 0) THEN
    ANS = 1
ELSE
    ANS = 0
ENDIF

END FUNCTION
!************************************
end module
