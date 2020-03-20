module pair_dist

use mod_angles


contains
!************************************
subroutine pair_distribution(cell,coor,numIons,atomtype,type1,type2,rmax,bins,pdf,angle_distr)
implicit NONE
real*8, dimension(:,:), intent(in)                 :: cell, coor
integer, dimension(:), intent(in)                  :: atomtype, numIons
real*8, intent(in)                                 :: Rmax
integer, intent(in)                                :: bins, type1, type2
real*8, dimension(:), allocatable, intent(out)     :: pdf, angle_distr
integer, dimension(:), allocatable                 :: N
integer, dimension(:,:), allocatable               :: dist_atoms
real*8                                             :: V
real*8, dimension(:,:), allocatable                :: dist_matrix

integer, allocatable :: neighbor_list(:,:)
integer              :: N_neighbor(size(atomtype))

allocate(pdf(bins), angle_distr(bins))

call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
call compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,numIons(type1)*numIons(type2),V,bins,rmax,pdf)

call get_neighbor_list( dist_matrix, dist_atoms, size(atomtype), &
                        N_neighbor, neighbor_list )
call compute_angle_distr( dist_matrix, dist_atoms, neighbor_list, N_neighbor, &
                          atomtype, type1, type2, bins, angle_distr)


end subroutine
!************************************
!************************************
subroutine compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,norm,V,bins,rcut,gdr2)
implicit none
real*8, dimension(:,:), intent(in)      :: dist_matrix
integer, dimension(:,:), intent(in)     :: dist_atoms
integer, dimension(:), intent(in)       :: atomType
integer, intent(in)                     :: bins, type1, type2, norm
real*8, intent(in)                      :: rcut, V
real*8, dimension(bins), intent(out)    :: gdr2
real*8, dimension(bins)                 :: gdr
real*8                                  :: r
integer                                 :: i, j, k
integer                                 :: len
real*8                                  :: x, c


c=0.2d0
gdr=0.0d0
len = size(dist_matrix,1)
k=0
do i = 1, len
    if (atomType(dist_atoms(i,1))==Type1) then
        if (atomType(dist_atoms(i,2))==Type2) then
            r=norm2(dist_matrix(i,:))
            k=k+1
            do j=1, bins
                x=rcut/real(bins,8)*j
                gdr(j)=gdr(j)+exp(-(x-r)**2/(2.0d0*c**2))/(c*sqrt(2.0d0*3.1416))
            enddo
        endif
    endif
enddo



do i = 1, bins
    if (i==1) then
        gdr2(i) = gdr(i)*rcut/real(bins,8)/(4.0*3.1416*((rcut/real(bins,8))**3*(i**3-(i-1)**3))/3.0d0)*V/(real(4)**2)
    else
        gdr2(i) = ((gdr(i)+gdr(i-1))/2.0d0)*(rcut/real(bins,8))/ &
        (4.0*3.1416*((rcut/real(bins,8))**3*(i**3-(i-1)**3))/3.0d0)*V/(real(4)**2)
    endif
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
