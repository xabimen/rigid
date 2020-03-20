! gfortran read_write.f90 check_min.f90 angles.f90 pair_dist.f90 main.f90

program main
use read_write
use pair_dist
use mod_check_min
use mod_angles

implicit none

integer                                 :: natoms, io, bins, i
integer                                 :: type1, type2
integer, dimension(:), allocatable      :: numIons, atomtype
real*8, dimension(:), allocatable       :: pdf, angle_distr
real*8, dimension(3,3)                  :: cell
real*8, dimension(:,:), allocatable     :: coor
real*8, parameter                       :: pi = acos(-1.0d0)
real*8                                  :: rmax, V
integer                                 :: N_species, min_pdf

integer, dimension(:), allocatable   :: N
integer, dimension(:,:), allocatable :: dist_atoms
real*8, dimension(:,:), allocatable  :: dist_matrix
integer, allocatable                 :: neighbor_list(:,:), N_neighbor(:)
real*8, allocatable                  :: mat_rcut(:,:), mat_theta(:,:,:)

rmax = 8.0
bins = 100
type1 = 3
type2 = 3

call atom_number("input",natoms)
allocate(coor(natoms,3),atomtype(natoms), N_neighbor(natoms),&
         pdf(bins), angle_distr(bins))


open(unit=123,file="input",status='old',action='read')
open(unit=124,file="kaka_POSCAR",status='replace',action='write')

call read_trj(123,cell,coor,numIons,atomtype,io)
call write_vasp(124,cell,coor)

N_species = size(numIons)
allocate(mat_rcut(N_species,N_species), mat_theta(N_species,N_species,2))


call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
call get_neighbor_list( dist_matrix, dist_atoms, natoms, N_neighbor, neighbor_list )


do type1 = 1, N_species

    do type2 = type1, N_species

        ! Compute radial distribution function
        call compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,numIons(type1)*numIons(type2),V,bins,rmax,pdf)
        
        ! Find the first minimum
        call check_min(pdf, int(bins*0.5/rmax), min_pdf)
        mat_rcut(type1,type2) = rmax/bins*min_pdf
        mat_rcut(type2,type1) = rmax/bins*min_pdf

        ! Compute angle distribution function
        ! Angles 1-2-1
        call compute_angle_distr( dist_matrix, dist_atoms, neighbor_list, N_neighbor, &
                                  atomtype, type1, type2, bins, rmax/bins*min_pdf, angle_distr)

        ! Find angle range
        ! call check_angle_range()
        ! mat_theta(type1,type2,:) = (/theta_min, theta_max/)


        ! Angles 2-1-2
        !call compute_angle_distr( dist_matrix, dist_atoms, neighbor_list, N_neighbor, &
        !                          atomtype, type2, type1, bins, rmax/bins*min_pdf, angle_distr)

        print*, type1, type2, mat_rcut(type1,type2)

    enddo
enddo


! type1 = 1
! type2 = 3
! call compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,numIons(type1)*numIons(type2),V,bins,rmax,pdf)
! call compute_angle_distr( dist_matrix, dist_atoms, neighbor_list, N_neighbor, &
!                 atomtype, type1, type2, bins, mat_rcut(type1,type2), angle_distr)

! print*, "# Minimum: ", rmax/bins*min_pdf
! do i = 1, bins
!    print*, rmax/real(bins)*i, pdf(i)
! enddo

! ! do i = 1, bins
! !     if (angle_distr(i) > 1.0D-8) then
! !         print*, pi/bins*i, angle_distr(i)
! ! 	endif
! ! enddo



end program
