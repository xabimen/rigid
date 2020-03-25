! gfortran read_write.f90 check_min.f90 angles.f90 pair_dist.f90 main.f90

program main
use read_write
use pair_dist
use mod_check_min
use mod_angles

implicit none

integer                                 :: natoms, io, bins, i
integer                                 :: type1, type2
integer, dimension(:), allocatable      :: numIons, atomtype, aux_numions
real*8, dimension(:), allocatable       :: pdf, angle_distr
real*8, dimension(3,3)                  :: cell
real*8, dimension(:,:), allocatable     :: coor
real*8, parameter                       :: pi = acos(-1.0d0)
real*8                                  :: rmax, V
integer                                 :: N_species, min_pdf, N_file, ifile

integer, dimension(:), allocatable   :: N
integer, dimension(:,:), allocatable :: dist_atoms
real*8, dimension(:,:), allocatable  :: dist_matrix
integer, allocatable                 :: neighbor_list(:,:), N_neighbor(:)
real*8, allocatable                  :: mat_rcut(:,:), mat_theta(:,:,:), mat_pdf(:,:,:), mat_adf(:,:,:)

character(len=2), allocatable   :: species(:)
character(len=100), allocatable :: file_list(:)
character(1) :: caux1, caux2

rmax = 5.0
bins = 1000
type1 = 3
type2 = 3

!*****************************************************************
!call atom_number("input",natoms)
!allocate(coor(natoms,3),atomtype(natoms), N_neighbor(natoms),&
!         pdf(bins), angle_distr(bins))

!open(unit=123,file="input",status='old',action='read')
!open(unit=124,file="kaka_POSCAR",status='replace',action='write')

!call read_trj(123,cell,coor,numIons,atomtype,io)
!call write_vasp(124,cell,coor)
!*****************************************************************


!*****************************************************************
allocate(species(2))
species = (/ "Si", "O " /)
!*****************************************************************


N_species = size(species)
allocate(aux_numions(N_species))
allocate(mat_rcut(N_species,N_species), mat_theta(N_species,N_species,2))
allocate(mat_pdf(N_species,N_species,bins), mat_adf(N_species,N_species,bins))
allocate(pdf(bins), angle_distr(bins))
mat_adf = 0.0d0
mat_pdf = 0.0d0



call get_all_files( "trajectory", N_file, file_list )



do ifile = 1, N_file

    call read_xsf(file_list(ifile), species, natoms, cell, coor, numions, atomtype, io)
    allocate(N_neighbor(natoms))

    aux_numions = numions

    call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
    call get_neighbor_list( dist_matrix, dist_atoms, natoms, N_neighbor, neighbor_list )


    do type1 = 1, N_species

        do type2 = 1, N_species

            ! Compute radial distribution function
            call compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,numIons,V,bins,rmax,pdf)
            
            mat_pdf(type1,type2,:) = mat_pdf(type1,type2,:) + pdf(:)/N_file
        enddo
    enddo

    deallocate( numions, atomtype, dist_matrix, dist_atoms, N_neighbor, neighbor_list )

enddo


do type1 = 1, N_species

    do type2 = 1, N_species

        ! Find the first minimum
        call check_min(mat_pdf(type1,type2,:), int(bins*0.5/rmax), min_pdf)
        mat_rcut(type1,type2) = rmax/bins*min_pdf

    enddo
enddo


do ifile = 1, N_file

    call read_xsf(file_list(ifile), species, natoms, cell, coor, numions, atomtype, io)
    allocate(N_neighbor(natoms))

    call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
    call get_neighbor_list( dist_matrix, dist_atoms, natoms, N_neighbor, neighbor_list )

    do type1 = 1, N_species

        do type2 = 1, N_species
            ! Compute angle distribution function
            ! Angles 1-2-1
            call compute_angle_distr( dist_matrix, dist_atoms, neighbor_list, N_neighbor, &
                                      atomtype, type1, type2, bins, rmax/bins*min_pdf, angle_distr)

            mat_adf(type1,type2,:) = mat_adf(type1,type2,:) + angle_distr(:)/N_file


            ! Find angle range
            ! call check_angle_range()
            ! mat_theta(type1,type2,:) = (/theta_min, theta_max/)


        enddo
    enddo

    deallocate( numions, atomtype, dist_matrix, dist_atoms, N_neighbor, neighbor_list )

enddo


do type1 = 1, N_species

    do type2 = 1, N_species

        min_pdf = nint(bins*mat_rcut(type1,type2)/rmax)

        write(*,"(2i5,2f12.6)") type1, type2, mat_rcut(type1,type2), &
         integrate(mat_pdf(type1,type2,:),min_pdf,aux_numions(type2)/V,rmax/real(bins)) 

        write(caux1,"(i1)") type1
        write(caux2,"(i1)") type2
        call write_gdr( "rdf"//caux1//caux2//".dat", mat_pdf(type1,type2,:), rmax )

        call write_angle_distr( "adf"//caux1//caux2//".dat", mat_adf(type1,type2,:) )

    enddo
enddo


end program
