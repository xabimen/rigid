! gfortran read_write.f90 check_min.f90 angles.f90 pair_dist.f90 main.f90

program main
use read_write
use pair_dist
use mod_check_min
use mod_angles

implicit none

integer                                 :: natoms, io, bins, i, j, iat, iesp
integer                                 :: type1, type2
integer, dimension(:), allocatable      :: numIons, atomtype, aux_numions, aux_atomtype
real*8, dimension(:), allocatable       :: pdf, angle_distr
real*8, dimension(3,3)                  :: cell
real*8, dimension(:,:), allocatable     :: coor
real*8, parameter                       :: pi = acos(-1.0d0)
real*8                                  :: rmax, V
integer                                 :: N_species, min_pdf, N_file, ifile, max_pair, N_pair

integer, dimension(:), allocatable   :: N
integer, dimension(:,:), allocatable :: dist_atoms
real*8, dimension(:,:), allocatable  :: dist_matrix
integer, allocatable                 :: neighbor_list(:,:,:,:), N_neighbor(:,:), mat_neighbor(:,:), &
                                        neighbor_order_list(:,:,:)
real*8, allocatable                  :: mat_rcut(:,:), mat_theta(:,:,:), mat_pdf(:,:,:), mat_adf(:,:,:,:), &
                                        mean(:,:,:), sigma(:,:,:)

character(len=2), allocatable   :: species(:)
character(len=100), allocatable :: file_list(:)
character(1) :: caux1, caux2

rmax = 5.0
bins = 200


allocate(species(2))
species = (/ "Si", "O " /)

N_species = size(species)
allocate(aux_numions(N_species))
allocate(mat_rcut(N_species,N_species), mat_neighbor(N_species, N_species), mat_theta(N_species,N_species,2))
allocate(mat_pdf(N_species,N_species,bins))
allocate(pdf(bins), angle_distr(bins))
mat_pdf = 0.0d0



! 1. READ TRAJECTORY AND COMPUTE RADIAL DISTRIBUTION FUNCTIONS 
! ------------------------------------------------------------
call get_all_files( "trajectory", N_file, file_list )
do ifile = 1, N_file

    call read_xsf(file_list(ifile), species, natoms, cell, coor, numions, atomtype, io)

    aux_numions = numions

    call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)

    do type1 = 1, N_species

        do type2 = 1, N_species

            if (type1 == type2) cycle

            ! Compute radial distribution function
            call compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,numIons,V,bins,rmax,pdf)
            
            mat_pdf(type1,type2,:) = mat_pdf(type1,type2,:) + pdf(:)/N_file
        enddo
    enddo

    deallocate( numions, atomtype, dist_matrix, dist_atoms )
enddo 
! ------------------------------------------------------------


! 2. FIND THE MINIMUM OF THE RDF FOR EACH PAIR OF SPECIES
! ------------------------------------------------------------
do type1 = 1, N_species

    do type2 = 1, N_species

        if (type1 == type2) cycle

        ! Find the first minimum
        call check_min(mat_pdf(type1,type2,:), int(bins*0.5/rmax), min_pdf)
        mat_rcut(type1,type2) = rmax/bins*min_pdf
        mat_neighbor(type1,type2) = ceiling( integrate(mat_pdf(type1,type2,:),min_pdf,aux_numions(type2)/V,rmax/real(bins)) )

    enddo
enddo 
! ------------------------------------------------------------


! 3. READ FIRST CONFIGURATION AND GET THE NEIGHBOR TAGS
! ------------------------------------------------------------
call read_xsf(file_list(1), species, natoms, cell, coor, numions, atomtype, io)
allocate(N_neighbor(natoms,N_species), aux_atomtype(natoms))
call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
call get_neighbor_list2( dist_matrix, dist_atoms, natoms, N_species, atomtype,  mat_neighbor,N_neighbor, neighbor_list )
aux_atomtype = atomtype

allocate(neighbor_order_list(natoms,N_species,maxval(mat_neighbor)))
do i = 1, natoms
    do type1 = 1, N_species
        do j = 1, mat_neighbor(atomtype(i),type1)
            neighbor_order_list(i,type1,j) = neighbor_list(i,type1,j,2)    
        enddo
    enddo
enddo
deallocate( numions, atomtype, dist_matrix, dist_atoms, N_neighbor, neighbor_list )


max_pair = maxval( mat_neighbor ) * ( maxval(mat_neighbor)-1 )/2
allocate( mat_adf(natoms,N_species,max_pair,bins) )
mat_adf = 0.0d0
! ------------------------------------------------------------


! 4. READ TRAJECTORY AND COMPUTE THE DISTRIBUTION OF EACH ANGLE
! ------------------------------------------------------------
do ifile = 1, N_file

    call read_xsf(file_list(ifile), species, natoms, cell, coor, numions, atomtype, io)
    allocate(N_neighbor(natoms,N_species))

    call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
    call get_neighbor_list2( dist_matrix, dist_atoms, natoms, N_species, atomtype,  mat_neighbor,N_neighbor, neighbor_list )

    call update_angle_distr(neighbor_order_list, atomtype, mat_neighbor, dist_matrix, &
                              neighbor_list, mat_adf )
    
    deallocate( numions, atomtype, dist_matrix, dist_atoms, N_neighbor, neighbor_list )

enddo
! ------------------------------------------------------------


! 5. COMPUTE THE STANDAR DEVIATION OF EACH ANGLE DISTRIBUTION
! ------------------------------------------------------------
allocate( mean(natoms,N_species,max_pair), sigma(natoms,N_species,max_pair) )

do iat = 1, natoms
    do iesp = 1, N_species

        N_pair = mat_neighbor(aux_atomtype(iat),iesp)*(mat_neighbor(aux_atomtype(iat),iesp)-1)/2

        if ( aux_atomtype(iat) == iesp ) cycle


        call get_mean_sigma( N_pair, mat_adf(iat,iesp,:,:), mean(iat,iesp,:), sigma(iat,iesp,:) )
        write(*,"(2i5,100f10.3)") iat,iesp, mean(iat,iesp,:N_pair)
        write(*,"(2i5,100f10.3)") iat,iesp, sigma(iat,iesp,:N_pair)
        write(*,*)

        write(caux1,"(i1)") iat
        write(caux2,"(i1)") iesp

        call write_angle_distr_full( "adf_at"//caux1//"_esp"//caux2//".dat", &
                                      N_pair, mat_adf(iat,iesp,:,:) )

    enddo
enddo
! ------------------------------------------------------------


end program
