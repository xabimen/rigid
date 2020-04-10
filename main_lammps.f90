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
integer, allocatable                 :: neighbor_list(:,:,:,:), N_neighbor(:,:), &
                                        neighbor_order_list(:,:,:)
real*8, allocatable                  :: mat_rcut(:,:), mat_theta(:,:,:), mat_pdf(:,:,:), mat_adf(:,:,:,:), &
                                        mean(:,:,:), sigma(:,:,:), mat_neighbor(:,:)

character(len=2), allocatable   :: species(:)
character(len=100), allocatable :: file_list(:)
character(1) :: caux1, caux2
character(13) :: infile

infile = "qe1.lammpstrj"
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
print*, "1. READ TRAJECTORY AND COMPUTE RADIAL DISTRIBUTION FUNCTIONS"
call atom_number(infile,natoms)
allocate(coor(natoms,3),atomtype(natoms), numions(N_species),N_neighbor(natoms,N_species))
open(unit=123,file=infile,status='old',action='read')
call read_trj(123,cell,coor,numIons,atomtype,io)

print*, numIons(:)
open(unit=124,file="kaka_POSCAR",status='replace',action='write')
call write_vasp(124,cell,coor)


N_file = 0
do
    call read_trj(123,cell,coor,numIons,atomtype,io)
    if (io < 0) exit
    if(N_file == 10) exit

    N_file = N_file + 1
    print*, "READING FILE ", N_file

    call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)

    do type1 = 1, N_species

        do type2 = 1, N_species

            if (type1 == type2) cycle

            ! Compute radial distribution function
            call compute_gdr(dist_matrix,dist_atoms,atomType,type1,type2,numIons,V,bins,rmax,pdf)
            
            mat_pdf(type1,type2,:) = mat_pdf(type1,type2,:) + pdf(:)
        enddo
    enddo

    deallocate(dist_matrix, dist_atoms )

enddo 

do type1 = 1, N_species
    do type2 = 1, N_species
        if (type1==type2) cycle
        mat_pdf(type1,type2,:) = mat_pdf(type1,type2,:)/real(N_file,8)
    enddo
enddo
rewind(unit=123)
!print*, N_file
!print*, "kaka"
! ------------------------------------------------------------


! 2. FIND THE MINIMUM OF THE RDF FOR EACH PAIR OF SPECIES
! ------------------------------------------------------------
print*, "2. FIND THE MINIMUM OF THE RDF FOR EACH PAIR OF SPECIES"
do type1 = 1, N_species

    do type2 = 1, N_species
        if (type1 == type2) cycle

        write(caux1,"(i1)") type1
        write(caux2,"(i1)") type2
        call write_gdr( "gdr"//caux1//caux2//".dat", mat_pdf(type1,type2,:), rmax )

        ! Find the first minimum
        call check_min(mat_pdf(type1,type2,:), int(bins*0.5/rmax), min_pdf)
        mat_rcut(type1,type2) = rmax/bins*min_pdf
        mat_neighbor(type1,type2) = integrate(mat_pdf(type1,type2,:),min_pdf,numions(type2)/V,rmax/real(bins))
        write(*,"(2i5,2f10.6)") type1, type2, mat_rcut(type1,type2), mat_neighbor(type1,type2)
    enddo
enddo 
! ------------------------------------------------------------

! 3. READ FIRST CONFIGURATION AND GET THE NEIGHBOR TAGS
! ------------------------------------------------------------
print*, "3. READ FIRST CONFIGURATION AND GET THE NEIGHBOR TAGS"
call read_trj(123,cell,coor,numIons,atomtype,io)
call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
call get_neighbor_list2( dist_matrix, dist_atoms, natoms, N_species, atomtype,  ceiling(mat_neighbor),N_neighbor, neighbor_list )

allocate(neighbor_order_list(natoms,N_species,maxval(ceiling(mat_neighbor))))
do i = 1, natoms
    do type1 = 1, N_species
        do j = 1, mat_neighbor(atomtype(i),type1)
            neighbor_order_list(i,type1,j) = neighbor_list(i,type1,j,2)    
        enddo
    enddo
enddo
deallocate( dist_matrix, dist_atoms, neighbor_list )


max_pair = maxval( mat_neighbor ) * ( maxval(mat_neighbor)-1 )/2
allocate( mat_adf(natoms,N_species,max_pair,bins) )
mat_adf = 0.0d0
! ------------------------------------------------------------


! 4. READ TRAJECTORY AND COMPUTE THE DISTRIBUTION OF EACH ANGLE
! ------------------------------------------------------------
print*, "4. READ TRAJECTORY AND COMPUTE THE DISTRIBUTION OF EACH ANGLE"
rewind(unit = 123)
do ifile = 1, 100

    call read_trj(123,cell,coor,numIons,atomtype,io)
    if (io < 0) exit

    call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
    call get_neighbor_list2( dist_matrix, dist_atoms, natoms, N_species, atomtype, ceiling(mat_neighbor),N_neighbor, neighbor_list )

    call update_angle_distr(neighbor_order_list, atomtype, ceiling(mat_neighbor), dist_matrix, &
                              neighbor_list, mat_adf )
    
    deallocate( dist_matrix, dist_atoms, neighbor_list )

enddo
close(unit = 123)
! ------------------------------------------------------------


! 5. COMPUTE THE STANDAR DEVIATION OF EACH ANGLE DISTRIBUTION
! ------------------------------------------------------------
print*, "5. COMPUTE THE STANDAR DEVIATION OF EACH ANGLE DISTRIBUTION"
allocate( mean(natoms,N_species,max_pair), sigma(natoms,N_species,max_pair) )

do iat = 1, natoms
    do iesp = 1, N_species

        N_pair = mat_neighbor(atomtype(iat),iesp)*(mat_neighbor(atomtype(iat),iesp)-1)/2

        if ( atomtype(iat) == iesp ) cycle


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
