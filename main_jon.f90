program main
use read_write
use pair_dist
use mod_check_min
use mod_angles
use options_main

implicit none

integer                                 :: natoms, bins
integer, dimension(:), allocatable      :: numIons, atomtype
real*8, dimension(:,:), allocatable     :: coor
real*8                                  :: rmax, V_mean
integer                                 :: N_species, N_file

integer, allocatable                 :: neighbor_list(:,:,:,:), N_neighbor(:,:), &
                                        neighbor_order_list(:,:,:), &
                                        cont_pdf(:,:,:), cont_adf(:,:,:), &
                                        mat_neighbor(:,:), mat_pairs(:,:)
real*8, allocatable                  :: mat_rcut(:,:), mat_pdf(:,:,:,:), mat_adf(:,:,:,:), &
                                        tot_pdf(:,:,:), tot_pdf_compute_neigh(:,:,:), &
                                        contribution_pdf(:,:,:,:), contribution_adf(:,:,:,:), tot_adf(:,:,:), &
                                        mean_pdf(:,:,:), sigma_pdf(:,:,:), mean_adf(:,:,:), sigma_adf(:,:,:)

character(len=2), allocatable   :: species(:)
character(13) :: trjfile

logical :: compute_first_neighbor, plot_results, compute_deviations, check_constraints
integer :: ext
real*8  :: max_sigma_angle

integer :: inunit, outunit


open(unit = 1, action = "read", status = "old", file = "input")
read(1,*) trjfile
read(1,*) N_species
allocate(species(N_species))
read(1,*) species(:)
read(1,*) rmax
read(1,*) bins
read(1,*) ext
read(1,*) 
read(1,*) compute_first_neighbor
read(1,*) plot_results
read(1,*) compute_deviations
read(1,*) 
read(1,*) check_constraints
read(1,*) max_sigma_angle
close(unit = 1)

species = (/ "Si", "O " /)
inunit = 123
outunit = 234


call atom_number(trjfile,natoms)
allocate(coor(natoms,3),atomtype(natoms), numions(N_species), N_neighbor(natoms,N_species))

allocate(mat_neighbor(N_species, N_species),mat_rcut(N_species, N_species))


if (.not. check_constraints) then
    ! 1. READ TRAJECTORY AND COMPUTE RADIAL DISTRIBUTION FUNCTIONS 
    ! ------------------------------------------------------------
    open(unit=inunit,file=trjfile,status='old',action='read')
    open(unit=outunit,file="output",status='replace',action='write')

    !print*, "0"

    if (compute_first_neighbor) then
        mat_neighbor = 0
        mat_rcut = 0.0
        allocate(tot_pdf_compute_neigh(N_species,N_species,bins))
        call get_mat_rcut_neighbor(inunit, outunit, natoms, N_species, bins, Rmax, &
                                    N_file, mat_neighbor, mat_rcut, tot_pdf_compute_neigh)
    else
        call get_N_file(inunit, outunit, natoms, N_file)
        mat_neighbor = 6
        !mat_pairs = 15
    endif

    !print*, "1"


    ! ------------------------------------------------------------



    ! 3. READ FIRST CONFIGURATION AND GET THE NEIGHBOR TAGS
    ! ------------------------------------------------------------
    call get_neighbor_tags ( inunit, outunit, natoms, N_species, bins, ext, rmax, mat_neighbor, &
                             neighbor_order_list, mat_pdf, tot_pdf, contribution_pdf, cont_pdf, &
                             mat_adf, tot_adf, contribution_adf, cont_adf, &
                             numIons, atomtype  )

    !print*, "2"

    ! ------------------------------------------------------------



    ! 4. READ TRAJECTORY AND COMPUTE THE DISTRIBUTION OF EACH ANGLE
    ! ------------------------------------------------------------

    call get_distributions_dist_angles ( inunit, outunit, N_file, natoms, N_species, bins, ext, rmax, mat_neighbor, &
                                         neighbor_order_list, mat_pdf, tot_pdf, contribution_pdf, cont_pdf, &
                                         mat_adf, tot_adf, contribution_adf, cont_adf, V_mean )
    ! ------------------------------------------------------------

    !print*, "3"




    ! 5. COMPUTE THE STANDAR DEVIATION OF EACH RADIAL AND ANGLE DISTRIBUTION
    ! ------------------------------------------------------------

    if (compute_deviations) then
        call get_deviation_each_dist_angle(outunit, natoms, N_species, atomtype, numions, rmax, V_mean, ext, mat_neighbor, &
                                         mat_pdf, mat_adf, sigma_pdf, sigma_adf)
    endif

    ! ------------------------------------------------------------


    !print*, "4"


  
    call total_pdf(ext, V_mean, N_file, rmax, mat_pdf, mat_neighbor, neighbor_list, &
                    atomtype, numIons, contribution_pdf, tot_pdf, cont_pdf)

    call total_adf(ext, mat_adf, mat_neighbor, neighbor_list, atomtype, numIons, contribution_adf, tot_adf, cont_adf)


    call get_deviation_contribution(outunit, natoms, N_species, numions, rmax, V_mean, ext, mat_neighbor, &
                                      contribution_pdf, contribution_adf)


    call dist_distr2pdf(ext, natoms, N_species, mat_neighbor, rmax, V_mean, numions, atomtype, contribution_pdf, tot_pdf)


    if (compute_first_neighbor) then
        call write_plot_contribution_total(outunit, N_species, mat_neighbor, ext, rmax, V_mean, numions, plot_results, &
                                              contribution_pdf, tot_pdf_compute_neigh, contribution_adf, tot_adf)
    else
        call write_plot_contribution_total(outunit, N_species, mat_neighbor, ext, rmax, V_mean, numions, plot_results, &
                                              contribution_pdf, tot_pdf, contribution_adf, tot_adf)
    endif

    close(unit=outunit)
else
    call compute_number_constraints(max_sigma_angle)
endif

end program
