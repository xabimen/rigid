module options_main

use read_write
use pair_dist
use mod_angles
use mod_check_min

contains

subroutine get_mat_rcut_neighbor(inunit, outunit, natoms, N_species, bins, Rmax, N_file, mat_neighbor, mat_rcut)
    implicit none
    integer, intent(in)    :: inunit, outunit, natoms, N_species, bins
    real*8, intent(in)     :: Rmax
    real*8, intent(inout)  :: mat_rcut(:,:)
    integer, intent(inout) ::mat_neighbor(:,:)
    integer, intent(out)   :: N_file

    real*8, allocatable  :: dist_matrix(:,:)
    integer, allocatable :: dist_atoms(:,:), N(:)

    integer :: numIons(N_species), atomtype(natoms)
    real*8  :: mat_pdf(N_species,N_species,bins), pdf(bins), &
               cell(3,3), coor(natoms,3)
    real*8  :: V, aux
    integer :: io, type1, type2, min_pdf
    character(4) :: caux1, caux2

    write(outunit,*) "1. READ TRAJECTORY AND COMPUTE RADIAL DISTRIBUTION FUNCTIONS"

    N_file = 0
    do
        call read_trj(inunit,cell,coor,numIons,atomtype,io)
        if (io < 0) exit
        if(N_file == 10) exit

        N_file = N_file + 1
        write(outunit,*) "READING FILE ", N_file

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
    rewind(unit=inunit)
    !print*, N_file
    !print*, "kaka"
    ! ------------------------------------------------------------


    ! 2. FIND THE MINIMUM OF THE RDF FOR EACH PAIR OF SPECIES
    ! ------------------------------------------------------------
    write(unit=outunit,fmt=*) "2. FIND THE MINIMUM OF THE RDF FOR EACH PAIR OF SPECIES"
    do type1 = 1, N_species

        do type2 = 1, N_species
            if (type1 == type2) cycle

            write(caux1,"(i1)") type1
            write(caux2,"(i1)") type2
            call write_gdr( "gdr"//trim(adjustl(caux1))//trim(adjustl(caux2))//".dat", mat_pdf(type1,type2,:), rmax )

            ! Find the first minimum
            call check_min(mat_pdf(type1,type2,:), int(bins*0.5/rmax), min_pdf)
            mat_rcut(type1,type2) = rmax/bins*min_pdf
            aux = integrate(mat_pdf(type1,type2,:),min_pdf,numions(type2)/V,rmax/real(bins))
            mat_neighbor(type1,type2) = nint(aux)
            write(outunit,"(2i5,2f10.6)") type1, type2, mat_rcut(type1,type2), aux
        enddo
    enddo 

end subroutine get_mat_rcut_neighbor



subroutine get_N_file(inunit, outunit, natoms, N_file)
    implicit none
    integer, intent(in)  :: inunit, outunit, natoms
    integer, intent(out) :: N_file
    integer              :: io, i
    N_file = 0
    do 
        read(inunit, fmt=*, iostat=io)

        if (io < 0) exit
        do i = 1, 4
            read(inunit,fmt=*)
        enddo

        do i = 1, 3
            read(inunit,fmt=*) 
        enddo

        read(inunit,fmt=*)

        do i = 1, natoms
            read(inunit,fmt=*) 
        enddo
        N_file = N_file + 1
    enddo

    rewind(unit=inunit)

end subroutine



subroutine get_neighbor_tags(inunit, outunit, natoms, N_species, bins, ext, rmax, mat_neighbor, &
                             neighbor_order_list, mat_pdf, tot_pdf, contribution_pdf, cont_pdf, &
                             mat_adf, tot_adf, contribution_adf, cont_adf, &
                             numIons, atomtype )
    implicit none
    integer, intent(in)               :: inunit, outunit, natoms, N_species, bins, ext, mat_neighbor(:,:)
    real*8, intent(in)                :: rmax
    real*8, intent(out), allocatable  :: mat_pdf(:,:,:,:), tot_pdf(:,:,:), contribution_pdf(:,:,:,:), &
                                         mat_adf(:,:,:,:), tot_adf(:,:,:), contribution_adf(:,:,:,:)
    integer, intent(out), allocatable :: neighbor_order_list(:,:,:), cont_pdf(:,:,:), cont_adf(:,:,:)
    integer, intent(out)              :: numIons(N_species), atomtype(natoms)

    real*8, allocatable  :: dist_matrix(:,:)
    integer, allocatable :: dist_atoms(:,:), N(:), neighbor_list(:,:,:,:)

    integer :: N_neighbor(natoms, N_species), io
    real*8  :: pdf(bins), cell(3,3), coor(natoms,3)
    real*8  :: V

    integer :: max_neigh, max_pair, i, j, type1, type2


    write(outunit,*)  "3. READ FIRST CONFIGURATION AND GET THE NEIGHBOR TAGS"

    call read_trj(inunit,cell,coor,numIons,atomtype,io)
    call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
    call get_neighbor_list2( ext, dist_matrix, dist_atoms, natoms, N_species, atomtype, &
                             mat_neighbor,N_neighbor, neighbor_list )


    allocate(neighbor_order_list(natoms,N_species,maxval(mat_neighbor)+ext))
    do i = 1, natoms
        do type1 = 1, N_species
            do j = 1, mat_neighbor(atomtype(i),type1)+ext
                neighbor_order_list(i,type1,j) = neighbor_list(i,type1,j,2)    
            enddo
        enddo
    enddo
    deallocate( dist_matrix, dist_atoms, neighbor_list )


    max_neigh = maxval(mat_neighbor)!+ext
    max_pair = max_neigh*(max_neigh-1)/2! + ext
    max_neigh = max_neigh + ext

    allocate( mat_adf(natoms,N_species,max_pair,bins), mat_pdf(natoms,N_species,max_neigh,bins), &
              tot_pdf(N_species, N_species, bins), contribution_pdf(N_species,N_species,max_neigh,bins), &
              cont_pdf(natoms,N_species,max_neigh), cont_adf(natoms,N_species,max_pair), &
              tot_adf(N_species, N_species, bins), contribution_adf(N_species,N_species,max_pair,bins) )
    mat_adf = 0.0d0
    cont_pdf = 0
    cont_adf = 0
    contribution_pdf = 0.0d0
    contribution_adf = 0.0d0

end subroutine get_neighbor_tags



subroutine get_distributions_dist_angles(inunit, outunit, N_file, natoms, N_species, bins, ext, rmax, mat_neighbor, &
                                         neighbor_order_list, mat_pdf, tot_pdf, contribution_pdf, cont_pdf, &
                                         mat_adf, tot_adf, contribution_adf, cont_adf )
    integer, intent(in)    :: inunit, outunit, natoms, N_species, bins, ext, N_file, mat_neighbor(:,:)
    real*8, intent(in)     :: rmax
    real*8, intent(inout)  :: mat_pdf(:,:,:,:), tot_pdf(:,:,:), contribution_pdf(:,:,:,:), &
                              mat_adf(:,:,:,:), tot_adf(:,:,:), contribution_adf(:,:,:,:)
    integer, intent(inout) :: neighbor_order_list(:,:,:), cont_pdf(:,:,:), cont_adf(:,:,:)

    real*8, allocatable  :: dist_matrix(:,:)
    integer, allocatable :: dist_atoms(:,:), N(:), neighbor_list(:,:,:,:)

    integer :: numIons(N_species), atomtype(natoms), N_neighbor(natoms, N_species), io
    real*8  :: cell(3,3), coor(natoms,3)
    real*8  :: V

    integer :: ifile


    write(outunit,*)  "4. READ TRAJECTORY AND COMPUTE THE DISTRIBUTION OF EACH ANGLE AND EACH BOND"
    rewind(unit = inunit)
    do ifile = 1, N_file
        call read_trj(inunit,cell,coor,numIons,atomtype,io)
        if (io < 0) exit

        if(ifile == 10) exit

        call makeMatrices(cell,coor,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)
        call get_neighbor_list2( ext, dist_matrix, dist_atoms, natoms, N_species, atomtype, &
                                 mat_neighbor, N_neighbor, neighbor_list )


        call update_angle_distr ( ext, neighbor_order_list, atomtype, mat_neighbor, dist_matrix, &
                                  neighbor_list, mat_adf, cont_adf )
        
        call update_dist_distr ( ext, V, neighbor_order_list, atomtype, numIons, mat_neighbor, &
                                  dist_matrix, neighbor_list, rmax, mat_pdf, cont_pdf )

        deallocate( dist_matrix, dist_atoms, neighbor_list )

    enddo
    close(unit = inunit)
end subroutine get_distributions_dist_angles



subroutine get_deviation_each_dist_angle(outunit, natoms, N_species, atomtype, ext, mat_neighbor, &
                                         mat_pdf, mat_adf, sigma_pdf, sigma_adf)
    implicit none
    integer, intent(in)              :: outunit, natoms, N_species, atomtype(:), ext, mat_neighbor(:,:)
    real*8, intent(inout)            :: mat_adf(:,:,:,:), mat_pdf(:,:,:,:)
    real*8, intent(out), allocatable :: sigma_pdf(:,:,:), sigma_adf(:,:,:)

    real*8, allocatable :: mean_pdf(:,:,:), mean_adf(:,:,:)
    integer :: iat, iesp, max_pair, max_neigh, N_neigh, N_pair, jesp

    character(4) :: caux1, caux2

    write(outunit,*)  "5. COMPUTE THE STANDAR DEVIATION OF EACH RADIAL AND ANGLE DISTRIBUTION"

    max_neigh = size(mat_pdf,3)
    max_pair  = size(mat_adf,3)

    allocate( mean_adf(natoms,N_species,max_pair), sigma_adf(natoms,N_species,max_pair), &
              mean_pdf(natoms,N_species,max_neigh), sigma_pdf(natoms,N_species,max_neigh) )

    open(unit = 1, action = "write", status = "replace", file = "deviations.log")
    write(1,"(2i5)") natoms, N_species
    write(1,"(1000i5)") atomType(:)
    do iesp = 1, N_species
        write(1,"(100i5)") mat_neighbor(iesp,:)
    enddo
    write(1,*)
    do iat = 1, natoms
        do iesp = 1, N_species

            N_neigh = mat_neighbor(atomtype(iat),iesp)
            N_pair = N_neigh*(N_neigh-1)/2

            if ( atomtype(iat) == iesp ) cycle

            ! N_pair + ext
            call get_mean_sigma_angle( N_pair, mat_adf(iat,iesp,:,:), mean_adf(iat,iesp,:), sigma_adf(iat,iesp,:), .true. )

            !if (atomType(iat)==2 .and. iesp==1) cycle

            write(1,"(2i5,100f10.3)") iat,iesp, mean_adf(iat,iesp,:N_pair)
            write(1,"(2i5,100f10.3)") iat,iesp, sigma_adf(iat,iesp,:N_pair)
            write(1,*)

            write(caux1,"(i0)") iat
            write(caux2,"(i0)") iesp

            !call write_angle_distr_full( "adf_at"//trim(adjustl(caux1))//"_esp"//trim(adjustl(caux2))//".dat", &
            !                              N_pair, mat_adf(iat,iesp,:,:) )


        enddo
    enddo
    close(unit = 1)

end subroutine get_deviation_each_dist_angle



subroutine write_plot_contribution_total (outunit, N_species, mat_neighbor, ext, rmax, plot_results, &
                                          contribution_pdf, tot_pdf, contribution_adf, tot_adf)
    implicit none
    integer, intent(in) :: outunit, N_species, ext, mat_neighbor(:,:)
    real*8, intent(in)  :: rmax
    real*8, intent(inout)  :: tot_pdf(:,:,:), contribution_pdf(:,:,:,:), &
                           tot_adf(:,:,:), contribution_adf(:,:,:,:)
    logical, intent(in) :: plot_results

    integer :: type1, type2
    character(4) :: caux1, caux2, caux3
    character(19) :: datafile
    character(17) :: datafile2


    real*8, allocatable :: mean(:,:,:), sigma(:,:,:)
    integer :: N_pair, max_pair

     max_pair  = size(contribution_adf,3)

    allocate( mean(N_species,N_species,max_pair), sigma(N_species,N_species,max_pair))


    open(unit = 11, action = "write", status = "replace", file = "instructions.gnuplot")
    do type1 = 1, N_species
        do type2 = 1, N_species
            if (type1 == type2) cycle
            write(caux1,"(i1)") type1
            write(caux2,"(i1)") type2
            

            ! DISTANCES
            !-------------------------------------------------------------------------------
            datafile = "contr_tot_gdr"//trim(adjustl(caux1))//trim(adjustl(caux2))//".dat"
            call write_gdr_tot_contr( datafile, &
                                     mat_neighbor(type1,type2)+ext, contribution_pdf(type1,type2,:,:), &
                                     tot_pdf(type1,type2,:), rmax )

            write(caux3,"(i0)") mat_neighbor(type1,type2)+ext
            write(unit = 11,fmt=*) "gnuplot -e "//'"'//"datafile='"//datafile//"'; outfile='"//&
                              "fig_pdf"//trim(adjustl(caux1))//trim(adjustl(caux2))//".pdf'; nneigh="//&
                              trim(adjustl(caux3))//'"'//" plot_contr_pdf.gnuplot"


            if (plot_results) then
                call system( "gnuplot -e "//'"'//"datafile='"//datafile//"'; outfile='"//&
                              "fig_pdf"//trim(adjustl(caux1))//trim(adjustl(caux2))//".pdf'; nneigh="//&
                              trim(adjustl(caux3))//'"'//" plot_contr_pdf.gnuplot" )
            endif
            !-------------------------------------------------------------------------------


            
            ! ANGLES
            !-------------------------------------------------------------------------------
            datafile = "contr_tot_adf"//trim(adjustl(caux1))//trim(adjustl(caux2))//".dat"
            call write_angle_tot_contr( datafile, ext, &
                                     mat_neighbor(type1,type2), contribution_adf(type1,type2,:,:), &
                                     tot_adf(type1,type2,:) )

            write(caux3,"(i0)") mat_neighbor(type1,type2)*(mat_neighbor(type1,type2)-1)/2!+ext
            write(unit = 11,fmt=*) "gnuplot -e "//'"'//"datafile='"//datafile//"'; outfile='"//&
                              "fig_adf"//trim(adjustl(caux1))//trim(adjustl(caux2))//".pdf'; nneigh="//&
                              trim(adjustl(caux3))//'"'//" plot_contr_adf.gnuplot"


            ! Gnuplot
            if (plot_results) then
                call system( "gnuplot -e "//'"'//"datafile='"//datafile//"'; outfile='"//&
                              "fig_adf"//trim(adjustl(caux1))//trim(adjustl(caux2))//".pdf'; nneigh="//&
                              trim(adjustl(caux3))//'"'//" plot_contr_adf.gnuplot" )
            endif


            N_pair = mat_neighbor(type1,type2)*(mat_neighbor(type1,type2)-1)/2!+ext
            call get_mean_sigma_angle( N_pair, contribution_adf(type1,type2,:N_pair,:), &
                                       mean(type1,type2,:), sigma(type1,type2,:), .false. )

            datafile2 = "sigma_angle"//trim(adjustl(caux1))//trim(adjustl(caux2))//".dat"
            call write_angle_neigh_sigma( datafile2, ext, N_pair, mean(type1,type2,:N_pair), sigma(type1,type2,:N_pair) )

            write(unit = 11,fmt=*) "gnuplot -e "//'"'//"datafile='"//datafile2//"'; outfile='"//&
                              "fig_sigma_angle"//trim(adjustl(caux1))//trim(adjustl(caux2))//".pdf'; nneigh="//&
                              trim(adjustl(caux3))//'"'//" plot_sigma_neigh.gnuplot"

            if (plot_results) then
                call system( "gnuplot -e "//'"'//"datafile='"//datafile2//"'; outfile='"//&
                              "fig_sigma_angle"//trim(adjustl(caux1))//trim(adjustl(caux2))//".pdf'; nneigh="//&
                              trim(adjustl(caux3))//'"'//" plot_sigma_neigh.gnuplot" )
            endif


            ! write(*,"(2i5,100f10.3)") type1, type2, mean(type1,type2,:N_pair)
            ! write(*,"(2i5,100f10.3)") type1, type2, sigma(type1,type2,:N_pair)
            ! write(*,*)
            !-------------------------------------------------------------------------------

        enddo
    enddo
    close (unit = 11)

end subroutine write_plot_contribution_total


subroutine compute_number_constraints(max_sigma_angle)
    implicit none
    real*8, intent(in) :: max_sigma_angle
    integer              :: iat, iesp, jesp, natoms, N_species, max_neigh, &
                            max_pair, N_neigh, N_pair, d1, d2, i
    integer, allocatable :: mat_neighbor(:,:), atomtype(:), numions(:), N_const_angle(:,:)
    real*8, allocatable  :: sigma_angle(:)

    open(unit = 1, action = "read", status = "old", file = "deviations.log")

    read(1,*) natoms, N_species

    allocate(mat_neighbor(N_species, N_species), atomtype(natoms), numions(N_species), &
             N_const_angle(N_species, N_species))
    read(1,*) atomtype(:)
    do iesp = 1, N_species
        read(1,*) mat_neighbor(iesp,:)
    enddo
    read(1,*)

    max_neigh = maxval(mat_neighbor)
    max_pair  = max_neigh*(max_neigh-1)/2
    allocate(sigma_angle(max_pair))

    numions = 0
    N_const_angle = 0
    do iat = 1, natoms
        numions(atomtype(iat)) = numions(atomtype(iat)) + 1
        do iesp = 1, N_species

            if ( atomtype(iat) == iesp ) cycle

            N_neigh = mat_neighbor(atomtype(iat),iesp)
            N_pair = N_neigh*(N_neigh-1)/2

            read(1,*) 
            read(1,*) d1, d2, sigma_angle(:)
            read(1,*) 

            do i = 1, N_pair
                if (sigma_angle(i) < max_sigma_angle) then
                    N_const_angle(atomtype(iat),iesp) = N_const_angle(atomtype(iat),iesp) + 1
                !else
                !    print*, iat, iesp, sigma_angle(i)
                endif
            enddo

        enddo
    enddo

    do iesp = 1, N_species
        do jesp = 1, N_species
            if (iesp == jesp) cycle
            N_neigh = mat_neighbor(iesp,jesp)
            N_pair = N_neigh*(N_neigh-1)/2
            print*, N_const_angle(iesp,jesp), real(numions(jesp),8), real(N_pair,8)
        enddo
    enddo

    print*, N_const_angle!, numions

    close(unit = 1)

end subroutine compute_number_constraints


end module options_main
