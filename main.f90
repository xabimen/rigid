program main
use read_write
use pair_dist
use mod_check_min
implicit none

integer                                 :: natoms, io, bins, i
integer                                 :: type1, type2
integer, dimension(:), allocatable      :: numIons, atomtype
real*8, dimension(:), allocatable       :: pdf, angle_distr
real*8, dimension(3,3)                  :: cell
real*8, dimension(:,:), allocatable     :: coor
real*8                                  :: rmax
real*8, parameter                       :: pi = acos(-1.0d0)
integer                                 :: min_pdf

rmax = 8.0
bins = 1000
type1 = 1
type2 = 2

call atom_number("input",natoms)
allocate(coor(natoms,3),atomtype(natoms))


open(unit=123,file="input",status='old',action='read')
open(unit=124,file="kaka_POSCAR",status='replace',action='write')

call read_trj(123,cell,coor,numIons,atomtype,io)
call write_vasp(124,cell,coor)

call pair_distribution(cell,coor,numIons,atomtype,type1,type2,rmax,bins,pdf,angle_distr)

call check_min(pdf, min_pdf)

print*, "# Minimum: ", rmax/bins*min_pdf
do i = 1, bins
    print*, rmax/real(bins)*i, pdf(i)
enddo
!do i = 1, bins
!    if (angle_distr(i) > 1.0D-8) then
!        print*, pi/bins*i, angle_distr(i)
!	endif
!enddo



end program
