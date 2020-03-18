program main
use read_write
use pair_dist
implicit none

integer                                 :: natoms, io, bins, i
integer                                 :: type1, type2
integer, dimension(:), allocatable      :: numIons, atomtype
real*8, dimension(:), allocatable       :: pdf
real*8, dimension(3,3)                  :: cell
real*8, dimension(:,:), allocatable     :: coor
real*8                                  :: rmax

rmax =10.0
bins =1000
type1 = 1
type2 = 2

call atom_number("atoms.lammpstrj",natoms)
allocate(coor(natoms,3),atomtype(natoms))


open(unit=123,file="atoms.lammpstrj",status='old',action='read')
open(unit=124,file="kaka_POSCAR",status='replace',action='write')

call read_trj(123,cell,coor,numIons,atomtype,io)
call write_vasp(124,cell,coor)

call pair_distribution(cell,coor,numIons,atomtype,type1,type2,rmax,bins,pdf)

do i = 1, bins
    print*, rmax/real(bins)*i, pdf(i)
enddo


end program
