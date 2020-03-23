module read_write


contains
!*********************************************
subroutine atom_number(filename,natoms)
implicit none
character(len=*), intent(in)          :: filename
integer, intent(out)                     :: natoms
integer                                 :: i

open(unit=123,file=filename,status='old',action='read')

do i =1, 3
    read(unit=123,fmt=*)
enddo

read(unit=123,fmt=*) natoms
close(unit=123)
end subroutine
!*********************************************
!*********************************************
subroutine read_trj(in,cell,coor,numIons,atomtype,io)
implicit none
integer, intent(in)                     :: in
real*8, dimension(:,:), intent(out)     :: cell
real*8, dimension(:,:), intent(out)     :: coor
integer, intent(out)                    :: io
integer, dimension(size(coor,1))        :: atomtype
integer, dimension(:), allocatable      :: numIons
real*8, dimension(3)                    :: r
real*8, dimension(3,3)                  :: box
integer                                 :: i, n

read(in, fmt=*, iostat=io)

if (io < 0) return
do i = 1, 4
    read(in,fmt=*)
enddo

do i = 1, 3
    read(in,fmt=*) box(i,:)
enddo

call box2cell(box,cell)
read(in,fmt=*)

do i = 1, size(coor,1)
    read(in,fmt=*) n, atomtype(n), r
    coor(n,:) = r
enddo

allocate(numIons(maxval(atomType)))

numIons=0
do i = 1 , size(coor,1)
    numIons(atomType(i)) = numIons(atomType(i)) + 1
enddo

!call sort(coor,atom_type,atomType)


end subroutine
!*********************************************
!*********************************************
subroutine sort(coor,atom_type,atomType)
implicit NONE
real*8, dimension(:,:), intent(inout)   :: coor
real*8, dimension(size(coor,1),3)       :: coor2
integer, dimension(:), intent(in)       :: atom_type
integer, dimension(:), allocatable      :: atomType
INTEGER                                 :: i, j, k, l


allocate(atomType(maxval(atom_type)))

coor2(:,:)=coor(:,:)

k = 1
do i = 1, maxval(atom_type)
    l = 0
    do j = 1, size(coor,1)
        if (atom_type(j)==i) THEN
            coor(k,:) = coor2(j,:)
            l = l + 1
            k = k + 1
        endif
    enddo
    atomType(i) = l
enddo


end subroutine
!*********************************************
!*********************************************
subroutine box2cell(box,cell)
real*8, dimension(:,:), intent(in)      :: box
real*8, dimension(:,:), intent(out)     :: cell
real*8                                  :: xlo
real*8                                  :: xhi
real*8                                  :: ylo
real*8                                  :: yhi
real*8                                  :: zlo
real*8                                  :: zhi
real*8                                  :: xy
real*8                                  :: xz
real*8                                  :: yz

xy = box(1,3)
xz = box(2,3)
yz = box(3,3)

xlo = box(1,1) - minval((/0.0d0,xy,xz,xy+xz/))
xhi = box(1,2) - maxval((/0.0d0,xy,xz,xy+xz/))
ylo = box(2,1) - minval((/0.0d0,yz/))
yhi = box(2,2) - maxval((/0.0d0,yz/))
zlo = box(3,1)
zhi = box(3,2)

cell = 0.0

cell(1,1) = xhi - xlo
cell(2,2) = yhi - ylo
cell(3,3) = zhi - zlo
cell(2,1) = xy
cell(3,1) = xz
cell(3,2) = yz

end subroutine
!*********************************************
!*********************************************
subroutine write_vasp(fu,cell,coor)
integer, intent(in)                     :: fu
real*8, dimension(:,:), intent(in)      :: cell
real*8, dimension(:,:), intent(in)      :: coor
integer                                 :: i

write(fu,fmt='(a10,f20.10)') "Energy = "
write(fu,fmt='(f5.3)') 1.0

do i = 1, 3
    write(fu,fmt='(3f15.8)') cell(i,:)
enddo

write(fu,fmt='(a)') "Ca Si O H"
write(fu,fmt='(a)') "4 4 24 24"
write(fu,fmt='(a)') "Cartesian"

do i = 1, 56
    write(fu,fmt='(3f15.8)') matmul(coor(i,:),cell)
enddo

end subroutine
!*********************************************
!*********************************************
subroutine read_xsf(filename, species, natoms, cell, coor, numions, atomtype, io)
implicit none
character(len=*), intent(in)      :: filename
character(len=2), intent(in)      :: species(:)
real*8, intent(out)               :: cell(:,:)
real*8, allocatable, intent(out)  :: coor(:,:)
integer, allocatable, intent(out) :: numions(:), atomtype(:)
integer, intent(out)              :: natoms, io
integer                           :: nspecies, i, ind
real*8                            :: cartesian(3), fractional(3), inv_cell(3,3)
character(len=2)                  :: at_type

open(unit = 1, action = "read", status = "old", file = filename)

do i = 1, 5
    read(1,*)
enddo

do i = 1, 3
    read(1,*) cell(i,:)
enddo
call M33INV(cell, inv_cell)

read(1,*)
read(1,*) natoms

nspecies = size(species)
allocate( coor(natoms,3), atomtype(natoms), numions(nspecies) )
numions = 0

do i = 1, natoms
    read(1,*) at_type, cartesian(:)
    fractional = matmul(cartesian, inv_cell)
    coor(i,:)  = fractional(:)
    do ind = 1, nspecies
        if ( trim(adjustl(at_type)) == trim(adjustl(species(ind))) ) then
            atomtype(i) = ind
            numions(ind) = numions(ind) + 1
            exit
        endif
    enddo
enddo


close(unit = 1)

end subroutine read_xsf
!*********************************************
!*********************************************

SUBROUTINE M33INV (A, AINV)

IMPLICIT NONE

real*8, DIMENSION(3,3), INTENT(IN)  :: A
real*8, DIMENSION(3,3), INTENT(OUT) :: AINV

real*8, PARAMETER :: EPS = 1.0D-10
real*8 :: DET
real*8, DIMENSION(3,3) :: COFACTOR


DET =   A(1,1)*A(2,2)*A(3,3)  &
    - A(1,1)*A(2,3)*A(3,2)  &
    - A(1,2)*A(2,1)*A(3,3)  &
    + A(1,2)*A(2,3)*A(3,1)  &
    + A(1,3)*A(2,1)*A(3,2)  &
    - A(1,3)*A(2,2)*A(3,1)

COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

AINV = TRANSPOSE(COFACTOR) / DET

END SUBROUTINE M33INV
!*********************************************
end module
