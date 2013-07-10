program parsedoscar

use class_species

implicit none

integer, parameter :: natommax=400
integer, parameter :: nmaxspec=5

integer i, j, k, l, nspec, unitnum
character(len=50) :: filename
character(len=50) :: system_title
character(len=20) :: dummy(5)
character(len=1)  :: maxl
real :: fermi_energy_shift
real, allocatable :: energy(:)
logical :: spinpol, noncoll, dproj
integer :: nedos

type sumdos
    type(spin) up, down
end type sumdos

type(sumdos) :: summeddos
type(species), allocatable :: spec(:)

!**************************************************

! Read from the command line the filename, number of species and the number of
! atoms of each species

write(6,*) "input DOS filename:"
read(5,*) filename

write(6,*) "input number of species:"
read(5,*) nspec
allocate( spec( nspec ) )

if (nspec > nmaxspec) then
    write(6,*) "nspec > nmaxspec"
    stop
endif

do i=1,nspec
    write(6,"(A,I2,A)") "number of atoms of species ",i,":"
    read(5,*) spec(i)%numatoms
enddo

write(6,*) "non-collinear calculation? (T/F):"
read(5,*) noncoll
write(6,*) "noncoll: ", noncoll

if ( noncoll ) then
    spinpol = .false.
else
    write(6,*) "spin polarised calculation? (T/F):"
    read(5,*) spinpol
    write(6,*) "spinpol: ", spinpol
endif

write(6,*) "max angular momentum: (d/f):"
read(5,*) maxl

if ((maxl /= "d").and.(maxl /= "f")) then
    write(6,*)"maxl = (d|f)"
    stop
endif

write(6,*) "output separate d orbitals? (T/F):"
read(5,*) dproj

! read header of DOS file
write(6,*)"opening DOS file"
open(unit=10, file=filename, status="old", action="read")

write(6,*)"reading DOS file"

do i=1, 4
    read(10,*) 
enddo

read(10,*) system_title

read(10,*) dummy(1), dummy(1), nedos, fermi_energy_shift, dummy(1)

do i=1, nspec
    call spec(i)%init( nedos )
end do

allocate( energy( nedos ) )

! read in summed DOS

allocate( summeddos%up%bandno( nedos ) )
allocate( summeddos%down%bandno( nedos ) )

if (spinpol) then
    do i=1, nedos
        read(10,*) energy(i), summeddos%up%bandno(i), summeddos%down%bandno(i), dummy(1), dummy(1)
    enddo
else
    do i=1, nedos
        read(10,*) energy(i), summeddos%up%bandno(i), dummy(1)
    enddo
endif

write(6,*) energy(nedos), summeddos%up%bandno(nedos)

write(6,*)" reading atom and angular momenta decomposed dos values"
! read in atom and angular momenta decomposed dos values
if (spinpol) then
    if (maxl=="f") then
        write(6,*)"spin-polarised / s,p,d,f" 
        do i=1,nspec
            do j=1,spec(i)%numatoms
                call spec(i)%atomno(j)%read_spin_pol_f
            enddo
        enddo
    else
        write(6,*)"spin-polarised / s,p,d" 
        do i=1,nspec
            do j=1,spec(i)%numatoms
                call spec(i)%atomno(j)%read_spin_pol_d
            enddo
        enddo
    endif
else if (noncoll) then
    if (maxl == "f") then
        write(6,*) "non-collinear / s,p,d,f"
        do i=1, nspec
            do j=1, spec(i)%numatoms
                call spec(i)%atomno(j)%read_noncoll_f
            end do
        end do
        ! up to l = f
    else
        write(6,*) "non-collinear / s,p,d"
        do i=1, nspec
            do j=1, spec(i)%numatoms
                call spec(i)%atomno(j)%read_noncoll_d
            end do
        end do
    endif
else
    if (maxl=="f") then
        write(6,*) "spin-unpolarised / s,p,d,f"
        do i=1, nspec
            do j=1,spec(i)%numatoms
                call spec(i)%atomno(j)%read_nonspin_pol_f
            enddo
        enddo
    else
        write(6,*)"spin-unpolarised / s,p,d"
        do i=1,nspec
            do j=1,spec(i)%numatoms
                call spec(i)%atomno(j)%read_nonspin_pol_d
            enddo
        enddo
    endif
endif


if( spinpol )then
! multiply all down spin values by -1
    do i=1, nspec
        do j=1, spec(i)%numatoms
            call spec(i)%atomno(j)%invert_down_bands
        enddo
    enddo
endif

! average over angular momenta for each atom
do i=1, nspec
    do j=1, spec(i)%numatoms
        call spec(i)%atomno(i)%totals
    enddo
enddo
        
! average over angular momenta for each species

do i=1, nspec
    call spec(i)%totals
enddo

! shift energy values to give correct Fermi level
energy = energy - fermi_energy_shift

if (dproj) then ! output d orbital projections
    do i=1, nspec
        do j=1, spec(i)%numatoms
            call spec(i)%atomno(j)%write_proj_dos_d( "fortd."//char(48+i)//"."//char(48+j), energy(1:nedos), energy(1), energy(nedos) )
        end do
    end do
end if

! write angular momenta and species decomposed DOS

if (spinpol) then
    do i=1, nspec
        call spec(i)%writedos( "fort."//char(48+i), energy(1:nedos), energy(1), energy(nedos) )
    enddo
else
    do i=1, nspec
        call spec(i)%writedosnospin("fort."//char(48+i), energy(1:nedos), energy(1), energy(nedos))
    enddo
endif

end program parsedoscar
