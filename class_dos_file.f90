module class_dos_file

    implicit none

    type dos_file
        character(len=20) :: filename
        integer :: unit
        character(len=1024) :: format_string
    contains
        procedure :: open
        procedure :: close
        procedure :: write
    end type dos_file

contains

    subroutine open(this, filename, number_of_entries, header )
        implicit none
        class(dos_file) :: this
        character(len=*) :: filename, header
        integer :: number_of_entries

        this%filename = filename
        open( file=this%filename, action="write", newunit=this%unit )
        write( this%filename, * ) header
        this%format_string = trim('(f12.6,'//char(48+number_of_entries)//'(x,e11.4))')
    end subroutine open

    subroutine write( this, energy, band_energies )
        implicit none
        class(dos_file) :: this
        real, intent(in) :: energy
        real, intent(in), dimension(:) :: band_energies

        write( this%unit, this%format_string ) energy, band_energies
    end subroutine write

    subroutine close( this )
        implicit none
        class(dos_file) :: this
        close( this%unit )
    end subroutine close

end module class_dos_file