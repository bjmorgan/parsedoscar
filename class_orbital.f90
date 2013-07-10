module class_orbital

    use class_spin
    implicit none

    type orbital
        type(spin) up, down
    contains
        procedure :: init => orbital_init
    end type orbital

contains

    subroutine orbital_init( this, nedos )
        implicit none
        class(orbital) :: this
        integer, intent(in) :: nedos

        call this%up%init( nedos )
        call this%down%init( nedos )
    end subroutine orbital_init

end module class_orbital
