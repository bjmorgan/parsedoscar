module class_spin

    implicit none

    type spin
        real, allocatable :: bandno(:)
    contains
        procedure :: init => spin_init
    end type spin

contains

    subroutine spin_init( this, nedos )
        implicit none
        class(spin) :: this
        integer, intent(in) :: nedos

        allocate( this%bandno( nedos ) )
        this%bandno = 0.0d0
    end subroutine spin_init

end module class_spin
