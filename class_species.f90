module class_species

    use class_atom
    use class_dos_file

    implicit none

    type species
        integer :: numatoms
        type(atom), allocatable :: atomno( : )
        type(orbital) :: d(5)
        type(orbital) stot, ptot, dtot, ftot
    contains
        procedure :: init => species_init
        procedure :: totals => species_totals
        procedure :: writedos
        procedure :: writedosnospin
    end type species

contains

    subroutine species_init( this, nedos )
        implicit none
        class(species) :: this
        integer, intent(in) :: nedos
        integer :: i

        allocate( this%atomno( this%numatoms ) )

        do i=1, this%numatoms
            call this%atomno(i)%init( nedos )
        enddo

        call this%stot%init( nedos )
        call this%ptot%init( nedos )
        call this%dtot%init( nedos )
        call this%ftot%init( nedos )

    end subroutine species_init

    subroutine species_totals( this )
        implicit none
        class(species) :: this
        integer :: j

        do j=1, this%numatoms
            this%stot%up%bandno   = this%stot%up%bandno   + this%atomno(j)%stot%up%bandno
            this%ptot%up%bandno   = this%ptot%up%bandno   + this%atomno(j)%ptot%up%bandno
            this%dtot%up%bandno   = this%dtot%up%bandno   + this%atomno(j)%dtot%up%bandno
            this%ftot%up%bandno   = this%ftot%up%bandno   + this%atomno(j)%ftot%up%bandno
            this%stot%down%bandno = this%stot%down%bandno + this%atomno(j)%stot%down%bandno
            this%ptot%down%bandno = this%ptot%down%bandno + this%atomno(j)%ptot%down%bandno
            this%dtot%down%bandno = this%dtot%down%bandno + this%atomno(j)%dtot%down%bandno
            this%ftot%down%bandno = this%ftot%down%bandno + this%atomno(j)%ftot%down%bandno
        enddo

    end subroutine species_totals 

    subroutine writedos( this, filename, energy, e_min, e_max )
        implicit none
        class( species ) :: this
        character(len=*), intent(in) :: filename
        real, dimension(:), intent(in) :: energy
        real, intent(in) :: e_min, e_max
        integer :: j

        open(unit=11, file=filename, action="write")
            write(11,*)"# Energy/eV s(up) s(down) p(up) p(down) d(up) d(down) f(up) f(down)"
            do j=1, size(energy)
                if ( ( energy(j) .ge. e_min ) .and. ( energy(j) .le. e_max ) ) then
                    write(11,'(f12.6,8(x,e11.4))') energy(j), this%stot%up%bandno(j), &
                                                              this%stot%down%bandno(j), &
                                                              this%ptot%up%bandno(j), &
                                                              this%ptot%down%bandno(j), &
                                                              this%dtot%up%bandno(j), &
                                                              this%dtot%down%bandno(j), &
                                                              this%ftot%up%bandno(j), &
                                                              this%ftot%down%bandno(j)
                end if 
            enddo
        close(11)

    end subroutine writedos

    subroutine writedosnospin( this, filename, energy, e_min, e_max )
        implicit none
        class( species ) :: this
        character(len=*), intent(in) :: filename
        real, dimension(:), intent(in) :: energy
        real, intent(in) :: e_min, e_max
        integer :: j, i

        open(unit=11, file=filename, action="write")
            write(11,*)"# Energy/eV s p d f"
            do j=1, size(energy)
                if ( ( energy(j) .ge. e_min ) .and. ( energy(j) .le. e_max ) ) then
                    write(11,'(f12.6,4(x,e11.4))') energy(j), this%stot%up%bandno(j), &
                                                              this%ptot%up%bandno(j), &
                                                              this%dtot%up%bandno(j), &
                                                              this%ftot%up%bandno(j)
                end if
            end do
        close(11)

end subroutine writedosnospin

end module class_species
