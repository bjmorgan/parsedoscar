module class_atom

    use class_orbital
    use class_dos_file
    implicit none

    type atom
        character(len=2) :: name
        type(orbital) :: s(1), p(3), d(5), f(7)
        type(orbital) :: stot, ptot, dtot, ftot
        integer :: nedos
    contains
        procedure :: read_noncoll_d
        procedure :: read_noncoll_f
        procedure :: read_spin_pol_d
        procedure :: read_spin_pol_f
        procedure :: read_nonspin_pol_d
        procedure :: read_nonspin_pol_f
        procedure :: invert_down_bands
        procedure :: init => atom_init
        procedure :: totals => atom_totals
        ! procedure :: write_proj_dos
        procedure :: write_proj_dos_p
        procedure :: write_proj_dos_d
    end type atom

contains

    subroutine atom_init( this, nedos )
        implicit none
        class( atom ) :: this
        integer, intent(in) :: nedos
        integer k

        this%nedos = nedos

        k=1
        call this%s(k)%init( this%nedos )
        call this%stot%init( this%nedos )        

        do k=1, 3
            call this%p(k)%init( this%nedos )            
        enddo
        call this%ptot%init( this%nedos )       

        do k=1, 5
            call this%d(k)%init( this%nedos )            
        enddo
        call this%dtot%init( this%nedos )        

        do k=1, 7
            call this%f(k)%init( this%nedos )
            
        enddo
        call this%ftot%init( this%nedos )        

    end subroutine atom_init

    subroutine atom_totals( this )

        implicit none
        class( atom ) :: this
        integer :: l

        l=1        
        this%stot%up%bandno = this%stot%up%bandno + this%s(l)%up%bandno
        this%stot%down%bandno = this%stot%down%bandno + this%s(l)%down%bandno

        do l=1,3
            this%ptot%up%bandno = this%ptot%up%bandno + this%p(l)%up%bandno
            this%ptot%down%bandno = this%ptot%down%bandno + this%p(l)%down%bandno
        enddo

        do l=1,5
            this%dtot%up%bandno = this%dtot%up%bandno + this%d(l)%up%bandno
            this%dtot%down%bandno = this%dtot%down%bandno + this%d(l)%down%bandno
        enddo
    
        do l=1,7
            this%ftot%up%bandno = this%ftot%up%bandno + this%f(l)%up%bandno
            this%ftot%down%bandno = this%ftot%down%bandno + this%f(l)%down%bandno
        enddo

    end subroutine atom_totals

    subroutine read_noncoll_d( this )
        implicit none
        class( atom ) :: this
        integer :: k
        double precision :: dummy(3)

        read(10,*) dummy
        do k=1, this%nedos
            read(10,*) dummy(1), &
                this%s(1)%up%bandno( k ), dummy(1:3), &

                this%p(1)%up%bandno( k ), dummy(1:3), &
                this%p(2)%up%bandno( k ), dummy(1:3), &
                this%p(3)%up%bandno( k ), dummy(1:3), &

                this%d(1)%up%bandno( k ), dummy(1:3), &
                this%d(2)%up%bandno( k ), dummy(1:3), &
                this%d(3)%up%bandno( k ), dummy(1:3), &
                this%d(4)%up%bandno( k ), dummy(1:3), &
                this%d(5)%up%bandno( k ), dummy(1:3) 
        end do

    end subroutine read_noncoll_d

    subroutine read_noncoll_f( this )
        implicit none
        class( atom ) :: this
        integer :: k
        double precision :: dummy(3)

        read(10,*) dummy
        do k=1, this%nedos
            read(10,*) dummy(1), &
                this%s(1)%up%bandno( k ), dummy(1:3), &

                this%p(1)%up%bandno( k ), dummy(1:3), &
                this%p(2)%up%bandno( k ), dummy(1:3), &
                this%p(3)%up%bandno( k ), dummy(1:3), &

                this%d(1)%up%bandno( k ), dummy(1:3), &
                this%d(2)%up%bandno( k ), dummy(1:3), &
                this%d(3)%up%bandno( k ), dummy(1:3), &
                this%d(4)%up%bandno( k ), dummy(1:3), &
                this%d(5)%up%bandno( k ), dummy(1:3), & 

                this%f(1)%up%bandno( k ), dummy(1:3), &
                this%f(2)%up%bandno( k ), dummy(1:3), &
                this%f(3)%up%bandno( k ), dummy(1:3), &
                this%f(4)%up%bandno( k ), dummy(1:3), &
                this%f(5)%up%bandno( k ), dummy(1:3), &
                this%f(6)%up%bandno( k ), dummy(1:3), &
                this%f(7)%up%bandno( k ), dummy(1:3) 
        end do
        
    end subroutine read_noncoll_f

    subroutine read_spin_pol_d( this )
        implicit none
        class( atom ) :: this
        integer :: k
        double precision :: dummy

        read(10,*) dummy

        do k=1, this%nedos
            read(10,*) dummy, &
                this%s(1)%up%bandno(k), this%s(1)%down%bandno(k), &

                this%p(1)%up%bandno(k), this%p(1)%down%bandno(k), &
                this%p(2)%up%bandno(k), this%p(2)%down%bandno(k), &
                this%p(3)%up%bandno(k), this%p(3)%down%bandno(k), &

                this%d(1)%up%bandno(k), this%d(1)%down%bandno(k), &
                this%d(2)%up%bandno(k), this%d(2)%down%bandno(k), &
                this%d(3)%up%bandno(k), this%d(3)%down%bandno(k), &
                this%d(4)%up%bandno(k), this%d(4)%down%bandno(k), &
                this%d(5)%up%bandno(k), this%d(5)%down%bandno(k)
        end do

    end subroutine read_spin_pol_d

    subroutine read_nonspin_pol_d( this )
        implicit none
        class( atom ) :: this
        integer :: k
        double precision :: dummy

        read(10,*) dummy

        do k=1, this%nedos
            read(10,*) dummy, &
                this%s(1)%up%bandno(k), &

                this%p(1)%up%bandno(k), &
                this%p(2)%up%bandno(k), &
                this%p(3)%up%bandno(k), &

                this%d(1)%up%bandno(k), &
                this%d(2)%up%bandno(k), &
                this%d(3)%up%bandno(k), &
                this%d(4)%up%bandno(k), &
                this%d(5)%up%bandno(k)
        end do

    end subroutine read_nonspin_pol_d

    subroutine read_spin_pol_f( this )
        implicit none
        class( atom ) :: this
        integer :: k
        double precision :: dummy

        read(10,*) dummy

        do k=1, this%nedos
            read(10,*) dummy, &
                this%s(1)%up%bandno(k), this%s(1)%down%bandno(k), &

                this%p(1)%up%bandno(k), this%p(1)%down%bandno(k), &
                this%p(2)%up%bandno(k), this%p(2)%down%bandno(k), &
                this%p(3)%up%bandno(k), this%p(3)%down%bandno(k), &

                this%d(1)%up%bandno(k), this%d(1)%down%bandno(k), &
                this%d(2)%up%bandno(k), this%d(2)%down%bandno(k), &
                this%d(3)%up%bandno(k), this%d(3)%down%bandno(k), &
                this%d(4)%up%bandno(k), this%d(4)%down%bandno(k), &
                this%d(5)%up%bandno(k), this%d(5)%down%bandno(k), &

                this%f(1)%up%bandno(k), this%f(1)%down%bandno(k), &
                this%f(2)%up%bandno(k), this%f(2)%down%bandno(k), &
                this%f(3)%up%bandno(k), this%f(3)%down%bandno(k), &
                this%f(4)%up%bandno(k), this%f(4)%down%bandno(k), &
                this%f(5)%up%bandno(k), this%f(5)%down%bandno(k), &
                this%f(6)%up%bandno(k), this%f(6)%down%bandno(k), &
                this%f(7)%up%bandno(k), this%f(7)%down%bandno(k)
        end do

    end subroutine read_spin_pol_f

    subroutine read_nonspin_pol_f( this )
        implicit none
        class( atom ) :: this
        integer :: k
        double precision :: dummy

        read(10,*) dummy

        do k=1, this%nedos
            read(10,*) dummy, &
                this%s(1)%up%bandno(k), &

                this%p(1)%up%bandno(k), &
                this%p(2)%up%bandno(k), &
                this%p(3)%up%bandno(k), &

                this%d(1)%up%bandno(k), &
                this%d(2)%up%bandno(k), &
                this%d(3)%up%bandno(k), &
                this%d(4)%up%bandno(k), &
                this%d(5)%up%bandno(k), &

                this%f(1)%up%bandno(k), &
                this%f(2)%up%bandno(k), &
                this%f(3)%up%bandno(k), &
                this%f(4)%up%bandno(k), &
                this%f(5)%up%bandno(k), &
                this%f(6)%up%bandno(k), &
                this%f(7)%up%bandno(k)
        end do

    end subroutine read_nonspin_pol_f

    subroutine invert_down_bands( this )
        implicit none
        class( atom ) :: this
        integer :: i

        this%s(1)%down%bandno = -this%s(1)%down%bandno

        forall( i=1:3 )
            this%p(i)%down%bandno = -this%p(i)%down%bandno
        end forall

        forall( i=1:5 )
            this%d(i)%down%bandno = -this%d(i)%down%bandno
        end forall

        forall( i=1:7 )
            this%f(i)%down%bandno = -this%f(i)%down%bandno
        end forall

    end subroutine invert_down_bands

    subroutine write_proj_dos_d( this, filename, energy, e_min, e_max )
        implicit none
        class( atom ) :: this
        character(len=*), intent(in) :: filename
        real, dimension(:), intent(in) :: energy
        real, intent(in) :: e_min, e_max
        integer :: j, i
        real :: output(5)
        type( dos_file ) :: dos
        character(len=20) :: header

        header = "# Energy/eV d_xy d_yz d_z2-r2 d_xz d_x2-y2"

        call dos%open( filename, 5, header ) 
            do j=1, size(energy)
                if ( ( energy(j) .ge. e_min ) .and. ( energy(j) .le. e_max ) ) then
                        forall (i=1:5)
                            output(i) = this%d(i)%up%bandno(j)
                        end forall
                    call dos%write( energy(j), output ) 
                end if
            end do
        call dos%close
    end subroutine write_proj_dos_d

    subroutine write_proj_dos_p( this, filename, energy, e_min, e_max )
        implicit none
        class( atom ) :: this
        character(len=*), intent(in) :: filename
        real, dimension(:), intent(in) :: energy
        real, intent(in) :: e_min, e_max
        integer :: j, i
        real :: output(3)
        type( dos_file ) :: dos
        character(len=20) :: header

        header = "# Energy/eV p_x p_y p_z"

        call dos%open( filename, 3, header ) 
            do j=1, size(energy)
                if ( ( energy(j) .ge. e_min ) .and. ( energy(j) .le. e_max ) ) then
                        forall (i=1:3)
                            output(i) = this%p(i)%up%bandno(j)
                        end forall
                    call dos%write( energy(j), output ) 
                end if
            end do
        call dos%close
    end subroutine write_proj_dos_p

end module class_atom
