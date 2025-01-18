module mod_writer

    use mod_kinds, only: ik, rk
    use mod_config, only: config => main_config
    use mod_io, only: write_time_slice, advance_time
    use mod_fields, only: isd, ied, jsd, jed, &
                          h, ud, vd

    implicit none

    private
    public :: write_output, allocate_writer

    real(rk), allocatable :: output_field(:,:)
    real(rk), allocatable :: gather_coarray(:,:)[:]
    
contains

    subroutine write_output(time)
        real(rk), intent(in) :: time
        integer :: i, j

        if (this_image() == 1) call advance_time(time)

        call write(h(isd:ied, jsd:jed), 'h', time)

        do concurrent (i=isd:ied, j=jsd:jed)
            output_field(i,j) = .5 * (ud(i,j) + ud(i,j+1))
        end do

        call write(output_field, 'u', time)

        do concurrent (i=isd:ied, j=jsd:jed)
            output_field(i,j) = .5 * (vd(i,j) + vd(i+1,j))
        end do

        call write(output_field, 'v', time)

    end subroutine write_output

    subroutine write(field, name, time)
        real(rk), intent(in) :: field(isd:ied, jsd:jed)
        character(len=*), intent(in) :: name
        real(rk), intent(in) :: time

        gather_coarray(isd:ied, jsd:jed)[1] = field(:,:)
        sync all
        if (this_image() == 1) then
            call write_time_slice(gather_coarray(:,:), name)
        end if
    end subroutine write

    subroutine allocate_writer()
    
        if (.not. allocated(output_field))   allocate(output_field(isd:ied, jsd:jed))
        if (.not. allocated(gather_coarray)) allocate(gather_coarray(config % nx, config % ny)[*])
        
    end subroutine allocate_writer
    
end module mod_writer