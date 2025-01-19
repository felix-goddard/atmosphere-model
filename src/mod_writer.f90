module mod_writer

    use mod_kinds, only: ik, rk
    use mod_config, only: config => main_config
    use mod_io, only: write_time_slice, advance_time
    use mod_tiles, only: isd, ied, jsd, jed
    use mod_fields, only: h, ud, vd

    implicit none

    private
    public :: allocate_writer, accumulate_output, write_output

    integer(ik), parameter :: n_output_fields = 3 ! number of outputs

    real(rk), allocatable :: output_field(:,:,:)
    real(rk)              :: accumulation_time
    real(rk), allocatable :: gather_coarray(:,:)[:]
    real(rk), allocatable :: gather(:,:)

    ! these indices into `output_field` determine which variable we are
    ! accumulating or outputting; make sure these are unique
    integer(ik), parameter :: H_IDX = 1
    integer(ik), parameter :: U_IDX = 2
    integer(ik), parameter :: V_IDX = 3
    
contains

    subroutine accumulate_output(dt)
        real(rk), intent(in) :: dt
        integer :: i, j

        output_field(:,:,H_IDX) = output_field(:,:,H_IDX) &
            + dt * h(isd:ied, jsd:jed)

        do concurrent (i=isd:ied, j=jsd:jed)
            output_field(i,j,U_IDX) = output_field(i,j,U_IDX) &
                + dt * .5 * (ud(i,j) + ud(i,j+1))

            output_field(i,j,V_IDX) = output_field(i,j,V_IDX) &
                + dt * .5 * (vd(i,j) + vd(i+1,j))
        end do
        
        accumulation_time = accumulation_time + dt

    end subroutine accumulate_output

    subroutine write_output(previous_time, time)
        real(rk), intent(in) :: previous_time, time
        integer :: i, j, k

        if (this_image() == 1) call advance_time(.5 * (previous_time + time))

        output_field(:,:,:) = output_field(:,:,:) / accumulation_time

        ! Output height field
        call write(output_field(:,:,H_IDX), 'h', time)

        ! Output u wind
        call write(output_field(:,:,U_IDX), 'u', time)

        ! Output v wind
        call write(output_field(:,:,V_IDX), 'v', time)

        output_field(:,:,:) = 0.
        accumulation_time = 0.

    end subroutine write_output

    subroutine write(field, name, time)
        real(rk), intent(in) :: field(isd:ied, jsd:jed)
        character(len=*), intent(in) :: name
        real(rk), intent(in) :: time

        sync all
        gather_coarray(isd:ied, jsd:jed)[1] = field(:,:)
        sync all

        if (this_image() == 1) then
            gather(:,:) = gather_coarray(:,:)
            call write_time_slice(gather, name)
        end if
    end subroutine write

    subroutine allocate_writer()
    
        if (.not. allocated(output_field)) &
            allocate(output_field(isd:ied, jsd:jed, n_output_fields))

        accumulation_time = 0.

        if (.not. allocated(gather_coarray)) &
            allocate(gather_coarray(1:config % nx, 1:config % ny)[*])

        if (this_image() == 1 .and. .not. allocated(gather)) &
            allocate(gather(1:config % nx, 1:config % ny))
        
    end subroutine allocate_writer
    
end module mod_writer