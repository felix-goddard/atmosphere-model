module mod_writer

    use mod_kinds, only: ik, rk
    use mod_config, only: config => main_config
    use mod_io, only: write_time_slice, advance_time
    use mod_tiles, only: isd, ied, jsd, jed
    use mod_fields, only: h, ud, vd

    implicit none

    private
    public :: allocate_writer, accumulate_output, clear_accumulator, write_output

    integer(ik), parameter :: n_output_fields = 3 ! number of outputs

    real(rk), allocatable :: output_field(:,:,:)
    real(rk), allocatable :: compensation(:,:,:)
    real(rk)              :: accumulation_time
    real(rk), allocatable :: gather_coarray(:,:)[:]
    real(rk), allocatable :: gather(:,:)

    ! these indices into `output_field` determine which variable we are
    ! accumulating or outputting; make sure these are unique
    integer(ik), parameter :: H_IDX = 1
    integer(ik), parameter :: U_IDX = 2
    integer(ik), parameter :: V_IDX = 3
    
contains

    subroutine compensated_sum(val, i, j, idx)
        integer(ik), intent(in) :: i, j, idx
        real(rk), intent(in) :: val
        real(rk) :: y, sum, c

        y = val - compensation(i,j,idx)
        sum = output_field(i,j,idx) + y
        c = (sum - output_field(i,j,idx)) - y
        output_field(i,j,idx) = sum

    end subroutine compensated_sum

    subroutine accumulate_output(dt)
        real(rk), intent(in) :: dt
        integer(ik) :: i, j

        do i = isd, ied
            do j = jsd, jed
                ! height
                call compensated_sum(dt * h(i,j), i, j, H_IDX)

                ! u wind
                call compensated_sum(dt * .5 * (ud(i,j) + ud(i,j+1)), i, j, U_IDX)

                ! v wind
                call compensated_sum(dt * .5 * (vd(i,j) + vd(i+1,j)), i, j, V_IDX)
            end do
        end do
        
        accumulation_time = accumulation_time + dt

    end subroutine accumulate_output

    subroutine clear_accumulator()
        
        output_field(:,:,:) = 0.
        compensation(:,:,:) = 0.
        accumulation_time = 0.
        
    end subroutine clear_accumulator

    subroutine write_output(previous_time, time)
        real(rk), intent(in) :: previous_time, time

        if (this_image() == 1) call advance_time(.5 * (previous_time + time))

        output_field(:,:,:) = output_field(:,:,:) / accumulation_time

        ! Output height field
        call write(output_field(:,:,H_IDX), 'h')

        ! Output u wind
        call write(output_field(:,:,U_IDX), 'u')

        ! Output v wind
        call write(output_field(:,:,V_IDX), 'v')

        call clear_accumulator()

    end subroutine write_output

    subroutine write(field, name)
        real(rk), intent(in) :: field(isd:ied, jsd:jed)
        character(len=*), intent(in) :: name

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

        if (.not. allocated(compensation)) &
            allocate(compensation(isd:ied, jsd:jed, n_output_fields))

        accumulation_time = 0.

        if (.not. allocated(gather_coarray)) &
            allocate(gather_coarray(1:config % nx, 1:config % ny)[*])

        if (this_image() == 1 .and. .not. allocated(gather)) &
            allocate(gather(1:config % nx, 1:config % ny))
        
    end subroutine allocate_writer
    
end module mod_writer