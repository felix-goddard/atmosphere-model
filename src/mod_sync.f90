module mod_sync

    use mod_kinds, only: ik, rk
    use mod_config, only: config => main_config
    use mod_tiles, only: tile_neighbors_2d
    use mod_fields, only: is, ie, js, je, isd, ied, jsd, jed, halo_width, &
                          h, ud, vd

    implicit none

    private
    public :: sync_halos, init_halo_sync

    integer, parameter :: n_prognostics = 3
    integer(ik) :: neighbours(8)
    real(rk), allocatable :: edge_buffer(:,:,:,:)[:]
    real(rk), allocatable :: corner_buffer(:,:,:,:)[:]

    ! The edge_buffer has 4 dimensions:
    ! - along-boundary direction (is/js to ie/je)
    ! - across-boundary direction (1:halo_width)
    ! - neighbour tile (1:4)
    ! - variable (1:n_prognostics) so we can sync all prognostics at once
    ! While we could do each neighbouring tile or each prognostic individually,
    ! doing it this way minimises the number of sync calls.

    ! Similarly, the corner_buffer has:
    ! - along-boundary direction (1:halo_width)
    ! - across-boundary direction (1:halo_width)
    ! - neighbour tile (1:4)
    ! - variable (1:n_prognostics)
    
contains

    subroutine sync_halos()

        call copy_to_buffer(h,  1)
        call copy_to_buffer(ud, 2)
        call copy_to_buffer(vd, 3)

        sync all

        call copy_from_buffer(h,  1)
        call copy_from_buffer(ud, 2)
        call copy_from_buffer(vd, 3)
        
    end subroutine sync_halos


    subroutine copy_to_buffer(q, idx)
        real(rk), intent(in) :: q(is:ie,js:je)
        integer(ik), intent(in) :: idx
        integer :: i

        do i = 0, halo_width-1
            edge_buffer(jsd:jed, i+1, 1, idx)[neighbours(1)] = q(isd+i, jsd:jed) ! left neighbour
            edge_buffer(jsd:jed, i+1, 2, idx)[neighbours(2)] = q(ied-i, jsd:jed) ! right neighbour
            edge_buffer(isd:ied, i+1, 3, idx)[neighbours(3)] = q(isd:ied, jsd+i) ! bottom neighbour
            edge_buffer(isd:ied, i+1, 4, idx)[neighbours(4)] = q(isd:ied, jed-i) ! top neighbour

            corner_buffer(i+1, :, 1, idx)[neighbours(5)] = q(isd+i, jed-halo_width:jed) ! upper left neighbour
            corner_buffer(i+1, :, 2, idx)[neighbours(6)] = q(ied-i, jed-halo_width:jed) ! upper right neighbour
            corner_buffer(i+1, :, 3, idx)[neighbours(7)] = q(isd+i, jsd:jsd+halo_width) ! lower left neighbour
            corner_buffer(i+1, :, 4, idx)[neighbours(8)] = q(ied-i, jsd:jsd+halo_width) ! lower right neighbour
        end do
    
    end subroutine copy_to_buffer


    subroutine copy_from_buffer(q, idx)
        real(rk), intent(inout) :: q(is:ie,js:je)
        integer(ik), intent(in) :: idx
        integer :: i
        
        do i = 0, halo_width-1
            q(is+i, jsd:jed) = edge_buffer(jsd:jed, i+1, 2, idx) ! left neighbour
            q(ie-i, jsd:jed) = edge_buffer(jsd:jed, i+1, 1, idx) ! right neighbour
            q(isd:ied, js+i) = edge_buffer(isd:ied, i+1, 4, idx) ! bottom neighbour
            q(isd:ied, je-i) = edge_buffer(isd:ied, i+1, 3, idx) ! top neighbour

            q(is+i, jed+1:je) = corner_buffer(i+1, :, 4, idx) ! upper left neighbour
            q(ie-i, jed+1:je) = corner_buffer(i+1, :, 3, idx) ! upper right neighbour
            q(is+i, js:jsd-1) = corner_buffer(i+1, :, 2, idx) ! lower left neighbour
            q(ie-i, js:jsd-1) = corner_buffer(i+1, :, 1, idx) ! lower right neighbour
        end do
    
    end subroutine copy_from_buffer


    subroutine init_halo_sync()
    
        if (.not. allocated(edge_buffer)) &
            allocate(edge_buffer(                             &
                1:max(config % nx, config % ny)+2*halo_width, &
                halo_width, 4, n_prognostics)[*]              )
    
        if (.not. allocated(corner_buffer)) &
            allocate(corner_buffer(halo_width, halo_width, 4, n_prognostics)[*])

        neighbours = tile_neighbors_2d(periodic=.true.)
        
    end subroutine init_halo_sync

end module mod_sync