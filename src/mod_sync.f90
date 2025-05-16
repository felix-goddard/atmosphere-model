module mod_sync

   use mod_kinds, only: ik, rk
   use mod_config, only: config => main_config
   use mod_tiles, only: neighbours, &
                        is, ie, js, je, isd, ied, jsd, jed, halo_width
   use mod_util, only: set

   implicit none

   private
   public :: halo_exchange, allocate_sync_buffers

   interface halo_exchange
      module procedure :: halo_exchange_1arg, halo_exchange_4arg, &
         halo_exchange_5arg
   end interface

   integer, parameter :: max_n_args = 5
   real(rk), allocatable :: edge_buffer(:, :, :, :, :) [:]
   real(rk), allocatable :: corner_buffer(:, :, :, :, :) [:]

   ! The edge_buffer has 5 dimensions:
   ! - along-boundary direction (is/js to ie/je)
   ! - across-boundary direction (1:halo_width)
   ! - vertical direction (1:nlev)
   ! - neighbour tile (1:4)
   ! - variable (1:max_n_args) so we can multiple variables at once
   ! While we could do each neighbouring tile or each prognostic individually,
   ! doing it this way minimises the number of sync calls.

   ! Similarly, the corner_buffer has:
   ! - along-boundary direction (1:halo_width)
   ! - across-boundary direction (1:halo_width)
   ! - vertical direction (1:nlev)
   ! - neighbour tile (1:4)
   ! - variable (1:max_n_args)

contains

   subroutine halo_exchange_1arg(a)
      real(rk), dimension(is:ie, js:je, 1:config%nlev), intent(inout) :: a

      sync images(set(neighbours))

      call copy_to_buffer(a, 1)

      sync images(set(neighbours))

      call copy_from_buffer(a, 1)

   end subroutine halo_exchange_1arg

   subroutine halo_exchange_4arg(a, b, c, d)
      real(rk), dimension(is:ie, js:je, 1:config%nlev), intent(inout) :: &
         a, b, c, d

      sync images(set(neighbours))

      call copy_to_buffer(a, 1)
      call copy_to_buffer(b, 2)
      call copy_to_buffer(c, 3)
      call copy_to_buffer(d, 4)

      sync images(set(neighbours))

      call copy_from_buffer(a, 1)
      call copy_from_buffer(b, 2)
      call copy_from_buffer(c, 3)
      call copy_from_buffer(d, 4)

   end subroutine halo_exchange_4arg

   subroutine halo_exchange_5arg(a, b, c, d, e)
      real(rk), dimension(is:ie, js:je, 1:config%nlev), intent(inout) :: &
         a, b, c, d, e

      sync images(set(neighbours))

      call copy_to_buffer(a, 1)
      call copy_to_buffer(b, 2)
      call copy_to_buffer(c, 3)
      call copy_to_buffer(d, 4)
      call copy_to_buffer(e, 5)

      sync images(set(neighbours))

      call copy_from_buffer(a, 1)
      call copy_from_buffer(b, 2)
      call copy_from_buffer(c, 3)
      call copy_from_buffer(d, 4)
      call copy_from_buffer(e, 5)

   end subroutine halo_exchange_5arg

   subroutine copy_to_buffer(q, idx)
      real(rk), intent(in) :: q(is:ie, js:je, 1:config%nlev)
      integer(ik), intent(in) :: idx
      integer(ik) :: i

      do i = 0, halo_width - 1
         edge_buffer(jsd:jed, i + 1, :, 1, idx) [neighbours(1)] = &
            q(isd + i, jsd:jed, :) ! left neighbour
         edge_buffer(jsd:jed, i + 1, :, 2, idx) [neighbours(2)] = &
            q(ied - i, jsd:jed, :) ! right neighbour
         edge_buffer(isd:ied, i + 1, :, 3, idx) [neighbours(3)] = &
            q(isd:ied, jsd + i, :) ! bottom neighbour
         edge_buffer(isd:ied, i + 1, :, 4, idx) [neighbours(4)] = &
            q(isd:ied, jed - i, :) ! top neighbour

         corner_buffer(i + 1, :, :, 1, idx) [neighbours(5)] = &
            q(isd + i, jed - halo_width + 1:jed, :) ! upper left neighbour
         corner_buffer(i + 1, :, :, 2, idx) [neighbours(6)] = &
            q(ied - i, jed - halo_width + 1:jed, :) ! upper right neighbour
         corner_buffer(i + 1, :, :, 3, idx) [neighbours(7)] = &
            q(isd + i, jsd:jsd + halo_width - 1, :) ! lower left neighbour
         corner_buffer(i + 1, :, :, 4, idx) [neighbours(8)] = &
            q(ied - i, jsd:jsd + halo_width - 1, :) ! lower right neighbour
      end do

   end subroutine copy_to_buffer

   subroutine copy_from_buffer(q, idx)
      real(rk), intent(inout) :: q(is:ie, js:je, 1:config%nlev)
      integer(ik), intent(in) :: idx
      integer(ik) :: i

      do i = 1, halo_width
         q(isd - i, jsd:jed, :) = &
            edge_buffer(jsd:jed, i, :, 2, idx) ! left neighbour
         q(ied + i, jsd:jed, :) = &
            edge_buffer(jsd:jed, i, :, 1, idx) ! right neighbour
         q(isd:ied, jsd - i, :) = &
            edge_buffer(isd:ied, i, :, 4, idx) ! bottom neighbour
         q(isd:ied, jed + i, :) = &
            edge_buffer(isd:ied, i, :, 3, idx) ! top neighbour

         q(isd - i, jed + 1:je, :) = &
            corner_buffer(i, :, :, 4, idx) ! upper left neighbour
         q(ied + i, jed + 1:je, :) = &
            corner_buffer(i, :, :, 3, idx) ! upper right neighbour
         q(isd - i, js:jsd - 1, :) = &
            corner_buffer(i, :, :, 2, idx) ! lower left neighbour
         q(ied + i, js:jsd - 1, :) = &
            corner_buffer(i, :, :, 1, idx) ! lower right neighbour
      end do

   end subroutine copy_from_buffer

   subroutine allocate_sync_buffers()

      if (.not. allocated(edge_buffer)) &
         allocate (edge_buffer( &
                   1:max(config%nx, config%ny), &
                   halo_width, config%nlev, 4, max_n_args) [*])

      if (.not. allocated(corner_buffer)) &
         allocate (corner_buffer(halo_width, halo_width, config%nlev, &
                                 4, max_n_args) [*])

   end subroutine allocate_sync_buffers

end module mod_sync
