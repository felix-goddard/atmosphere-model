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
      module procedure :: halo_exchange_1arg_2d, halo_exchange_1arg_levels, &
         halo_exchange_4arg_layers, halo_exchange_5arg_layers
   end interface

   integer, parameter :: max_n_args = 5
   real(rk), allocatable :: edge_buffer(:, :, :, :, :) [:]
   real(rk), allocatable :: corner_buffer(:, :, :, :, :) [:]

   ! The edge_buffer has 5 dimensions:
   ! - along-boundary direction (is/js to ie/je)
   ! - across-boundary direction (1:halo_width)
   ! - vertical direction (1:nlev+1)
   ! - neighbour tile (1:4)
   ! - variable (1:max_n_args) so we can multiple variables at once
   ! While we could do each neighbouring tile or each prognostic individually,
   ! doing it this way minimises the number of sync calls.

   ! Similarly, the corner_buffer has:
   ! - along-boundary direction (1:halo_width)
   ! - across-boundary direction (1:halo_width)
   ! - vertical direction (1:nlev+1)
   ! - neighbour tile (1:4)
   ! - variable (1:max_n_args)

contains

   subroutine halo_exchange_1arg_2d(a)
      real(rk), dimension(is:ie, js:je), intent(inout) :: a

      sync images(set(neighbours))

      call copy_to_buffer_2d(a, 1)

      sync images(set(neighbours))

      call copy_from_buffer_2d(a, 1)

   end subroutine halo_exchange_1arg_2d

   subroutine halo_exchange_1arg_levels(a)
      real(rk), dimension(is:ie, js:je, 1:config%nlev + 1), intent(inout) :: a

      sync images(set(neighbours))

      call copy_to_buffer_3d(a, 1, config%nlev + 1)

      sync images(set(neighbours))

      call copy_from_buffer_3d(a, 1, config%nlev + 1)

   end subroutine halo_exchange_1arg_levels

   subroutine halo_exchange_4arg_layers(a, b, c, d)
      real(rk), dimension(is:ie, js:je, 1:config%nlev), intent(inout) :: &
         a, b, c, d

      sync images(set(neighbours))

      call copy_to_buffer_3d(a, 1, config%nlev)
      call copy_to_buffer_3d(b, 2, config%nlev)
      call copy_to_buffer_3d(c, 3, config%nlev)
      call copy_to_buffer_3d(d, 4, config%nlev)

      sync images(set(neighbours))

      call copy_from_buffer_3d(a, 1, config%nlev)
      call copy_from_buffer_3d(b, 2, config%nlev)
      call copy_from_buffer_3d(c, 3, config%nlev)
      call copy_from_buffer_3d(d, 4, config%nlev)

   end subroutine halo_exchange_4arg_layers

   subroutine halo_exchange_5arg_layers(a, b, c, d, e)
      real(rk), dimension(is:ie, js:je, 1:config%nlev), intent(inout) :: &
         a, b, c, d, e

      sync images(set(neighbours))

      call copy_to_buffer_3d(a, 1, config%nlev)
      call copy_to_buffer_3d(b, 2, config%nlev)
      call copy_to_buffer_3d(c, 3, config%nlev)
      call copy_to_buffer_3d(d, 4, config%nlev)
      call copy_to_buffer_3d(e, 5, config%nlev)

      sync images(set(neighbours))

      call copy_from_buffer_3d(a, 1, config%nlev)
      call copy_from_buffer_3d(b, 2, config%nlev)
      call copy_from_buffer_3d(c, 3, config%nlev)
      call copy_from_buffer_3d(d, 4, config%nlev)
      call copy_from_buffer_3d(e, 5, config%nlev)

   end subroutine halo_exchange_5arg_layers

   subroutine copy_to_buffer_3d(q, idx, k_max)
      integer(ik), intent(in) :: k_max
      real(rk), intent(in) :: q(is:ie, js:je, 1:k_max)
      integer(ik), intent(in) :: idx
      integer(ik) :: i

      do i = 0, halo_width - 1
         edge_buffer(jsd:jed, i + 1, 1:k_max, 1, idx) [neighbours(1)] = &
            q(isd + i, jsd:jed, 1:k_max) ! left neighbour
         edge_buffer(jsd:jed, i + 1, 1:k_max, 2, idx) [neighbours(2)] = &
            q(ied - i, jsd:jed, 1:k_max) ! right neighbour
         edge_buffer(isd:ied, i + 1, 1:k_max, 3, idx) [neighbours(3)] = &
            q(isd:ied, jsd + i, 1:k_max) ! bottom neighbour
         edge_buffer(isd:ied, i + 1, 1:k_max, 4, idx) [neighbours(4)] = &
            q(isd:ied, jed - i, 1:k_max) ! top neighbour

         corner_buffer(i + 1, :, 1:k_max, 1, idx) [neighbours(5)] = &
            q(isd + i, jed - halo_width + 1:jed, 1:k_max) ! upper left neighbour
         corner_buffer(i + 1, :, 1:k_max, 2, idx) [neighbours(6)] = &
            q(ied - i, jed - halo_width + 1:jed, 1:k_max) ! upper right neighbour
         corner_buffer(i + 1, :, 1:k_max, 3, idx) [neighbours(7)] = &
            q(isd + i, jsd:jsd + halo_width - 1, 1:k_max) ! lower left neighbour
         corner_buffer(i + 1, :, 1:k_max, 4, idx) [neighbours(8)] = &
            q(ied - i, jsd:jsd + halo_width - 1, 1:k_max) ! lower right neighbour
      end do

   end subroutine copy_to_buffer_3d

   subroutine copy_to_buffer_2d(q, idx)
      real(rk), intent(in) :: q(is:ie, js:je)
      integer(ik), intent(in) :: idx
      integer(ik) :: i

      do i = 0, halo_width - 1
         edge_buffer(jsd:jed, i + 1, 1, 1, idx) [neighbours(1)] = &
            q(isd + i, jsd:jed) ! left neighbour
         edge_buffer(jsd:jed, i + 1, 1, 2, idx) [neighbours(2)] = &
            q(ied - i, jsd:jed) ! right neighbour
         edge_buffer(isd:ied, i + 1, 1, 3, idx) [neighbours(3)] = &
            q(isd:ied, jsd + i) ! bottom neighbour
         edge_buffer(isd:ied, i + 1, 1, 4, idx) [neighbours(4)] = &
            q(isd:ied, jed - i) ! top neighbour

         corner_buffer(i + 1, :, 1, 1, idx) [neighbours(5)] = &
            q(isd + i, jed - halo_width + 1:jed) ! upper left neighbour
         corner_buffer(i + 1, :, 1, 2, idx) [neighbours(6)] = &
            q(ied - i, jed - halo_width + 1:jed) ! upper right neighbour
         corner_buffer(i + 1, :, 1, 3, idx) [neighbours(7)] = &
            q(isd + i, jsd:jsd + halo_width - 1) ! lower left neighbour
         corner_buffer(i + 1, :, 1, 4, idx) [neighbours(8)] = &
            q(ied - i, jsd:jsd + halo_width - 1) ! lower right neighbour
      end do

   end subroutine copy_to_buffer_2d

   subroutine copy_from_buffer_3d(q, idx, k_max)
      integer(ik), intent(in) :: k_max
      real(rk), intent(inout) :: q(is:ie, js:je, 1:k_max)
      integer(ik), intent(in) :: idx
      integer(ik) :: i

      do i = 1, halo_width
         q(isd - i, jsd:jed, 1:k_max) = &
            edge_buffer(jsd:jed, i, 1:k_max, 2, idx) ! left neighbour
         q(ied + i, jsd:jed, 1:k_max) = &
            edge_buffer(jsd:jed, i, 1:k_max, 1, idx) ! right neighbour
         q(isd:ied, jsd - i, 1:k_max) = &
            edge_buffer(isd:ied, i, 1:k_max, 4, idx) ! bottom neighbour
         q(isd:ied, jed + i, 1:k_max) = &
            edge_buffer(isd:ied, i, 1:k_max, 3, idx) ! top neighbour

         q(isd - i, jed + 1:je, 1:k_max) = &
            corner_buffer(i, :, 1:k_max, 4, idx) ! upper left neighbour
         q(ied + i, jed + 1:je, 1:k_max) = &
            corner_buffer(i, :, 1:k_max, 3, idx) ! upper right neighbour
         q(isd - i, js:jsd - 1, 1:k_max) = &
            corner_buffer(i, :, 1:k_max, 2, idx) ! lower left neighbour
         q(ied + i, js:jsd - 1, 1:k_max) = &
            corner_buffer(i, :, 1:k_max, 1, idx) ! lower right neighbour
      end do

   end subroutine copy_from_buffer_3d

      subroutine copy_from_buffer_2d(q, idx)
      real(rk), intent(inout) :: q(is:ie, js:je)
      integer(ik), intent(in) :: idx
      integer(ik) :: i

      do i = 1, halo_width
         q(isd - i, jsd:jed) = &
            edge_buffer(jsd:jed, i, 1, 2, idx) ! left neighbour
         q(ied + i, jsd:jed) = &
            edge_buffer(jsd:jed, i, 1, 1, idx) ! right neighbour
         q(isd:ied, jsd - i) = &
            edge_buffer(isd:ied, i, 1, 4, idx) ! bottom neighbour
         q(isd:ied, jed + i) = &
            edge_buffer(isd:ied, i, 1, 3, idx) ! top neighbour

         q(isd - i, jed + 1:je) = &
            corner_buffer(i, :, 1, 4, idx) ! upper left neighbour
         q(ied + i, jed + 1:je) = &
            corner_buffer(i, :, 1, 3, idx) ! upper right neighbour
         q(isd - i, js:jsd - 1) = &
            corner_buffer(i, :, 1, 2, idx) ! lower left neighbour
         q(ied + i, js:jsd - 1) = &
            corner_buffer(i, :, 1, 1, idx) ! lower right neighbour
      end do

   end subroutine copy_from_buffer_2d

   subroutine allocate_sync_buffers()

      if (.not. allocated(edge_buffer)) &
         allocate (edge_buffer( &
                   1:max(config%nx, config%ny), &
                   halo_width, config%nlev + 1, 4, max_n_args) [*])

      if (.not. allocated(corner_buffer)) &
         allocate (corner_buffer( &
                   halo_width, halo_width, config%nlev + 1, 4, max_n_args) [*])

   end subroutine allocate_sync_buffers

end module mod_sync
