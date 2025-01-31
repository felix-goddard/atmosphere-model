module mod_tiles

   use mod_kinds, only: ik, rk
   use mod_config, only: config => main_config

   implicit none

   private
   public :: init_tiles, neighbours, &
             is, ie, js, je, isd, ied, jsd, jed, halo_width

   integer(ik), parameter :: halo_width = 9
   integer(ik) :: is, ie, js, je ! data bounds (including halo)
   integer(ik) :: isd, ied, jsd, jed ! domain bounds (not including halo)
   integer(ik) :: neighbours(8)

   interface tile_indices
      module procedure :: tile_indices_1d, tile_indices_2d
   end interface tile_indices

contains

   pure function divisors(n)
      ! Returns all integer divisors of n.

      integer(ik), intent(in) :: n
      integer(ik), allocatable :: divisors(:)
      integer(ik) :: i

      divisors = [integer(ik) ::]
      do i = 1, n
         if (mod(n, i) == 0) divisors = [divisors, i]
      end do

   end function divisors

   pure function num_tiles(n)
      ! Returns the optimal number of tiles in 2 dimensions
      ! given total number of tiles n.

      integer(ik), intent(in) :: n
      integer(ik) :: num_tiles(2)
      integer(ik), allocatable :: divs(:)
      integer(ik) :: dim1, dim2, i

      divs = divisors(n)

      num_tiles = [1, n]
      do i = 1, size(divs)
         dim1 = n/divs(i)
         dim2 = divs(i)
         if (dim1 + dim2 < num_tiles(1) + num_tiles(2)) then
            num_tiles = [dim1, dim2]
         end if
      end do

   end function num_tiles

   pure function tile_indices_1d(dims, i, n) result(indices)
      ! Given input global array size, return start and end index
      ! of a parallel 1-d tile that correspond to this image.

      integer(ik), intent(in) :: dims, i, n
      integer(ik) :: indices(2)
      integer(ik) :: offset, tile_size

      tile_size = dims/n

      ! start and end indices assuming equal tile sizes
      indices(1) = (i - 1)*tile_size + 1
      indices(2) = indices(1) + tile_size - 1

      ! if we have any remainder, distribute it to the tiles at the end
      offset = n - mod(dims, n)
      if (i > offset) then
         indices(1) = indices(1) + i - offset - 1
         indices(2) = indices(2) + i - offset
      end if

   end function tile_indices_1d

   pure function tile_indices_2d(dims) result(indices)
      ! Given an input x- and y- dimensions of the total computational domain [im, jm].
      ! returns an array of start and end indices in x- and y-, [is, ie, js, je].

      integer(ik), intent(in) :: dims(2)
      integer(ik) :: indices(4)
      integer(ik) :: tiles(2), tiles_ij(2)

      tiles = num_tiles(num_images())
      tiles_ij = tile_n2ij(this_image())
      indices(1:2) = tile_indices_1d(dims(1), tiles_ij(1), tiles(1))
      indices(3:4) = tile_indices_1d(dims(2), tiles_ij(2), tiles(2))

   end function tile_indices_2d

   pure function tile_n2ij(n) result(ij)
      ! Given tile index in a 1-d layout, returns the
      ! corresponding tile indices in a 2-d layout.
      !
      !    +---+---+---+
      !  2 | 4 | 5 | 6 |
      !    +---+---+---+
      !  1 | 1 | 2 | 3 |
      !  j +---+---+---+
      !    i 1   2   3
      !
      ! Examples:
      !   * tile_n2ij(2) = [2, 1]
      !   * tile_n2ij(4) = [1, 2]
      !   * tile_n2ij(6) = [3, 2]
      !

      integer(ik), intent(in) :: n
      integer(ik) :: ij(2), i, j, tiles(2)

      if (n == 0) then
         ij = 0
      else
         tiles = num_tiles(num_images())
         j = (n - 1)/tiles(1) + 1
         i = n - (j - 1)*tiles(1)
         ij = [i, j]
      end if

   end function tile_n2ij

   pure function tile_ij2n(ij) result(n)
      ! Given tile indices in a 2-d layout, returns the
      ! corresponding tile index in a 1-d layout:
      !
      !    +---+---+---+
      !  2 | 4 | 5 | 6 |
      !    +---+---+---+
      !  1 | 1 | 2 | 3 |
      !  j +---+---+---+
      !    i 1   2   3
      !
      ! Examples:
      !   * tile_ij2n([2, 1]) = 2
      !   * tile_ij2n([1, 2]) = 4
      !   * tile_ij2n([3, 2]) = 6
      !

      integer(ik), intent(in) :: ij(2)
      integer(ik) :: n, tiles(2)

      if (any(ij == 0)) then
         n = 0
      else
         tiles = num_tiles(num_images())
         n = (ij(2) - 1)*tiles(1) + ij(1)
      end if

   end function tile_ij2n

   pure function tile_neighbors_2d(periodic) result(neighbors)
      ! Returns the neighbor image indices given.

      logical, intent(in) :: periodic
      integer(ik) :: neighbors(8)
      integer(ik) :: tiles(2), tiles_ij(2), itile, jtile
      integer(ik) :: ij_L(2), ij_R(2), ij_D(2), ij_U(2), &
                     ij_UL(2), ij_UR(2), ij_DL(2), ij_DR(2)

      tiles = num_tiles(num_images())
      tiles_ij = tile_n2ij(this_image())
      itile = tiles_ij(1)
      jtile = tiles_ij(2)

      ! i, j tile indices for each of the neighbors
      ij_L = [itile - 1, jtile]
      ij_R = [itile + 1, jtile]
      ij_D = [itile, jtile - 1]
      ij_U = [itile, jtile + 1]
      ij_UL = [itile - 1, jtile + 1]
      ij_UR = [itile + 1, jtile + 1]
      ij_DL = [itile - 1, jtile - 1]
      ij_DR = [itile + 1, jtile - 1]

      if (periodic) then
         ! set neighbor to wrap around the edge
         if (ij_L(1) < 1) ij_L(1) = tiles(1)
         if (ij_R(1) > tiles(1)) ij_R(1) = 1
         if (ij_UL(1) < 1) ij_UL(1) = tiles(1)
         if (ij_UR(1) > tiles(1)) ij_UR(1) = 1
         if (ij_DL(1) < 1) ij_DL(1) = tiles(1)
         if (ij_DR(1) > tiles(1)) ij_DR(1) = 1

         if (ij_D(2) < 1) ij_D(2) = tiles(2)
         if (ij_U(2) > tiles(2)) ij_U(2) = 1
         if (ij_DL(2) < 1) ij_DL(2) = tiles(2)
         if (ij_UL(2) > tiles(2)) ij_UL(2) = 1
         if (ij_DR(2) < 1) ij_DR(2) = tiles(2)
         if (ij_UR(2) > tiles(2)) ij_UR(2) = 1
      else
         ! set neighbor to 0 -- no neighbor
         if (ij_L(1) < 1) ij_L(1) = 0
         if (ij_R(1) > tiles(1)) ij_R(1) = 0
         if (ij_UL(1) < 1) ij_UL(1) = 0
         if (ij_UR(1) > tiles(1)) ij_UR(1) = 0
         if (ij_DL(1) < 1) ij_DL(1) = 0
         if (ij_DR(1) > tiles(1)) ij_DR(1) = 0

         if (ij_D(2) < 1) ij_D(2) = 0
         if (ij_U(2) > tiles(2)) ij_U(2) = 0
         if (ij_DL(2) < 1) ij_DL(2) = 0
         if (ij_UL(2) > tiles(2)) ij_UL(2) = 0
         if (ij_DR(2) < 1) ij_DR(2) = 0
         if (ij_UR(2) > tiles(2)) ij_UR(2) = 0
      end if

      neighbors = [ &
                  tile_ij2n(ij_L), tile_ij2n(ij_R), &
                  tile_ij2n(ij_D), tile_ij2n(ij_U), &
                  tile_ij2n(ij_UL), tile_ij2n(ij_UR), &
                  tile_ij2n(ij_DL), tile_ij2n(ij_DR) &
                  ]

   end function tile_neighbors_2d

   subroutine init_tiles()
      integer(ik) :: indices(4), upper(2), lower(2)

      indices = tile_indices([config%nx, config%ny])
      lower = indices([1, 3])
      upper = indices([2, 4])

      isd = lower(1)
      ied = upper(1)
      jsd = lower(2)
      jed = upper(2)

      is = isd - halo_width
      ie = ied + halo_width
      js = jsd - halo_width
      je = jed + halo_width

      neighbours = tile_neighbors_2d(periodic=.true.)

   end subroutine init_tiles

end module mod_tiles
