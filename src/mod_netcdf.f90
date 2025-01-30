module mod_netcdf

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_config, only: config => main_config
   use mod_util, only: abort_now
   use netcdf, only: &
      nf90_open, nf90_inquire, nf90_inquire_dimension, nf90_inquire_variable, &
      nf90_create, nf90_inq_dimid, nf90_def_dim, nf90_inq_varid, nf90_def_var, &
      nf90_put_var, nf90_get_var, nf90_enddef, nf90_redef, nf90_close, &
      nf90_float, nf90_nowrite, nf90_clobber, nf90_strerror, nf90_noerr, &
      nf90_unlimited

   implicit none

   private
   public :: netcdf_file, open_netcdf, create_netcdf

   integer(ik), parameter :: max_name_length = 10

   type :: netcdf_axis
      character(len=:), allocatable :: name
      integer(ik) :: size, dimid, varid
   end type

   type :: netcdf_var
      character(len=:), allocatable :: name
      integer(ik) :: varid
      character(len=max_name_length), allocatable :: dims(:)
   end type

   type :: netcdf_file
      integer(ik) :: ncid
      type(netcdf_axis) :: t_axis ! time axis
      type(netcdf_axis), allocatable :: axes(:)
      type(netcdf_var), allocatable :: vars(:)
   contains
      procedure :: close => close_netcdf
      procedure :: create_time_axis, advance_time
      procedure :: create_axis, get_axis_index, read_axis
      procedure :: create_variable, get_variable_index
      procedure :: write_variable, read_variable
   end type

contains

   function open_netcdf(filename) result(out)
      character(len=*), intent(in) :: filename
      type(netcdf_file) :: out
      type(netcdf_axis), allocatable :: axes(:)
      type(netcdf_var), allocatable :: vars(:)
      character(len=max_name_length) :: name
      integer(ik) :: status, i, j, length
      integer(ik) :: ncid, varid, ndims, nvars, tdim
      integer(ik), allocatable :: dimids(:)
      character(len=max_name_length), allocatable :: dimnames(:)

      status = nf90_open(trim(filename), nf90_nowrite, ncid)
      call handle_netcdf_error(status, 'open_netcdf')

      status = nf90_inquire(ncid, nDimensions=ndims, nVariables=nvars, &
                            unlimitedDimId=tdim)
      call handle_netcdf_error(status, 'open_netcdf')

      out = netcdf_file(ncid)
      axes = [netcdf_axis ::]
      vars = [netcdf_var ::]

      ! Get the unlimited (time) dimension, if it exists
      if (tdim /= -1) then
         status = nf90_inquire_dimension(ncid, tdim, name=name, len=length)
         call handle_netcdf_error(status, 'open_netcdf')

         varid = get_variable_matching_dimension(ncid, name)

         out%t_axis = netcdf_axis(name, length, tdim, varid)
      else
         out%t_axis = netcdf_axis('', 0, -1, -1)
      end if

      ! Get the other dimensions
      do i = 1, ndims
         status = nf90_inquire_dimension(ncid, i, name=name, len=length)
         call handle_netcdf_error(status, 'open_netcdf')

         if (trim(name) == trim(out%t_axis%name)) cycle

         varid = get_variable_matching_dimension(ncid, name)
         axes = [axes, netcdf_axis(name, length, i, varid)]
      end do

      ! Get the variables
      var_loop: do i = 1, nvars
         status = nf90_inquire_variable(ncid, i, name=name, ndims=ndims)
         call handle_netcdf_error(status, 'open_netcdf')

         do j = 1, size(axes)
            if (trim(name) == trim(axes(j)%name)) cycle var_loop
         end do

         allocate (dimids(1:ndims))

         status = nf90_inquire_variable(ncid, i, dimids=dimids)
         call handle_netcdf_error(status, 'open_netcdf')

         allocate (dimnames(1:ndims))
         do j = 1, ndims
            status = nf90_inquire_dimension(ncid, dimids(j), name=dimnames(j))
            call handle_netcdf_error(status, 'open_netcdf')
         end do

         vars = [vars, netcdf_var(name, i, dimnames)]

         deallocate (dimids)
         deallocate (dimnames)
      end do var_loop

      out%axes = axes
      out%vars = vars

   end function

   function get_variable_matching_dimension(ncid, name) result(varid)
      integer(ik), intent(in) :: ncid
      character(len=*), intent(in) :: name
      character(len=max_name_length) :: other_name
      integer(ik) :: status, i, nvars, varid

      status = nf90_inquire(ncid, nVariables=nvars)
      call handle_netcdf_error(status, 'get_variable_matching_dimension')

      varid = -1
      do i = 1, nvars
         status = nf90_inquire_variable(ncid, i, name=other_name)
         call handle_netcdf_error(status, 'get_variable_matching_dimension')

         if (trim(name) == trim(other_name)) then
            varid = i
            exit
         end if
      end do

      if (varid == -1) then
         write (log_str, '(a)') &
            'Could not find variable matching dimension `'//trim(name)//'`.'
         call logger%fatal('get_variable_matching_dimension', log_str)
         call abort_now()
      end if

   end function get_variable_matching_dimension

   function create_netcdf(filename) result(out)
      character(len=*), intent(in) :: filename
      type(netcdf_file) :: out
      integer(ik) :: status, ncid

      status = nf90_create(path=filename, cmode=nf90_clobber, ncid=ncid)
      call handle_netcdf_error(status, 'create_netcdf')

      status = nf90_enddef(ncid)
      call handle_netcdf_error(status, 'create_netcdf')

      out = netcdf_file(ncid)

      out%axes = [netcdf_axis ::]
      out%vars = [netcdf_var ::]

   end function

   subroutine close_netcdf(self)
      class(netcdf_file), intent(in) :: self
      integer(ik) :: status

      status = nf90_close(self%ncid)
      call handle_netcdf_error(status, 'close_netcdf')

   end subroutine close_netcdf

   subroutine create_time_axis(self)
      class(netcdf_file), intent(inout) :: self
      integer(ik) :: status, dimid, varid

      status = nf90_redef(self%ncid)
      call handle_netcdf_error(status, 'create_time_axis')

      status = nf90_def_dim(self%ncid, 't', nf90_unlimited, dimid)
      call handle_netcdf_error(status, 'create_time_axis')

      status = nf90_def_var(self%ncid, 't', nf90_float, dimid, varid)
      call handle_netcdf_error(status, 'create_time_axis')

      self%t_axis = netcdf_axis('t', 0, dimid, varid)

      status = nf90_enddef(self%ncid)
      call handle_netcdf_error(status, 'create_time_axis')

   end subroutine create_time_axis

   subroutine advance_time(self, time)
      class(netcdf_file), intent(inout) :: self
      real(rk), intent(in) :: time
      integer(ik) :: status

      self%t_axis%size = self%t_axis%size + 1

      status = nf90_put_var(self%ncid, self%t_axis%varid, &
                            time, start=[self%t_axis%size])
      call handle_netcdf_error(status, 'advance_time')

   end subroutine advance_time

   subroutine create_axis(self, name, points)
      class(netcdf_file), intent(inout) :: self
      character(len=*), intent(in) :: name
      real(rk), intent(in) :: points(:)
      integer(ik) :: status, dimid, varid, i
      type(netcdf_axis), allocatable :: axes(:)
      type(netcdf_axis) :: axis

      status = nf90_redef(self%ncid)
      call handle_netcdf_error(status, 'create_axis')

      do i = 1, size(self%axes)
         if (trim(self%axes(i)%name) == trim(name)) then
            write (log_str, '(a)') &
               'Axis by name `'//trim(name)//'` already exists.'
            call logger%fatal('create_axis', log_str)
            call abort_now()
         end if
      end do

      status = nf90_def_dim(self%ncid, trim(name), size(points), dimid)
      call handle_netcdf_error(status, 'create_axis')

      status = nf90_def_var(self%ncid, trim(name), nf90_float, dimid, varid)
      call handle_netcdf_error(status, 'create_axis')

      axis = netcdf_axis(name, size(points), dimid, varid)

      status = nf90_enddef(self%ncid)
      call handle_netcdf_error(status, 'create_axis')

      status = nf90_put_var(self%ncid, varid, points)
      call handle_netcdf_error(status, 'create_axis')

      axes = [self%axes, axis]
      self%axes = axes

   end subroutine create_axis

   function get_axis_index(self, name) result(idx)
      class(netcdf_file), intent(in) :: self
      character(len=*), intent(in) :: name
      integer(ik) :: i, idx

      idx = -1
      do i = 1, size(self%axes)
         if (trim(self%axes(i)%name) == trim(name)) then
            idx = i
         end if
      end do

   end function get_axis_index

   function read_axis(self, name) result(values)
      class(netcdf_file), intent(in) :: self
      character(len=*), intent(in) :: name
      integer(ik), allocatable :: values(:)
      integer(ik) :: status, idx, size, varid

      if (trim(self%t_axis%name) == trim(name)) then
         size = self%t_axis%size
         varid = self%t_axis%varid
      else
         idx = self%get_axis_index(name)
         size = self%axes(idx)%size
         varid = self%axes(idx)%varid
      end if

      allocate (values(1:size))

      status = nf90_get_var(self%ncid, varid, values)
      call handle_netcdf_error(status, 'read_axis')

   end function read_axis

   subroutine create_variable(self, name, dims)
      class(netcdf_file), intent(inout) :: self
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: dims(:)
      type(netcdf_var) :: var
      type(netcdf_var), allocatable :: vars(:)
      integer(ik) :: dimids(size(dims))
      character(len=max_name_length) :: dimnames(size(dims))
      integer(ik) :: status, varid, i, j, n

      do i = 1, size(self%vars)
         if (self%vars(i)%name == trim(name)) then
            write (log_str, '(a)') &
               'Variable by name `'//trim(name)//'` already exists.'
            call logger%fatal('create_variable', log_str)
            call abort_now()
         end if
      end do

      ! Translate dimension names into dimension IDs; do time first
      ! since we have to put it as the last dimension if it's present

      dimids(:) = -1
      dimnames(:) = ''

      do i = 1, size(dims)
         if (trim(dims(i)) == self%t_axis%name) then
            dimids(size(dims)) = self%t_axis%dimid
            dimnames(size(dims)) = self%t_axis%name
         end if
      end do

      n = 1
      do i = 1, size(dims)
         if (trim(dims(i)) == 't') cycle
         do j = 1, size(self%axes)
            if (trim(dims(i)) == trim(self%axes(j)%name)) then
               dimids(n) = self%axes(j)%dimid
               dimnames(n) = self%axes(j)%name
            end if
         end do
         if (dimids(n) == -1) then
            write (log_str, '(a)') 'Unknown dimension `'//trim(dims(i))//'`.'
            call logger%fatal('create_variable', log_str)
            call abort_now()
         end if
         n = n + 1
      end do

      status = nf90_redef(self%ncid)
      call handle_netcdf_error(status, 'create_variable')

      status = nf90_def_var(self%ncid, name, nf90_float, dimids, varid)
      call handle_netcdf_error(status, 'create_variable')

      status = nf90_enddef(self%ncid)
      call handle_netcdf_error(status, 'create_variable')

      var = netcdf_var(name, varid, dimnames)

      vars = [self%vars, var]
      self%vars = vars

   end subroutine create_variable

   function get_variable_index(self, name) result(idx)
      class(netcdf_file), intent(in) :: self
      character(len=*), intent(in) :: name
      integer(ik) :: i, idx

      idx = -1
      do i = 1, size(self%vars)
         if (trim(self%vars(i)%name) == trim(name)) then
            idx = i
         end if
      end do

   end function get_variable_index

   subroutine write_variable(self, name, data)
      class(netcdf_file), intent(inout) :: self
      character(len=*), intent(in) :: name
      real(rk), intent(in) :: data(:, :)
      type(netcdf_var) :: var
      integer(ik) :: status, idx

      idx = self%get_variable_index(name)

      if (idx == -1) then
         write (log_str, '(a)') 'Unknown variable `'//trim(name)//'`.'
         call logger%fatal('write_variable', log_str)
         call abort_now()
      end if

      var = self%vars(idx)

      if (any(var%dims == self%t_axis%name)) then
         status = nf90_put_var( &
                  self%ncid, var%varid, data, &
                  start=[1, 1, self%t_axis%size], &
                  count=[size(data, 1), size(data, 2), 1])
      else
         status = nf90_put_var( &
                  self%ncid, var%varid, data, &
                  start=[1, 1], &
                  count=[size(data, 1), size(data, 2)])
      end if

      call handle_netcdf_error(status, 'write_variable')

   end subroutine write_variable

   subroutine read_variable(self, name, data)
      class(netcdf_file), intent(in) :: self
      character(len=*), intent(in) :: name
      real(rk), intent(inout) :: data(:, :)
      integer(ik), allocatable :: start(:), count(:), dim_map(:)
      integer(ik) :: status, idx, i, j

      idx = self%get_variable_index(name)

      if (idx == -1) then
         write (log_str, '(a)') 'Unknown variable `'//trim(name)//'`.'
         call logger%fatal('read_variable', log_str)
         call abort_now()
      end if

      allocate (start(1:size(self%vars(idx)%dims)))
      allocate (count(1:size(self%vars(idx)%dims)))
      allocate (dim_map(1:size(self%vars(idx)%dims)))

      do i = 1, size(self%vars(idx)%dims)
         if (self%vars(idx)%dims(i) == self%t_axis%name) then
            start(i) = self%t_axis%size
            count(i) = 1
            dim_map(i) = 0
         else
            start(i) = 1

            j = self%get_axis_index(self%vars(idx)%dims(i))
            count(i) = self%axes(j)%size

            if (any(self%vars(idx)%dims(i) == ['x ', 'xc', 'xf'])) then
               dim_map(i) = 1
            else if (any(self%vars(idx)%dims(i) == ['y ', 'yc', 'yf'])) then
               dim_map(i) = config%nx
            end if
         end if
      end do

      status = nf90_get_var( &
               self%ncid, self%vars(idx)%varid, data, &
               start=start, count=count, map=dim_map)

      call handle_netcdf_error(status, 'read_variable')

   end subroutine read_variable

   subroutine handle_netcdf_error(errcode, source)
      integer(ik), intent(in) :: errcode
      character(len=*), intent(in) :: source

      if (errcode /= nf90_noerr) then
         write (log_str, *) trim(nf90_strerror(errcode))
         call logger%fatal(source, log_str)
         call abort_now()
      end if

   end subroutine handle_netcdf_error

end module mod_netcdf
