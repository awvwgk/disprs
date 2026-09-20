module multicharge
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   implicit none
   private
   integer, parameter :: message_capacity = 512
   public :: mchrg_model_type, mchrg_model, new_eeq2019_model, new_eeqbc2025_model
   public :: get_charges, get_eeq_charges, get_eeqbc_charges, get_multicharge_version

   type :: mchrg_model_type
      private
      integer(c_int) :: kind = 0
   end type
   type :: charge_model_enum
      integer :: eeq2019 = 1, eeqbc2025 = 2
   end type
   type(charge_model_enum), parameter :: mchrg_model = charge_model_enum()
   interface
      function c_new_error() bind(C, name='disprs_d4_new_error') result(handle)
         import c_ptr
         type(c_ptr) :: handle
      end function
      function c_check(error) bind(C, name='disprs_d4_check_error') result(status)
         import c_ptr, c_int
         type(c_ptr), value :: error
         integer(c_int) :: status
      end function
      subroutine c_get_error(error, message, length) bind(C, name='disprs_d4_get_error')
         import c_ptr, c_char, c_int
         type(c_ptr), value :: error
         integer(c_int), intent(in) :: length
         character(kind=c_char), intent(out) :: message(length)
      end subroutine
      subroutine c_delete_error(error) bind(C, name='disprs_d4_delete_error')
         import c_ptr
         type(c_ptr) :: error
      end subroutine
      function c_new_structure(error, nat, numbers, xyz, charge, lattice, periodic) &
         bind(C, name='disprs_d4_new_structure') result(handle)
         import c_ptr, c_int
         type(c_ptr), value :: error, numbers, xyz, charge, lattice, periodic
         integer(c_int), value :: nat
         type(c_ptr) :: handle
      end function
      subroutine c_delete_structure(handle) bind(C, name='disprs_d4_delete_structure')
         import c_ptr
         type(c_ptr) :: handle
      end subroutine
      subroutine c_charges(error, structure, kind, charges, cartesian, strain) bind(C, name='disprs_get_charges')
         import c_ptr, c_int
         type(c_ptr), value :: error, structure, charges, cartesian, strain
         integer(c_int), value :: kind
      end subroutine
   end interface
contains
   subroutine new_eeq2019_model(mol, model, error)
      type(structure_type), intent(in) :: mol
      class(mchrg_model_type), allocatable, intent(out) :: model
      type(error_type), allocatable, intent(out) :: error
        if (mol%nat <= 0 .or. any(mol%num < 1) .or. any(mol%num > 118) &
           .or. any(mol%num >= 104 .and. mol%num <= 111)) then
         call fatal_error(error, 'Unsupported EEQ structure')
         return
      end if
      allocate(model)
   end subroutine

   subroutine new_eeqbc2025_model(mol, model, error)
      type(structure_type), intent(in) :: mol
      class(mchrg_model_type), allocatable, intent(out) :: model
      type(error_type), allocatable, intent(out) :: error
      if (any(mol%num > 103)) then
         call fatal_error(error, 'EEQBC supports atomic numbers 1 through 103')
         return
      end if
      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return
      model%kind = 1
   end subroutine

   subroutine get_charges(mchrg_model, mol, error, qvec, dqdr, dqdL)
      class(mchrg_model_type), intent(in) :: mchrg_model
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error
      real(wp), contiguous, target, intent(out) :: qvec(:)
      real(wp), contiguous, target, optional, intent(out) :: dqdr(:, :, :), dqdL(:, :, :)
      integer(c_int), target :: numbers(mol%nat)
      real(c_double), target :: positions(3, mol%nat), charge, lattice(3, 3)
      logical(c_bool), target :: periodic(3)
      type(c_ptr) :: failure, structure, lattice_ptr, cartesian_ptr, strain_ptr
      character(kind=c_char) :: buffer(message_capacity)
      character(len=message_capacity) :: message
      integer :: index
      if (size(qvec) /= mol%nat) then
         call fatal_error(error, 'Invalid charge array shape')
         return
      end if
      if (present(dqdr)) then
         if (any(shape(dqdr) /= [3, mol%nat, mol%nat])) then
            call fatal_error(error, 'Invalid Cartesian charge derivative shape')
            return
         end if
      end if
      if (present(dqdL)) then
         if (any(shape(dqdL) /= [3, 3, mol%nat])) then
            call fatal_error(error, 'Invalid strain charge derivative shape')
            return
         end if
      end if
      numbers = mol%num(mol%id)
      positions = mol%xyz
      charge = mol%charge
      periodic = .false.
      if (size(mol%periodic) == 1) periodic = mol%periodic(1)
      if (size(mol%periodic) == 3) periodic = mol%periodic
      lattice_ptr = c_null_ptr
      if (size(mol%lattice) == 9) then
         lattice = mol%lattice
         lattice_ptr = c_loc(lattice)
      end if
      cartesian_ptr = c_null_ptr
      strain_ptr = c_null_ptr
      if (present(dqdr)) cartesian_ptr = c_loc(dqdr)
      if (present(dqdL)) strain_ptr = c_loc(dqdL)
      failure = c_new_error()
      structure = c_new_structure(failure, int(mol%nat, c_int), c_loc(numbers), c_loc(positions), &
         c_loc(charge), lattice_ptr, c_loc(periodic))
      if (c_check(failure) == 0) then
         call c_charges(failure, structure, mchrg_model%kind, c_loc(qvec), cartesian_ptr, strain_ptr)
      end if
      if (c_check(failure) /= 0) then
         call c_get_error(failure, buffer, int(size(buffer), c_int))
         message = ''
         read_message: do index = 1, size(buffer)
            if (buffer(index) == c_null_char) exit read_message
            message(index:index) = buffer(index)
         end do read_message
         call fatal_error(error, trim(message))
      end if
      call c_delete_structure(structure)
      call c_delete_error(failure)
   end subroutine

   subroutine get_eeq_charges(mol, error, qvec, dqdr, dqdL)
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error
      real(wp), contiguous, intent(out) :: qvec(:)
      real(wp), contiguous, optional, intent(out) :: dqdr(:, :, :), dqdL(:, :, :)
      class(mchrg_model_type), allocatable :: model
      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return
      call get_charges(model, mol, error, qvec, dqdr, dqdL)
   end subroutine

   subroutine get_eeqbc_charges(mol, error, qvec, dqdr, dqdL)
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error
      real(wp), contiguous, intent(out) :: qvec(:)
      real(wp), contiguous, optional, intent(out) :: dqdr(:, :, :), dqdL(:, :, :)
      class(mchrg_model_type), allocatable :: model
      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call get_charges(model, mol, error, qvec, dqdr, dqdL)
   end subroutine

   subroutine get_multicharge_version(major, minor, patch, string)
      integer, optional, intent(out) :: major, minor, patch
      character(len=:), allocatable, optional, intent(out) :: string
      if (present(major)) major = 0
      if (present(minor)) minor = 5
      if (present(patch)) patch = 0
      if (present(string)) string = '0.5.0'
   end subroutine
end module multicharge