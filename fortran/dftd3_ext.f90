module dftd3
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_cutoff, only : get_lattice_points
   implicit none
   private
   integer, parameter :: damping_parameter_count = 9, gcp_scalar_count = 7

   type, public :: realspace_cutoff
      real(wp) :: cn = 40.0_wp
      real(wp) :: disp2 = 60.0_wp
      real(wp) :: disp3 = 40.0_wp
      real(wp) :: width2 = 0.0_wp
      real(wp) :: width3 = 0.0_wp
      real(wp) :: gcp = 60.0_wp
      real(wp) :: srb = 60.0_wp
   end type

   type, public :: d3_param
      real(wp) :: s6 = 1.0_wp, s8 = 1.0_wp, s9 = 0.0_wp
      real(wp) :: rs6 = 1.0_wp, rs8 = 1.0_wp
      real(wp) :: a1 = 0.4_wp, a2 = 5.0_wp, alp = 14.0_wp, bet = 0.0_wp
   end type

   type, public :: d3_lowrank_config
      integer :: rank = 0
      real(wp) :: tolerance = 1.0e-4_wp, kcut = 0.0_wp
      integer :: mesh = 0
   end type

   type, public :: work_partition
      integer :: part = 0, nparts = 1
   end type
   type(work_partition), public, parameter :: serial_work_partition = work_partition()

   type, abstract, public :: damping_param
      type(c_ptr) :: handle = c_null_ptr
   contains
      procedure :: close => close_damping_param
   end type
   type, extends(damping_param), public :: zero_damping_param; end type
   type, extends(damping_param), public :: rational_damping_param; end type
   type, extends(damping_param), public :: mzero_damping_param; end type
   type, extends(damping_param), public :: optimizedpower_damping_param; end type
   type, extends(damping_param), public :: cso_damping_param; end type
   type, extends(damping_param), public :: z_damping_param; end type

   type, public :: d3_model
      type(c_ptr) :: error = c_null_ptr
      type(c_ptr) :: structure = c_null_ptr
      type(c_ptr) :: handle = c_null_ptr
      integer :: nat = 0
   contains
      procedure :: close => close_d3_model
   end type

   type, public :: gcp_param
      type(c_ptr) :: error = c_null_ptr
      type(c_ptr) :: structure = c_null_ptr
      type(c_ptr) :: handle = c_null_ptr
      integer :: nat = 0
      real(c_double) :: eta = 0.0_c_double
      logical :: base = .false., srb = .false., damp = .false.
      real(wp), allocatable :: emiss(:), xv(:), slater(:), rvdw(:, :), rvdw_srb(:, :)
      integer, allocatable :: zeff(:)
      real(wp) :: sigma = 0.0_wp, alpha = 0.0_wp, beta = 0.0_wp
      real(wp) :: dmp_scal = 4.0_wp, dmp_exp = 6.0_wp, rscal = 0.0_wp, qscal = 0.0_wp
      real(wp), private :: reference_eta = 0.0_wp
   contains
      final :: delete_gcp_param
   end type

   public :: new_d3_model, new_zero_damping, new_rational_damping
   public :: get_lattice_points
   public :: new_mzero_damping, new_optimizedpower_damping, new_cso_damping, new_z_damping
   public :: get_dispersion, new_work_partition, get_dftd3_version
   public :: get_pairwise_dispersion, get_properties
   public :: get_zero_damping, get_rational_damping, get_mzero_damping, get_mrational_damping
   public :: get_optimizedpower_damping, get_cso_damping, get_z_damping
   public :: get_gcp_param, get_geometric_counterpoise, get_geometric_counterpoise_hessian

   interface get_dispersion
      module procedure get_dispersion_scalar
      module procedure get_dispersion_scalar_error
      module procedure get_dispersion_atomic
      module procedure get_dispersion_atomic_error
   end interface

   interface get_pairwise_dispersion
      module procedure get_pairwise_dispersion_legacy
      module procedure get_pairwise_dispersion_error
   end interface

   interface get_properties
      module procedure get_properties_legacy
      module procedure get_properties_error
   end interface

   interface
      subroutine c_get_error(error, message, size) bind(C, name='disprs_d3_get_error')
         import c_ptr, c_char, c_int
         type(c_ptr), value :: error
         integer(c_int), intent(in) :: size
         character(kind=c_char), intent(out) :: message(size)
      end subroutine
      subroutine c_properties(error, structure, model, cn, c6) bind(C, name='disprs_d3_get_properties')
         import c_ptr
         type(c_ptr), value :: error, structure, model, cn, c6
      end subroutine
      subroutine c_property_response(error, structure, model, cn, c6, dcndr, dcndL, dc6dr, dc6dL) &
         bind(C, name='disprs_d3_get_property_response')
         import c_ptr
         type(c_ptr), value :: error, structure, model, cn, c6, dcndr, dcndL, dc6dr, dc6dL
      end subroutine
      subroutine c_get_gcp_parameters(error, handle, nat, scalars, flags, zeff, emiss, xv, slater, rvdw, rvdw_srb) &
         bind(C, name='disprs_d3_get_gcp_parameters')
         import c_ptr, c_int
         type(c_ptr), value :: error, handle, scalars, flags, zeff, emiss, xv, slater, rvdw, rvdw_srb
         integer(c_int), value :: nat
      end subroutine
      subroutine c_set_gcp_parameters(error, handle, nat, scalars, flags, zeff, emiss, xv, slater, rvdw, rvdw_srb) &
         bind(C, name='disprs_d3_set_gcp_parameters')
         import c_ptr, c_int
         type(c_ptr), value :: error, handle, scalars, flags, zeff, emiss, xv, slater, rvdw, rvdw_srb
         integer(c_int), value :: nat
      end subroutine
      subroutine c_get_gcp_controls(error, handle, eta, base, srb) bind(C, name='disprs_d3_get_gcp_controls')
         import c_ptr, c_double, c_bool
         type(c_ptr), value :: error, handle
         real(c_double), intent(out) :: eta
         logical(c_bool), intent(out) :: base, srb
      end subroutine
      subroutine c_set_gcp_controls(error, handle, eta, base, srb) bind(C, name='disprs_d3_set_gcp_controls')
         import c_ptr, c_double, c_bool
         type(c_ptr), value :: error, handle
         real(c_double), intent(in) :: eta
         logical(c_bool), intent(in) :: base, srb
      end subroutine
      subroutine c_named_parameters(error, damping, method, values) bind(C, name='disprs_d3_get_named_parameters')
         import c_ptr, c_int, c_char, c_double, damping_parameter_count
         type(c_ptr), value :: error
         integer(c_int), value :: damping
         character(c_char), intent(in) :: method(*)
         real(c_double), intent(out) :: values(damping_parameter_count)
      end subroutine
      function c_new_error() bind(C, name='disprs_d3_new_error') result(ptr)
         import c_ptr; type(c_ptr) :: ptr
      end function
      function c_new_structure(error, nat, num, xyz, lattice, periodic) bind(C, name='disprs_d3_new_structure') result(ptr)
         import c_ptr, c_int
         type(c_ptr), value :: error, num, xyz, lattice, periodic
         integer(c_int), value :: nat
         type(c_ptr) :: ptr
      end function
      function c_new_model(error, structure, kind) bind(C, name='disprs_d3_new_model_kind') result(ptr)
         import c_ptr, c_int
         type(c_ptr), value :: error, structure
         integer(c_int), value :: kind
         type(c_ptr) :: ptr
      end function
      function c_load_gcp(error, structure, method, basis) bind(C, name='disprs_d3_load_gcp') result(ptr)
         import c_ptr
         type(c_ptr), value :: error, structure, method, basis
         type(c_ptr) :: ptr
      end function
      subroutine c_set_gcp_cutoff(error, gcp, bas, srb) bind(C, name='disprs_d3_set_gcp_realspace_cutoff')
         import c_ptr, c_double
         type(c_ptr), value :: error, gcp
         real(c_double), value :: bas, srb
      end subroutine
      subroutine c_set_gcp_partition(error, gcp, part, nparts) bind(C, name='disprs_d3_set_gcp_work_partition')
         import c_ptr, c_int
         type(c_ptr), value :: error, gcp
         integer(c_int), value :: part, nparts
      end subroutine
      subroutine c_eval_gcp(error, structure, gcp, energy, gradient, sigma) bind(C, name='disprs_d3_get_counterpoise')
         import c_ptr
         type(c_ptr), value :: error, structure, gcp, energy, gradient, sigma
      end subroutine
      subroutine c_hessian_gcp(error, structure, gcp, energy, hessian) bind(C, name='disprs_d3_get_counterpoise_hessian')
         import c_ptr
         type(c_ptr), value :: error, structure, gcp, energy, hessian
      end subroutine
      subroutine c_delete_gcp(gcp) bind(C, name='disprs_d3_delete_gcp')
         import c_ptr
         type(c_ptr) :: gcp
      end subroutine
      subroutine c_delete_structure(structure) bind(C, name='disprs_d3_delete_structure')
         import c_ptr
         type(c_ptr) :: structure
      end subroutine
      subroutine c_delete_model(model) bind(C, name='disprs_d3_delete_model')
         import c_ptr
         type(c_ptr) :: model
      end subroutine
      subroutine c_delete_param(param) bind(C, name='disprs_d3_delete_param')
         import c_ptr
         type(c_ptr) :: param
      end subroutine
      subroutine c_delete_error(error) bind(C, name='disprs_d3_delete_error')
         import c_ptr
         type(c_ptr) :: error
      end subroutine
      subroutine c_update(error, structure, xyz, lattice) bind(C, name='disprs_d3_update_structure')
         import c_ptr; type(c_ptr), value :: error, structure, xyz, lattice
      end subroutine
      function c_new_param(error, s6, s8, s9, p1, p2, alp) bind(C, name='disprs_d3_new_rational_damping') result(ptr)
         import c_ptr, c_double
         type(c_ptr), value :: error
         real(c_double), value :: s6, s8, s9, p1, p2, alp
         type(c_ptr) :: ptr
      end function
      function c_new_zero(error, s6, s8, s9, rs6, rs8, alp) bind(C, name='disprs_d3_new_zero_damping') result(ptr)
         import c_ptr, c_double
         type(c_ptr), value :: error; real(c_double), value :: s6, s8, s9, rs6, rs8, alp; type(c_ptr) :: ptr
      end function
      function c_new_mzero(error, s6, s8, s9, rs6, rs8, alp, bet) bind(C, name='disprs_d3_new_mzero_damping') result(ptr)
         import c_ptr, c_double
         type(c_ptr), value :: error; real(c_double), value :: s6, s8, s9, rs6, rs8, alp, bet; type(c_ptr) :: ptr
      end function
      function c_new_op(error, s6, s8, s9, a1, a2, alp, bet) bind(C, name='disprs_d3_new_optimizedpower_damping') result(ptr)
         import c_ptr, c_double
         type(c_ptr), value :: error; real(c_double), value :: s6, s8, s9, a1, a2, alp, bet; type(c_ptr) :: ptr
      end function
      function c_new_cso(error, s6, s9, a1, a2, a3, a4, alp) bind(C, name='disprs_d3_new_cso_damping') result(ptr)
         import c_ptr, c_double
         type(c_ptr), value :: error; real(c_double), value :: s6, s9, a1, a2, a3, a4, alp; type(c_ptr) :: ptr
      end function
      function c_new_z(error, s6, s8, s9, a1, alp) bind(C, name='disprs_d3_new_z_damping') result(ptr)
         import c_ptr, c_double
         type(c_ptr), value :: error; real(c_double), value :: s6, s8, s9, a1, alp; type(c_ptr) :: ptr
      end function
      subroutine c_set_cutoff(error, model, disp2, disp3, cn, width2, width3) &
         bind(C, name='disprs_d3_set_model_realspace_cutoff_smooth')
         import c_ptr, c_double
         type(c_ptr), value :: error, model
         real(c_double), value :: disp2, disp3, cn, width2, width3
      end subroutine
      subroutine c_set_ewald(error, model, rank, tolerance, kcut, mesh) bind(C, name='disprs_d3_set_model_ewald')
         import c_ptr, c_int, c_double
         type(c_ptr), value :: error, model
         integer(c_int), value :: rank, mesh
         real(c_double), value :: tolerance, kcut
      end subroutine
      subroutine c_set_ghost(error, model, ghost, count) bind(C, name='disprs_d3_set_model_ghost_index')
         import c_ptr, c_int
         type(c_ptr), value :: error, model, ghost
         integer(c_int), value :: count
      end subroutine
      subroutine c_set_partition(error, model, part, parts) bind(C, name='disprs_d3_set_model_work_partition')
         import c_ptr, c_int
         type(c_ptr), value :: error, model
         integer(c_int), value :: part, parts
      end subroutine
      subroutine c_eval(error, structure, model, param, energy, gradient, sigma) bind(C, name='disprs_d3_get_dispersion')
         import c_ptr; type(c_ptr), value :: error, structure, model, param, energy, gradient, sigma
      end subroutine
      subroutine c_hessian(error, structure, model, param, energy, hessian) bind(C, name='disprs_d3_get_dispersion_hessian')
         import c_ptr; type(c_ptr), value :: error, structure, model, param, energy, hessian
      end subroutine
      subroutine c_pairwise(error, structure, model, param, pair2, pair3) bind(C, name='disprs_d3_get_pairwise_dispersion')
         import c_ptr
         type(c_ptr), value :: error, structure, model, param, pair2, pair3
      end subroutine
      function c_check(error) bind(C, name='disprs_d3_check_error') result(status)
         import c_ptr, c_int; type(c_ptr), value :: error; integer(c_int) :: status
      end function
   end interface

contains

   subroutine get_named_damping(param, method, error, damping, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      integer(c_int), intent(in) :: damping
      real(wp), intent(in), optional :: s9
      type(c_ptr) :: handle
      real(c_double) :: values(damping_parameter_count)
      handle = c_new_error()
      call c_named_parameters(handle, damping, trim(method)//c_null_char, values)
      if (c_check(handle) /= 0) then
         call fatal_error(error, "No entry for '"//method//"' present")
      else
         param = d3_param(s6=values(1), s8=values(2), s9=values(3), rs6=values(4), &
            rs8=values(5), a1=values(6), a2=values(7), alp=values(8), bet=values(9))
         if (present(s9)) param%s9 = s9
      end if
      call c_delete_error(handle)
   end subroutine

   subroutine get_zero_damping(param, method, error, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: s9
      call get_named_damping(param, method, error, 0_c_int, s9)
   end subroutine

   subroutine get_rational_damping(param, method, error, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: s9
      call get_named_damping(param, method, error, 1_c_int, s9)
   end subroutine

   subroutine get_mzero_damping(param, method, error, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: s9
      call get_named_damping(param, method, error, 2_c_int, s9)
   end subroutine

   subroutine get_mrational_damping(param, method, error, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: s9
      call get_named_damping(param, method, error, 3_c_int, s9)
   end subroutine

   subroutine get_optimizedpower_damping(param, method, error, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: s9
      call get_named_damping(param, method, error, 4_c_int, s9)
   end subroutine

   subroutine get_cso_damping(param, method, error, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: s9
      call get_named_damping(param, method, error, 5_c_int, s9)
   end subroutine

   subroutine get_z_damping(param, method, error, s9)
      type(d3_param), intent(out) :: param
      character(len=*), intent(in) :: method
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: s9
      call get_named_damping(param, method, error, 6_c_int, s9)
   end subroutine

   subroutine get_gcp_param(self, mol, method, basis, eta, error)
      type(gcp_param), intent(out) :: self
      type(structure_type), intent(in) :: mol
      character(len=*), intent(in), optional :: method, basis
      real(wp), intent(in), optional :: eta
      type(error_type), allocatable, intent(out), optional :: error
      type(error_type), allocatable :: failure
      integer(c_int), allocatable, target :: numbers(:)
      real(c_double), allocatable, target :: coordinates(:, :), lattice(:, :)
      logical(c_bool), target :: periodic(3)
      character(kind=c_char, len=:), allocatable, target :: c_method, c_basis
      type(c_ptr) :: method_ptr, basis_ptr, lattice_ptr
      integer :: atom, other, species
      real(c_double), target :: scalars(gcp_scalar_count), emiss(mol%nat), xv(mol%nat), slater(mol%nat)
      real(c_double), target :: rvdw(mol%nat,mol%nat), rvdw_srb(mol%nat,mol%nat)
      integer(c_int), target :: zeff(mol%nat)
      logical(c_bool), target :: flags(3)
      if (.not.present(method) .and. .not.present(basis) .and. .not.present(eta)) return
      allocate(numbers(mol%nat))
      coordinates = mol%xyz
      do atom = 1, mol%nat; numbers(atom) = mol%num(mol%id(atom)); end do
      periodic = .false._c_bool
      if (size(mol%periodic) == 1) periodic = logical(mol%periodic(1), c_bool)
      if (size(mol%periodic) == 3) periodic = logical(mol%periodic, c_bool)
      lattice_ptr = c_null_ptr
      if (size(mol%lattice) == 9) then
         lattice = mol%lattice
         lattice_ptr = c_loc(lattice)
      end if
      method_ptr = c_null_ptr
      if (present(method)) then
         c_method = method // c_null_char
         method_ptr = c_loc(c_method)
      end if
      basis_ptr = c_null_ptr
      if (present(basis)) then
         c_basis = basis // c_null_char
         basis_ptr = c_loc(c_basis)
      end if
      self%error = c_new_error()
      self%structure = c_new_structure(self%error, mol%nat, c_loc(numbers), c_loc(coordinates), &
         lattice_ptr, c_loc(periodic))
      setup: block
         if (gcp_failed(self%error, failure)) exit setup
         self%handle = c_load_gcp(self%error, self%structure, method_ptr, basis_ptr)
         if (gcp_failed(self%error, failure)) exit setup
         call c_get_gcp_controls(self%error, self%handle, self%eta, flags(2), flags(3))
         if (gcp_failed(self%error, failure)) exit setup
         self%base = flags(2); self%srb = flags(3)
         if (present(eta)) self%eta = eta
         self%reference_eta = self%eta
         call c_set_gcp_controls(self%error, self%handle, self%eta, logical(self%base,c_bool), logical(self%srb,c_bool))
         if (gcp_failed(self%error, failure)) exit setup
         call c_get_gcp_parameters(self%error, self%handle, mol%nat, c_loc(scalars), c_loc(flags), &
            c_loc(zeff), c_loc(emiss), c_loc(xv), c_loc(slater), c_loc(rvdw), c_loc(rvdw_srb))
         if (gcp_failed(self%error, failure)) exit setup
      end block setup
      if (allocated(failure)) then
         call delete_gcp_param(self)
         call gcp_error(error, failure%message)
         return
      end if
      self%nat = mol%nat
      self%sigma = scalars(1); self%alpha = scalars(2); self%beta = scalars(3)
      self%dmp_scal = scalars(4); self%dmp_exp = scalars(5); self%rscal = scalars(6); self%qscal = scalars(7)
      self%damp = flags(1); self%base = flags(2); self%srb = flags(3)
      species = size(mol%num)
      allocate(self%zeff(species), self%rvdw(species,species), self%rvdw_srb(species,species))
      if (any(abs(emiss) > tiny(0.0_wp))) allocate(self%emiss(species), self%xv(species))
      if (self%eta > 0.0_wp) allocate(self%slater(species))
      do atom = 1, mol%nat
         self%zeff(mol%id(atom)) = zeff(atom)
         if (allocated(self%emiss)) self%emiss(mol%id(atom)) = emiss(atom)
         if (allocated(self%xv)) self%xv(mol%id(atom)) = xv(atom)
         if (allocated(self%slater)) self%slater(mol%id(atom)) = slater(atom)
         do other = 1, mol%nat
            self%rvdw(mol%id(atom),mol%id(other)) = rvdw(atom,other)
            self%rvdw_srb(mol%id(atom),mol%id(other)) = rvdw_srb(atom,other)
         end do
      end do
   end subroutine

   subroutine gcp_error(error, message)
      type(error_type), allocatable, intent(out), optional :: error
      character(len=*), intent(in) :: message
      if (present(error)) then
         call fatal_error(error, message)
      else
         error stop message
      end if
   end subroutine

   logical function gcp_failed(handle, error) result(failed)
      type(c_ptr), intent(in) :: handle
      type(error_type), allocatable, intent(out), optional :: error
      character(kind=c_char, len=512) :: message
      failed = c_check(handle) /= 0
      if (.not. failed) return
      call c_get_error(handle, message, len(message, kind=c_int))
      call gcp_error(error, message(:index(message, c_null_char)-1))
   end function

   logical function sync_gcp_parameters(param, mol, error) result(valid)
      type(gcp_param), intent(in) :: param
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out), optional :: error
      real(c_double), target :: scalars(gcp_scalar_count), emiss(mol%nat), xv(mol%nat), slater(mol%nat)
      real(c_double), target :: rvdw(mol%nat,mol%nat), rvdw_srb(mol%nat,mol%nat)
      integer(c_int), target :: zeff(mol%nat)
      logical(c_bool), target :: flags(3)
      integer :: species
      valid = .false.
      species = size(mol%num)
      if (.not. allocated(param%zeff)) then
         call gcp_error(error, 'gCP parameters are missing')
         return
      end if
      if (size(param%zeff) /= species) then
         call gcp_error(error, 'invalid gCP zeff shape')
         return
      end if
      zeff = param%zeff(mol%id)
      emiss = 0.0_wp; xv = 0.0_wp; slater = 0.0_wp
      if (allocated(param%emiss)) then
         if (size(param%emiss) /= species) then
            call gcp_error(error, 'invalid gCP emiss shape')
            return
         end if
         emiss = param%emiss(mol%id)
      end if
      if (allocated(param%xv)) then
         if (size(param%xv) /= species) then
            call gcp_error(error, 'invalid gCP xv shape')
            return
         end if
         xv = param%xv(mol%id)
      end if
      if (allocated(param%slater)) then
         if (size(param%slater) /= species) then
            call gcp_error(error, 'invalid gCP slater shape')
            return
         end if
         slater = param%slater(mol%id)
         if (param%reference_eta > 0.0_wp) slater = slater * (param%eta / param%reference_eta)
      end if
      if (.not. allocated(param%rvdw) .or. .not. allocated(param%rvdw_srb)) then
         call gcp_error(error, 'gCP radii are missing')
         return
      end if
      if (any(shape(param%rvdw) /= [species,species]) .or. &
          any(shape(param%rvdw_srb) /= [species,species])) then
         call gcp_error(error, 'invalid gCP radii shape')
         return
      end if
      rvdw = param%rvdw(mol%id,mol%id); rvdw_srb = param%rvdw_srb(mol%id,mol%id)
      scalars = [param%sigma, param%alpha, param%beta, param%dmp_scal, param%dmp_exp, param%rscal, param%qscal]
      flags = [param%damp, param%base, param%srb]
      call c_set_gcp_controls(param%error, param%handle, param%eta, logical(param%base,c_bool), logical(param%srb,c_bool))
      if (gcp_failed(param%error, error)) return
      call c_set_gcp_parameters(param%error, param%handle, mol%nat, c_loc(scalars), c_loc(flags), &
         c_loc(zeff), c_loc(emiss), c_loc(xv), c_loc(slater), c_loc(rvdw), c_loc(rvdw_srb))
      if (gcp_failed(param%error, error)) return
      valid = .true.
   end function

   logical function prepare_gcp(mol, param, cutoff, partition, error) result(valid)
      type(structure_type), intent(in) :: mol
      type(gcp_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      type(work_partition), intent(in), optional :: partition
      type(error_type), allocatable, intent(out), optional :: error
      real(c_double), allocatable, target :: coordinates(:, :), lattice(:, :)
      type(c_ptr) :: lattice_ptr
      integer :: part, nparts
      valid = .false.
      if (param%nat /= mol%nat) then
         call gcp_error(error, 'gCP atom count changed')
         return
      end if
      coordinates = mol%xyz
      lattice_ptr = c_null_ptr
      if (size(mol%lattice) == 9) then
         lattice = mol%lattice
         lattice_ptr = c_loc(lattice)
      end if
      call c_update(param%error, param%structure, c_loc(coordinates), lattice_ptr)
      if (gcp_failed(param%error, error)) return
      if (.not. sync_gcp_parameters(param, mol, error)) return
      call c_set_gcp_cutoff(param%error, param%handle, cutoff%gcp, cutoff%srb)
      if (gcp_failed(param%error, error)) return
      part = 0; nparts = 1
      if (present(partition)) then
         part = partition%part; nparts = partition%nparts
      end if
      call c_set_gcp_partition(param%error, param%handle, part, nparts)
      if (gcp_failed(param%error, error)) return
      valid = .true.
   end function

   subroutine get_geometric_counterpoise(mol, param, cutoff, energy, gradient, sigma, partition, error)
      type(structure_type), intent(in) :: mol
      type(gcp_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(inout), target :: energy
      real(wp), intent(inout), target, contiguous, optional :: gradient(:, :), sigma(:, :)
      type(work_partition), intent(in), optional :: partition
      type(error_type), allocatable, intent(out), optional :: error
      type(c_ptr) :: gradient_ptr, sigma_ptr
      if (present(gradient)) then
         if (any(shape(gradient) /= [3, mol%nat])) then
            call gcp_error(error, 'invalid gCP gradient shape')
            return
         end if
      end if
      if (present(sigma)) then
         if (any(shape(sigma) /= [3, 3])) then
            call gcp_error(error, 'invalid gCP virial shape')
            return
         end if
      end if
      if (.not. prepare_gcp(mol, param, cutoff, partition, error)) return
      gradient_ptr = c_null_ptr; sigma_ptr = c_null_ptr
      if (present(gradient)) gradient_ptr = c_loc(gradient)
      if (present(sigma)) sigma_ptr = c_loc(sigma)
      call c_eval_gcp(param%error, param%structure, param%handle, c_loc(energy), gradient_ptr, sigma_ptr)
      if (gcp_failed(param%error, error)) return
   end subroutine

   subroutine get_geometric_counterpoise_hessian(mol, param, cutoff, hessian, partition, error)
      type(structure_type), intent(in) :: mol
      type(gcp_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(inout), target, contiguous :: hessian(:, :)
      type(work_partition), intent(in), optional :: partition
      type(error_type), allocatable, intent(out), optional :: error
      real(wp), target :: energy
      if (any(shape(hessian) /= [3*mol%nat, 3*mol%nat])) then
         call gcp_error(error, 'invalid gCP Hessian shape')
         return
      end if
      if (.not. prepare_gcp(mol, param, cutoff, partition, error)) return
      call c_hessian_gcp(param%error, param%structure, param%handle, c_loc(energy), c_loc(hessian))
      if (gcp_failed(param%error, error)) return
   end subroutine

   subroutine delete_gcp_param(self)
      type(gcp_param), intent(inout) :: self
      if (c_associated(self%handle)) call c_delete_gcp(self%handle)
      if (c_associated(self%structure)) call c_delete_structure(self%structure)
      if (c_associated(self%error)) call c_delete_error(self%error)
      self%nat = 0
   end subroutine

   subroutine close_d3_model(self)
      class(d3_model), intent(inout) :: self
      call c_delete_model(self%handle)
      call c_delete_structure(self%structure)
      call c_delete_error(self%error)
      self%nat = 0
   end subroutine

   subroutine close_damping_param(self)
      class(damping_param), intent(inout) :: self
      call c_delete_param(self%handle)
   end subroutine

   subroutine new_d3_model(self, mol, lowrank, ghost, error, d3s, model_kind)
      type(d3_model), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      type(d3_lowrank_config), intent(in), optional :: lowrank
      integer, intent(in), optional :: ghost(:)
      type(error_type), allocatable, intent(out), optional :: error
      logical, intent(in), optional :: d3s
      integer, intent(in), optional :: model_kind
      integer(c_int), allocatable, target :: numbers(:)
      integer(c_int), allocatable, target :: cghost(:)
      real(c_double), allocatable, target :: coordinates(:, :), lattice(:, :)
      logical(c_bool), target :: periodic(3)
      type(c_ptr) :: lattice_ptr
      character(kind=c_char, len=512) :: message
      integer :: atom
      integer(c_int) :: selected_model
      call self%close()
      selected_model = 0_c_int
      if (present(d3s)) selected_model = merge(1_c_int, 0_c_int, d3s)
      if (present(model_kind)) then
         if (present(d3s)) then
            if (model_kind /= selected_model) then
               if (present(error)) then
                  call fatal_error(error, 'model_kind and d3s select different D3 models')
                  return
               else
                  error stop 'model_kind and d3s select different D3 models'
               end if
            end if
         end if
         selected_model = int(model_kind, c_int)
      end if
      allocate(numbers(mol%nat))
      coordinates = mol%xyz
      do atom = 1, mol%nat; numbers(atom) = mol%num(mol%id(atom)); end do
      periodic = .false._c_bool
      if (size(mol%periodic) == 1) periodic = logical(mol%periodic(1), c_bool)
      if (size(mol%periodic) == 3) periodic = logical(mol%periodic, c_bool)
      self%error = c_new_error()
      lattice_ptr = c_null_ptr
      if (size(mol%lattice) == 9) then
         lattice = mol%lattice
         lattice_ptr = c_loc(lattice)
      end if
      self%structure = c_new_structure(self%error, mol%nat, c_loc(numbers), c_loc(coordinates), &
         lattice_ptr, c_loc(periodic))
      setup: block
         if (c_check(self%error) /= 0) exit setup
         self%handle = c_new_model(self%error, self%structure, selected_model)
         if (c_check(self%error) /= 0) exit setup
         if (present(ghost)) then
            if (size(ghost) > 0) then
               cghost = int(ghost - 1, c_int)
               call c_set_ghost(self%error, self%handle, c_loc(cghost), int(size(cghost), c_int))
               if (c_check(self%error) /= 0) exit setup
            end if
         end if
         if (present(lowrank)) call c_set_ewald(self%error, self%handle, lowrank%rank, &
            lowrank%tolerance, lowrank%kcut, lowrank%mesh)
      end block setup
      if (c_check(self%error) /= 0) then
         call c_get_error(self%error, message, len(message, kind=c_int))
         call self%close()
         if (present(error)) then
            call fatal_error(error, message(:index(message, c_null_char)-1))
         else
            error stop message(:index(message, c_null_char)-1)
         end if
         return
      end if
      self%nat = mol%nat
   end subroutine

   subroutine finish_damping_construction(self, handle, error)
      class(damping_param), intent(inout) :: self
      type(c_ptr), intent(inout) :: handle
      type(error_type), allocatable, intent(out), optional :: error
      character(kind=c_char, len=512) :: message
      if (c_check(handle) /= 0) then
         call c_get_error(handle, message, len(message, kind=c_int))
         call self%close()
         call c_delete_error(handle)
         if (present(error)) then
            call fatal_error(error, message(:index(message, c_null_char)-1))
         else
            error stop message(:index(message, c_null_char)-1)
         end if
      else
         call c_delete_error(handle)
      end if
   end subroutine

   subroutine new_rational_damping(self, param, error)
      type(rational_damping_param), intent(inout) :: self; type(d3_param), intent(in) :: param
      type(error_type), allocatable, intent(out), optional :: error
      type(c_ptr) :: handle
      call self%close()
      handle = c_new_error()
      self%handle = c_new_param(handle, param%s6, param%s8, param%s9, param%a1, param%a2, param%alp)
      call finish_damping_construction(self, handle, error)
   end subroutine

   subroutine new_zero_damping(self, param, error)
      type(zero_damping_param), intent(inout) :: self; type(d3_param), intent(in) :: param
      type(error_type), allocatable, intent(out), optional :: error
      type(c_ptr) :: handle
      call self%close()
      handle = c_new_error()
      self%handle = c_new_zero(handle, param%s6, param%s8, param%s9, param%rs6, param%rs8, param%alp)
      call finish_damping_construction(self, handle, error)
   end subroutine

   subroutine new_mzero_damping(self, param, error)
      type(mzero_damping_param), intent(inout) :: self; type(d3_param), intent(in) :: param
      type(error_type), allocatable, intent(out), optional :: error
      type(c_ptr) :: handle
      call self%close()
      handle = c_new_error()
      self%handle = c_new_mzero(handle, param%s6, param%s8, param%s9, param%rs6, param%rs8, param%alp, param%bet)
      call finish_damping_construction(self, handle, error)
   end subroutine
   subroutine new_optimizedpower_damping(self, param, error)
      type(optimizedpower_damping_param), intent(inout) :: self; type(d3_param), intent(in) :: param
      type(error_type), allocatable, intent(out), optional :: error
      type(c_ptr) :: handle
      call self%close()
      handle = c_new_error()
      self%handle = c_new_op(handle, param%s6, param%s8, param%s9, param%a1, param%a2, param%alp, param%bet)
      call finish_damping_construction(self, handle, error)
   end subroutine
   subroutine new_cso_damping(self, param, error)
      type(cso_damping_param), intent(inout) :: self; type(d3_param), intent(in) :: param
      type(error_type), allocatable, intent(out), optional :: error
      type(c_ptr) :: handle
      call self%close()
      handle = c_new_error()
      self%handle = c_new_cso(handle, param%s6, param%s9, param%a1, param%a2, param%rs6, param%rs8, param%alp)
      call finish_damping_construction(self, handle, error)
   end subroutine
   subroutine new_z_damping(self, param, error)
      type(z_damping_param), intent(inout) :: self; type(d3_param), intent(in) :: param
      type(error_type), allocatable, intent(out), optional :: error
      type(c_ptr) :: handle
      call self%close()
      handle = c_new_error()
      self%handle = c_new_z(handle, param%s6, param%s8, param%s9, param%a1, param%alp)
      call finish_damping_construction(self, handle, error)
   end subroutine

   subroutine get_properties_legacy(mol, disp, cutoff, cn, c6, dcndr, dcndL, dc6dr, dc6dL)
      class(structure_type), intent(in) :: mol
      class(d3_model), intent(in) :: disp
      type(realspace_cutoff), intent(in), optional :: cutoff
      real(wp), intent(inout), target, contiguous, optional :: cn(:), c6(:, :)
      real(wp), intent(inout), target, contiguous, optional :: dcndr(:, :, :), dcndL(:, :, :)
      real(wp), intent(inout), target, contiguous, optional :: dc6dr(:, :, :, :), dc6dL(:, :, :, :)
      type(error_type), allocatable :: error
      call get_properties_error(error, mol, disp, cutoff, cn, c6, dcndr, dcndL, dc6dr, dc6dL)
      if (allocated(error)) error stop error%message
   end subroutine

   subroutine get_properties_error(error, mol, disp, cutoff, cn, c6, dcndr, dcndL, dc6dr, dc6dL)
      type(error_type), allocatable, intent(out) :: error
      class(structure_type), intent(in) :: mol
      class(d3_model), intent(in) :: disp
      type(realspace_cutoff), intent(in), optional :: cutoff
      real(wp), intent(inout), target, contiguous, optional :: cn(:), c6(:, :)
      real(wp), intent(inout), target, contiguous, optional :: dcndr(:, :, :), dcndL(:, :, :)
      real(wp), intent(inout), target, contiguous, optional :: dc6dr(:, :, :, :), dc6dL(:, :, :, :)
      real(c_double), allocatable, target :: coordinates(:, :), lattice(:, :)
      type(c_ptr) :: lattice_ptr, cn_ptr, c6_ptr
      type(c_ptr) :: response_ptr(4)
      if (mol%nat /= disp%nat) call fatal_error(error, 'D3 model atom count mismatch')
      if (present(cn)) then
         if (size(cn) /= mol%nat) call fatal_error(error, 'Invalid D3 CN shape')
      end if
      if (present(c6)) then
         if (any(shape(c6) /= [mol%nat, mol%nat])) call fatal_error(error, 'Invalid D3 C6 shape')
      end if
      if (present(dcndr)) then
         if (any(shape(dcndr) /= [3, mol%nat, mol%nat])) call fatal_error(error, 'Invalid D3 dcndr shape')
      end if
      if (present(dcndL)) then
         if (any(shape(dcndL) /= [3, 3, mol%nat])) call fatal_error(error, 'Invalid D3 dcndL shape')
      end if
      if (present(dc6dr)) then
         if (any(shape(dc6dr) /= [3, mol%nat, mol%nat, mol%nat])) call fatal_error(error, 'Invalid D3 dc6dr shape')
      end if
      if (present(dc6dL)) then
         if (any(shape(dc6dL) /= [3, 3, mol%nat, mol%nat])) call fatal_error(error, 'Invalid D3 dc6dL shape')
      end if
      if (allocated(error)) return
      coordinates = mol%xyz
      lattice_ptr = c_null_ptr
      if (size(mol%lattice) == 9) then
         lattice = mol%lattice
         lattice_ptr = c_loc(lattice)
      end if
      call c_update(disp%error, disp%structure, c_loc(coordinates), lattice_ptr)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 structure update failed')
         return
      end if
      if (present(cutoff)) call c_set_cutoff(disp%error, disp%handle, cutoff%disp2, &
         cutoff%disp3, cutoff%cn, cutoff%width2, cutoff%width3)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 cutoff failed')
         return
      end if
      cn_ptr = c_null_ptr
      c6_ptr = c_null_ptr
      if (present(cn)) cn_ptr = c_loc(cn)
      if (present(c6)) c6_ptr = c_loc(c6)
      if (present(dcndr) .or. present(dcndL) .or. present(dc6dr) .or. present(dc6dL)) then
         response_ptr = c_null_ptr
         if (present(dcndr)) response_ptr(1) = c_loc(dcndr)
         if (present(dcndL)) response_ptr(2) = c_loc(dcndL)
         if (present(dc6dr)) response_ptr(3) = c_loc(dc6dr)
         if (present(dc6dL)) response_ptr(4) = c_loc(dc6dL)
         call c_property_response(disp%error, disp%structure, disp%handle, cn_ptr, c6_ptr, &
            response_ptr(1), response_ptr(2), response_ptr(3), response_ptr(4))
      else
         call c_properties(disp%error, disp%structure, disp%handle, cn_ptr, c6_ptr)
      end if
      if (c_check(disp%error) /= 0) call fatal_error(error, 'disprs D3 property evaluation failed')
   end subroutine

   subroutine get_dispersion_scalar(mol, model, param, cutoff, energy, gradient, sigma, hessian, partition)
      class(structure_type), intent(in), target :: mol
      class(d3_model), intent(in) :: model
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in), optional :: cutoff
      real(wp), intent(out), target :: energy
      real(wp), intent(out), target, contiguous, optional :: gradient(:, :), sigma(:, :), hessian(:, :)
      type(work_partition), intent(in), optional :: partition
      type(error_type), allocatable :: error
      call get_dispersion_scalar_error(error, mol, model, param, cutoff, energy, gradient, sigma, hessian, partition)
      if (allocated(error)) error stop error%message
   end subroutine

   subroutine get_dispersion_scalar_error(error, mol, disp, param, cutoff, energy, gradient, sigma, hessian, partition)
      type(error_type), allocatable, intent(out) :: error
      class(structure_type), intent(in), target :: mol
      class(d3_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in), optional :: cutoff
      real(wp), intent(out), target :: energy
      real(wp), intent(out), target, contiguous, optional :: gradient(:, :), sigma(:, :), hessian(:, :)
      type(work_partition), intent(in), optional :: partition
      type(work_partition) :: selected_partition
      real(c_double), allocatable, target :: coordinates(:, :), lattice(:, :)
      type(c_ptr) :: gradient_ptr, sigma_ptr, lattice_ptr
      if (mol%nat /= disp%nat) call fatal_error(error, 'D3 model atom count mismatch')
      if (present(gradient)) then
         if (any(shape(gradient) /= [3, mol%nat])) call fatal_error(error, 'Invalid D3 gradient shape')
      end if
      if (present(sigma)) then
         if (any(shape(sigma) /= [3, 3])) call fatal_error(error, 'Invalid D3 virial shape')
      end if
      if (present(hessian)) then
         if (any(shape(hessian) /= [3*mol%nat, 3*mol%nat])) call fatal_error(error, 'Invalid D3 Hessian shape')
      end if
      if (allocated(error)) return
      coordinates = mol%xyz
      lattice_ptr = c_null_ptr
      if (size(mol%lattice) == 9) then
         lattice = mol%lattice
         lattice_ptr = c_loc(lattice)
      end if
      call c_update(disp%error, disp%structure, c_loc(coordinates), lattice_ptr)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 structure update failed')
         return
      end if
      selected_partition = serial_work_partition
      if (present(partition)) selected_partition = partition
      call c_set_partition(disp%error, disp%handle, selected_partition%part, selected_partition%nparts)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 partition failed')
         return
      end if
      if (present(cutoff)) call c_set_cutoff(disp%error, disp%handle, cutoff%disp2, &
         cutoff%disp3, cutoff%cn, cutoff%width2, cutoff%width3)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 cutoff failed')
         return
      end if
      if (present(hessian)) then
         call c_hessian(disp%error, disp%structure, disp%handle, param%handle, c_loc(energy), c_loc(hessian))
         if (c_check(disp%error) /= 0) then
            call fatal_error(error, 'disprs D3 Hessian failed')
            return
         end if
      end if
      if (.not.present(hessian) .or. present(gradient) .or. present(sigma)) then
         gradient_ptr = c_null_ptr
         sigma_ptr = c_null_ptr
         if (present(gradient)) gradient_ptr = c_loc(gradient)
         if (present(sigma)) sigma_ptr = c_loc(sigma)
         call c_eval(disp%error, disp%structure, disp%handle, param%handle, c_loc(energy), &
            gradient_ptr, sigma_ptr)
         if (c_check(disp%error) /= 0) call fatal_error(error, 'disprs D3 evaluation failed')
      end if
   end subroutine

   subroutine get_dispersion_atomic(mol, disp, param, cutoff, energies, gradient, sigma, hessian, partition)
      class(structure_type), intent(in) :: mol
      class(d3_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out) :: energies(:)
      real(wp), intent(out), contiguous, optional :: gradient(:, :), sigma(:, :), hessian(:, :)
      type(work_partition), intent(in), optional :: partition
      type(error_type), allocatable :: error
      call get_dispersion_atomic_error(error, mol, disp, param, cutoff, energies, gradient, sigma, hessian, partition)
      if (allocated(error)) error stop error%message
   end subroutine

   subroutine get_dispersion_atomic_error(error, mol, disp, param, cutoff, energies, gradient, sigma, hessian, partition)
      type(error_type), allocatable, intent(out) :: error
      class(structure_type), intent(in) :: mol
      class(d3_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out) :: energies(:)
      real(wp), intent(out), contiguous, optional :: gradient(:, :), sigma(:, :), hessian(:, :)
      type(work_partition), intent(in), optional :: partition
      real(wp) :: energy
      real(wp), allocatable :: energy2(:, :), energy3(:, :)
      if (size(energies) /= mol%nat) then
         call fatal_error(error, 'Invalid D3 atomic energy shape')
         return
      end if
      if (present(gradient) .or. present(sigma) .or. present(hessian)) then
         call get_dispersion_scalar_error(error, mol, disp, param, cutoff, energy, gradient, sigma, hessian, partition)
         if (allocated(error)) return
      end if
      ! ponytail: O(N**2) pair storage; expose native atomic vectors if memory becomes limiting.
      allocate(energy2(mol%nat, mol%nat), energy3(mol%nat, mol%nat))
      call get_pairwise_dispersion_error(error, mol, disp, param, cutoff, energy2, energy3, partition)
      if (allocated(error)) return
      energies = sum(energy2, dim=1) + sum(energy3, dim=1)
   end subroutine

   subroutine get_pairwise_dispersion_legacy(mol, model, param, cutoff, pair2, pair3, partition)
      type(structure_type), intent(in) :: mol
      type(d3_model), intent(inout) :: model
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out), target, contiguous :: pair2(:, :), pair3(:, :)
      type(work_partition), intent(in), optional :: partition
      type(error_type), allocatable :: error
      call get_pairwise_dispersion_error(error, mol, model, param, cutoff, pair2, pair3, partition)
      if (allocated(error)) error stop error%message
   end subroutine

   subroutine get_pairwise_dispersion_error(error, mol, disp, param, cutoff, energy2, energy3, partition)
      type(error_type), allocatable, intent(out) :: error
      class(structure_type), intent(in) :: mol
      class(d3_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out), target, contiguous :: energy2(:, :), energy3(:, :)
      type(work_partition), intent(in), optional :: partition
      type(work_partition) :: selected_partition
      real(c_double), allocatable, target :: coordinates(:, :), lattice(:, :)
      type(c_ptr) :: lattice_ptr
      if (any(shape(energy2) /= [mol%nat, mol%nat]) .or. &
          any(shape(energy3) /= [mol%nat, mol%nat])) then
         call fatal_error(error, 'Invalid D3 pairwise shape')
         return
      end if
      if (mol%nat /= disp%nat) then
         call fatal_error(error, 'D3 model atom count mismatch')
         return
      end if
      coordinates = mol%xyz
      lattice_ptr = c_null_ptr
      if (size(mol%lattice) == 9) then
         lattice = mol%lattice
         lattice_ptr = c_loc(lattice)
      end if
      call c_update(disp%error, disp%structure, c_loc(coordinates), lattice_ptr)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 structure update failed')
         return
      end if
      selected_partition = serial_work_partition
      if (present(partition)) selected_partition = partition
      call c_set_partition(disp%error, disp%handle, selected_partition%part, selected_partition%nparts)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 partition failed')
         return
      end if
      call c_set_cutoff(disp%error, disp%handle, cutoff%disp2, cutoff%disp3, &
         cutoff%cn, cutoff%width2, cutoff%width3)
      if (c_check(disp%error) /= 0) then
         call fatal_error(error, 'disprs D3 cutoff failed')
         return
      end if
      call c_pairwise(disp%error, disp%structure, disp%handle, param%handle, c_loc(energy2), c_loc(energy3))
      if (c_check(disp%error) /= 0) call fatal_error(error, 'disprs D3 pairwise evaluation failed')
   end subroutine

   subroutine new_work_partition(error, partition, part, nparts)
      type(error_type), allocatable, intent(out) :: error
      type(work_partition), intent(out) :: partition
      integer, intent(in) :: part, nparts
      if (nparts <= 0 .or. part < 0 .or. part >= nparts) then
         call fatal_error(error, 'Invalid dispersion work partition')
      else
         partition = work_partition(part, nparts)
      end if
   end subroutine

   subroutine get_dftd3_version(major, minor, patch, string)
      integer, optional, intent(out) :: major, minor, patch
      character(len=:), allocatable, optional, intent(out) :: string
      if (present(major)) major = 1
      if (present(minor)) minor = 6
      if (present(patch)) patch = 0
      if (present(string)) string = '1.6.0'
   end subroutine
end module