module dftd4
   use, intrinsic :: iso_c_binding
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, new
   use mctc_cutoff, only : get_lattice_points
   use disprs, only : disprs_d4_energy, d4_lowrank_config
   implicit none
   private
   integer, parameter :: damping_parameter_count = 6
   public :: structure_type, new, realspace_cutoff, dispersion_model, d4_model, d4s_model
   public :: get_lattice_points
   public :: damping_param, rational_damping_param, d4_qmod
   public :: new_dispersion_model, new_d4_model, new_d4s_model, get_rational_damping
   public :: set_fixed_charges
   public :: d4_lowrank_config
   public :: get_dispersion, get_properties, get_pairwise_dispersion, get_dispersion_hessian, get_dftd4_version

   type :: realspace_cutoff
      sequence
      real(wp) :: cn = 30.0_wp, disp2 = 60.0_wp, disp3 = 40.0_wp
      real(wp) :: width2 = 0.0_wp, width3 = 0.0_wp
   end type
   type :: enum_qmod
      integer :: eeq = 1, gfn2 = 2, eeqbc = 3
   end type
   type(enum_qmod), parameter :: d4_qmod = enum_qmod()
   type, abstract :: dispersion_model
      real(wp) :: ga = 3.0_wp, gc = 2.0_wp
      integer, private :: kind = 0, qmod = 1
      integer, private :: nat = 0
      real(wp), allocatable, private :: fixed_charges(:)
      type(d4_lowrank_config), allocatable, private :: lowrank
   end type
   type, extends(dispersion_model) :: d4_model
      real(wp) :: wf = 6.0_wp
   end type
   type, extends(dispersion_model) :: d4s_model
   end type
   type, abstract :: damping_param
   end type
   type, extends(damping_param) :: rational_damping_param
      real(wp) :: s6 = 1.0_wp, s8 = 0.0_wp, s9 = 1.0_wp
      real(wp) :: a1 = 0.0_wp, a2 = 0.0_wp, alp = 16.0_wp
   end type
   interface get_rational_damping
      module procedure get_rational_damping_name, get_rational_damping_id
   end interface
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
      subroutine c_delete_error(error) bind(C, name='disprs_d4_delete_error')
         import c_ptr
         type(c_ptr) :: error
      end subroutine
      subroutine c_named(error, method, s9, values) bind(C, name='disprs_d4_get_named_parameters_s9')
         import c_ptr, c_char, c_double, damping_parameter_count
         type(c_ptr), value :: error, s9
         character(c_char), intent(in) :: method(*)
         real(c_double), intent(out) :: values(damping_parameter_count)
      end subroutine
   end interface
contains
   subroutine new_d4_model(error, d4, mol, ga, gc, wf, qmod, fixed_charges, lowrank)
      type(error_type), allocatable, intent(out) :: error
      type(d4_model), intent(out) :: d4
      class(structure_type), intent(in) :: mol
      real(wp), optional, intent(in) :: ga, gc, wf
      real(wp), optional, intent(in) :: fixed_charges(:)
      integer, optional, intent(in) :: qmod
      type(d4_lowrank_config), optional, intent(in) :: lowrank
      if (present(wf)) d4%wf = wf
      call configure(error, d4, mol, ga, gc, qmod, fixed_charges, lowrank)
      if (allocated(error)) return
      if (.not.ieee_is_finite(d4%wf) .or. d4%wf <= 0) call fatal_error(error, 'Invalid D4 weighting factor')
   end subroutine

   subroutine new_d4s_model(error, d4, mol, ga, gc, qmod, fixed_charges, lowrank)
      type(error_type), allocatable, intent(out) :: error
      type(d4s_model), intent(out) :: d4
      class(structure_type), intent(in) :: mol
      real(wp), optional, intent(in) :: ga, gc
      real(wp), optional, intent(in) :: fixed_charges(:)
      integer, optional, intent(in) :: qmod
      type(d4_lowrank_config), optional, intent(in) :: lowrank
      d4%kind = 1
      call configure(error, d4, mol, ga, gc, qmod, fixed_charges, lowrank)
   end subroutine

   subroutine configure(error, self, mol, ga, gc, qmod, fixed_charges, lowrank)
      type(error_type), allocatable, intent(out) :: error
      class(dispersion_model), intent(inout) :: self
      class(structure_type), intent(in) :: mol
      real(wp), optional, intent(in) :: ga, gc
      real(wp), optional, intent(in) :: fixed_charges(:)
      integer, optional, intent(in) :: qmod
      type(d4_lowrank_config), optional, intent(in) :: lowrank
      if (present(ga)) self%ga = ga
      if (present(gc)) self%gc = gc
      if (present(qmod)) self%qmod = qmod
      if (mol%nat <= 0 .or. any(mol%num < 1) .or. any(mol%num > 118) &
          .or. any(mol%num >= 104 .and. mol%num <= 111)) then
         call fatal_error(error, 'Unsupported D4 structure')
      else if (self%qmod == d4_qmod%eeqbc .and. any(mol%num > 103)) then
         call fatal_error(error, 'EEQBC supports atomic numbers 1 through 103')
      else if (.not.all(ieee_is_finite([self%ga, self%gc])) .or. min(self%ga, self%gc) <= 0) then
         call fatal_error(error, 'Invalid D4 charge scaling')
      else if (self%qmod /= d4_qmod%eeq .and. self%qmod /= d4_qmod%eeqbc) then
         call fatal_error(error, 'Only EEQ and EEQBC reference charge models are supported')
      end if
      if (allocated(error)) return
      self%nat = mol%nat
      if (present(lowrank)) self%lowrank = lowrank
      call set_fixed_charges(error, self, fixed_charges)
   end subroutine

   subroutine set_fixed_charges(error, disp, charges)
      type(error_type), allocatable, intent(out) :: error
      class(dispersion_model), intent(inout) :: disp
      real(wp), optional, intent(in) :: charges(:)
      if (.not.present(charges)) then
         if (allocated(disp%fixed_charges)) deallocate(disp%fixed_charges)
         return
      end if
      if (disp%nat <= 0 .or. size(charges) /= disp%nat) then
         call fatal_error(error, 'Fixed D4 charges require one value per atom')
      else if (.not.all(ieee_is_finite(charges))) then
         call fatal_error(error, 'Fixed D4 charges must be finite')
      else
         disp%fixed_charges = charges
      end if
   end subroutine

   subroutine new_dispersion_model(error, d4, mol, model, ga, gc, wf, qmod, fixed_charges, lowrank)
      type(error_type), allocatable, intent(out) :: error
      class(dispersion_model), allocatable, intent(out) :: d4
      type(structure_type), intent(in) :: mol
      character(len=*), optional, intent(in) :: model
      real(wp), optional, intent(in) :: ga, gc, wf
      real(wp), optional, intent(in) :: fixed_charges(:)
      integer, optional, intent(in) :: qmod
      type(d4_lowrank_config), optional, intent(in) :: lowrank
      character(len=:), allocatable :: selected
      integer :: index, code
      selected = 'd4'
      if (present(model)) selected = trim(model)
      do index = 1, len(selected)
         code = iachar(selected(index:index))
         if (code >= iachar('A') .and. code <= iachar('Z')) selected(index:index) = achar(code + 32)
      end do
      select case(selected)
      case('d4')
         allocate(d4_model :: d4)
         select type(d4)
         type is(d4_model)
            call new_d4_model(error, d4, mol, ga, gc, wf, qmod, fixed_charges, lowrank)
         end select
      case('d4s')
         allocate(d4s_model :: d4)
         select type(d4)
         type is(d4s_model)
            call new_d4s_model(error, d4, mol, ga, gc, qmod, fixed_charges, lowrank)
         end select
      case default
         call fatal_error(error, 'Unknown D4 dispersion model')
      end select
   end subroutine

   subroutine get_rational_damping_name(functional, param, s9)
      character(len=*), intent(in) :: functional
      class(damping_param), allocatable, intent(out) :: param
      real(wp), optional, target, intent(in) :: s9
      type(c_ptr) :: error, s9_ptr
      real(c_double) :: values(damping_parameter_count)
      integer :: slash
      slash = index(functional, '/')
      if (slash == 0) slash = len_trim(functional) + 1
      error = c_new_error()
      s9_ptr = c_null_ptr
      if (present(s9)) s9_ptr = c_loc(s9)
      call c_named(error, functional(:slash-1)//c_null_char, s9_ptr, values)
      if (c_check(error) == 0) then
         allocate(param, source=rational_damping_param(values(1), values(2), values(3), &
            values(4), values(5), values(6)))
      end if
      call c_delete_error(error)
   end subroutine

   subroutine get_rational_damping_id(id, param, s9)
      integer, intent(in) :: id
      class(damping_param), allocatable, intent(out) :: param
      real(wp), optional, intent(in) :: s9
      character(len=20), parameter :: names(*) = [character(len=20) :: &
         'hf', 'blyp', 'bpbe', 'bp', 'bpw', 'lb94', 'mpwlyp', 'mpwpw', &
         'olyp', 'opbe', 'pbe', 'rpbe', 'revpbe', 'pw86pbe', 'rpw86pbe', 'pw91', &
         'pwp', 'xlyp', 'b97', 'tpss', 'revtpss', 'scan', 'rscan', 'r2scan', &
         'b1lyp', 'b3lyp', 'bhlyp', 'b1p', 'b3p', 'b1pw', 'b3pw', 'o3lyp', &
         'revpbe0', 'revpbe38', 'pbe0', 'pwp1', 'pw1pw', 'mpw1pw', 'mpw1lyp', 'pw6b95', &
         'tpssh', 'tpss0', 'x3lyp', 'm06l', 'm06', 'b97d', 'wb97', 'wb97x_2008', &
         'b97m', 'wb97m', 'camb3lyp', 'lcblyp', 'lh07tsvwn', 'lh07ssvwn', &
         'lh12ctssirpw92', 'lh12ctssifpw92', 'lh14tcalpbe', 'lh20t', 'b2plyp', 'b2gpplyp', &
         'mpw2plyp', 'pwpb95', 'dsdblyp', 'dsdpbe', 'dsdpbeb95', 'dsdpbep86', 'dsdsvwn', &
         'dodblyp', 'dodpbe', 'dodpbeb95', 'dodpbep86', 'dodsvwn', 'pbe0_2', 'pbe0_dh', &
         'hsesol', 'dftb_3ob', 'dftb_mio', 'dftb_ob2', 'dftb_matsci', 'dftb_pbc', &
         'b1b95', 'pbesol', 'hse06', 'mpwb1k', 'hse03', 'revtpssh', 'mn12sx', 'glyp', &
         'mpw1b95', 'revpbe0dh', 'revtpss0', 'revdsdpbep86', 'revdsdpbe', 'revdsdblyp', &
         'revdodpbep86', 'am05', 'hse12', 'hse12s', 'r2scanh', 'r2scan0', 'r2scan50', &
         'r2scan_3c', 'camqtp01', 'lcwpbe', 'lcwpbeh', 'wb97x_rev', 'wb97m_rev', &
         'wb97x_3c', 'wr2scan', 'r2scan0_dh', 'r2scan_cidh', 'r2scan_qidh', 'r2scan0_2', &
         'pr2scan50', 'pr2scan69', 'kpr2scan50', 'wpr2scan50', 'wb97x']
      if (id < 1 .or. id > size(names)) return
      call get_rational_damping_name(trim(names(id)), param, s9)
   end subroutine

   subroutine evaluate(mol, disp, param, cutoff, energy, gradient, sigma, cn, q, c6, alpha, energy2, energy3, hessian, &
      dcndr, dcndL, dqdr, dqdL, dc6dr, dc6dL, dalphadr, dalphadL)
      class(structure_type), intent(in) :: mol
      class(dispersion_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out) :: energy
      real(wp), optional, intent(out) :: gradient(:, :), sigma(:, :), cn(:), q(:), c6(:, :), alpha(:)
      real(wp), optional, intent(out) :: energy2(:, :), energy3(:, :)
      real(wp), optional, intent(out) :: hessian(:, :, :, :)
      real(wp), optional, intent(inout) :: dcndr(:, :, :), dcndL(:, :, :), dqdr(:, :, :), dqdL(:, :, :), &
         dc6dr(:, :, :, :), dc6dL(:, :, :, :), dalphadr(:, :, :), dalphadL(:, :, :)
      integer(c_int) :: status
      integer(c_int) :: numbers(mol%nat)
      real(c_double) :: grad(3, mol%nat), virial(3, 3), damping(damping_parameter_count), wf
      logical(c_bool) :: periodic(3)
      real(wp) :: lattice(3, 3)
      select type(param)
      type is(rational_damping_param)
         damping = [param%s6, param%s8, param%s9, param%a1, param%a2, param%alp]
      class default
         error stop 'Unsupported D4 damping type'
      end select
      wf = 6.0_wp
      select type(disp)
      type is(d4_model)
         wf = disp%wf
      end select
      numbers = mol%num(mol%id)
      periodic = .false.
      if (size(mol%periodic) == 1) periodic = mol%periodic(1)
      if (size(mol%periodic) == 3) periodic = mol%periodic
      lattice = 0.0_wp
      if (size(mol%lattice) == 9) lattice = mol%lattice
      call disprs_d4_energy(numbers, mol%xyz, '', mol%charge, .true._c_bool, energy, grad, status, &
         lattice=lattice, periodic=periodic, virial=virial, model_kind=int(disp%kind, c_int), damping=damping, &
         coordination=cn, charges=q, c6=c6, polarizabilities=alpha, pair2=energy2, pair3=energy3, &
         ga=disp%ga, gc=disp%gc, wf=wf, &
         cutoff=[cutoff%cn, cutoff%disp2, cutoff%disp3, cutoff%width2, cutoff%width3], &
         charge_model=int(merge(1, 0, disp%qmod == d4_qmod%eeqbc), c_int), hessian=hessian, &
         dcndr=dcndr, dcndL=dcndL, dqdr=dqdr, dqdL=dqdL, dc6dr=dc6dr, dc6dL=dc6dL, &
         dalphadr=dalphadr, dalphadL=dalphadL, fixed_charges=disp%fixed_charges, lowrank=disp%lowrank)
      if (status /= 0) error stop 'Native D4 evaluation failed'
      if (present(gradient)) gradient = grad
      if (present(sigma)) sigma = virial
   end subroutine

   subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma)
      class(structure_type), intent(in) :: mol
      class(dispersion_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out) :: energy
      real(wp), optional, contiguous, intent(out) :: gradient(:, :), sigma(:, :)
      call evaluate(mol, disp, param, cutoff, energy, gradient, sigma)
   end subroutine

   subroutine get_properties(mol, disp, cutoff, cn, q, c6, alpha, dcndr, dcndL, dqdr, dqdL, dc6dr, dc6dL, dalphadr, dalphadL)
      class(structure_type), intent(in) :: mol
      class(dispersion_model), intent(in) :: disp
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out) :: cn(:), c6(:, :), alpha(:)
      real(wp), contiguous, intent(out) :: q(:)
      real(wp), optional, intent(inout) :: dcndr(:, :, :), dcndL(:, :, :), dqdr(:, :, :), dqdL(:, :, :), &
         dc6dr(:, :, :, :), dc6dL(:, :, :, :), dalphadr(:, :, :), dalphadL(:, :, :)
      real(wp) :: energy
      call evaluate(mol, disp, rational_damping_param(s6=0.0_wp, s8=0.0_wp, s9=0.0_wp, &
         a1=0.4_wp, a2=5.0_wp), cutoff, energy, cn=cn, q=q, c6=c6, alpha=alpha, &
         dcndr=dcndr, dcndL=dcndL, dqdr=dqdr, dqdL=dqdL, dc6dr=dc6dr, dc6dL=dc6dL, &
         dalphadr=dalphadr, dalphadL=dalphadL)
   end subroutine

   subroutine get_pairwise_dispersion(mol, disp, param, cutoff, energy2, energy3)
      class(structure_type), intent(in) :: mol
      class(dispersion_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out) :: energy2(:, :), energy3(:, :)
      real(wp) :: energy
      call evaluate(mol, disp, param, cutoff, energy, energy2=energy2, energy3=energy3)
   end subroutine

   subroutine get_dispersion_hessian(mol, disp, param, cutoff, hessian)
      class(structure_type), intent(in) :: mol
      class(dispersion_model), intent(in) :: disp
      class(damping_param), intent(in) :: param
      type(realspace_cutoff), intent(in) :: cutoff
      real(wp), intent(out) :: hessian(:, :, :, :)
      real(wp) :: energy
      call evaluate(mol, disp, param, cutoff, energy, hessian=hessian)
   end subroutine

   subroutine get_dftd4_version(major, minor, patch, string)
      integer, optional, intent(out) :: major, minor, patch
      character(len=:), allocatable, optional, intent(out) :: string
      if (present(major)) major = 4
      if (present(minor)) minor = 2
      if (present(patch)) patch = 0
      if (present(string)) string = '4.2.0'
   end subroutine
end module dftd4