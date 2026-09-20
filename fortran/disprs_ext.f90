module disprs
   use, intrinsic :: iso_c_binding
   implicit none
   private
   public :: disprs_d3_energy, disprs_d4_energy
   integer, parameter :: d4_parameter_count = 6, cutoff_count = 5, response_count = 8

   type, public :: d4_lowrank_config
      integer(c_int) :: rank = 0, mesh = 0
      real(c_double) :: tolerance = 1.0e-4_c_double, kcut = 0.0_c_double
   end type

   interface
      function c_d3_new_error() bind(C, name="disprs_d3_new_error") result(handle)
         import c_ptr
         type(c_ptr) :: handle
      end function
      function c_d3_check_error(error) bind(C, name="disprs_d3_check_error") result(status)
         import c_int, c_ptr
         type(c_ptr), value :: error
         integer(c_int) :: status
      end function
      subroutine c_d3_delete_error(error) bind(C, name="disprs_d3_delete_error")
         import c_ptr
         type(c_ptr) :: error
      end subroutine
      function c_d3_new_structure(error, natoms, numbers, positions, lattice, periodic) &
            bind(C, name="disprs_d3_new_structure") result(handle)
         import c_int, c_ptr
         type(c_ptr), value :: error, numbers, positions, lattice, periodic
         integer(c_int), value :: natoms
         type(c_ptr) :: handle
      end function
      subroutine c_d3_delete_structure(structure) bind(C, name="disprs_d3_delete_structure")
         import c_ptr
         type(c_ptr) :: structure
      end subroutine
      function c_d3_new_model(error, structure, kind) bind(C, name="disprs_d3_new_model_kind") result(handle)
         import c_ptr, c_int
         type(c_ptr), value :: error, structure
         integer(c_int), value :: kind
         type(c_ptr) :: handle
      end function
      subroutine c_d3_delete_model(model) bind(C, name="disprs_d3_delete_model")
         import c_ptr
         type(c_ptr) :: model
      end subroutine
      function c_d3_load_param(error, damping, method, atm) &
            bind(C, name="disprs_d3_load_param") result(handle)
         import c_bool, c_char, c_int, c_ptr
         type(c_ptr), value :: error
         integer(c_int), value :: damping
         character(kind=c_char), intent(in) :: method(*)
         logical(c_bool), value :: atm
         type(c_ptr) :: handle
      end function
      subroutine c_d3_delete_param(param) bind(C, name="disprs_d3_delete_param")
         import c_ptr
         type(c_ptr) :: param
      end subroutine
      subroutine c_d3_get_dispersion(error, structure, model, param, energy, gradient, virial) &
            bind(C, name="disprs_d3_get_dispersion")
         import c_ptr
         type(c_ptr), value :: error, structure, model, param, energy, gradient, virial
      end subroutine

      function c_d4_new_error() bind(C, name="disprs_d4_new_error") result(handle)
         import c_ptr
         type(c_ptr) :: handle
      end function
      function c_d4_check_error(error) bind(C, name="disprs_d4_check_error") result(status)
         import c_int, c_ptr
         type(c_ptr), value :: error
         integer(c_int) :: status
      end function
      subroutine c_d4_delete_error(error) bind(C, name="disprs_d4_delete_error")
         import c_ptr
         type(c_ptr) :: error
      end subroutine
      function c_d4_new_structure(error, natoms, numbers, positions, charge, lattice, periodic) &
            bind(C, name="disprs_d4_new_structure") result(handle)
         import c_int, c_ptr
         type(c_ptr), value :: error, numbers, positions, charge, lattice, periodic
         integer(c_int), value :: natoms
         type(c_ptr) :: handle
      end function
      subroutine c_d4_delete_structure(structure) bind(C, name="disprs_d4_delete_structure")
         import c_ptr
         type(c_ptr) :: structure
      end subroutine
      function c_d4_new_model(error, structure, model, ga, gc, wf) &
            bind(C, name="disprs_d4_new_custom_model") result(handle)
         import c_int, c_ptr, c_double
         type(c_ptr), value :: error, structure
         integer(c_int), value :: model
         real(c_double), value :: ga, gc, wf
         type(c_ptr) :: handle
      end function
      subroutine c_d4_set_cutoff(error, model, cn, disp2, disp3, width2, width3) &
            bind(C, name="disprs_d4_set_realspace_cutoff")
         import c_ptr, c_double
         type(c_ptr), value :: error, model
         real(c_double), value :: cn, disp2, disp3, width2, width3
      end subroutine
      subroutine c_d4_set_charge_model(error, model, charge_model) bind(C, name="disprs_d4_set_charge_model")
         import c_ptr, c_int
         type(c_ptr), value :: error, model
         integer(c_int), value :: charge_model
      end subroutine
      subroutine c_d4_set_ewald(error, model, rank, tolerance, kcut, mesh) bind(C, name="disprs_d4_set_model_ewald")
         import c_ptr, c_int, c_double
         type(c_ptr), value :: error, model
         integer(c_int), value :: rank, mesh
         real(c_double), value :: tolerance, kcut
      end subroutine
      subroutine c_d4_set_fixed_charges(error, model, charges, count) bind(C, name="disprs_d4_set_fixed_charges")
         import c_ptr, c_int
         type(c_ptr), value :: error, model, charges
         integer(c_int), value :: count
      end subroutine
      subroutine c_d4_set_charge_cutoff(error, model, cutoff) bind(C, name="disprs_d4_set_charge_cutoff")
         import c_ptr, c_double
         type(c_ptr), value :: error, model
         real(c_double), value :: cutoff
      end subroutine
      subroutine c_d4_delete_model(model) bind(C, name="disprs_d4_delete_model")
         import c_ptr
         type(c_ptr) :: model
      end subroutine
      subroutine c_d4_set_ghost(error, model, ghost, count) bind(C, name="disprs_d4_set_model_ghost_index")
         import c_ptr, c_int
         type(c_ptr), value :: error, model, ghost
         integer(c_int), value :: count
      end subroutine
      subroutine c_d4_set_partition(error, model, part, parts) bind(C, name="disprs_d4_set_model_work_partition")
         import c_ptr, c_int
         type(c_ptr), value :: error, model
         integer(c_int), value :: part, parts
      end subroutine
      function c_d4_load_param(error, method, atm) &
            bind(C, name="disprs_d4_load_param") result(handle)
         import c_bool, c_char, c_ptr
         type(c_ptr), value :: error
         character(kind=c_char), intent(in) :: method(*)
         logical(c_bool), value :: atm
         type(c_ptr) :: handle
      end function
      subroutine c_d4_delete_param(param) bind(C, name="disprs_d4_delete_param")
         import c_ptr
         type(c_ptr) :: param
      end subroutine
      function c_d4_new_rational_damping(error, s6, s8, s9, a1, a2, alp) &
         bind(C, name="disprs_d4_new_rational_damping") result(handle)
         import c_ptr, c_double
         type(c_ptr), value :: error
         real(c_double), value :: s6, s8, s9, a1, a2, alp
         type(c_ptr) :: handle
      end function
      subroutine c_d4_get_dispersion(error, structure, model, param, energy, gradient, virial) &
            bind(C, name="disprs_d4_get_dispersion")
         import c_ptr
         type(c_ptr), value :: error, structure, model, param, energy, gradient, virial
      end subroutine
      subroutine c_d4_get_properties(error, structure, model, cn, charges, c6, alpha) &
            bind(C, name="disprs_d4_get_properties")
         import c_ptr
         type(c_ptr), value :: error, structure, model, cn, charges, c6, alpha
      end subroutine
      subroutine c_d4_property_response(error, structure, model, cn, charges, c6, alpha, &
            dcndr, dcndL, dqdr, dqdL, dc6dr, dc6dL, dalphadr, dalphadL) &
            bind(C, name="disprs_d4_get_property_response")
         import c_ptr
         type(c_ptr), value :: error, structure, model, cn, charges, c6, alpha
         type(c_ptr), value :: dcndr, dcndL, dqdr, dqdL, dc6dr, dc6dL, dalphadr, dalphadL
      end subroutine
      subroutine c_d4_get_pairwise(error, structure, model, param, pair2, pair3) &
            bind(C, name="disprs_d4_get_pairwise_dispersion")
         import c_ptr
         type(c_ptr), value :: error, structure, model, param, pair2, pair3
      end subroutine
      subroutine c_d4_hessian(error, structure, model, param, hessian) &
         bind(C, name="disprs_d4_get_dispersion_hessian")
         import c_ptr
         type(c_ptr), value :: error, structure, model, param, hessian
      end subroutine
   end interface

contains

   subroutine disprs_d3_energy(numbers, positions, method, damping, atm, energy, gradient, status, d3s, model_kind)
      integer(c_int), intent(in), target, contiguous :: numbers(:)
      real(c_double), intent(in), target, contiguous :: positions(:, :)
      character(len=*), intent(in) :: method
      integer(c_int), intent(in) :: damping
      logical(c_bool), intent(in) :: atm
      real(c_double), intent(out), target :: energy
      real(c_double), intent(out), target, contiguous :: gradient(:, :)
      integer(c_int), intent(out) :: status
      logical(c_bool), intent(in), optional :: d3s
      integer(c_int), intent(in), optional :: model_kind
      integer(c_int) :: selected_model
      character(kind=c_char, len=:), allocatable :: c_method
      type(c_ptr) :: error, structure, model, param

      if (size(positions, 1) /= 3 .or. size(positions, 2) /= size(numbers) .or. &
            any(shape(gradient) /= shape(positions))) error stop "invalid disprs D3 array shape"
      selected_model = 0_c_int
      if (present(d3s)) selected_model = merge(1_c_int, 0_c_int, d3s)
      if (present(model_kind)) then
         status = 1_c_int
         if (model_kind < 0 .or. model_kind > 1) return
         if (present(d3s)) then
            if (model_kind /= selected_model) return
         end if
         selected_model = model_kind
      end if
      c_method = trim(method) // c_null_char
      error = c_d3_new_error()
      structure = c_d3_new_structure(error, size(numbers), c_loc(numbers), c_loc(positions), &
         c_null_ptr, c_null_ptr)
      model = c_d3_new_model(error, structure, selected_model)
      param = c_d3_load_param(error, damping, c_method, atm)
      call c_d3_get_dispersion(error, structure, model, param, c_loc(energy), &
         c_loc(gradient), c_null_ptr)
      status = c_d3_check_error(error)
      call c_d3_delete_param(param)
      call c_d3_delete_model(model)
      call c_d3_delete_structure(structure)
      call c_d3_delete_error(error)
   end subroutine

      subroutine disprs_d4_energy(numbers, positions, method, charge, atm, energy, gradient, status, &
         lattice, periodic, virial, model_kind, damping, coordination, charges, c6, polarizabilities, pair2, pair3, &
         ga, gc, wf, cutoff, charge_model, charge_cutoff, ghost, partition, hessian, &
         dcndr, dcndL, dqdr, dqdL, dc6dr, dc6dL, dalphadr, dalphadL, fixed_charges, lowrank)
      integer(c_int), intent(in), target, contiguous :: numbers(:)
      real(c_double), intent(in), target, contiguous :: positions(:, :)
      character(len=*), intent(in) :: method
      real(c_double), intent(in), target :: charge
      logical(c_bool), intent(in) :: atm
      real(c_double), intent(out), target :: energy
      real(c_double), intent(out), target, contiguous :: gradient(:, :)
      integer(c_int), intent(out) :: status
      real(c_double), intent(in), target, contiguous, optional :: lattice(:, :)
      logical(c_bool), intent(in), target, contiguous, optional :: periodic(:)
      real(c_double), intent(out), target, contiguous, optional :: virial(:, :)
      integer(c_int), intent(in), optional :: model_kind
      integer(c_int), intent(in), optional :: charge_model
      integer(c_int), intent(in), optional :: ghost(:), partition(2)
      integer(c_int), allocatable, target :: cghost(:)
      real(c_double), intent(in), optional :: damping(d4_parameter_count)
      real(c_double), intent(in), optional :: ga, gc, wf, cutoff(cutoff_count), charge_cutoff
      real(c_double), intent(in), target, contiguous, optional :: fixed_charges(:)
      type(d4_lowrank_config), intent(in), optional :: lowrank
      real(c_double) :: selected_ga, selected_gc, selected_wf
      real(c_double), intent(out), target, contiguous, optional :: coordination(:), charges(:), &
         c6(:, :), polarizabilities(:), pair2(:, :), pair3(:, :)
      real(c_double), intent(out), target, contiguous, optional :: hessian(:, :, :, :)
      real(c_double), intent(inout), target, contiguous, optional :: dcndr(:, :, :), dcndL(:, :, :), &
         dqdr(:, :, :), dqdL(:, :, :), dalphadr(:, :, :), dalphadL(:, :, :)
      real(c_double), intent(inout), target, contiguous, optional :: dc6dr(:, :, :, :), dc6dL(:, :, :, :)
      character(kind=c_char, len=:), allocatable :: c_method
      type(c_ptr) :: error, structure, model, param, lattice_ptr, periodic_ptr, virial_ptr
      integer(c_int) :: selected_model
      type(c_ptr) :: cn_ptr, charges_ptr, c6_ptr, alpha_ptr
      type(c_ptr) :: response_ptr(response_count)
      integer :: nat

      if (size(positions, 1) /= 3 .or. size(positions, 2) /= size(numbers) .or. &
            any(shape(gradient) /= shape(positions))) error stop "invalid disprs D4 array shape"
      nat = size(numbers)
      status = 1_c_int
      if (present(fixed_charges)) then
         if (size(fixed_charges) /= nat) return
      end if
      response_ptr = c_null_ptr
      if (present(dcndr)) then
         if (any(shape(dcndr) /= [3, nat, nat])) return
         response_ptr(1) = c_loc(dcndr)
      end if
      if (present(dcndL)) then
         if (any(shape(dcndL) /= [3, 3, nat])) return
         response_ptr(2) = c_loc(dcndL)
      end if
      if (present(dqdr)) then
         if (any(shape(dqdr) /= [3, nat, nat])) return
         response_ptr(3) = c_loc(dqdr)
      end if
      if (present(dqdL)) then
         if (any(shape(dqdL) /= [3, 3, nat])) return
         response_ptr(4) = c_loc(dqdL)
      end if
      if (present(dc6dr)) then
         if (any(shape(dc6dr) /= [3, nat, nat, nat])) return
         response_ptr(5) = c_loc(dc6dr)
      end if
      if (present(dc6dL)) then
         if (any(shape(dc6dL) /= [3, 3, nat, nat])) return
         response_ptr(6) = c_loc(dc6dL)
      end if
      if (present(dalphadr)) then
         if (any(shape(dalphadr) /= [3, nat, nat])) return
         response_ptr(7) = c_loc(dalphadr)
      end if
      if (present(dalphadL)) then
         if (any(shape(dalphadL) /= [3, 3, nat])) return
         response_ptr(8) = c_loc(dalphadL)
      end if
      if (present(ghost)) then
         if (any(ghost < 1_c_int) .or. any(ghost > size(numbers))) then
            status = 1_c_int
            return
         end if
      end if
      c_method = trim(method) // c_null_char
      lattice_ptr = c_null_ptr
      periodic_ptr = c_null_ptr
      virial_ptr = c_null_ptr
      selected_model = 0_c_int
      if (present(lattice)) then
         if (any(shape(lattice) /= [3, 3])) error stop "invalid D4 lattice shape"
         lattice_ptr = c_loc(lattice)
      end if
      if (present(periodic)) then
         if (size(periodic) /= 3) error stop "invalid D4 periodic shape"
         periodic_ptr = c_loc(periodic)
      end if
      if (present(virial)) then
         if (any(shape(virial) /= [3, 3])) error stop "invalid D4 virial shape"
         virial_ptr = c_loc(virial)
      end if
      if (present(model_kind)) selected_model = model_kind
      cn_ptr = c_null_ptr
      charges_ptr = c_null_ptr
      c6_ptr = c_null_ptr
      alpha_ptr = c_null_ptr
      if (present(coordination)) then
         if (size(coordination) /= size(numbers)) error stop "invalid D4 coordination shape"
         cn_ptr = c_loc(coordination)
      end if
      if (present(charges)) then
         if (size(charges) /= size(numbers)) error stop "invalid D4 charges shape"
         charges_ptr = c_loc(charges)
      end if
      if (present(c6)) then
         if (any(shape(c6) /= [size(numbers), size(numbers)])) error stop "invalid D4 C6 shape"
         c6_ptr = c_loc(c6)
      end if
      if (present(polarizabilities)) then
         if (size(polarizabilities) /= size(numbers)) error stop "invalid D4 polarizabilities shape"
         alpha_ptr = c_loc(polarizabilities)
      end if
      if (present(pair2) .neqv. present(pair3)) error stop "both D4 pairwise outputs are required"
      if (present(hessian)) then
         if (any(shape(hessian) /= [3, size(numbers), 3, size(numbers)])) error stop "invalid D4 Hessian shape"
      end if
      if (present(pair2) .and. present(pair3)) then
         if (any(shape(pair2) /= [size(numbers), size(numbers)]) .or. &
             any(shape(pair3) /= [size(numbers), size(numbers)])) error stop "invalid D4 pairwise shape"
      end if
      error = c_d4_new_error()
      model = c_null_ptr
      param = c_null_ptr
      structure = c_d4_new_structure(error, size(numbers), c_loc(numbers), c_loc(positions), &
         c_loc(charge), lattice_ptr, periodic_ptr)
      if (c_d4_check_error(error) /= 0) go to 900
      selected_ga = 3.0_c_double; selected_gc = 2.0_c_double; selected_wf = 6.0_c_double
      if (present(ga)) selected_ga = ga
      if (present(gc)) selected_gc = gc
      if (present(wf)) selected_wf = wf
      model = c_d4_new_model(error, structure, selected_model, selected_ga, selected_gc, selected_wf)
      if (c_d4_check_error(error) /= 0) go to 900
      if (present(lowrank)) call c_d4_set_ewald(error, model, lowrank%rank, &
         lowrank%tolerance, lowrank%kcut, lowrank%mesh)
      if (c_d4_check_error(error) /= 0) go to 900
      if (present(ghost)) then
         if (size(ghost) > 0) then
            cghost = ghost - 1_c_int
            call c_d4_set_ghost(error, model, c_loc(cghost), int(size(cghost), c_int))
            if (c_d4_check_error(error) /= 0) go to 900
         end if
      end if
      if (present(partition)) call c_d4_set_partition(error, model, partition(1), partition(2))
      if (c_d4_check_error(error) /= 0) go to 900
      if (present(charge_model)) call c_d4_set_charge_model(error, model, charge_model)
      if (c_d4_check_error(error) /= 0) go to 900
      if (present(fixed_charges)) then
         call c_d4_set_fixed_charges(error, model, c_loc(fixed_charges), int(size(fixed_charges), c_int))
         if (c_d4_check_error(error) /= 0) go to 900
      end if
      if (present(charge_cutoff)) call c_d4_set_charge_cutoff(error, model, charge_cutoff)
      if (c_d4_check_error(error) /= 0) go to 900
      if (present(cutoff)) then
         call c_d4_set_cutoff(error, model, cutoff(1), cutoff(2), cutoff(3), cutoff(4), cutoff(5))
         if (c_d4_check_error(error) /= 0) go to 900
      end if
      if (present(damping)) then
         param = c_d4_new_rational_damping(error, damping(1), damping(2), &
            merge(damping(3), 0.0_c_double, atm), damping(4), damping(5), damping(6))
      else
         param = c_d4_load_param(error, c_method, atm)
      end if
      if (c_d4_check_error(error) /= 0) go to 900
      call c_d4_get_dispersion(error, structure, model, param, c_loc(energy), &
         c_loc(gradient), virial_ptr)
         if (c_d4_check_error(error) == 0 .and. &
               (present(dcndr) .or. present(dcndL) .or. present(dqdr) .or. present(dqdL) .or. &
                present(dc6dr) .or. present(dc6dL) .or. present(dalphadr) .or. present(dalphadL))) then
             call c_d4_property_response(error, structure, model, cn_ptr, charges_ptr, c6_ptr, alpha_ptr, &
                  response_ptr(1), response_ptr(2), response_ptr(3), response_ptr(4), &
                  response_ptr(5), response_ptr(6), response_ptr(7), response_ptr(8))
         else if (c_d4_check_error(error) == 0 .and. &
          (present(coordination) .or. present(charges) .or. present(c6) .or. present(polarizabilities))) then
         call c_d4_get_properties(error, structure, model, cn_ptr, charges_ptr, c6_ptr, alpha_ptr)
      end if
      if (c_d4_check_error(error) == 0 .and. present(pair2) .and. present(pair3)) then
         call c_d4_get_pairwise(error, structure, model, param, c_loc(pair2), c_loc(pair3))
      end if
      if (c_d4_check_error(error) == 0 .and. present(hessian)) then
         call c_d4_hessian(error, structure, model, param, c_loc(hessian))
      end if
   900   status = c_d4_check_error(error)
      call c_d4_delete_param(param)
      call c_d4_delete_model(model)
      call c_d4_delete_structure(structure)
      call c_d4_delete_error(error)
   end subroutine

end module