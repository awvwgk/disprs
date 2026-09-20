program test_dftd3
   use, intrinsic :: iso_c_binding, only : c_associated
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use dftd3, only : d3_model, d3_param, d3_lowrank_config, rational_damping_param, &
      & realspace_cutoff, gcp_param, new_d3_model, new_rational_damping, get_dispersion, &
      & get_pairwise_dispersion, get_properties, work_partition, &
      & get_rational_damping, get_z_damping, &
      & get_dftd3_version, get_lattice_points, &
      & get_gcp_param, get_geometric_counterpoise, get_geometric_counterpoise_hessian
   implicit none
   type(structure_type) :: mol, periodic_mol
   type(d3_model) :: model
   type(rational_damping_param) :: damping
   type(d3_param) :: named
   type(error_type), allocatable :: error
   type(gcp_param) :: gcp
   type(realspace_cutoff) :: cutoff
   integer, parameter :: numbers(2) = [6, 6]
   integer, parameter :: ndim = 3 * size(numbers)
   real(wp), parameter :: positions(3, 2) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, 0.0_wp, 0.0_wp], [3, 2])
   real(wp), parameter :: periodic_positions(3, 2) = reshape([ &
      & 0.4_wp, 0.8_wp, 1.2_wp, 4.4_wp, 2.8_wp, 3.6_wp], [3, 2])
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 8.0_wp], [3, 3])
   real(wp) :: energy, gradient(3, 2), sigma(3, 3), hessian(ndim, ndim)
   real(wp) :: pair2(2, 2), pair3(2, 2)
   real(wp) :: serial_energy, total_energy, serial_gradient(3, 2), total_gradient(3, 2)
   integer :: part
   character(len=32) :: mode
   character(len=:), allocatable :: version
   real(wp), allocatable :: translations(:, :)
   call get_dftd3_version(string=version)
   if (version /= '1.6.0') error stop 'D3 version failed'
   call get_lattice_points([.false., .false., .false.], lattice, 10._wp, translations)
   if (size(translations, 2) /= 1) error stop 'D3 lattice helper failed'

   block
      type(structure_type) :: smooth_mol
      type(d3_model) :: smooth_model
      type(rational_damping_param) :: smooth_damping
      real(wp) :: ordinary
      call new(smooth_mol, [1, 28], positions)
      call new_rational_damping(smooth_damping, d3_param(s6=1._wp, s8=0.7875_wp, a1=0.4289_wp, a2=4.4407_wp))
      call new_d3_model(smooth_model, smooth_mol)
      call get_dispersion(smooth_mol, smooth_model, smooth_damping, cutoff, ordinary, gradient, sigma)
      call new_d3_model(smooth_model, smooth_mol, d3s=.true., error=error)
      if (allocated(error)) error stop 'D3S construction failed'
      call get_dispersion(error, smooth_mol, smooth_model, smooth_damping, cutoff, energy, gradient, sigma, hessian)
      if (allocated(error)) error stop 'D3S derivatives failed'
      if (abs(energy-ordinary) < 1e-10_wp) error stop 'D3S model selection ignored'
      ordinary = energy
      call new_d3_model(smooth_model, smooth_mol, model_kind=1, error=error)
      if (allocated(error)) error stop 'D3S kind construction failed'
      call get_dispersion(smooth_mol, smooth_model, smooth_damping, cutoff, energy, gradient, sigma)
      if (abs(energy-ordinary) > 1e-14_wp) error stop 'D3S kind differs from alias'
      call new_d3_model(smooth_model, smooth_mol, model_kind=2, error=error)
      if (.not.allocated(error)) error stop 'Invalid D3 kind accepted'
      call new_d3_model(smooth_model, smooth_mol, model_kind=0, d3s=.true., error=error)
      if (.not.allocated(error)) error stop 'Conflicting D3 selectors accepted'
      call new_d3_model(smooth_model, smooth_mol, d3s=.true., lowrank=d3_lowrank_config(), error=error)
      if (.not.allocated(error)) error stop 'D3S Fourier accepted'
      deallocate(error)
   end block

   block
      use dftd3, only : damping_param, zero_damping_param, mzero_damping_param, &
         optimizedpower_damping_param, cso_damping_param, z_damping_param, &
         new_zero_damping, new_mzero_damping, new_optimizedpower_damping, new_cso_damping, new_z_damping
      class(damping_param), allocatable :: parameter
      type(d3_param) :: values
      integer :: family, attempt
      do family = 1, 6
         select case(family)
         case(1)
            allocate(rational_damping_param :: parameter)
         case(2)
            allocate(zero_damping_param :: parameter)
         case(3)
            allocate(mzero_damping_param :: parameter)
         case(4)
            allocate(optimizedpower_damping_param :: parameter)
         case(5)
            allocate(cso_damping_param :: parameter)
         case(6)
            allocate(z_damping_param :: parameter)
         case default
            error stop 'Unknown damping test family'
         end select
         do attempt = 1, 3
            values = d3_param()
            if (attempt == 2) values%alp = ieee_value(0._wp, ieee_quiet_nan)
            select type(parameter)
            type is(rational_damping_param)
               call new_rational_damping(parameter, values, error=error)
            type is(zero_damping_param)
               call new_zero_damping(parameter, values, error=error)
            type is(mzero_damping_param)
               call new_mzero_damping(parameter, values, error=error)
            type is(optimizedpower_damping_param)
               call new_optimizedpower_damping(parameter, values, error=error)
            type is(cso_damping_param)
               call new_cso_damping(parameter, values, error=error)
            type is(z_damping_param)
               call new_z_damping(parameter, values, error=error)
            end select
            if (attempt == 2) then
               if (.not.allocated(error)) error stop 'Invalid damping constructor accepted'
               if (index(error%message, 'finite') == 0) error stop 'Native damping error lost'
               if (c_associated(parameter%handle)) error stop 'Failed damping constructor retained handle'
            else
               if (allocated(error)) error stop 'Damping constructor recovery failed'
               if (.not.c_associated(parameter%handle)) error stop 'Damping constructor returned no handle'
            end if
         end do
         call parameter%close()
         deallocate(parameter)
      end do
   end block

   call new(mol, numbers, positions)
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case(mode)
      case('--invalid-model')
         call new_d3_model(model, mol, lowrank=d3_lowrank_config(), ghost=[3])
      case('--invalid-damping')
         call new_rational_damping(damping, d3_param(alp=ieee_value(0._wp, ieee_quiet_nan)))
      case('--invalid-gcp-cutoff')
         call get_gcp_param(gcp, mol, 'pbeh3c')
         cutoff%gcp = -1._wp
         energy = 42._wp
         call get_geometric_counterpoise(mol, gcp, cutoff, energy)
      case('--invalid-gcp-constructor')
         call get_gcp_param(gcp, mol, 'unknown')
      case('--invalid-gcp-hessian')
         call get_gcp_param(gcp, mol, 'pbeh3c')
         call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian, work_partition(2, 2))
      case default
         error stop 'Unknown test mode'
      end select
      stop 0
   end if
   call get_rational_damping(named, 'PBE', error, s9=1.0_wp)
   if (allocated(error)) error stop 'named D3 loading failed'
   if (abs(named%a1-0.4289_wp) > 1e-15_wp .or. abs(named%s9-1.0_wp) > 1e-15_wp) error stop 'named D3 values failed'
   call get_z_damping(named, 'unknown', error)
   if (.not. allocated(error)) error stop 'unknown Z accepted'
   deallocate(error)
   call new_d3_model(model, mol, lowrank=d3_lowrank_config(), ghost=[3], error=error)
   if (.not.allocated(error)) error stop 'Invalid D3 ghost was masked by lowrank setup'
   if (index(error%message, 'ghost index') == 0) error stop 'Original D3 ghost error lost'
   if (c_associated(model%handle) .or. c_associated(model%structure) .or. c_associated(model%error)) then
     error stop 'Failed D3 construction leaked handles'
   end if
   if (model%nat /= 0) error stop 'Failed D3 model retained atom count'
   mol%xyz(1, 1) = ieee_value(0._wp, ieee_quiet_nan)
   call new_d3_model(model, mol, error=error)
   if (.not.allocated(error)) error stop 'Nonfinite D3 construction accepted'
   if (index(error%message, 'nonfinite geometry') == 0) error stop 'Original D3 structure error lost'
   mol%xyz = positions
   call new_d3_model(model, mol, lowrank=d3_lowrank_config(tolerance=ieee_value(0._wp, ieee_quiet_nan)), error=error)
   if (.not.allocated(error)) error stop 'Nonfinite D3 lowrank setup accepted'
   if (index(error%message, 'Ewald controls') == 0) error stop 'Original D3 lowrank error lost'
   call new_d3_model(model, mol, error=error)
   if (allocated(error)) error stop 'D3 construction recovery failed'
   call new_d3_model(model, mol)
   block
      real(wp) :: cn(2), c6(2, 2), reference_cn(2), reference_c6(2, 2)
      real(wp) :: dcndr(3, 2, 2), dcndL(3, 3, 2), dc6dr(3, 2, 2, 2), dc6dL(3, 3, 2, 2)
      real(wp) :: plus_cn(2), minus_cn(2), plus_c6(2, 2), minus_c6(2, 2)
      real(wp), parameter :: response_step = 1e-5_wp
      call get_properties(mol, model, cutoff, reference_cn, reference_c6)
      call get_properties(error, mol, model, cutoff, dcndr=dcndr, dcndL=dcndL, dc6dr=dc6dr, dc6dL=dc6dL)
      if (allocated(error)) error stop 'D3 response query failed'
      mol%xyz(1, 2) = positions(1, 2) + response_step
      call get_properties(mol, model, cutoff, plus_cn, plus_c6)
      mol%xyz(1, 2) = positions(1, 2) - response_step
      call get_properties(mol, model, cutoff, minus_cn, minus_c6)
      mol%xyz = positions
      if (maxval(abs((plus_cn-minus_cn)/(2*response_step)-dcndr(1, 2, :))) > 1e-8_wp) then
        error stop 'D3 CN response layout failed'
      end if
      if (maxval(abs((plus_c6-minus_c6)/(2*response_step)-dc6dr(1, 2, :, :))) > 1e-7_wp) then
        error stop 'D3 C6 response layout failed'
      end if
      if (maxval(abs(dcndL(1, 1, :)-6*dcndr(1, 2, :))) > 1e-12_wp .or. &
          maxval(abs(dc6dL(1, 1, :, :)-6*dc6dr(1, 2, :, :))) > 1e-12_wp) error stop 'D3 strain layout failed'
      dcndr = 42._wp
      call get_properties(error, mol, model, dcndr=dcndr(:, :1, :))
      if (.not.allocated(error) .or. any(abs(dcndr-42._wp) > 0._wp)) error stop 'D3 response shape guard failed'
      call get_properties(error, mol, model, dcndr=dcndr)
      if (allocated(error)) error stop 'D3 response recovery failed'
      if (any(reference_cn <= 0._wp) .or. any(reference_c6 <= 0._wp)) error stop 'D3 properties missing'
      if (maxval(abs(reference_c6-transpose(reference_c6))) > 0._wp) error stop 'D3 C6 symmetry failed'
      call get_properties(error, mol, model, cn=cn)
      if (allocated(error) .or. maxval(abs(cn-reference_cn)) > 0._wp) error stop 'D3 CN-only query failed'
      call get_properties(error, mol, model, c6=c6)
      if (allocated(error) .or. maxval(abs(c6-reference_c6)) > 0._wp) error stop 'D3 C6-only query failed'
      call get_properties(error, mol, model, cutoff, cn(:1), c6)
      if (.not.allocated(error)) error stop 'Invalid D3 CN shape accepted'
      call get_properties(error, mol, model, cutoff, cn, c6(:1, :))
      if (.not.allocated(error)) error stop 'Invalid D3 C6 shape accepted'
      cutoff%cn = -1._wp
      cn = 42._wp; c6 = 42._wp
      call get_properties(error, mol, model, cutoff, cn, c6)
      if (.not.allocated(error)) error stop 'Invalid D3 property cutoff accepted'
      if (maxval(abs(cn-42._wp)) > 0._wp .or. maxval(abs(c6-42._wp)) > 0._wp) then
        error stop 'Failed D3 properties modified outputs'
      end if
      cutoff%cn = 0._wp
      call get_properties(error, mol, model, cutoff, cn, c6)
      if (allocated(error) .or. maxval(abs(cn)) > 0._wp) error stop 'D3 property cutoff recovery failed'
      cutoff = realspace_cutoff()
      mol%xyz(1, 2) = 4._wp
      call get_properties(error, mol, model, cutoff, cn, c6)
      if (allocated(error) .or. maxval(abs(cn-reference_cn)) < 1e-12_wp) then
        error stop 'D3 property geometry update failed'
      end if
      mol%xyz = positions
   end block
   call new_rational_damping(damping, d3_param( &
      & s6=1.0_wp, s8=0.7875_wp, a1=0.4289_wp, a2=4.4407_wp))
   call get_dispersion(mol, model, damping, cutoff, energy, gradient, sigma)
   if (abs(energy + 0.0005341413931338267_wp) > 1.0e-14_wp) error stop
   call get_pairwise_dispersion(mol, model, damping, cutoff, pair2, pair3)
   if (abs(sum(pair2)+sum(pair3)-energy) > 1.0e-14_wp) error stop 'D3 pairwise sum failed'
   if (maxval(abs(pair2-transpose(pair2))) > 1.0e-14_wp) error stop 'D3 pairwise symmetry failed'
   call get_pairwise_dispersion(error=error, mol=mol, disp=model, param=damping, cutoff=cutoff, &
      energy2=pair2, energy3=pair3)
   if (allocated(error)) error stop 'D3 pairwise error overload failed'
   if (abs(sum(pair2)+sum(pair3)-energy) > 1e-14_wp) error stop 'D3 pairwise error overload sum failed'
   call get_pairwise_dispersion(error, mol, model, damping, cutoff, pair2(:1, :), pair3)
   if (.not.allocated(error)) error stop 'invalid D3 pairwise shape accepted'
   call get_pairwise_dispersion(error, mol, model, damping, cutoff, pair2, pair3, work_partition(2, 2))
   if (.not.allocated(error)) error stop 'invalid D3 pairwise partition accepted'
   cutoff%cn = -1.0_wp
   call get_pairwise_dispersion(error, mol, model, damping, cutoff, pair2, pair3)
   if (.not.allocated(error)) error stop 'invalid D3 pairwise cutoff accepted'
   cutoff = realspace_cutoff()
   call get_pairwise_dispersion(error, mol, model, damping, cutoff, pair2, pair3)
   if (allocated(error)) error stop 'D3 pairwise error recovery failed'
   serial_energy = energy; serial_gradient = gradient
   block
      real(wp) :: reference_sigma(3, 3), atomic(2), reference_hessian(ndim, ndim)
      reference_sigma = sigma
      gradient = 42.0_wp; sigma = 42.0_wp
      call get_dispersion(error, mol, model, damping, cutoff, energy, gradient, sigma, hessian)
      if (allocated(error)) error stop 'combined D3 derivatives failed'
      if (maxval(abs(gradient-serial_gradient)) > 1e-14_wp) error stop 'combined D3 gradient missing'
      if (maxval(abs(sigma-reference_sigma)) > 1e-14_wp) error stop 'combined D3 virial missing'
      reference_hessian = hessian
      call get_dispersion(error=error, mol=mol, disp=model, param=damping, cutoff=cutoff, &
         energies=atomic, gradient=gradient, sigma=sigma, hessian=hessian)
      if (allocated(error)) error stop 'D3 atomic overload failed'
      if (abs(sum(atomic)-serial_energy) > 1e-14_wp) error stop 'D3 atomic sum failed'
      if (maxval(abs(gradient-serial_gradient)) > 1e-14_wp) error stop 'D3 atomic gradient failed'
      if (maxval(abs(sigma-reference_sigma)) > 1e-14_wp) error stop 'D3 atomic virial failed'
      if (maxval(abs(hessian-reference_hessian)) > 1e-14_wp) error stop 'D3 atomic Hessian failed'
      call get_dispersion(mol, model, damping, cutoff, atomic)
      if (abs(sum(atomic)-serial_energy) > 1e-14_wp) error stop 'D3 legacy atomic sum failed'
      call get_dispersion(error, mol, model, damping, cutoff, atomic(:1))
      if (.not.allocated(error)) error stop 'invalid D3 atomic shape accepted'
      call get_dispersion(error, mol, model, damping, cutoff, atomic, gradient=gradient(:2, :))
      if (.not.allocated(error)) error stop 'invalid D3 atomic gradient shape accepted'
      call get_dispersion(error, mol, model, damping, cutoff, energy, sigma=sigma(:2, :))
      if (.not.allocated(error)) error stop 'invalid D3 virial shape accepted'
      call get_dispersion(error, mol, model, damping, cutoff, energy, hessian=hessian(:2, :))
      if (.not.allocated(error)) error stop 'invalid D3 Hessian shape accepted'
      call get_dispersion(error, mol, model, damping, cutoff, energy, partition=work_partition(2, 2))
      if (.not.allocated(error)) error stop 'invalid D3 partition accepted'
      cutoff%cn = -1.0_wp
      call get_dispersion(error, mol, model, damping, cutoff, energy)
      if (.not.allocated(error)) error stop 'invalid D3 cutoff accepted'
      cutoff = realspace_cutoff()
      call get_dispersion(error, mol, model, damping, cutoff, energy)
      if (allocated(error)) error stop 'D3 error recovery failed'
   end block
   total_energy = 0.0_wp; total_gradient = 0.0_wp
   do part = 0, 2
      call get_dispersion(mol, model, damping, cutoff, energy, gradient, sigma, partition=work_partition(part, 3))
      call get_pairwise_dispersion(mol, model, damping, cutoff, pair2, pair3, partition=work_partition(part, 3))
      if (abs(sum(pair2)+sum(pair3)-energy) > 1e-14_wp) error stop 'D3 partition pairs failed'
      total_energy = total_energy + energy; total_gradient = total_gradient + gradient
   end do
   if (abs(total_energy-serial_energy) > 1e-14_wp) error stop 'D3 partition sum failed'
   if (maxval(abs(total_gradient-serial_gradient)) > 1e-14_wp) error stop 'D3 partition gradient failed'
   call get_dispersion(mol, model, damping, cutoff, energy)
   if (abs(energy-serial_energy) > 1e-14_wp) error stop 'D3 serial reset failed'
   call new_d3_model(model, mol, ghost=[1])
   call get_dispersion(mol, model, damping, cutoff, energy, gradient, sigma)
   if (abs(energy) > tiny(energy) .or. maxval(abs(gradient)) > tiny(energy)) error stop
   call get_gcp_param(gcp, mol, 'pbeh3c')
   energy = 0.0_wp; gradient = 0.0_wp; sigma = 0.0_wp
   call get_geometric_counterpoise(mol, gcp, cutoff, energy, gradient, sigma)
   if (abs(energy - 2.581420786404301e-5_wp) > 1.0e-15_wp) error stop
   if (abs(gradient(1, 1) - 5.324722896549186e-5_wp) > 1.0e-14_wp) error stop
   if (abs(sigma(1, 1) + 3.1948337379295114e-4_wp) > 1.0e-13_wp) error stop
   call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian)
   if (abs(hessian(1, 1) - 9.51416946673085e-5_wp) > 1.0e-13_wp) error stop
   block
      type(structure_type) :: smaller
      real(wp) :: reference_energy, reference_gradient(3, 2), reference_sigma(3, 3)
      real(wp) :: reference_hessian(ndim, ndim), strided_gradient(ndim, 2)
      real(wp) :: reference_eta, reference_alpha
      integer :: failure
      reference_energy = energy; reference_gradient = gradient; reference_sigma = sigma
      reference_hessian = hessian; reference_eta = gcp%eta; reference_alpha = gcp%alpha
      do failure = 1, 4
         select case(failure)
         case(1)
            cutoff%gcp = -1._wp
         case(2)
            mol%xyz(1, 1) = ieee_value(0._wp, ieee_quiet_nan)
         case(3)
            gcp%eta = ieee_value(0._wp, ieee_quiet_nan)
         case(4)
            gcp%alpha = ieee_value(0._wp, ieee_quiet_nan)
         case default
            error stop 'Unknown gCP failure case'
         end select
         energy = 42._wp; gradient = 42._wp; sigma = 42._wp; hessian = 42._wp
         call get_geometric_counterpoise(mol, gcp, cutoff, energy, gradient, sigma, error=error)
         if (.not.allocated(error)) error stop 'Invalid gCP setup accepted'
         if (failure == 1) then
            if (index(error%message, 'cutoffs') == 0) error stop 'Original gCP cutoff error lost'
         else if (failure == 2) then
            if (index(error%message, 'nonfinite geometry') == 0) error stop 'Original gCP update error lost'
         end if
          if (abs(energy-42._wp) > 0._wp .or. any(abs(gradient-42._wp) > 0._wp) .or. &
             any(abs(sigma-42._wp) > 0._wp)) then
            error stop 'Failed gCP modified outputs'
          end if
         call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian, error=error)
         if (.not.allocated(error)) error stop 'Invalid gCP Hessian setup accepted'
         if (any(abs(hessian-42._wp) > 0._wp)) error stop 'Failed gCP modified Hessian'
         cutoff = realspace_cutoff(); mol%xyz = positions; gcp%eta = reference_eta; gcp%alpha = reference_alpha
         call get_geometric_counterpoise(mol, gcp, cutoff, energy, gradient, sigma, error=error)
         if (allocated(error)) error stop 'gCP recovery failed'
         if (abs(energy-reference_energy) > 1e-15_wp) error stop 'gCP recovery changed energy'
         call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian, error=error)
         if (allocated(error)) error stop 'gCP Hessian recovery failed'
         if (maxval(abs(hessian-reference_hessian)) > 1e-13_wp) error stop 'gCP Hessian recovery changed result'
      end do
      call get_geometric_counterpoise(mol, gcp, cutoff, energy, partition=work_partition(2, 2), error=error)
      if (.not.allocated(error)) error stop 'Invalid gCP partition accepted'
      call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian, work_partition(2, 2), error)
      if (.not.allocated(error)) error stop 'Invalid gCP Hessian partition accepted'
      call get_geometric_counterpoise(mol, gcp, cutoff, energy, gradient(:2, :), error=error)
      if (.not.allocated(error)) error stop 'Invalid gCP gradient shape accepted'
      call get_geometric_counterpoise(mol, gcp, cutoff, energy, sigma=sigma(:2, :), error=error)
      if (.not.allocated(error)) error stop 'Invalid gCP virial shape accepted'
      call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian(:2, :), error=error)
      if (.not.allocated(error)) error stop 'Invalid gCP Hessian shape accepted'
      call new(smaller, [6], positions(:, :1))
      call get_geometric_counterpoise(smaller, gcp, cutoff, energy, error=error)
      if (.not.allocated(error)) error stop 'Changed gCP atom count accepted'
      call get_geometric_counterpoise_hessian(smaller, gcp, cutoff, hessian(:3, :3), error=error)
      if (.not.allocated(error)) error stop 'Changed gCP Hessian atom count accepted'
      gcp%zeff = [6, 6]
      call get_geometric_counterpoise(mol, gcp, cutoff, energy, error=error)
      if (.not.allocated(error)) error stop 'Invalid gCP parameter shape accepted'
      call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian, error=error)
      if (.not.allocated(error)) error stop 'Invalid gCP Hessian parameter shape accepted'
      gcp%zeff = [6]
      strided_gradient = 42._wp
      call get_geometric_counterpoise(mol, gcp, cutoff, energy, strided_gradient(::2, :), sigma, error=error)
      if (allocated(error)) error stop 'Strided gCP output failed'
      if (maxval(abs(strided_gradient(::2, :)-reference_gradient)) > 1e-14_wp) error stop 'Strided gCP gradient wrong'
      if (any(abs(strided_gradient(2::2, :)-42._wp) > 0._wp)) error stop 'Strided gCP overwrote unrelated elements'
      if (maxval(abs(sigma-reference_sigma)) > 1e-13_wp) error stop 'gCP recovery changed virial'
      do failure = 1, 3
         select case(failure)
         case(1)
            call get_gcp_param(gcp, mol, 'unknown', error=error)
         case(2)
            mol%xyz(1, 1) = ieee_value(0._wp, ieee_quiet_nan)
            call get_gcp_param(gcp, mol, 'pbeh3c', error=error)
         case(3)
            call get_gcp_param(gcp, mol, 'pbeh3c', eta=ieee_value(0._wp, ieee_quiet_nan), error=error)
         case default
            error stop 'Unknown gCP constructor failure case'
         end select
         if (.not.allocated(error)) error stop 'Invalid gCP constructor accepted'
         if (failure == 2) then
            if (index(error%message, 'nonfinite geometry') == 0) error stop 'Original gCP structure error lost'
         end if
         if (c_associated(gcp%handle) .or. c_associated(gcp%structure) .or. c_associated(gcp%error)) then
           error stop 'Failed gCP constructor retained handles'
         end if
         if (gcp%nat /= 0 .or. allocated(gcp%zeff)) error stop 'Failed gCP constructor retained parameters'
         mol%xyz = positions
         call get_gcp_param(gcp, mol, 'pbeh3c', error=error)
         if (allocated(error)) error stop 'gCP constructor recovery failed'
         call get_geometric_counterpoise(mol, gcp, cutoff, energy, error=error)
         if (allocated(error)) error stop 'Reconstructed gCP evaluation failed'
         if (abs(energy-reference_energy) > 1e-15_wp) error stop 'Reconstructed gCP changed energy'
      end do
      call get_gcp_param(gcp, mol, error=error)
      if (allocated(error)) error stop 'Empty gCP construction failed'
      if (allocated(gcp%zeff) .or. allocated(gcp%emiss) .or. allocated(gcp%xv) .or. allocated(gcp%slater)) then
        error stop 'Empty gCP construction allocated parameters'
      end if
      if (c_associated(gcp%handle) .or. c_associated(gcp%structure) .or. c_associated(gcp%error)) then
        error stop 'Empty gCP construction retained handles'
      end if
      call get_gcp_param(gcp, mol, 'pbeh3c', error=error)
      if (allocated(error)) error stop 'Reconstruction after empty gCP failed'
   end block
   gcp%eta = 0.0_wp; gcp%base = .false.; gcp%srb = .false.
   energy = 0.0_wp
   call get_geometric_counterpoise(mol, gcp, cutoff, energy)
   if (abs(energy) > tiny(energy)) error stop 'disabled gCP failed'
   call get_geometric_counterpoise_hessian(mol, gcp, cutoff, hessian)
   if (maxval(abs(hessian)) > tiny(energy)) error stop 'disabled gCP Hessian failed'
   call get_gcp_param(gcp, mol, 'b973c')
   if (.not. gcp%srb) error stop 'gCP defaults failed'
   gcp%srb = .false.
   energy = 0.0_wp
   call get_geometric_counterpoise(mol, gcp, cutoff, energy)
   if (abs(energy) > tiny(energy)) error stop 'disabled SRB failed'
   call new(periodic_mol, [6, 8], periodic_positions, lattice=lattice, periodic=[.true., .true., .true.])
   call get_gcp_param(gcp, periodic_mol, 'pbeh3c')
   cutoff%gcp = 12.0_wp; cutoff%srb = 12.0_wp; energy = 0.0_wp
   call get_geometric_counterpoise(periodic_mol, gcp, cutoff, energy)
   if (abs(energy - 7.503293440007774e-4_wp) > 1.0e-13_wp) error stop
   call get_geometric_counterpoise_hessian(periodic_mol, gcp, cutoff, hessian)
   if (abs(hessian(1, 1) - 6.505276577134137e-4_wp) > 1.0e-12_wp) error stop
   if (abs(hessian(4, 1) + 6.505276577134137e-4_wp) > 1.0e-12_wp) error stop
   call model%close()
   call model%close()
   call damping%close()
   call damping%close()
   call get_dispersion(error, mol, model, damping, cutoff, energy)
   if (.not.allocated(error)) error stop 'Closed D3 handles accepted'
   do part = 1, 3
      call new_d3_model(model, mol)
      call new_rational_damping(damping, d3_param(s6=1._wp, s8=0.7875_wp, a1=0.4289_wp, a2=4.4407_wp))
      call get_dispersion(error, mol, model, damping, cutoff, energy)
      if (allocated(error)) error stop 'D3 reconstruction failed'
      if (abs(energy + 0.0005341413931338267_wp) > 1e-14_wp) error stop 'D3 reconstruction changed energy'
   end do
   call damping%close()
   call model%close()
end program