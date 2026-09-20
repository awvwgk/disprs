program test_dftd4
   use mctc_env, only : wp, error_type
   use dftd4
   implicit none
   type(structure_type) :: mol
   class(dispersion_model), allocatable :: model
   class(damping_param), allocatable :: param
   type(d4_model) :: custom
   type(d4s_model) :: custom_s
   type(error_type), allocatable :: error
   type(realspace_cutoff) :: cutoff
   real(wp) :: energy, plus, minus, grad(3, 3), sigma(3, 3), check(3, 3), hessian(3, 3, 3, 3)
   real(wp) :: cn(3), charges(3), c6(3, 3), alpha(3), pair2(3, 3), pair3(3, 3), saved(3, 3)
   real(wp) :: dcndr(3, 3, 3), dcndL(3, 3, 3), dqdr(3, 3, 3), dqdL(3, 3, 3)
   real(wp) :: dc6dr(3, 3, 3, 3), dc6dL(3, 3, 3, 3), dalphadr(3, 3, 3), dalphadL(3, 3, 3)
   real(wp), parameter :: step = 2e-4_wp
   real(wp), parameter :: xyz(3, 3) = reshape([0._wp, 0._wp, 0._wp, 2._wp, 1._wp, 0._wp, &
      -1._wp, 2._wp, 0.5_wp], [3, 3])
   real(wp), parameter :: lattice(3, 3) = reshape([12._wp, 0._wp, 0._wp, 1._wp, 13._wp, 0._wp, &
      0.5_wp, 0.2_wp, 14._wp], [3, 3])
   integer :: kind, dims, charge_kind, atom, axis
   character(len=:), allocatable :: version
   real(wp), allocatable :: translations(:, :)
   call get_lattice_points([.false., .false., .false.], lattice, 10._wp, translations)
   if (size(translations, 2) /= 1) error stop 'D4 lattice helper failed'
   call get_dftd4_version(string=version)
   if (version /= '4.2.0') error stop 'D4 version failed'
   call get_rational_damping('unknown', param)
   if (allocated(param)) error stop 'Unknown D4 method accepted'
   call get_rational_damping(id=0, param=param)
   if (allocated(param)) error stop 'Invalid D4 ID accepted'
   call get_rational_damping(id=119, param=param)
   if (allocated(param)) error stop 'Out-of-range D4 ID accepted'
   call get_rational_damping(id=76, param=param, s9=0._wp)
   if (.not.allocated(param)) error stop 'D4 DFTB ID missing'
   select type(param)
   type is(rational_damping_param)
      if (maxval(abs([param%s8, param%a1, param%a2, param%s9] &
         - [0.4727337_wp, 0.5467502_wp, 4.4955068_wp, 0._wp])) > 1e-15_wp) then
        error stop 'D4 two-body-only fit failed'
      end if
   end select
   call get_rational_damping(id=11, param=param, s9=0.5_wp)
   if (.not.allocated(param)) error stop 'D4 PBE ID missing'
   select type(param)
   type is(rational_damping_param)
      if (maxval(abs([param%a1, param%a2, param%s9] - [0.38574991_wp, 4.80688534_wp, 0.5_wp])) &
         > 1e-15_wp) error stop 'D4 PBE ID values failed'
   end select
   call get_rational_damping('pbe/def2-TZVP', param)
   if (.not.allocated(param)) error stop 'D4 named damping missing'
   cutoff = realspace_cutoff(cn=7._wp, disp2=10._wp, disp3=8._wp, width2=1._wp, width3=1._wp)
   do kind = 0, 1
      do charge_kind = 1, 3, 2
         do dims = 0, 3
            call new(mol, [6, 8, 1], xyz, charge=0.25_wp, lattice=lattice, &
               periodic=[dims >= 1, dims >= 2, dims >= 3])
            call new_dispersion_model(error, model, mol, merge('d4 ', 'd4s', kind == 0), qmod=charge_kind)
            if (allocated(error)) error stop 'D4 model failed'
            call get_dispersion(mol, model, param, cutoff, energy, grad, sigma)
            call get_pairwise_dispersion(mol, model, param, cutoff, pair2, pair3)
            if (abs(sum(pair2)+sum(pair3)-energy) > 1e-13_wp) error stop 'D4 pair sums failed'
            call get_properties(mol, model, cutoff, cn, charges, c6, alpha, &
               dcndr, dcndL, dqdr, dqdL, dc6dr, dc6dL, dalphadr, dalphadL)
            if (maxval(abs(sum(dqdr, dim=3))) > 1e-12_wp .or. &
               maxval(abs(sum(dqdL, dim=3))) > 1e-12_wp) error stop 'D4 response charge conservation failed'
            block
               real(wp) :: plus_cn(3), minus_cn(3), plus_q(3), minus_q(3), plus_c6(3, 3), minus_c6(3, 3)
               real(wp) :: plus_alpha(3), minus_alpha(3)
               real(wp), parameter :: delta = 1e-5_wp
               mol%xyz(1, 2) = xyz(1, 2) + delta
               call get_properties(mol, model, cutoff, plus_cn, plus_q, plus_c6, plus_alpha)
               mol%xyz(1, 2) = xyz(1, 2) - delta
               call get_properties(mol, model, cutoff, minus_cn, minus_q, minus_c6, minus_alpha)
               mol%xyz = xyz
               if (maxval(abs((plus_cn-minus_cn)/(2*delta)-dcndr(1, 2, :))) > 1e-7_wp .or. &
                  maxval(abs((plus_q-minus_q)/(2*delta)-dqdr(1, 2, :))) > 1e-7_wp .or. &
                  maxval(abs((plus_c6-minus_c6)/(2*delta)-dc6dr(1, 2, :, :))) > 1e-6_wp .or. &
                  maxval(abs((plus_alpha-minus_alpha)/(2*delta)-dalphadr(1, 2, :))) > 1e-6_wp) then
                 error stop 'D4 Cartesian response layout failed'
               end if
               mol%xyz(1, :) = xyz(1, :) + delta * xyz(2, :)
               mol%lattice(1, :) = lattice(1, :) + delta * lattice(2, :)
               call get_properties(mol, model, cutoff, plus_cn, plus_q, plus_c6, plus_alpha)
               mol%xyz(1, :) = xyz(1, :) - delta * xyz(2, :)
               mol%lattice(1, :) = lattice(1, :) - delta * lattice(2, :)
               call get_properties(mol, model, cutoff, minus_cn, minus_q, minus_c6, minus_alpha)
               mol%xyz = xyz
               mol%lattice = lattice
               if (maxval(abs((plus_cn-minus_cn)/(2*delta)-dcndL(1, 2, :))) > 1e-7_wp .or. &
                  maxval(abs((plus_q-minus_q)/(2*delta)-dqdL(1, 2, :))) > 1e-7_wp .or. &
                  maxval(abs((plus_c6-minus_c6)/(2*delta)-dc6dL(1, 2, :, :))) > 1e-6_wp .or. &
                  maxval(abs((plus_alpha-minus_alpha)/(2*delta)-dalphadL(1, 2, :))) > 1e-6_wp) then
                 error stop 'D4 strain response layout failed'
               end if
            end block
            if (abs(sum(charges)-mol%charge) > 1e-12_wp) error stop 'D4 charge sum failed'
            if (any(c6 <= 0) .or. any(alpha <= 0)) error stop 'D4 properties failed'
            call get_dispersion(mol, model, param, cutoff, plus, sigma=check)
            if (maxval(abs(check-sigma)) > 1e-13_wp) error stop 'D4 virial only failed'
            saved = mol%xyz
            call get_dispersion_hessian(mol, model, param, cutoff, hessian)
            if (any(mol%xyz /= saved)) error stop 'D4 Hessian mutated structure'
            do atom = 1, 3
               do axis = 1, 3
                  mol%xyz(axis, atom) = saved(axis, atom) + step
                  call get_dispersion(mol, model, param, cutoff, plus, check)
                  mol%xyz(axis, atom) = saved(axis, atom) - step
                  block
                     real(wp) :: minus_grad(3, 3)
                     call get_dispersion(mol, model, param, cutoff, minus, minus_grad)
                     if (maxval(abs((check-minus_grad)/(2*step)-hessian(:, :, axis, atom))) > 2e-7_wp) then
                       error stop 'D4 Hessian differences failed'
                     end if
                  end block
                  if (abs((plus-minus)/(2*step)-grad(axis, atom)) > 2e-8_wp) error stop 'D4 gradient failed'
                  mol%xyz = saved
               end do
            end do
            call set_fixed_charges(error, model, charges)
            if (allocated(error)) error stop 'D4 fixed charge setup failed'
            call get_dispersion(mol, model, param, cutoff, plus)
            if (abs(plus-energy) > 1e-13_wp) error stop 'D4 fixed charge snapshot changed energy'
            call set_fixed_charges(error, model, [0.2_wp, -0.4_wp, 0.1_wp])
            if (allocated(error)) error stop 'D4 fixed charge update failed'
            call set_fixed_charges(error, model, [0.0_wp])
            if (.not.allocated(error)) error stop 'D4 invalid fixed charge count accepted'
            call get_properties(mol, model, cutoff, cn, charges, c6, alpha, dqdr=dqdr, dqdL=dqdL)
            if (maxval(abs(charges-[0.2_wp, -0.4_wp, 0.1_wp])) > 0._wp) error stop 'D4 fixed charges lost'
            if (maxval(abs(dqdr)) > 0._wp .or. maxval(abs(dqdL)) > 0._wp) error stop 'D4 fixed charge response nonzero'
            call get_dispersion(mol, model, param, cutoff, energy, grad)
            mol%xyz(1, 2) = xyz(1, 2) + step
            call get_dispersion(mol, model, param, cutoff, plus)
            mol%xyz(1, 2) = xyz(1, 2) - step
            call get_dispersion(mol, model, param, cutoff, minus)
            mol%xyz = xyz
            if (abs((plus-minus)/(2*step)-grad(1, 2)) > 2e-8_wp) error stop 'D4 fixed charge force failed'
            call set_fixed_charges(error, model)
            if (allocated(error)) error stop 'D4 fixed charge reset failed'
            call get_properties(mol, model, cutoff, cn, charges, c6, alpha)
            if (abs(sum(charges)-mol%charge) > 1e-12_wp) error stop 'D4 charge solver not restored'
         end do
      end do
   end do
   call new_d4_model(error=error, d4=custom, mol=mol, ga=2._wp, gc=1._wp, wf=4._wp)
   if (allocated(error)) error stop 'Custom D4 failed'
   call new_d4s_model(error=error, d4=custom_s, mol=mol, ga=2._wp, gc=1._wp)
   if (allocated(error)) error stop 'Custom D4S failed'
   call new_d4_model(error, custom, mol, fixed_charges=[0.2_wp, -0.4_wp, 0.1_wp])
   if (allocated(error)) error stop 'Fixed D4 constructor failed'
   call get_properties(mol, custom, cutoff, cn, charges, c6, alpha)
   if (maxval(abs(charges-[0.2_wp, -0.4_wp, 0.1_wp])) > 0._wp) error stop 'Fixed D4 constructor ignored charges'
   call new_d4s_model(error, custom_s, mol, fixed_charges=[0.2_wp, -0.4_wp, 0.1_wp])
   if (allocated(error)) error stop 'Fixed D4S constructor failed'
   call get_properties(mol, custom_s, cutoff, cn, charges, c6, alpha)
   if (maxval(abs(charges-[0.2_wp, -0.4_wp, 0.1_wp])) > 0._wp) error stop 'Fixed D4S constructor ignored charges'
   call new_d4_model(error, custom, mol, wf=-1._wp)
   if (.not.allocated(error)) error stop 'Invalid D4 weighting accepted'
   call new_dispersion_model(error, model, mol, 'unknown')
   if (.not.allocated(error)) error stop 'Unknown D4 model accepted'
   call new_d4_model(error, custom, mol, qmod=d4_qmod%gfn2)
   if (.not.allocated(error)) error stop 'Unsupported GFN2 accepted'
   call get_rational_damping('pbe', param, s9=0._wp)
   cutoff = realspace_cutoff()
   do kind = 0, 1
      call new_dispersion_model(error, model, mol, merge('d4 ', 'd4s', kind == 0), &
         lowrank=d4_lowrank_config(mesh=-1, tolerance=1e-10_wp))
      if (allocated(error)) error stop 'D4 direct Ewald constructor failed'
      call get_dispersion(mol, model, param, cutoff, energy, grad, sigma)
      call new_dispersion_model(error, model, mol, merge('d4 ', 'd4s', kind == 0), &
         lowrank=d4_lowrank_config(mesh=32, tolerance=1e-10_wp))
      if (allocated(error)) error stop 'D4 SPME constructor failed'
      call get_dispersion(mol, model, param, cutoff, plus, check)
      if (abs(energy-plus) > 5e-8_wp .or. maxval(abs(grad-check)) > 5e-8_wp) then
        error stop 'D4 direct/SPME mismatch'
      end if
   end do
end program test_dftd4