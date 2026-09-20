program test_multicharge
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use multicharge
   implicit none
   type(structure_type) :: mol, displaced
   class(mchrg_model_type), allocatable :: model, copied
   type(error_type), allocatable :: error
   real(wp) :: charges(3), plus(3), minus(3), cartesian(3, 3, 3), strain(3, 3, 3), single(3, 3, 3)
   real(wp), parameter :: step = 1e-5_wp
   real(wp), parameter :: xyz(3, 3) = reshape([0._wp, 0._wp, 0._wp, 2._wp, 1._wp, 0._wp, &
      -1._wp, 2._wp, 0.5_wp], [3, 3])
   real(wp), parameter :: lattice(3, 3) = reshape([12._wp, 0._wp, 0._wp, 1._wp, 13._wp, 0._wp, &
      0.5_wp, 0.2_wp, 14._wp], [3, 3])
   character(len=:), allocatable :: version
   integer :: kind, dims, atom, axis, column
   call get_multicharge_version(string=version)
   if (version /= '0.5.0') error stop 'Multicharge version failed'
   do kind = mchrg_model%eeq2019, mchrg_model%eeqbc2025
      do dims = 0, 3
         call new(mol, [6, 8, 1], xyz, charge=0.25_wp, lattice=lattice, &
            periodic=[dims >= 1, dims >= 2, dims >= 3])
         if (kind == mchrg_model%eeq2019) then
            call new_eeq2019_model(mol, model, error)
            if (allocated(error)) error stop 'EEQ model failed'
            call get_eeq_charges(mol, error, charges)
         else
            call new_eeqbc2025_model(mol, model, error)
            if (allocated(error)) error stop 'EEQBC model failed'
            call get_eeqbc_charges(mol, error, charges)
         end if
         if (allocated(error)) error stop 'Charge convenience call failed'
         copied = model
         call get_charges(copied, mol, error, plus, cartesian, strain)
         if (allocated(error)) error stop 'Charge derivatives failed'
         if (maxval(abs(charges-plus)) > 1e-13_wp) error stop 'Copied model charges differ'
         if (abs(sum(charges)-mol%charge) > 1e-12_wp) error stop 'Charge conservation failed'
         if (maxval(abs(sum(cartesian, dim=3))) > 1e-12_wp) error stop 'Cartesian charge conservation failed'
         if (maxval(abs(sum(strain, dim=3))) > 1e-12_wp) error stop 'Strain charge conservation failed'
         call get_charges(model, mol, error, plus, dqdr=single)
         if (allocated(error)) error stop 'Cartesian only failed'
         if (maxval(abs(single-cartesian)) > 1e-13_wp) error stop 'Cartesian only differs'
         call get_charges(model, mol, error, plus, dqdL=single)
         if (allocated(error)) error stop 'Strain only failed'
         if (maxval(abs(single-strain)) > 1e-13_wp) error stop 'Strain only differs'
         do atom = 1, 3
            do axis = 1, 3
               displaced = mol
               displaced%xyz(axis, atom) = mol%xyz(axis, atom) + step
               call get_charges(model, displaced, error, plus)
               if (allocated(error)) error stop 'Charge plus displacement failed'
               displaced%xyz(axis, atom) = mol%xyz(axis, atom) - step
               call get_charges(model, displaced, error, minus)
               if (allocated(error)) error stop 'Charge minus displacement failed'
               if (maxval(abs((plus-minus)/(2*step)-cartesian(axis, atom, :))) > 1e-7_wp) then
                 error stop 'Cartesian charge difference failed'
               end if
            end do
         end do
         do column = 1, 3
            do axis = 1, 3
               displaced = mol
               displaced%xyz(axis, :) = mol%xyz(axis, :) + step*mol%xyz(column, :)
               displaced%lattice(axis, :) = mol%lattice(axis, :) + step*mol%lattice(column, :)
               call get_charges(model, displaced, error, plus)
               if (allocated(error)) error stop 'Charge plus strain failed'
               displaced%xyz(axis, :) = mol%xyz(axis, :) - step*mol%xyz(column, :)
               displaced%lattice(axis, :) = mol%lattice(axis, :) - step*mol%lattice(column, :)
               call get_charges(model, displaced, error, minus)
               if (allocated(error)) error stop 'Charge minus strain failed'
               if (maxval(abs((plus-minus)/(2*step)-strain(axis, column, :))) > 1e-7_wp) then
                 error stop 'Strain charge difference failed'
               end if
            end do
         end do
      end do
   end do
   call get_charges(model, mol, error, charges(:2))
   if (.not.allocated(error)) error stop 'Invalid charge shape accepted'
   mol%xyz(:, 2) = mol%xyz(:, 1)
   call get_charges(model, mol, error, charges)
   if (.not.allocated(error)) error stop 'Coincident charge structure accepted'
   mol%xyz = xyz
   call get_charges(model, mol, error, charges)
   if (allocated(error)) error stop 'Charge error recovery failed'
end program test_multicharge