program fortran_smoke
   use, intrinsic :: iso_c_binding
   use disprs
   implicit none
   integer(c_int), parameter :: numbers(2) = [6, 6]
   real(c_double), parameter :: positions(3, 2) = reshape([ &
      0.0_c_double, 0.0_c_double, 0.0_c_double, &
      6.0_c_double, 0.0_c_double, 0.0_c_double], [3, 2])
   real(c_double) :: energy, gradient(3, 2), virial(3, 3)
   real(c_double) :: reference_energy, reference_gradient(3, 2)
   real(c_double) :: dc6dr(3, 2, 2, 2), dc6dL(3, 3, 2, 2)
   real(c_double) :: coordination(2), charges(2), c6(2, 2), polarizabilities(2), pair2(2, 2), pair3(2, 2)
   real(c_double), parameter :: lattice(3, 3) = reshape([ &
      20.0_c_double, 0.0_c_double, 0.0_c_double, &
      0.0_c_double, 20.0_c_double, 0.0_c_double, &
      0.0_c_double, 0.0_c_double, 20.0_c_double], [3, 3])
   logical(c_bool), parameter :: periodic(3) = .true._c_bool
   logical(c_bool) :: directional(3)
   real(c_double) :: partial_lattice(3, 3)
   integer :: dimensions, kind, charge_kind
   integer(c_int) :: status

   call disprs_d3_energy(numbers, positions, "pbe", 1_c_int, .false._c_bool, &
      energy, gradient, status)
   if (status /= 0 .or. energy >= 0.0_c_double) error stop "D3 failed"

   call disprs_d3_energy([1_c_int, 28_c_int], positions, "blyp", 1_c_int, .false._c_bool, &
      reference_energy, reference_gradient, status)
   if (status /= 0) error stop "D3 reference failed"
   call disprs_d3_energy([1_c_int, 28_c_int], positions, "blyp", 1_c_int, .false._c_bool, &
      energy, gradient, status, d3s=.true._c_bool)
   if (status /= 0 .or. energy >= 0.0_c_double) error stop "D3S failed"
   if (abs(energy-reference_energy) < 1e-10_c_double) error stop "D3S model selection ignored"
   reference_energy = energy
   call disprs_d3_energy([1_c_int, 28_c_int], positions, "blyp", 1_c_int, .false._c_bool, &
      energy, gradient, status, model_kind=1_c_int)
   if (status /= 0 .or. abs(energy-reference_energy) > 1e-14_c_double) error stop "D3S kind differs from alias"
   call disprs_d3_energy(numbers, positions, "pbe", 1_c_int, .false._c_bool, &
      energy, gradient, status, model_kind=2_c_int)
   if (status == 0) error stop "Invalid D3 model kind accepted"
   call disprs_d3_energy(numbers, positions, "pbe", 1_c_int, .false._c_bool, &
      energy, gradient, status, model_kind=0_c_int, d3s=.true._c_bool)
   if (status == 0) error stop "Conflicting D3 selectors accepted"
   call disprs_d3_energy([1_c_int, 95_c_int], positions, "blyp", 1_c_int, .false._c_bool, &
      energy, gradient, status, d3s=.true._c_bool)
   if (status == 0) error stop "Unsupported D3S element accepted"

   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .false._c_bool, &
      energy, gradient, status)
   if (status /= 0 .or. energy >= 0.0_c_double) error stop "D4 failed"
   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .false._c_bool, &
      energy, gradient, status, fixed_charges=[0.2_c_double])
   if (status == 0) error stop "Invalid fixed D4 charge count accepted"
   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .false._c_bool, &
      energy, gradient, status, fixed_charges=[0.2_c_double, -0.3_c_double], charges=charges)
   if (status /= 0 .or. maxval(abs(charges-[0.2_c_double, -0.3_c_double])) > 0) error stop "Fixed D4 charges failed"
   dc6dr = 42.0_c_double
   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .false._c_bool, &
      energy, gradient, status, dc6dr=dc6dr(:, :1, :, :))
   if (status == 0 .or. any(abs(dc6dr-42.0_c_double) > 0)) error stop "D4 response shape guard failed"
   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .false._c_bool, &
      energy, gradient, status, dc6dr=dc6dr, dc6dL=dc6dL)
   if (status /= 0) error stop "D4 response recovery failed"
   if (maxval(abs(dc6dL(1, 1, :, :)-6*dc6dr(1, 2, :, :))) > 1e-12_c_double) then
     error stop "D4 two-atom response layout failed"
   end if
   reference_energy = energy
   reference_gradient = gradient
   call disprs_d4_energy(numbers, positions, "", 0.0_c_double, .false._c_bool, &
      energy, gradient, status, damping=[1.0_c_double, 0.95948085_c_double, &
      1.0_c_double, 0.38574991_c_double, 4.80688534_c_double, 16.0_c_double])
   if (status /= 0 .or. abs(energy-reference_energy) > 1e-15_c_double) error stop "explicit D4 failed"
   if (maxval(abs(gradient-reference_gradient)) > 1e-15_c_double) error stop "explicit D4 gradient failed"
   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .true._c_bool, &
      energy, gradient, status, coordination=coordination, charges=charges, &
      c6=c6, polarizabilities=polarizabilities, pair2=pair2, pair3=pair3)
   if (status /= 0) error stop "D4 outputs failed"
   if (abs(sum(pair2)+sum(pair3)-energy) > 1e-14_c_double) error stop "D4 pairwise sum failed"
   if (abs(sum(charges)) > 1e-14_c_double) error stop "D4 total charge failed"
   if (any(c6 <= 0) .or. any(polarizabilities <= 0)) error stop "D4 properties failed"
   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .true._c_bool, &
      energy, gradient, status, ga=2.0_c_double, gc=1.0_c_double, wf=5.0_c_double, &
      cutoff=[30.0_c_double, 8.0_c_double, 8.0_c_double, 4.0_c_double, 4.0_c_double], &
      pair2=pair2, pair3=pair3)
   if (status /= 0 .or. energy >= 0.0_c_double) error stop "custom D4 failed"
   if (abs(sum(pair2)+sum(pair3)-energy) > 1e-14_c_double) error stop "custom D4 pairwise failed"

   call disprs_d4_energy(numbers, positions, "pbe", 0.5_c_double, .true._c_bool, &
      energy, gradient, status, charge_model=1_c_int, charges=charges, pair2=pair2, pair3=pair3)
   if (status /= 0 .or. energy >= 0.0_c_double) error stop "EEQBC D4 failed"
   if (abs(sum(charges)-0.5_c_double) > 1e-14_c_double) error stop "EEQBC total charge failed"
   if (abs(sum(pair2)+sum(pair3)-energy) > 1e-14_c_double) error stop "EEQBC pairwise failed"

   call disprs_d4_energy(numbers, positions, "pbe", 0.0_c_double, .false._c_bool, &
      energy, gradient, status, lattice, periodic, virial, 1_c_int, charge_model=1_c_int)
   if (status /= 0 .or. energy >= 0.0_c_double) error stop "periodic D4S failed"
   do dimensions = 1, 2
      directional = [.true._c_bool, .false._c_bool, logical(dimensions == 2, c_bool)]
      partial_lattice = lattice
      partial_lattice(:, 2) = 0.0_c_double
      if (dimensions == 1) partial_lattice(:, 3) = 0.0_c_double
      do kind = 0, 1
         do charge_kind = 0, 1
            call disprs_d4_energy(numbers, positions, "pbe", 0.5_c_double, .true._c_bool, &
               energy, gradient, status, lattice=partial_lattice, periodic=directional, virial=virial, &
               model_kind=int(kind, c_int), charge_model=int(charge_kind, c_int), &
               charge_cutoff=80.0_c_double, charges=charges, pair2=pair2, pair3=pair3)
            if (status /= 0 .or. energy >= 0.0_c_double) error stop "directional D4 failed"
            if (abs(sum(charges)-0.5_c_double) > 1e-13_c_double) error stop "directional D4 charges failed"
            if (abs(sum(pair2)+sum(pair3)-energy) > 1e-13_c_double) error stop "directional D4 pairs failed"
         end do
      end do
   end do
   block
      integer(c_int), parameter :: atoms(4) = [6, 8, 7, 1]
      real(c_double), parameter :: xyz(3, 4) = reshape([ &
         0.2_c_double, 0.3_c_double, 0.4_c_double, 2.6_c_double, 0.7_c_double, 0.8_c_double, &
         1.2_c_double, 2.8_c_double, 0.6_c_double, 0.7_c_double, 1.1_c_double, 2.3_c_double], [3, 4])
      real(c_double), parameter :: cell(3, 3) = reshape([ &
         5.3_c_double, 0.0_c_double, 0.0_c_double, 0.6_c_double, 5.7_c_double, 0.0_c_double, &
         0.3_c_double, 0.4_c_double, 6.2_c_double], [3, 3])
      real(c_double), parameter :: ranges(*) = [6.0_c_double, 7.0_c_double, 6.0_c_double, 1.0_c_double, 1.0_c_double]
      real(c_double) :: serial, total, grad(3, 4), serial_grad(3, 4), total_grad(3, 4)
      real(c_double) :: sigma(3, 3), serial_sigma(3, 3), total_sigma(3, 3), pairs2(4, 4), pairs3(4, 4)
      real(c_double) :: partial_charges(4), serial_charges(4)
      integer(c_int) :: part
      do dimensions = 0, 3
         directional = [dimensions >= 1, dimensions >= 2, dimensions >= 3]
         do kind = 0, 1
            do charge_kind = 0, 1
               call disprs_d4_energy(atoms, xyz, "pbe", 0.5_c_double, .true._c_bool, &
                  serial, serial_grad, status, lattice=cell, periodic=directional, virial=serial_sigma, &
                  model_kind=int(kind, c_int), charge_model=int(charge_kind, c_int), cutoff=ranges, &
                  charge_cutoff=10.0_c_double, ghost=[4_c_int], charges=serial_charges)
               if (status /= 0) error stop "D4 ghost evaluation failed"
               total = 0.0_c_double; total_grad = 0.0_c_double; total_sigma = 0.0_c_double
               do part = 0, 2
                  call disprs_d4_energy(atoms, xyz, "pbe", 0.5_c_double, .true._c_bool, &
                     energy, grad, status, lattice=cell, periodic=directional, virial=sigma, &
                     model_kind=int(kind, c_int), charge_model=int(charge_kind, c_int), cutoff=ranges, &
                     charge_cutoff=10.0_c_double, ghost=[4_c_int], partition=[part, 3_c_int], &
                     charges=partial_charges, pair2=pairs2, pair3=pairs3)
                  if (status /= 0) error stop "D4 partition evaluation failed"
                  if (abs(sum(pairs2)+sum(pairs3)-energy) > 1e-13_c_double) error stop "D4 partition pairs failed"
                    if (maxval(abs(pairs2(4, :))) > tiny(energy) .or. &
                       maxval(abs(pairs3(4, :))) > tiny(energy) .or. &
                       maxval(abs(pairs2(:, 4))) > tiny(energy) .or. &
                       maxval(abs(pairs3(:, 4))) > tiny(energy)) error stop "D4 ghost pairs failed"
                  if (maxval(abs(partial_charges-serial_charges)) > 1e-14_c_double) error stop "D4 ghost charges failed"
                  total = total + energy; total_grad = total_grad + grad; total_sigma = total_sigma + sigma
               end do
               if (abs(total-serial) > 1e-13_c_double) error stop "D4 partition energy failed"
               if (maxval(abs(total_grad-serial_grad)) > 1e-12_c_double) error stop "D4 partition gradient failed"
               if (maxval(abs(total_sigma-serial_sigma)) > 1e-12_c_double) error stop "D4 partition virial failed"
            end do
         end do
      end do
      call disprs_d4_energy(atoms, xyz, "pbe", 0.0_c_double, .true._c_bool, &
         energy, grad, status, ghost=[0_c_int])
      if (status == 0) error stop "invalid D4 ghost accepted"
      call disprs_d4_energy(atoms, xyz, "pbe", 0.0_c_double, .true._c_bool, &
         energy, grad, status, partition=[2_c_int, 2_c_int])
      if (status == 0) error stop "invalid D4 partition accepted"
   end block
end program