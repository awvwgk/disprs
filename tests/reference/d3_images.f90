program complete_images
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd3
   use dftd3_cutoff, only : get_lattice_points
   use dftd3_ncoord, only : get_coordination_number
   implicit none
   type(rational_damping_param) :: rational
   type(zero_damping_param) :: zero
   call new_rational_damping(rational, d3_param(s6=1.0_wp, s9=0.0_wp, alp=14.0_wp, &
      a1=0.4289_wp, s8=0.7875_wp, a2=4.4407_wp))
   call evaluate("acetic", rational)
   call new_rational_damping(rational, d3_param(s6=1.0_wp, s9=0.0_wp, alp=14.0_wp, &
      a1=0.4466_wp, s8=2.9491_wp, a2=6.1742_wp))
   call evaluate("adaman", rational)
   call new_zero_damping(zero, d3_param(s6=1.0_wp, s9=0.0_wp, alp=14.0_wp, &
      rs8=1.0_wp, rs6=1.581_wp, s8=0.0_wp))
   call evaluate("cyanamide", zero)
contains
   subroutine evaluate(name, param)
      character(len=*), intent(in) :: name
      class(damping_param), intent(in) :: param
      type(structure_type) :: mol
      type(d3_model) :: model
      type(realspace_cutoff) :: cutoff
      real(wp), allocatable :: trans(:, :), cn(:), weights(:, :), c6(:, :), energies(:)
      real(wp) :: original
      call get_structure(mol, "X23", name)
      call new_d3_model(model, mol)
      cutoff = realspace_cutoff(cn=30.0_wp, disp2=60.0_wp, disp3=15.0_wp)
      call get_dispersion(mol, model, param, cutoff, original)
      allocate(cn(mol%nat), weights(maxval(model%ref), mol%nat), c6(mol%nat, mol%nat), energies(mol%nat))
      call get_lattice_points(mol%periodic, mol%lattice, 120.0_wp, trans)
      call get_coordination_number(mol, trans, 30.0_wp, model%rcov, cn)
      call model%weight_references(mol, cn, weights)
      call model%get_atomic_c6(mol, weights, c6=c6)
      energies = 0.0_wp
      call param%get_dispersion2(mol, trans, 60.0_wp, 0.0_wp, model%rvdw, model%r4r2, c6, energy=energies)
      print '(a,2es26.17)', name, original, sum(energies)
   end subroutine
end program