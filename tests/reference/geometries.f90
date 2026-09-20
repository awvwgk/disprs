program export_geometries
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   implicit none
   type(structure_type) :: mol
   integer :: record, atom
   character(len=2) :: identifier
   character(len=10), parameter :: crystals(*) = &
      [character(len=10) :: 'ammonia', 'acetic', 'adaman', 'anthracene']

   do record = 1, 21
      if (record <= 16) then
         write(identifier, '(i2.2)') record
         call get_structure(mol, 'MB16-43', identifier)
      else if (record == 17) then
         call get_structure(mol, 'UPU23', '0a')
      else
         call get_structure(mol, 'X23', trim(crystals(record - 17)))
      end if
      write(*, '(i0, 1x, es26.17)') mol%nat, mol%charge
      do atom = 1, mol%nat
         write(*, '(i0, 3(1x, es26.17))') mol%num(mol%id(atom)), mol%xyz(:, atom)
      end do
      if (record >= 18) write(*, '(9(1x, es26.17))') mol%lattice
   end do
end program export_geometries