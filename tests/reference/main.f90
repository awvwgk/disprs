program reference_tests
   use, intrinsic :: iso_fortran_env, only : error_unit
   use testdrive, only : run_testsuite
   use test_dftd3, only : collect_dftd3
   use test_param, only : collect_param
   use test_gcp, only : collect_gcp
   use test_pairwise, only : collect_pairwise
   use test_regression, only : collect_regression
   use test_hessian, only : collect_hessian
   use test_gcp_hessian, only : collect_gcp_hessian
   use test_partition, only : collect_partition
   use test_fourier, only : collect_fourier
   use test_periodic_1d, only : collect_periodic_1d
   use test_periodic_2d, only : collect_periodic_2d
   use test_periodic_3d, only : collect_periodic_3d
   use test_periodic_atm, only : collect_periodic_atm
   implicit none
   integer :: failures
   character(len=32) :: suite

   failures = 0
   call get_command_argument(1, suite)
   if (trim(suite) == 'regression') then
      call run_testsuite(collect_regression, error_unit, failures)
      if (failures /= 0) error stop 1
      stop
   end if
   if (trim(suite) == 'hessians') then
      call run_testsuite(collect_hessian, error_unit, failures)
      call run_testsuite(collect_gcp_hessian, error_unit, failures)
      if (failures /= 0) error stop 1
      stop
   end if
   if (trim(suite) == 'partitions') then
      call run_testsuite(collect_partition, error_unit, failures)
      if (failures /= 0) error stop 1
      stop
   end if
   if (trim(suite) == 'fourier') then
      call run_testsuite(collect_fourier, error_unit, failures)
      if (failures /= 0) error stop 1
      stop
   end if
   if (len_trim(suite) /= 0) error stop 'Unknown reference suite'
   call run_testsuite(collect_dftd3, error_unit, failures)
   call run_testsuite(collect_param, error_unit, failures)
   call run_testsuite(collect_gcp, error_unit, failures)
   call run_testsuite(collect_pairwise, error_unit, failures)
   call run_testsuite(collect_periodic_1d, error_unit, failures)
   call run_testsuite(collect_periodic_2d, error_unit, failures)
   call run_testsuite(collect_periodic_3d, error_unit, failures)
   call run_testsuite(collect_periodic_atm, error_unit, failures)
   if (failures /= 0) error stop 1
end program reference_tests