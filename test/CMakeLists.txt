# Copyright (c) 2024
#
# This file is part of TeaLeaf.
#
# TeaLeaf is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TeaLeaf is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# TeaLeaf. If not, see http://www.gnu.org/licenses/.

add_subdirectory(data)

add_pfunit_ctest(
  field_summary_kernel_test
  TEST_SOURCES
  field_summary_kernel.pf
  REGISTRY
  field_summary_kernel.inc
  LINK_LIBRARIES
  tea_leaf_kernels_f)

add_pfunit_ctest(
  field_summary_kernel_cpp_test
  TEST_SOURCES
  field_summary_kernel_cpp.pf
  REGISTRY
  field_summary_kernel_cpp.inc
  LINK_LIBRARIES
  tea_leaf_kernels_cpp)

add_pfunit_ctest(
  set_field_kernel_test
  TEST_SOURCES
  set_field_kernel.pf
  REGISTRY
  set_field_kernel.inc
  LINK_LIBRARIES
  tea_leaf_kernels_f)

add_pfunit_ctest(
  set_field_kernel_cpp_test
  TEST_SOURCES
  set_field_kernel_cpp.pf
  REGISTRY
  set_field_kernel_cpp.inc
  LINK_LIBRARIES
  tea_leaf_kernels_cpp)
