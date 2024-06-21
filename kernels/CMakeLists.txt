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

# Fortran kernels
add_library(
  tea_leaf_kernels_f
  calc_dt_kernel.f90
  field_summary_kernel.f90
  generate_chunk_kernel.f90
  initialise_chunk_kernel.f90
  pack_kernel.f90
  set_field_kernel.f90
  tea_leaf_cg_kernel.f90
  tea_leaf_cheby_kernel.f90
  tea_leaf_common_kernel.f90
  tea_leaf_jacobi_kernel.f90
  tea_leaf_ppcg_kernel.f90
  update_halo_kernel.f90
  update_internal_halo_kernel.f90)

target_link_libraries(tea_leaf_kernels_f PRIVATE tea_leaf_core)

if(TEALEAF_ENABLE_OPENMP)
  target_compile_definitions(tea_leaf_kernels_f PRIVATE WITH_OMP4)
endif()

target_link_libraries(tea_leaf_kernels_f PRIVATE tea_leaf_core)

target_include_directories(tea_leaf_kernels_f
                           PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# C++ kernels
add_library(
  tea_leaf_kernels_cpp
  field_summary_kernel.cpp field_summary_kernel_cpp_module.F90
  set_field_kernel.cpp set_field_kernel_cpp_module.F90)

if(TEALEAF_ENABLE_OPENMP)
  target_compile_definitions(tea_leaf_kernels_cpp PRIVATE WITH_OMP4)
endif()

target_link_libraries(tea_leaf_kernels_cpp PRIVATE tea_leaf_core)

target_include_directories(tea_leaf_kernels_cpp
                           PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
