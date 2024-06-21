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

cmake_minimum_required(VERSION 3.20)

project(
  "tealeaf"
  DESCRIPTION
    "TeaLeaf is a mini-app that solves the linear heat conduction equation"
  HOMEPAGE_URL "https://github.com/UK-MAC/TeaLeaf"
  LANGUAGES Fortran CXX C)

# Options
option(TEALEAF_ENABLE_OPENMP "Enable OpenMP" ON)
option(TEALEAF_ENABLE_TEST "Enable unit tests" OFF)

find_package(MPI REQUIRED)

if(TEALEAF_ENABLE_TEST)
  enable_testing()
  find_package(PFUNIT REQUIRED)
endif()

if(TEALEAF_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Set compiler flags specific to the compiler
set(TEA_LEAF_FORTRAN_COMPILER_FLAGS)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND TEA_LEAF_FORTRAN_COMPILER_FLAGS -fallow-argument-mismatch)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  # TODO: Add Intel OneAPI compiler flags
endif()

add_library(tea_leaf_timer timer_c.c)

add_library(tea_leaf_core data.f90 definitions.f90)

add_subdirectory(kernels)

add_executable(
  tea_leaf
  build_field.f90
  calc_dt.f90
  diffuse.f90
  field_summary.f90
  generate_chunk.f90
  global_mpi.f90
  initialise_chunk.f90
  initialise.f90
  pack.f90
  parse.f90
  read_input.f90
  report.f90
  set_field.f90
  start.f90
  tea.f90
  tea_leaf_cg.f90
  tea_leaf_cheby.f90
  tea_leaf_common.f90
  tea_leaf.f90
  tea_leaf_jacobi.f90
  tea_leaf_ppcg.f90
  tea_solve.f90
  timer.f90
  timestep.f90
  update_halo.f90
  visit.f90)

if(TEALEAF_ENABLE_OPENMP)
  target_compile_definitions(tea_leaf PRIVATE WITH_OMP4)
endif()

target_compile_options(tea_leaf PRIVATE ${TEA_LEAF_FORTRAN_COMPILER_FLAGS})

target_link_libraries(
  tea_leaf PRIVATE MPI::MPI_Fortran tea_leaf_core tea_leaf_timer
                   tea_leaf_kernels_f tea_leaf_kernels_cpp)

target_include_directories(tea_leaf PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if(TEALEAF_ENABLE_TEST)
  add_subdirectory(test)
endif()
