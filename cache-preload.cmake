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

# Options
option(TEALEAF_ENABLE_OPENMP "Enable OpenMP" ON)
option(TEALEAF_ENABLE_TEST "Enable unit tests" ON)

# MPI
set(MPI_PREFIX_PATH
    ""
)

# Google Test
set(GTEST_PREFIX_PATH
    ""
)

# pFUnit
set(PFUNIT_PREFIX_PATH
    ""
)

set(CMAKE_PREFIX_PATH "${MPI_PREFIX_PATH}" "${GTEST_PREFIX_PATH}"
                      "${PFUNIT_PREFIX_PATH}")

set(CMAKE_PREFIX_PATH
    ${CMAKE_PREFIX_PATH}
    CACHE PATH "CMake prefix path")
