// Copyright (c) 2024
//
// This file is part of TeaLeaf.
//
// TeaLeaf is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// TeaLeaf is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// TeaLeaf. If not, see http://www.gnu.org/licenses/.

//>  @brief Fortran set field kernel.
//>  @author David Beckingsale, Wayne Gaudin
//>  @details Copies all of the final start of step filed data to the end of
//>  step data.

#include "set_field_kernel.h"

void set_field_kernel_cpp(int32_t x_min, int32_t x_max, int32_t y_min,
                          int32_t y_max, int32_t halo_exchange_depth,
                          double *energy0, double *energy1) {
  //! @param[in] x_min is the index of the first non-ghost cells in the x axis
  //! @param[in] x_max is the index of the last non-ghost cells in the x axis
  //! @param[in] y_min is the index of the first non-ghost cells in the y axis
  //! @param[in] y_max is the index of the last non-ghost cells in the y axis
  //! @param[in] halo_exchange_depth is the number of halo layers
  //! @param[in] energy0 is the energy array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)
  //! @param[out] energy1 is the energy array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)

  int n_x{x_max - x_min + 1};
  int n_y{y_max - y_min + 1};

#pragma omp parallel private(cell_idx)
#pragma omp for
  for (int k = 0; k < n_y; ++k) {
    for (int j = 0; j < n_x; ++j) {
      int cell_idx{(k + halo_exchange_depth) *
                       (x_max + 2 * halo_exchange_depth) +
                   j + halo_exchange_depth};
      energy1[cell_idx] = energy0[cell_idx];
    }
  }
}
