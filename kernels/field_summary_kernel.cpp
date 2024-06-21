// Copyright (c) 2024
//
//  This file is part of TeaLeaf.
//
//  TeaLeaf is free software: you can redistribute it and/or modify it under
//  the terms of the GNU General Public License as published by the
//  Free Software Foundation, either version 3 of the License, or (at your
//  option) any later version.
//
//  TeaLeaf is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
//  more details.
//
//  You should have received a copy of the GNU General Public License along with
//  TeaLeaf. If not, see http://www.gnu.org/licenses/.

//>  @brief Fortran field summary kernel
//>  @author David Beckingsale, Wayne Gaudin
//>  @details The total mass, internal energy, temperature are calculated

#include "field_summary_kernel.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

void field_summary_kernel_cpp(int32_t x_min, int32_t x_max, int32_t y_min,
                              int32_t y_max, int32_t halo_exchange_depth,
                              double *volume, double *density, double *energy1,
                              double *u, double *vol, double *mass, double *ie,
                              double *temp) {
  //! @param[in] x_min is the index of the first non-ghost cells in the x axis
  //! @param[in] x_max is the index of the last non-ghost cells in the x axis
  //! @param[in] y_min is the index of the first non-ghost cells in the y axis
  //! @param[in] y_max is the index of the last non-ghost cells in the y axis
  //! @param[in] halo_exchange_depth is the number of halo layers
  //! @param[in] volume is the cell volume array with dimension
  //! [x_min-2:x_max+2, y_min-2:y_max+2] (two ghost halos)
  //! @param[in] density is the density array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)
  //! @param[in] energy1 is the energy array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)
  //! @param[in] u is the temperature array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)
  //! @param[out] vol the total cell volume
  //! @param[out] mass the total cell mass
  //! @param[out] ie the total cell internal energy
  //! @param[out] temp the total cell temperature

  // For debugging
  // write_field_summary_kernel_input_to_file(
  //     "field_summary_kernel_input.txt", x_min, x_max, y_min, y_max,
  //     halo_exchange_depth, volume, density, energy1, u);

  constexpr int n_volume_halos{2};

  double volume_sum{0.0};
  double mass_sum{0.0};
  double internal_energy_sum{0.0};
  double temperature_sum{0.0};

#pragma omp parallel private(cell_volume, cell_mass, volume_index, data_index)
#pragma omp parallel for reduction(+ : volume_sum, mass_sum, \
                                       internal_energy_sum, temperature_sum)
  for (int k{0}; k < y_max; ++k) {
    for (int j{0}; j < x_max; ++j) {
      int volume_index{j + n_volume_halos +
                       (k + n_volume_halos) * (x_max + 2 * n_volume_halos)};
      int data_index{j + halo_exchange_depth +
                     (k + halo_exchange_depth) *
                         (x_max + 2 * halo_exchange_depth)};
      double cell_volume = volume[volume_index];
      double cell_mass = cell_volume * density[data_index];
      volume_sum += cell_volume;
      mass_sum += cell_mass;
      internal_energy_sum += cell_mass * energy1[data_index];
      temperature_sum += cell_mass * u[data_index];
    }
  }

  *vol = volume_sum;
  *mass = mass_sum;
  *ie = internal_energy_sum;
  *temp = temperature_sum;
}

void write_field_summary_kernel_input_to_file(std::string filename, int x_min,
                                              int x_max, int y_min, int y_max,
                                              int halo_exchange_depth,
                                              double *volume, double *density,
                                              double *energy1, double *u) {
  //! @param[in] filename is the name of the file
  //! @param[in] x_min is the index of the first non-ghost cells in the x axis
  //! @param[in] x_max is the index of the last non-ghost cells in the x axis
  //! @param[in] y_min is the index of the first non-ghost cells in the y axis
  //! @param[in] y_max is the index of the last non-ghost cells in the y axis
  //! @param[in] halo_exchange_depth is the number of halo layers
  //! @param[in] volume is the cell volume array with dimension
  //! [x_min-2:x_max+2, y_min-2:y_max+2] (two ghost halos)
  //! @param[in] density is the density array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)
  //! @param[in] energy1 is the energy array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)
  //! @param[in] u is the temperature array with dimensions
  //! [x_min-halo_exchange_depth:x_max+halo_exchange_depth,
  //! y_min-halo_exchange_depth:y_max+halo_exchange_depth] (halo_exchange_depth
  //! ghost halo)

  std::ofstream out_file;
  out_file.open(filename);
  out_file << x_min << std::endl;
  out_file << x_max << std::endl;
  out_file << y_min << std::endl;
  out_file << y_max << std::endl;
  out_file << halo_exchange_depth << std::endl;

  // Write volume
  constexpr int n_volume_halos{2};
  int n_volume_x{x_max + 2 * n_volume_halos};
  int n_volume_y{y_max + 2 * n_volume_halos};
  out_file << n_volume_x << ' ' << n_volume_y << std::endl;
  for (int k{0}; k < n_volume_y; ++k) {
    for (int j{0}; j < n_volume_x; ++j) {
      int volume_idx{j + n_volume_x * k};
      out_file << std::setprecision(std::numeric_limits<double>::max_digits10)
               << volume[volume_idx] << ' ';
    }
  }
  out_file << std::endl;

  // Sizes
  int n_x{x_max + 2 * halo_exchange_depth};
  int n_y{y_max + 2 * halo_exchange_depth};

  // Write density
  out_file << n_x << ' ' << n_y << std::endl;
  for (int k{0}; k < n_y; ++k) {
    for (int j{0}; j < n_x; ++j) {
      int density_idx{j + n_x * k};
      out_file << std::setprecision(std::numeric_limits<double>::max_digits10)
               << density[density_idx] << ' ';
    }
  }
  out_file << std::endl;

  // Write energy
  out_file << n_x << ' ' << n_y << std::endl;
  for (int k{0}; k < n_y; ++k) {
    for (int j{0}; j < n_x; ++j) {
      int energy_idx{j + n_x * k};
      out_file << std::setprecision(std::numeric_limits<double>::max_digits10)
               << energy1[energy_idx] << ' ';
    }
  }
  out_file << std::endl;

  // Write temperature
  out_file << n_x << ' ' << n_y << std::endl;
  for (int k{0}; k < n_y; ++k) {
    for (int j{0}; j < n_x; ++j) {
      int temperature_idx{j + n_x * k};
      out_file << std::setprecision(std::numeric_limits<double>::max_digits10)
               << u[temperature_idx] << ' ';
    }
  }
  out_file << std::endl;

  out_file.close();
}
