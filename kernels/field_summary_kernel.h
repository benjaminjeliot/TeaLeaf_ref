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

#ifndef FIELD_SUMMARY_KERNEL_H_
#define FIELD_SUMMARY_KERNEL_H_

#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

void field_summary_kernel_cpp(int x_min, int x_max, int y_min, int y_max,
                              int halo_exchange_depth, double *volume,
                              double *density, double *energy1, double *u,
                              double *vol, double *mass, double *ie,
                              double *temp);

#ifdef __cplusplus
}
#endif

void write_field_summary_kernel_input_to_file(std::string filename, int x_min,
                                              int x_max, int y_min, int y_max,
                                              int halo_exchange_depth,
                                              double *volume, double *density,
                                              double *energy1, double *u);

#endif  // FIELD_SUMMARY_KERNEL_H_
