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

#ifndef SET_FIELD_KERNEL_H_
#define SET_FIELD_KERNEL_H_

#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

void set_field_kernel_cpp(int32_t x_min, int32_t x_max, int32_t y_min,
                          int32_t y_max, int32_t halo_exchange_depth,
                          double *energy0, double *energy1);

#ifdef __cplusplus
}
#endif

#endif  // SET_FIELD_KERNEL_H_
