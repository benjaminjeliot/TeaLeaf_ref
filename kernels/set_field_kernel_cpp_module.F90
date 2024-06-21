! Copyright (c) 2024
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

module set_field_kernel_cpp_module

  interface

    subroutine set_field_kernel_cpp(x_min, x_max, y_min, y_max, halo_exchange_depth, energy0, energy1) &
      bind(C, name='set_field_kernel_cpp')
      use, intrinsic :: iso_c_binding, only: c_double, c_int32_t

      implicit none

      integer(c_int32_t), value :: x_min, x_max, y_min, y_max, halo_exchange_depth
      real(c_double), dimension( &
        x_min-halo_exchange_depth:x_max+halo_exchange_depth, &
        y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: energy0, energy1

    end subroutine set_field_kernel_cpp

  end interface

end module set_field_kernel_cpp_module
