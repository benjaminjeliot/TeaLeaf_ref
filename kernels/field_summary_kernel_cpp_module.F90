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


module field_summary_kernel_cpp_module

  interface

    subroutine field_summary_kernel_cpp(x_min, x_max, y_min, y_max, halo_exchange_depth, volume, density, energy1, u, vol, &
      mass, ie, temp) &
      bind(C, name='field_summary_kernel_cpp')
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_int32_t
      implicit none
      integer (c_int32_t), intent(in), value :: x_min, x_max, y_min, y_max, halo_exchange_depth
      real(c_double), intent(in), dimension(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
      real(c_double), intent(in), dimension( &
        x_min-halo_exchange_depth:x_max+halo_exchange_depth, &
        y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: &
        density,energy1,u
      REAL(c_double), intent(out) :: vol,mass,ie,temp

    end subroutine field_summary_kernel_cpp

  end interface

end module field_summary_kernel_cpp_module
