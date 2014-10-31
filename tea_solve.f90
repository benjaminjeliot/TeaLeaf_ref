!Crown Copyright 2014 AWE.
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

!>  @brief Driver for the heat conduction kernel
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the user specified kernel for the heat conduction

MODULE tea_leaf_module

  USE report_module
  USE data_module
  USE tea_leaf_kernel_module
  USE tea_leaf_kernel_cg_module
  USE tea_leaf_kernel_cheby_module
  USE update_halo_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf()

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_THREAD_NUM
  INTEGER :: c, n
  REAL(KIND=8) :: ry,rx, error

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: kernel_time,timer

  ! For CG solver
  REAL(KIND=8) :: rro, pw, rrn, alpha, beta

  ! For chebyshev solver
  REAL(KIND=8), DIMENSION(max_iters) :: cg_alphas, cg_betas
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta
  REAL(KIND=8) :: it_alpha, cn, gamm, bb
  INTEGER :: est_itc, cheby_calc_steps, max_cheby_iters, info, switch_step
  LOGICAL :: ch_switch_check

  INTEGER :: cg_calc_steps
  REAL(KIND=8) :: cg_time, ch_time, solve_timer
  cg_time = 0.0_8
  ch_time = 0.0_8
  cg_calc_steps = 0

  IF(coefficient .nE. RECIP_CONDUCTIVITY .and. coefficient .ne. conductivity) THEN
    CALL report_error('tea_leaf', 'unknown coefficient option')
  endif

  error = 1e10
  cheby_calc_steps = 0

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      fields=0
      fields(FIELD_ENERGY1) = 1
      fields(FIELD_DENSITY1) = 1
      CALL update_halo(fields,2)

      ! INIT
      IF(profiler_on) kernel_time=timer()

      if (use_fortran_kernels .or. use_c_kernels) then
        rx = dt/(chunks(c)%field%celldx(chunks(c)%field%x_min)**2)
        ry = dt/(chunks(c)%field%celldy(chunks(c)%field%y_min)**2)
      endif

      IF(tl_use_cg .or. tl_use_chebyshev) then
        IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_init_cg_fortran(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%u,                           &
              chunks(c)%field%work_array1,                 & ! p
              chunks(c)%field%work_array2,                 & ! r
              chunks(c)%field%work_array3,                 & ! Mi
              chunks(c)%field%work_array4,                 & ! w
              chunks(c)%field%work_array5,                 & ! z
              chunks(c)%field%work_array6,                 & ! Kx
              chunks(c)%field%work_array7,                 & ! Ky
              rx, ry, rro, coefficient)
        ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_init_cg_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%u,                           &
              chunks(c)%field%work_array1,                 & ! p
              chunks(c)%field%work_array2,                 & ! r
              chunks(c)%field%work_array3,                 & ! Mi
              chunks(c)%field%work_array4,                 & ! w
              chunks(c)%field%work_array5,                 & ! z
              chunks(c)%field%work_array6,                 & ! Kx
              chunks(c)%field%work_array7,                 & ! Ky
              rx, ry, rro, coefficient)
        ENDIF

        ! need to update p at this stage
        fields=0
        fields(FIELD_U) = 1
        fields(FIELD_P) = 1
        CALL update_halo(fields,1)

        ! and globally sum rro
        call tea_allsum(rro)
      ELSEIF(tl_use_jacobi) then
        IF (use_fortran_kernels) THEN
          CALL tea_leaf_kernel_init(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%celldx,                      &
              chunks(c)%field%celldy,                      &
              chunks(c)%field%volume,                      &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%work_array1,                 & !u0
              chunks(c)%field%u,                           & !u1
              chunks(c)%field%work_array2,                 & !un
              chunks(c)%field%work_array4,                 & !Kx temp
              chunks(c)%field%work_array5,                 & !Ky temp
              chunks(c)%field%work_array6,                 & !Kx
              chunks(c)%field%work_array7,                 & !Ky
              coefficient)
        ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_init_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                         &
              chunks(c)%field%y_min,                         &
              chunks(c)%field%y_max,                         &
              chunks(c)%field%celldx,                        &
              chunks(c)%field%celldy,                        &
              chunks(c)%field%volume,                        &
              chunks(c)%field%density1,                      &
              chunks(c)%field%energy1,                       &
              chunks(c)%field%work_array1,                   & !u0
              chunks(c)%field%u,                             & !u1
              chunks(c)%field%work_array2,                   & !un
              chunks(c)%field%work_array4,                   & !Kx temp
              chunks(c)%field%work_array5,                   & !Ky temp
              chunks(c)%field%work_array6,                   & !Kx
              chunks(c)%field%work_array7,                   & !Ky
              coefficient)
        ENDIF

      ENDIF

      fields=0
      fields(FIELD_U) = 1

      ! need the original value of u
      if(tl_use_chebyshev) then
        IF(use_fortran_kernels) then
          call tea_leaf_kernel_cheby_copy_u(chunks(c)%field%x_min,&
            chunks(c)%field%x_max,                                &
            chunks(c)%field%y_min,                                &
            chunks(c)%field%y_max,                                &
            chunks(c)%field%u0,                                   &
            chunks(c)%field%u)
        endif
      endif

      DO n=1,max_iters

        if (profile_solver) solve_timer=timer()

        IF (tl_ch_cg_errswitch) then
            ! either the error has got below tolerance, or it's already going
            ch_switch_check = (cheby_calc_steps .gt. 0) .or. (error .le. tl_ch_cg_epslim)
        ELSE
            ! enough steps have passed
            ch_switch_check = n .ge. tl_ch_cg_presteps
        ENDIF

        IF (tl_use_chebyshev .and. ch_switch_check) then
          ! don't need to update p any more
          fields(FIELD_P) = 0

          ! on the first chebyshev steps, find the eigenvalues, coefficients,
          ! and expected number of iterations
          IF (cheby_calc_steps .eq. 0) then
            ! maximum number of iterations in chebyshev solver
            max_cheby_iters = max_iters - n + 2
            rro = error

            ! calculate eigenvalues
            call tea_calc_eigenvalues(cg_alphas, cg_betas, eigmin, eigmax, &
                max_iters, n-1, info)

            ! maximum number of iterations in chebyshev solver
            max_cheby_iters = max_iters - n + 2

            ! calculate chebyshev coefficients
            call tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                theta, max_cheby_iters)

            ! calculate 2 norm of u0
            IF(use_fortran_kernels) THEN
              call tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,&
                    chunks(c)%field%x_max,                          &
                    chunks(c)%field%y_min,                          &
                    chunks(c)%field%y_max,                          &
                    chunks(c)%field%u0,                             &
                    bb)
            ENDIF

            call tea_allsum(bb)

            ! initialise 'p' array
            IF(use_fortran_kernels) THEN
              call tea_leaf_kernel_cheby_init(chunks(c)%field%x_min,&
                    chunks(c)%field%x_max,                          &
                    chunks(c)%field%y_min,                          &
                    chunks(c)%field%y_max,                          &
                    chunks(c)%field%u,                              &
                    chunks(c)%field%u0,                             &
                    chunks(c)%field%work_array1,                    & ! p
                    chunks(c)%field%work_array2,                    & ! r
                    chunks(c)%field%work_array3,                    & ! Mi
                    chunks(c)%field%work_array4,                    & ! w
                    chunks(c)%field%work_array5,                    & ! z
                    chunks(c)%field%work_array6,                    & ! Kx
                    chunks(c)%field%work_array7,                    & ! Ky
                    ch_alphas, ch_betas, max_cheby_iters,           &
                    rx, ry, theta, error)
            ENDIF

            CALL update_halo(fields,1)

            IF(use_fortran_kernels) THEN
                call tea_leaf_kernel_cheby_iterate(chunks(c)%field%x_min,&
                    chunks(c)%field%x_max,                               &
                    chunks(c)%field%y_min,                               &
                    chunks(c)%field%y_max,                               &
                    chunks(c)%field%u,                                   &
                    chunks(c)%field%u0,                                  &
                    chunks(c)%field%work_array1,                         & ! p
                    chunks(c)%field%work_array2,                         & ! r
                    chunks(c)%field%work_array3,                         & ! Mi
                    chunks(c)%field%work_array4,                         & ! w
                    chunks(c)%field%work_array5,                         & ! z
                    chunks(c)%field%work_array6,                         & ! Kx
                    chunks(c)%field%work_array7,                         & ! Ky
                    ch_alphas, ch_betas, max_cheby_iters,                &
                    rx, ry, cheby_calc_steps)
            ENDIF

            IF(use_fortran_kernels) THEN
              call tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,&
                    chunks(c)%field%x_max,                          &
                    chunks(c)%field%y_min,                          &
                    chunks(c)%field%y_max,                          &
                    chunks(c)%field%work_array2,                    &
                    error)
            ENDIF

            call tea_allsum(error)

            it_alpha = eps*bb/(4.0_8*error)
            cn = eigmax/eigmin
            gamm = (sqrt(cn) - 1.0_8)/(sqrt(cn) + 1.0_8)
            est_itc = nint(log(it_alpha)/(2.0_8*log(gamm)))

            ! FIXME still not giving correct answer, but multiply by 2.5 does
            ! an 'okay' job for now
            est_itc = est_itc * 2.5

            if (parallel%boss) then
              write(g_out,'(a,i3,a,e15.7)') "Switching after ",n," steps, error ",rro
              write(g_out,'(5a11)')"eigmin", "eigmax", "cn", "error", "est itc"
              write(g_out,'(2f11.8,2e11.4,11i11)')eigmin, eigmax, cn, error, est_itc
              write(0,'(a,i3,a,e15.7)') "Switching after ",n," steps, error ",rro
              write(0,'(5a11)')"eigmin", "eigmax", "cn", "error", "est itc"
              write(0,'(2f11.8,2e11.4,11i11)')eigmin, eigmax, cn, error, est_itc
            endif

            if (info .ne. 0) then
              CALL report_error('tea_leaf', 'Error in calculating eigenvalues')
            endif

            switch_step = n
            cheby_calc_steps = 2
          else
            IF(use_fortran_kernels) THEN
                call tea_leaf_kernel_cheby_iterate(chunks(c)%field%x_min,&
                    chunks(c)%field%x_max,                               &
                    chunks(c)%field%y_min,                               &
                    chunks(c)%field%y_max,                               &
                    chunks(c)%field%u,                                   &
                    chunks(c)%field%u0,                                  &
                    chunks(c)%field%work_array1,                         & ! p
                    chunks(c)%field%work_array2,                         & ! r
                    chunks(c)%field%work_array3,                         & ! Mi
                    chunks(c)%field%work_array4,                         & ! w
                    chunks(c)%field%work_array5,                         & ! z
                    chunks(c)%field%work_array6,                         & ! Kx
                    chunks(c)%field%work_array7,                         & ! Ky
                    ch_alphas, ch_betas, max_cheby_iters,                &
                    rx, ry, cheby_calc_steps)
            ENDIF

            ! after estimated number of iterations has passed, calc resid
            ! Leaving 10 iterations between each global reduction won't affect
            ! total time spent much if at all (number of steps spent in
            ! chebyshev is typically O(300+)) but will greatyl reduce global
            ! synchronisations needed
            if ((n-switch_step .ge. est_itc) .and. (mod(n, 10) .eq. 0)) then
              IF(use_fortran_kernels) THEN
                call tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,&
                      chunks(c)%field%x_max,                          &
                      chunks(c)%field%y_min,                          &
                      chunks(c)%field%y_max,                          &
                      chunks(c)%field%work_array2,                    &
                      error)
              ENDIF

              call tea_allsum(error)
            else
              ! dummy to make it go smaller every time but not reach tolerance
              error = 1.0_8/(cheby_calc_steps)
            endif
          endif

          cheby_calc_steps = cheby_calc_steps + 1

        ELSEIF(tl_use_cg .or. tl_use_chebyshev) then
          fields(FIELD_P) = 1
          cg_calc_steps = cg_calc_steps + 1

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_w(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                         &
                chunks(c)%field%y_min,                                         &
                chunks(c)%field%y_max,                                         &
                chunks(c)%field%work_array1,                                   & ! p
                chunks(c)%field%work_array4,                                   & ! w
                chunks(c)%field%work_array6,                                   & ! Kx
                chunks(c)%field%work_array7,                                   & ! Ky
                rx, ry, pw)
          ELSEIF(use_c_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_c_calc_w(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                   &
                chunks(c)%field%y_min,                                   &
                chunks(c)%field%y_max,                                   &
                chunks(c)%field%work_array1,                             &
                chunks(c)%field%work_array4,                             &
                chunks(c)%field%work_array6,                             &
                chunks(c)%field%work_array7,                             &
                rx, ry, pw)
          ENDIF

          CALL tea_allsum(pw)
          alpha = rro/pw
          if(tl_use_chebyshev) cg_alphas(n) = alpha

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_ur(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                          &
                chunks(c)%field%y_min,                                          &
                chunks(c)%field%y_max,                                          &
                chunks(c)%field%u,                                              &
                chunks(c)%field%work_array1,                                    & ! p
                chunks(c)%field%work_array2,                                    & ! r
                chunks(c)%field%work_array3,                                    & ! Mi
                chunks(c)%field%work_array4,                                    & ! w
                chunks(c)%field%work_array5,                                    & ! z
                alpha, rrn)
          ELSEIF(use_c_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_c_calc_ur(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                    &
                chunks(c)%field%y_min,                                    &
                chunks(c)%field%y_max,                                    &
                chunks(c)%field%u,                                        &
                chunks(c)%field%work_array1,                              & ! p
                chunks(c)%field%work_array2,                              & ! r
                chunks(c)%field%work_array3,                              & ! Mi
                chunks(c)%field%work_array4,                              & ! w
                chunks(c)%field%work_array5,                              & ! z
                alpha, rrn)
          ENDIF

          CALL tea_allsum(rrn)
          beta = rrn/rro
          if(tl_use_chebyshev) cg_betas(n) = beta

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_p(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                         &
                chunks(c)%field%y_min,                                         &
                chunks(c)%field%y_max,                                         &
                chunks(c)%field%work_array1,                                   & ! p
                chunks(c)%field%work_array2,                                   & ! r
                chunks(c)%field%work_array5,                                   & ! z
                beta)
          ELSEIF(use_c_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_c_calc_p(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                   &
                chunks(c)%field%y_min,                                   &
                chunks(c)%field%y_max,                                   &
                chunks(c)%field%work_array1,                             & ! p
                chunks(c)%field%work_array2,                             & ! r
                chunks(c)%field%work_array5,                             & ! z
                beta)
          ENDIF

          error = rrn
          rro = rrn

          CALL tea_allsum(error)

        ELSEIF(tl_use_jacobi) then
          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                rx,                                          &
                ry,                                          &
                chunks(c)%field%work_array6,                 & ! Kx
                chunks(c)%field%work_array7,                 & ! Ky
                error,                                       &
                chunks(c)%field%work_array1,                 & ! u0
                chunks(c)%field%u,                           & ! u1
                chunks(c)%field%work_array2)                   ! un
          ELSEIF(use_C_kernels) THEN
            CALL tea_leaf_kernel_solve_c(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                         &
                chunks(c)%field%y_min,                         &
                chunks(c)%field%y_max,                         &
                rx,                                            &
                ry,                                            &
                chunks(c)%field%work_array6,                   & ! Kx
                chunks(c)%field%work_array7,                   & ! Ky
                error,                                         &
                chunks(c)%field%work_array1,                   & ! u0
                chunks(c)%field%u,                             & ! u1  
                chunks(c)%field%work_array2)                     ! un
          ENDIF

          CALL tea_max(error)
        ENDIF

        ! updates u and possibly p
        CALL update_halo(fields,1)

        if (profile_solver) then
          IF (tl_use_chebyshev .and. ch_switch_check) then
              ch_time=ch_time+(timer()-solve_timer)
          else
              cg_time=cg_time+(timer()-solve_timer)
          endif
        endif

        IF (abs(error) .LT. eps) EXIT

      ENDDO

      IF (parallel%boss) THEN
!$      IF(OMP_GET_THREAD_NUM().EQ.0) THEN
          WRITE(g_out,"('Conduction error ',e14.7)") error
          WRITE(g_out,"('Iteration count ',i8)") n-1
          WRITE(0,"('Conduction error ',e14.7)") error
          WRITE(0,"('Iteration count ', i8)") n-1
!$      ENDIF
      ENDIF

      ! RESET
      IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_finalise(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density1,                        &
              chunks(c)%field%u)
      ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_finalise_c(chunks(c)%field%x_min,&
              chunks(c)%field%x_max,                            &
              chunks(c)%field%y_min,                            &
              chunks(c)%field%y_max,                            &
              chunks(c)%field%energy1,                          &
              chunks(c)%field%density1,                         &
              chunks(c)%field%u)
      ENDIF

      fields=0
      fields(FIELD_ENERGY1) = 1
      CALL update_halo(fields,1)

    ENDIF

  ENDDO
  IF(profile_solver) profiler%tea=profiler%tea+(timer()-kernel_time)

  IF (profile_solver .and. tl_use_chebyshev) THEN
    call tea_sum(ch_time)
    call tea_sum(cg_time)
  endif
  IF (profile_solver .and. parallel%boss .and. tl_use_chebyshev) THEN
    write(0, "(a3, a16, a7, a16, a7)") "", "Time", "Steps", "Per it", "Ratio"
    write(0, "(a3, f16.10, i7, f16.10, f7.2)") "CG", cg_time + 0.0_8, cg_calc_steps, &
        merge(cg_time/cg_calc_steps, 0.0_8, cg_calc_steps .gt. 0), 1.0_8
    write(0, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
        merge(ch_time/cheby_calc_steps, 0.0_8, cheby_calc_steps .gt. 0), &
        merge((ch_time/cheby_calc_steps)/(cg_time/cg_calc_steps), 0.0_8, cheby_calc_steps .gt. 0)
    write(0, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
        cheby_calc_steps, cheby_calc_steps-est_itc

    write(g_out, "(a3, a16, a7, a16, a7)") "", "Time", "Steps", "Per it", "Ratio"
    write(g_out, "(a3, f16.10, i7, f16.10, f7.2)") "CG", cg_time + 0.0_8, cg_calc_steps, &
        merge(cg_time/cg_calc_steps, 0.0_8, cg_calc_steps .gt. 0), 1.0_8
    write(g_out, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
        merge(ch_time/cheby_calc_steps, 0.0_8, cheby_calc_steps .gt. 0), &
        merge((ch_time/cheby_calc_steps)/(cg_time/cg_calc_steps), 0.0_8, cheby_calc_steps .gt. 0)
    write(g_out, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
        cheby_calc_steps, cheby_calc_steps-est_itc
  endif

END SUBROUTINE tea_leaf

END MODULE tea_leaf_module
