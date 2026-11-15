module output_screen
use parameters_module
use mpi_module
use storage_module
use profiling_module, only: StartProfRange, EndProfRange
implicit none

private
public :: screen_output, compute_average_iteration_time, screen_elapse_time

contains

subroutine screen_output

        use, intrinsic :: ISO_FORTRAN_ENV, only: output_unit
        implicit none
        real(rp), parameter :: DivCheckToll = 1.0E-14_rp
        real(rp)            :: lcl_max(4), gbl_max(4)
        real(rp)            :: lcl_min(3), gbl_min(3)
        real(rp)            :: r_mx, p_mx, T_mx, M2_mx
        real(rp)            :: r_mn, p_mn, T_mn
        integer             :: i, j, k, k_lo, k_hi, err

        call StartProfRange('screen_output')

        err = 0

        r_mx  = 0.0_rp
        p_mx  = 0.0_rp
        T_mx  = 0.0_rp
        M2_mx = 0.0_rp
        r_mn  = huge(0.0_rp)
        p_mn  = huge(0.0_rp)
        T_mn  = huge(0.0_rp)

        if(dims == 2) then
          k_lo = 1
          k_hi = 1
        else
          k_lo = sz
          k_hi = ez
        endif

        !$acc parallel default(present) copy(r_mx,p_mx,T_mx,M2_mx,r_mn,p_mn,T_mn)
        !$acc loop gang, vector collapse(3) &
        !$acc   reduction(max:r_mx,p_mx,T_mx,M2_mx) &
        !$acc   reduction(min:r_mn,p_mn,T_mn)
        do k       = k_lo, k_hi
           do j    = sy, ey
              do i = sx, ex
                 r_mx  = max(phi(i,j,k,1), r_mx)
                 p_mx  = max(P(i,j,k),     p_mx)
                 T_mx  = max(T(i,j,k),     T_mx)
                 M2_mx = max((U(i,j,k)*U(i,j,k) + V(i,j,k)*V(i,j,k) + W(i,j,k)*W(i,j,k)) &
                                 /(GAMMA_T(i,j,k)*T(i,j,k)), M2_mx)

                 r_mn = min(phi(i,j,k,1), r_mn)
                 p_mn = min(P(i,j,k),     p_mn)
                 T_mn = min(T(i,j,k),     T_mn)
              enddo
           enddo
        enddo
        !$acc end parallel

        lcl_max(1) = r_mx
        lcl_max(2) = p_mx
        lcl_max(3) = T_mx
        lcl_max(4) = M2_mx
        lcl_min(1) = r_mn
        lcl_min(2) = p_mn
        lcl_min(3) = T_mn

        call MPI_allreduce(lcl_max, gbl_max, 4, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_min, gbl_min, 3, MPI_RP, MPI_MIN, mpi_comm_cart, err)

        gbl_max(4) = sqrt(gbl_max(4))

        if(rank == root) then
          if(mod(it,10*itout)==0) then
            write(*,'(A)') ' '//repeat('=',142)
            write(*,'(1x,A8,1x,A7,9(1x,A13))') &
                    '      it', ' time %', &
                    '           dt', '        M max', &
                    '        r max', '        r min', &
                    '        T max', '        T min', &
                    '        p max', '        p min', &
                    '        titer'
            write(*,'(A)') ' '//repeat('-',142)
          endif

          write(*,10) it , time/tmax*100 , dt, gbl_max(4),           &
                                               gbl_max(1), gbl_min(1), &
                                               gbl_max(3), gbl_min(3), &
                                               gbl_max(2), gbl_min(2), &
                                               end_it_time

          flush(output_unit)
        endif

        if(gbl_min(1) < DivCheckToll .or. &
           gbl_min(3) < DivCheckToll .or. &
           dt         < DivCheckToll) then
           if(rank == root) then
             print*, ' Divergency has been detected!'
             flush(output_unit)
           endif
           istop = 1
        endif

        10 format(1x,I8,1x,F7.3,9(1x,ES13.5))

        call EndProfRange

        return
end subroutine screen_output

subroutine compute_average_iteration_time
        implicit none

        end_it_time = MPI_WTIME()
        end_it_time = (end_it_time - str_it_time)/real(itOut,rp)
        str_it_time = MPI_WTIME()

        return
end subroutine compute_average_iteration_time

subroutine screen_elapse_time
        implicit none

        if(rank == root) then
          print*, ' -----------------------------------------'
          write(*,'(A,f18.6,A)') " cpu time = ", e_cpu_time-s_cpu_time, ' s'
          print*, ' -----------------------------------------'
        endif

        return
end subroutine screen_elapse_time

end module output_screen
