module output_residuals
use parameters_module
use storage_module
use profiling_module
use output_probes     , only: flush_probes

implicit none
private
public dump_residuals

contains

subroutine dump_residuals(iter)
        use, intrinsic :: IEEE_Arithmetic

        implicit none
        integer, intent(in) :: iter
        integer  :: i,j,k
        real(rp) :: iNp
        real(rp) :: rhs_1, rhs_2, rhs_3, rhs_4, rhs_5
        real(rp) :: r_, ru, rv, rw, re
        real(rp) :: lcl_res_1, lcl_res_2, lcl_res_3, lcl_res_4, lcl_res_5
        real(rp) :: lcl_max_r_, lcl_max_ru, lcl_max_rv, lcl_max_rw, lcl_max_re
        real(rp) :: lcl_min_r_, lcl_min_ru, lcl_min_rv, lcl_min_rw, lcl_min_re
        real(rp), dimension(5) :: lcl_res, gbl_res
        real(rp), dimension(5) :: lcl_max, gbl_max
        logical                :: nan_lcl, nan_any
        real(rp), dimension(5) :: lcl_min, gbl_min

        call StartProfRange('dump_residuals')

        lcl_res_1 = 0.0_rp
        lcl_res_2 = 0.0_rp
        lcl_res_3 = 0.0_rp
        lcl_res_4 = 0.0_rp
        lcl_res_5 = 0.0_rp

        lcl_max_r_ = 0.0_rp
        lcl_max_ru = 0.0_rp
        lcl_max_rv = 0.0_rp
        lcl_max_rw = 0.0_rp
        lcl_max_re = 0.0_rp

        lcl_min_r_ = huge(0.0_rp)
        lcl_min_ru = huge(0.0_rp)
        lcl_min_rv = huge(0.0_rp)
        lcl_min_rw = huge(0.0_rp)
        lcl_min_re = huge(0.0_rp)

        !$acc parallel loop gang vector collapse(3) default(present)        &
        !$acc reduction(+:lcl_res_1,lcl_res_2,lcl_res_3,lcl_res_4,lcl_res_5) &
        !$acc reduction(max:lcl_max_r_,lcl_max_ru,lcl_max_rv,lcl_max_rw,lcl_max_re) &
        !$acc reduction(min:lcl_min_r_,lcl_min_ru,lcl_min_rv,lcl_min_rw,lcl_min_re)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 rhs_1 = RHS(i,j,k,1)
                 rhs_2 = RHS(i,j,k,2)
                 rhs_3 = RHS(i,j,k,3)
                 rhs_4 = RHS(i,j,k,4)
                 rhs_5 = RHS(i,j,k,5)

                 lcl_res_1 = lcl_res_1 + rhs_1*rhs_1
                 lcl_res_2 = lcl_res_2 + rhs_2*rhs_2
                 lcl_res_3 = lcl_res_3 + rhs_3*rhs_3
                 lcl_res_4 = lcl_res_4 + rhs_4*rhs_4
                 lcl_res_5 = lcl_res_5 + rhs_5*rhs_5

                 r_ = abs(phi(i,j,k,1))
                 ru = abs(phi(i,j,k,2))
                 rv = abs(phi(i,j,k,3))
                 rw = abs(phi(i,j,k,4))
                 re = abs(phi(i,j,k,5))

                 lcl_min_r_ = min(r_,lcl_min_r_)
                 lcl_min_ru = min(ru,lcl_min_ru)
                 lcl_min_rv = min(rv,lcl_min_rv)
                 lcl_min_rw = min(rw,lcl_min_rw)
                 lcl_min_re = min(re,lcl_min_re)

                 lcl_max_r_ = max(r_,lcl_max_r_)
                 lcl_max_ru = max(ru,lcl_max_ru)
                 lcl_max_rv = max(rv,lcl_max_rv)
                 lcl_max_rw = max(rw,lcl_max_rw)
                 lcl_max_re = max(re,lcl_max_re)

              enddo
           enddo
        enddo

        lcl_res = [ lcl_res_1, lcl_res_2, lcl_res_3, lcl_res_4, lcl_res_5 ]
        lcl_max = [ lcl_max_r_, lcl_max_ru, lcl_max_rv, lcl_max_rw, lcl_max_re ]
        lcl_min = [ lcl_min_r_, lcl_min_ru, lcl_min_rv, lcl_min_rw, lcl_min_re ]

        nan_lcl = any(ieee_is_nan(lcl_res))
        call MPI_allreduce(nan_lcl, nan_any, 1, MPI_LOGICAL, MPI_LOR, mpi_comm_cart, err)
        if(nan_any) then
                if(nan_lcl) call report_divergency(lcl_res, iter)
                call flush_probes
                call MPI_Barrier(mpi_comm_cart, err)
                call MPI_Abort(mpi_comm_cart, 1, err)
        endif

        call MPI_reduce(lcl_res, gbl_res, 5, MPI_RP, MPI_SUM, root, mpi_comm_cart, err)
        call MPI_reduce(lcl_max, gbl_max, 5, MPI_RP, MPI_MAX, root, mpi_comm_cart, err)
        call MPI_reduce(lcl_min, gbl_min, 5, MPI_RP, MPI_MIN, root, mpi_comm_cart, err)

        if(rank == root) then

          iNp = 1.0_rp/(real(nx*ny*nz,rp))

          gbl_res = sqrt(gbl_res*iNp)*Dt

          call write_residuals_line(iter, time*u_inf,                            &
                  gbl_res(1), gbl_res(2), gbl_res(3), gbl_res(4), gbl_res(5),  &
                  gbl_max(1), gbl_min(1), gbl_max(2), gbl_min(2),              &
                  gbl_max(3), gbl_min(3), gbl_max(4), gbl_min(4),              &
                  gbl_max(5), gbl_min(5))

        endif

        call EndProfRange

        return
end subroutine dump_residuals

subroutine report_divergency(lcl_res, iter)
        use, intrinsic :: IEEE_Arithmetic
        use, intrinsic :: ISO_FORTRAN_ENV, only: error_unit
        implicit none
        real(rp), dimension(5), intent(in) :: lcl_res
        integer               , intent(in) :: iter
        character(len=8), dimension(5), parameter :: var = &
                ['rho   ','rho*u ','rho*v ','rho*w ','rho*E ']
        integer :: n

        write(error_unit,'(A)') ''
        write(error_unit,'(A)') '================================================================'
        write(error_unit,'(A)') ' URANOS: SIMULATION DIVERGED -- NaN detected in residuals'
        write(error_unit,'(A)') '----------------------------------------------------------------'
        write(error_unit,'(A,I0)')      '   iteration : ', iter
        write(error_unit,'(A,ES12.5)')  '   time      : ', time
        write(error_unit,'(A,I0)')      '   mpi rank  : ', rank
        write(error_unit,'(A)')         '   NaN in    :'
        do n = 1,5
           if(ieee_is_nan(lcl_res(n))) write(error_unit,'(A,A,A)') &
                '                ', trim(var(n)), '  (RHS L2-norm)'
        enddo
        write(error_unit,'(A)') '----------------------------------------------------------------'
        write(error_unit,'(A)') ' Likely causes: CFL too high, bad initial/boundary state,'
        write(error_unit,'(A)') ' inflow turbulence too aggressive, or shock-capturing failure.'
        write(error_unit,'(A)') ' Aborting MPI on errorcode 1.'
        write(error_unit,'(A)') '================================================================'
        flush(error_unit)

end subroutine report_divergency

subroutine write_residuals_line(it_in, t_phys,                  &
                res1, res2, res3, res4, res5,                   &
                max_r, min_r, max_ru, min_ru, max_rv, min_rv,   &
                max_rw, min_rw, max_re, min_re)
        implicit none
        integer , intent(in) :: it_in
        real(rp), intent(in) :: t_phys
        real(rp), intent(in) :: res1, res2, res3, res4, res5
        real(rp), intent(in) :: max_r, min_r, max_ru, min_ru
        real(rp), intent(in) :: max_rv, min_rv, max_rw, min_rw
        real(rp), intent(in) :: max_re, min_re

        character(len=dl) :: dirpath, filename
        logical           :: dir_exist, file_exist
        integer           :: unit, ios

        dirpath  = 'DATA/'//trim(data_dir)
        filename = trim(dirpath)//'/RESIDUALS_MONITOR.txt'

        inquire(file=trim(dirpath), exist=dir_exist)
        if(.not. dir_exist) call execute_command_line('mkdir -p '//trim(dirpath))

        inquire(file=trim(filename), exist=file_exist)
        if(it_in == 0) then
          if(file_exist) call execute_command_line('rm -f '//trim(filename))
          open(newunit=unit, file=trim(filename), status='new', iostat=ios)
        else
          open(newunit=unit, file=trim(filename), status='unknown', &
               position='append', iostat=ios)
        endif

        if(ios /= 0) then
          print*, ' FATAL ERROR! File ', trim(filename), ' has not been opened'
          stop
        endif

        if(it_in == 0 .or. .not. file_exist) then
          write(unit,'(A)') '# RESIDUALS MONITOR: one row every rsOut iterations. Non-dimensional.'
          write(unit,'(A)') '# q = rho, rho u, rho v, rho w, rho E (conservative variables).'
          write(unit,'(A)') '# RES(q) = Dt sqrt( sum of RHS(q)^2 over the interior nodes / (nx ny nz) ):'
          write(unit,'(A)') '#   RMS of the RHS of the last Runge-Kutta stage times the time step,'
          write(unit,'(A)') '#   i.e. the RMS change of q per step.'
          write(unit,'(A)') '# MAX/MIN = extrema of |q| over the interior nodes. Absolute values: a'
          write(unit,'(A)') '#   negative density shows up as its modulus, not as a negative minimum.'
          write(unit,'(A)') '# Column 1  : iteration'
          write(unit,'(A)') '# Column 2  : time u_inf (reference length 1)'
          write(unit,'(A)') '# Column 3  : RES(rho)'
          write(unit,'(A)') '# Column 4  : RES(rho u)'
          write(unit,'(A)') '# Column 5  : RES(rho v)'
          write(unit,'(A)') '# Column 6  : RES(rho w)'
          write(unit,'(A)') '# Column 7  : RES(rho E)'
          write(unit,'(A)') '# Column 8  : max |rho|'
          write(unit,'(A)') '# Column 9  : min |rho|'
          write(unit,'(A)') '# Column 10 : max |rho u|'
          write(unit,'(A)') '# Column 11 : min |rho u|'
          write(unit,'(A)') '# Column 12 : max |rho v|'
          write(unit,'(A)') '# Column 13 : min |rho v|'
          write(unit,'(A)') '# Column 14 : max |rho w|'
          write(unit,'(A)') '# Column 15 : min |rho w|'
          write(unit,'(A)') '# Column 16 : max |rho E|'
          write(unit,'(A)') '# Column 17 : min |rho E|'
        endif

        write(unit,10) it_in, t_phys,                                       &
                res1, res2, res3, res4, res5,                               &
                max_r, min_r, max_ru, min_ru, max_rv, min_rv,               &
                max_rw, min_rw, max_re, min_re

        flush(unit)
        close(unit)

        10 format(I7,16e18.9)

        return
end subroutine write_residuals_line

end module output_residuals
