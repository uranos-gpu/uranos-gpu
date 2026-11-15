module output_stats_tgv
use parameters_module
use storage_module
use mpi_module
use onlineStats, only: compute_3DStats
use FileModule

implicit none
private
public stat_taylor_green_vortex

real(rp), save :: Ek0    = 0.0_rp
real(rp), save :: u0_tgv = 1.0_rp
real(rp), save :: Ek_old = 0.0_rp
real(rp), save :: t_old  = 0.0_rp
logical , save :: first  = .true.

contains

subroutine stat_taylor_green_vortex
        implicit none
        real(rp) :: vol, Ek, Ew, Mmax, rhoS, rhoE, eps_k, eps_w

        call compute_3DStats(phi,U,V,W,T,sx,ex,sy,ey,sz,ez,mpi_comm_world,mpi_flag, &
                             vol,Ek,Ew,Mmax,rhoS,rhoE)

        if(first) then
          Ek0    = Ek
          u0_tgv = sqrt(8.0_rp*Ek0)
          eps_k  = 0.0_rp
          first  = .false.
        else
          eps_k = - (Ek - Ek_old)/max(time - t_old, tiny(1.0_rp))
        endif
        Ek_old = Ek
        t_old  = time

        eps_w = 2.0_rp*mu_inf*Ew

        if(rank == root) call print_tgv_monitor(Ek,Ew,eps_k,eps_w,Mmax,rhoS,rhoE)

        return
end subroutine stat_taylor_green_vortex

subroutine print_tgv_monitor(Ek, Ew, eps_k, eps_w, Mmax, rhoS, rhoE)
        implicit none
        real(rp), intent(in) :: Ek, Ew, eps_k, eps_w, Mmax, rhoS, rhoE
        logical, save      :: header_written = .false.
        character(len=600) :: fpath
        logical            :: dir_exist, file_exist
        integer            :: u, err = 0
        real(rp)           :: tstar

        inquire(file = 'DATA/'//trim(data_dir), exist = dir_exist)
        if(.not. dir_exist) call execute_command_line('mkdir -p DATA/'//trim(data_dir))

        fpath = 'DATA/'//trim(data_dir)//'/TGV_MONITOR.txt'

        inquire(file = trim(fpath), exist = file_exist)

        if(.not. header_written .and. file_exist .and. it > 0) call recover_reference_state(fpath)

        tstar = time*u0_tgv

        if(header_written .or. (file_exist .and. it > 0)) then
          open(newunit = u, file = trim(fpath), status = 'old', position = 'append', &
               iostat = err)
          header_written = .true.
        else
          open(newunit = u, file = trim(fpath), status = 'replace', iostat = err)
          write(u,'(A)')        '# Taylor-Green-Vortex monitor (volume-averaged, GCL-consistent)'
          write(u,'(A,es16.8)') '# Mach        = ', Mach
          write(u,'(A,es16.8)') '# Reynolds_nl = ', Reynolds
          write(u,'(A,es16.8)') '# mu_inf      = ', mu_inf
          write(u,'(A,es26.16e3)') '# U0          = ', u0_tgv
          write(u,'(A,es16.8)')    '# Re_flow     = ', u0_tgv/mu_inf
          write(u,'(A,es26.16e3)') '# Ek0         = ', Ek0
          write(u,'(A)')        '#'
          write(u,'(A)')        '# col 1 : it          iteration'
          write(u,'(A)')        '# col 2 : t           solver time'
          write(u,'(A)')        '# col 3 : t*          convective time  t*u0_tgv'
          write(u,'(A)')        '# col 4 : Ek          <1/2 rho |u|^2>'
          write(u,'(A)')        '# col 5 : Ek/Ek0      normalized kinetic energy'
          write(u,'(A)')        '# col 6 : Ew          <1/2 rho |omega|^2>   (enstrophy)'
          write(u,'(A)')        '# col 7 : eps_k       -dEk/dt               (KE dissipation)'
          write(u,'(A)')        '# col 8 : eps_w       2 (mu/rho0) Ew        (enstrophy dissipation)'
          write(u,'(A)')        '# col 9 : Mach_max    peak local Mach'
          write(u,'(A)')        '# col10 : rhoS        <rho s>               (mean entropy density)'
          write(u,'(A)')        '# col11 : rhoE        <rho E>               (mean total-energy density)'
          header_written = .true.
        endif

        if(err /= 0) then
          write(*,'(A,A)') ' WARNING: could not open TGV monitor file ', trim(fpath)
          return
        endif

        write(u,'(I9,10es26.16e3)') it, time, tstar, Ek, Ek/Ek0, Ew, eps_k, eps_w, Mmax, rhoS, rhoE
        flush(u)
        close(u)

        return
end subroutine print_tgv_monitor

subroutine recover_reference_state(fpath)
        implicit none
        character(len=*), intent(in) :: fpath
        character(len=200) :: line
        integer            :: u, err, ieq
        real(rp)           :: val
        logical            :: got_u0, got_ek0

        got_u0  = .false.
        got_ek0 = .false.

        open(newunit = u, file = trim(fpath), status = 'old', action = 'read', iostat = err)
        if(err /= 0) return

        do
          read(u,'(A)',iostat = err) line
          if(err /= 0)          exit
          if(line(1:1) /= '#')  exit

          ieq = index(line,'=')
          if(ieq == 0) cycle
          read(line(ieq+1:),*,iostat = err) val
          if(err /= 0) cycle

          if(index(line,' U0 ') > 0) then
            u0_tgv = val
            got_u0 = .true.
          elseif(index(line,' Ek0 ') > 0) then
            Ek0     = val
            got_ek0 = .true.
          endif
        enddo
        close(u)

        if(.not.(got_u0 .and. got_ek0)) write(*,'(A)') &
          ' WARNING: TGV monitor header incomplete, t* and Ek/Ek0 are rescaled at this restart'

        return
end subroutine recover_reference_state

end module output_stats_tgv
