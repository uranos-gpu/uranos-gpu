module input_module
use parameters_module
use mpi_module

implicit none
interface open_unit
  module procedure open_unit
endinterface


contains
subroutine read_input_data
! ------------------------------------------------------------------------------
!
!       This subroutine reads input data from a file and sets some input values
!
! ------------------------------------------------------------------------------
        implicit none
        integer       :: io_err = 0
        character(dl) :: input_file

        ! --- getting the argument --- !

        call get_command_argument(1,input_file)
        call get_command_argument(2,restart_file)

        if(len_trim(input_file) == 0) then
                
          if(rank == root) then
            print*, ' Input data file missed!. Program will be stopped.' 
          endif
          call secure_stop

        else
          
          if(rank == root) then
            print*, ' Reading input data from file ', "'"//trim(input_file)//"'", '.'
          endif

        endif
        
        ! --- reading the fine in argument --- !

        open(unit = 1, file = input_file, IOSTAT = io_err)
        if(io_err .ne. 0) then
                
          if(rank == root) print*, ' Error in reading file ', trim(input_file), '.'
          call secure_stop

        endif

        read(1, nml = input_list)
        close(1)
        !
        ! ==== adimensional grups
        !
        mu_inf = sqrt(gamma0)*Mach/Reynolds
         u_inf = sqrt(gamma0)*Mach
         q_inf = 0.5_rp*gamma0*Mach**2

        k_inf = gamma0/((gamma0-1._rp) * Prandtl)
        weno_num = (weno_order+1)/2

        !
        ! ==== set hybrid_weno falg
        !
        hybrid_weno = .false.

        if    (scheme == 'hybrid_weno5') then
          hybrid_weno = .true.
        elseif(scheme == 'hybrid_wenoEP') then
          hybrid_weno = .true.
        elseif(scheme == 'hybrid_tenoEP') then
          hybrid_weno = .true.
        elseif(scheme == 'hybrid_tenoaEP') then
          hybrid_weno = .true.
        endif

        ! Set inflow flag
        inflow = .false.
        if(trim(bc(1)) == 'supersonic_inflow'     .or. &
           trim(bc(1)) == 'subsonic_inflow'       .or. &
           trim(bc(1)) == 'nscbc_inflow'          .or. &
           trim(bc(1)) == 'nscbc_inflow_relaxed'  .or. &
           trim(bc(1)) == 'turb_supersonic_inflow'.or. &
           trim(bc(1)) == 'pressure_inflow'       .or. &
           trim(bc(1)) == 'shock_inflow')         then

           inflow = .true.

        endif

        wmles = .false.
        if(trim(bc(3)) == 'dws_isothermal'          .or. &
           trim(bc(3)) == 'dws_adiabatic'           .or. &
           trim(bc(3)) == 'idws_adiabatic'          .or. &
           trim(bc(3)) == 'idws_isothermal'         .or. &
           trim(bc(3)) == 'static_dws_adiabatic'    .or. &
           trim(bc(3)) == 'static_dws_isothermal'   .or. &
           trim(bc(3)) == 'istatic_dws_adiabatic'   .or. &
           trim(bc(3)) == 'vorticity_dws_adiabatic' .or. &
           trim(bc(3)) == 'aws_adiabatic'           .or. &
           trim(bc(3)) == 'aws_isothermal'          ) then
           wmles = .true.
        endif

        ! forcing viscous flux to be staggered if wmles is active
        if(wmles)     diffusion_scheme = 'staggered'

        return
end subroutine read_input_data

subroutine get_unit(iunit)

        implicit none
        integer, intent(out) :: iunit
        integer              :: i, ios = 0
        logical              :: isopen
        
        iunit = 0
        
        do i = 0, 1000
        
          if(i /= 5 .and. i /= 6 .and. i /= 9) then
        
            inquire(unit = i, opened = isopen, iostat = ios)
        
            if(ios == 0) then
              if(.not. isopen) then
                iunit = i
                return
              end if
            end if
          end if
        end do
        
        return
end subroutine get_unit

subroutine open_unit(iunit,filename,status,access,position)

        implicit none
        integer                     , intent(in) :: iunit
        character(len = *)          , intent(in) :: filename
        character(len = *), optional, intent(in) :: status
        character(len = *), optional, intent(in) :: access
        character(len = *), optional, intent(in) :: position
        integer                                  :: ios = 0
       
        if    (present(status)  ) then
          open(unit = iunit, file = filename, status = status, iostat = ios)

        elseif(present(access)  ) then
          open(unit = iunit, file = filename, access = access,iostat = ios)

        elseif(present(position)) then
          open(unit = iunit, file = filename, position = position,iostat = ios)

        elseif(present(status).and.present(access)  ) then
          open(unit = iunit, file = filename, status = status, access = access,iostat = ios)

        elseif(present(status).and.present(position)) then
          open(unit = iunit, file = filename, status = status, position = position,iostat = ios)

        elseif(present(access).and.present(position)) then
          open(unit = iunit, file = filename, access = access, position = position,iostat = ios)

        elseif(present(access).and.present(status).and.present(position)) then
          open(unit = iunit, file = filename, status = status, access = access, position = position,iostat = ios)

        else
          open(unit = iunit, file = filename)
        endif

        if(ios /= 0) then

          write(*,'(a)') ''
          write(*,'(a)') ' FATAL ERROR!'
          write(*,'(a)') ' Could not open file "' //trim(filename)// '"'

          stop

        endif

        return
end subroutine open_unit








subroutine make_saving_directories
! ----------------------------------------------------------------------
!
!       This subroutine create the saving directories
!
! ----------------------------------------------------------------------
        use FileModule

        implicit none
        logical :: dir_exist

        if(rank == root) then

          ! --- removing temporary data from previous simulations
          if(data_dir.eq.'temp') then
            inquire(file = 'DATA/'//trim(data_dir), exist= dir_exist)
            if(dir_exist) then
              call execute_command_line('rm -r DATA/'//trim(data_dir))
              write(*,*) " Previous DATA/"//trim(data_dir)//" has been removed."
            endif
          endif

          call mkdir('DATA/'//trim(data_dir))
          call mkdir('DATA/'//trim(data_dir)//'/BINARY')
          if(StFlg) call mkdir('DATA/'//trim(data_dir)//'/BINARY_STAT')
          if(wmles) call mkdir('DATA/'//trim(data_dir)//'/BINARY_WMLES_STAT')

        endif
        return
end subroutine make_saving_directories


subroutine initProbes(x,y,z,Probe)
        implicit none
        real(rp) , dimension(:), allocatable, intent(in)    :: x,y,z
        type(prb)                           , intent(inout) :: Probe

        selectcase(ic)
        case('supersonic_intake')
          call initProbesSupersonicIntake(x,y,z,Probe)

        case('turbulent_BL')
          call initProbesTurbulentBoundaryLayer(x,y,z,Probe)

        case('supersonic_ramp')
          call initProbesSupersonicRamp(x,y,z,Probe)

        case('blades')
          call initProbesBlades(x,y,z,Probe)

        endselect

        return
end subroutine initProbes


subroutine initProbesSupersonicRamp(x,y,z,Probe)

        use math_tools_module
        use real_to_integer_module
                
        implicit none
        real(rp) , dimension(:), allocatable, intent(in)    :: x,y,z
        type(prb)                           , intent(inout) :: Probe

        !local declarations
        real(rp) :: PidStr(3), PidEnd(3)
        real(rp) :: xp, yp, zp
        integer  :: n, err = 0

        Probe%n = nint(Lx/2.0_rp)
        allocate(Probe%i(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'
        allocate(Probe%x(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'
        allocate(Probe%rank(Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'

        do n = 1,Probe%n

           xp = xmin + (n-1)*Lx/(real(Probe%n))
           yp = 0.0_rp
           if(xp > 0.0_rp) yp = tan(DegToRad(24.0_rp))*xp
           zp = 0.5_rp*(zmax-zmin)

           Probe%x(1,n) = xp
           Probe%x(2,n) = yp
           Probe%x(3,n) = zp

        enddo

        ! find probes location within the procs
        PidStr = (/x(sx-1), y(sy-1), z(sz-1)/)
        PidEnd = (/x(ex+1), y(ey+1), z(ez+1)/)

        do n = 1,Probe%n

           if(Probe%x(1,n) >= PidStr(1) .and. Probe%x(1,n) <= PidEnd(1) .and. &
              Probe%x(2,n) >= PidStr(2) .and. Probe%x(2,n) <= PidEnd(2)) then

              Probe%i(1,n) = locNearest(x,sx,ex,Probe%x(1,n))
              Probe%i(2,n) = locNearest(y,sy,ey,Probe%x(2,n))

              Probe%rank(n) = rank

           endif
        enddo

        
        return
end subroutine initProbesSupersonicRamp


subroutine initProbesSupersonicIntake(x,y,z,Probe)

        use math_tools_module    , only: degToRad
        use real_to_integer_module, only: nearest_integer_opt

        implicit none
        real(rp) , dimension(:), allocatable, intent(in)    :: x,y,z
        type(prb)                           , intent(inout) :: Probe

        ! local declarations
        real(rp), parameter :: L0     = 150.00_rp
        real(rp), parameter :: r1     =  52.80_rp/L0
        real(rp), parameter :: r2     =  32.36_rp/L0
        real(rp), parameter :: theta1 =  10.00_rp
        real(rp), parameter :: theta2 =  22.00_rp

        real(rp) :: PidStr(3), PidEnd(3)
        real(rp) :: t1, t2, a1, a2, b1, b2
        real(rp) :: x1, x2, x3, x4, x5, x6
        real(rp) :: y1, y2, y3, y4, y5, y6
        real(rp) :: z1, z2, z3, z4, z5, z6
        real(rp) :: hChann
        integer  :: n, err = 0

        Probe%n = 6

        allocate(Probe%i(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'
        allocate(Probe%x(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'

        t1 = degToRad(theta1)
        t2 = degToRad(theta2)

        a1 = r1*cos(t1)
        b1 = r1*sin(t1)

        a2 = r2*cos(t2)
        b2 = r2*sin(t2)

        ! x Probes locations
        x1 = (L0 - 113.00_rp)/L0
        x2 = (L0 -  75.00_rp)/L0
        x3 = (L0 -  53.00_rp)/L0

        ! y Probes locations
        y1 = x1*tan(t1)
        y2 = b1 + (x2 - a1)*tan(t2)
        y3 = b1 + b2

        hChann = 0.24_rp - y3
       
        ! inner channel probes
        x4 = x3
        y4 = y3 + 1.0_rp/3.0_rp*hChann

        x5 = x3
        y5 = y3 + 1.0_rp/2.0_rp*hChann

        x6 = x3
        y6 = y3 + 2.0_rp/3.0_rp*hChann

        ! z Probes locations
        if(dims == 2) then
          z1 = 0.0_rp
          z2 = 0.0_rp
          z3 = 0.0_rp
          z4 = 0.0_rp
          z5 = 0.0_rp
          z6 = 0.0_rp
        else
          z1 = 0.5_rp*(zmax-zmin)
          z2 = 0.5_rp*(zmax-zmin)
          z3 = 0.5_rp*(zmax-zmin)
          z4 = 0.5_rp*(zmax-zmin)
          z5 = 0.5_rp*(zmax-zmin)
          z6 = 0.5_rp*(zmax-zmin)
        endif

        Probe%x(:,1) = (/x1, y1, z1/)
        Probe%x(:,2) = (/x2, y2, z2/)
        Probe%x(:,3) = (/x3, y3, z3/)
        Probe%x(:,4) = (/x4, y4, z4/)
        Probe%x(:,5) = (/x5, y5, z5/)
        Probe%x(:,6) = (/x6, y6, z6/)
        
        ! find probes location within the procs
        PidStr = (/x(sx)  , y(sy)  , z(sz)  /)
        PidEnd = (/x(ex+1), y(ey+1), z(ez+1)/)

        do n = 1,Probe%n

           if(Probe%x(1,n) >= PidStr(1) .and. Probe%x(1,n) <= PidEnd(1) .and. &
              Probe%x(2,n) >= PidStr(2) .and. Probe%x(2,n) <= PidEnd(2) .and. &
              Probe%x(3,n) >= PidStr(3) .and. Probe%x(3,n) <= PidEnd(3))then

              Probe%i(1,n) = nearest_integer_opt(x,sx,ex,Probe%x(1,n))
              Probe%i(2,n) = nearest_integer_opt(y,sy,ey,Probe%x(2,n))
              Probe%i(3,n) = nearest_integer_opt(z,sz,ez,Probe%x(3,n))

           endif
        enddo

        return
end subroutine initProbesSupersonicIntake




subroutine initProbesTurbulentBoundaryLayer(x,y,z,Probe)

        use real_to_integer_module, only: nearest_integer_opt

        implicit none
        real(rp) , dimension(:), allocatable, intent(in)    :: x,y,z
        type(prb)                           , intent(inout) :: Probe

        !local declarations
        real(rp) :: PidStr(3), PidEnd(3)
        real(rp) :: xp, yp,zp
        integer  :: n, err = 0


        Probe%n = 36
        allocate(Probe%i(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'
        allocate(Probe%x(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'

        do n = 1,Probe%n
           
           xp = xmin + n*Lx/real(Probe%n)
           yp = y(8)
           zp = 0.5_rp*(zmax-zmin)

           Probe%x(1,n) = xp
           Probe%x(2,n) = yp
           Probe%x(3,n) = zp

        enddo

        ! find probes location within the procs
        PidStr = (/x(sx)  , y(sy)  , z(sz)  /)
        PidEnd = (/x(ex+1), y(ey+1), z(ez+1)/)

        do n = 1,Probe%n

           if(Probe%x(1,n) >= PidStr(1) .and. Probe%x(1,n) <= PidEnd(1) .and. &
              Probe%x(2,n) >= PidStr(2) .and. Probe%x(2,n) <= PidEnd(2) .and. &
              Probe%x(3,n) >= PidStr(3) .and. Probe%x(3,n) <= PidEnd(3))then

              Probe%i(1,n) = nearest_integer_opt(x,sx,ex,Probe%x(1,n))
              Probe%i(2,n) = nearest_integer_opt(y,sy,ey,Probe%x(2,n))
              Probe%i(3,n) = nearest_integer_opt(z,sz,ez,Probe%x(3,n))

           endif
        enddo


        return
end subroutine initProbesTurbulentBoundaryLayer

subroutine initProbesBlades(x,y,z,Probe)

        use real_to_integer_module, only: nearest_integer_opt
        use allocate_module

        implicit none
        real(rp) , dimension(:), allocatable, intent(in)    :: x,y,z
        type(prb)                           , intent(inout) :: Probe

        real(rp), dimension(:), allocatable :: distx, disty
        real(rp), dimension(3) :: PidStr, PidEnd
        integer :: i, j, n, err 

        err = 0

        Probe%n = 14
        allocate(Probe%i(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'
        allocate(Probe%x(3,Probe%n), stat = err)
        if(err .ne. 0) stop ' Allocation error in initProbes'

        Probe%i    = 0
        Probe%x    = 0.0_rp

        ! frontal probes
        Probe%x(:, 1) = (/-0.3_rp, 0.0_rp, 0.0_rp/)

        ! suction side body 2 probes
        Probe%x(:, 2) = (/0.0850000739131077_rp, 1.8300575043978300_rp, 0.0_rp/)
        Probe%x(:, 3) = (/0.1946723431933419_rp, 1.8390009904356440_rp, 0.0_rp/)
        Probe%x(:, 4) = (/0.3043446124735760_rp, 1.7283382855115530_rp, 0.0_rp/)
        Probe%x(:, 5) = (/0.3591881384244681_rp, 1.6172468845625080_rp, 0.0_rp/)
        Probe%x(:, 6) = (/0.4140316643753603_rp, 1.4604919656451880_rp, 0.0_rp/)
        Probe%x(:, 7) = (/0.4414386447292562_rp, 1.3649666651884040_rp, 0.0_rp/)
        Probe%x(:, 8) = (/0.4688751903262524_rp, 1.2591378775111980_rp, 0.0_rp/)
        Probe%x(:, 9) = (/0.4962821706801484_rp, 1.1466569101364440_rp, 0.0_rp/)
        Probe%x(:,10) = (/0.5237039336555944_rp, 1.0332594202255830_rp, 0.0_rp/)
        
        ! wake points
        Probe%x(:,11) = (/0.5652920962199310_rp, 0.42068965517241380_rp, 0.0_rp/)
        Probe%x(:,12) = (/0.6013745704467353_rp, 0.32586206896551717_rp, 0.0_rp/)
        Probe%x(:,13) = (/0.6391752577319589_rp, 0.22758620689655162_rp, 0.0_rp/)
        Probe%x(:,14) = (/0.6838487972508591_rp, 0.12413793103448267_rp, 0.0_rp/)

        ! move the probes according the domain periodicity
        do i = 1,Probe%n
           if(Probe%x(2,i) > ymax) Probe%x(2,i) = Probe%x(2,i) - 0.85_rp
           if(Probe%x(2,i) < ymin) Probe%x(2,i) = Probe%x(2,i) + 0.85_rp
        enddo

        ! find probes location within the procs
        PidStr = (/x(sx)  , y(sy)  , z(sz)  /)
        PidEnd = (/x(ex+1), y(ey+1), z(ez+1)/)
                
        call AllocateReal(distx,sx,ex)
        call AllocateReal(disty,sy,ey)

        do n = 1,Probe%n

           if(Probe%x(1,n) >= PidStr(1) .and. Probe%x(1,n) <= PidEnd(1) .and. &
              Probe%x(2,n) >= PidStr(2) .and. Probe%x(2,n) <= PidEnd(2)) then
                
              do i = sx,ex
                 distx(i) = abs(x(i) - probe%x(1,n))
              enddo
              do j = sy,ey
                 disty(j) = abs(y(j) - probe%x(2,n))
              enddo
        
              Probe%i(1,n) = minloc(distx,1)
              Probe%i(2,n) = minloc(disty,1)

           endif
        enddo
        
        call DeallocateReal(distx)
        call DeallocateReal(disty)

        return
end subroutine initProbesBlades




end module

