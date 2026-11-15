module output_probes
use parameters_module
use mpi_module
use storage_module        , only: phi, U, V, W, P, T, time
use mesh_module           , only: xx, yy, nx_eta, ny_eta
use fluid_functions_module, only: laminar_viscosity
use FileModule            , only: prepare_row_file, prepare_record_file
use profiling_module      , only: StartProfRange, EndProfRange

implicit none
private

integer, parameter :: pb_bfr_sz = 10000

logical, public :: probes_on = .false.

integer :: nwl = 0, nfl = 0
integer , allocatable, dimension(:)   :: wl_n, wl_i
integer , allocatable, dimension(:)   :: fl_n, fl_i, fl_j
real(rp), allocatable, dimension(:,:) :: wl_g
real(rp), allocatable, dimension(:,:) :: wl_s
real(rp), allocatable, dimension(:,:) :: fl_s
real(rp), allocatable, dimension(:)   :: s_lcl, s_gbl
integer , allocatable, dimension(:)   :: bfr_it
real(rp), allocatable, dimension(:,:) :: bfr
integer            :: ncol = 0, n_bfr = 0
character(len=dl)  :: pb_file
character(len=32)  :: row_fmt

integer :: wm_comm = MPI_COMM_NULL, wm_fh, wm_type
integer(kind=MPI_OFFSET_KIND) :: wm_off, wm_rec
real(4), allocatable, dimension(:,:) :: wm_p

public :: init_probes, sample_probes, flush_probes, write_wall_map, end_probes

contains

subroutine init_probes
        implicit none
        real(rp), dimension(4,prb_max)    :: wloc, floc
        character(len=:), allocatable     :: hdr
        character(len=6), dimension(6)    :: fname = ['rho_f','u_f  ','v_f  ','w_f  ','p_f  ','T_f  ']
        character(len=25)                 :: cname
        integer :: n, c, iv, unit
        logical :: ok

        call init_wall_map
        if(itPrb <= 0) return

        if(prb_nwall == 0 .and. prb_nfield == 0) then
          if(rank == root) print*, ' ERROR: itPrb > 0 needs probes: prb_nwall, prb_xc, prb_side', &
                                   ' and/or prb_nfield, prb_x, prb_y.'
          call secure_stop
        endif

        ok = prb_nwall >= 0 .and. prb_nwall <= prb_max .and. prb_nfield >= 0 .and. prb_nfield <= prb_max
        do n = 1,min(prb_nwall,prb_max)
           if(abs(prb_side(n)) /= 1 .or. prb_xc(n) <= 0.0_rp .or. prb_xc(n) >= 1.0_rp) ok = .false.
        enddo
        do n = 1,min(prb_nfield,prb_max)
           if(prb_x(n) == huge(1.0_rp) .or. prb_y(n) == huge(1.0_rp)) ok = .false.
        enddo
        if(.not. ok) then
          if(rank == root) print*, ' ERROR: probes need 0 <= prb_nwall, prb_nfield <= ', prb_max, &
                                   ', 0 < prb_xc < 1, prb_side = +1/-1 and prb_x, prb_y for every field probe.'
          call secure_stop
        endif
        probes_on = .true.

        allocate(wl_n(max(prb_nwall,1)), wl_i(max(prb_nwall,1)), wl_g(4,max(prb_nwall,1)), &
                 wl_s(2,max(prb_nwall,1)))
        allocate(fl_n(max(prb_nfield,1)), fl_i(max(prb_nfield,1)), fl_j(max(prb_nfield,1)), &
                 fl_s(6,max(prb_nfield,1)))
        wl_n = 0; wl_i = 1; wl_g = 0.0_rp; wl_s = 0.0_rp
        fl_n = 0; fl_i = 1; fl_j = 1; fl_s = 0.0_rp
        if(prb_nwall  > 0) call locate_wall_probes(wloc)
        if(prb_nfield > 0) call locate_field_probes(floc)
        !$acc enter data copyin(wl_i, wl_g, fl_i, fl_j) create(wl_s, fl_s)

        ncol = 1 + 2*prb_nwall + 6*prb_nfield
        allocate(s_lcl(ncol-1), s_gbl(ncol-1))

        if(rank == root) then
          allocate(bfr_it(pb_bfr_sz), bfr(ncol,pb_bfr_sz))
          write(row_fmt,'(A,I0,A)') '(i12,', ncol, 'es25.16e3)'

          allocate(character(len=12+25*ncol) :: hdr)
          hdr(1:12) = '#         it'
          write(hdr(13:37),'(A25)') 't*u_inf/c'
          c = 37
          do n = 1,prb_nwall
             write(cname,'(A,I2.2,A,I0)') 'Cp_w', n, '_i', nint(wloc(1,n))
             write(hdr(c+1:c+25),'(A25)') trim(cname)
             write(cname,'(A,I2.2,A,I0)') 'Cf_w', n, '_i', nint(wloc(1,n))
             write(hdr(c+26:c+50),'(A25)') trim(cname)
             c = c + 50
          enddo
          do n = 1,prb_nfield
             do iv = 1,6
                write(cname,'(A,I2.2,A,I0,A,I0)') trim(fname(iv)), n, '_i', nint(floc(1,n)), '_j', nint(floc(2,n))
                write(hdr(c+1:c+25),'(A25)') trim(cname)
                c = c + 25
             enddo
          enddo
          pb_file = 'DATA/'//trim(data_dir)//'/PROBES.txt'
          call prepare_row_file(pb_file, hdr, it)

          open(newunit=unit, file='DATA/'//trim(data_dir)//'/PROBES_LOCATION.txt', status='replace')
          write(unit,'(A)') '# wall probes:  n, side, i, x, y, x/c'
          do n = 1,prb_nwall
             write(unit,'(2I4,I8,3es16.7)') n, prb_side(n), nint(wloc(1,n)), wloc(2:4,n)
          enddo
          write(unit,'(A)') '# field probes: n, i, j, x, y'
          do n = 1,prb_nfield
             write(unit,'(I4,2I8,2es16.7)') n, nint(floc(1,n)), nint(floc(2,n)), floc(3:4,n)
          enddo
          close(unit)

          write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' Probes: ', prb_nwall, ' wall + ', prb_nfield, ' field, every ', &
                                            itPrb, ' it, buffer ', pb_bfr_sz, ' rows -> '//trim(pb_file)
        endif

        return
end subroutine init_probes

subroutine locate_wall_probes(wloc)
        implicit none
        real(rp), dimension(4,prb_max), intent(out) :: wloc
        real(rp), allocatable, dimension(:) :: xw, yw, nxw, nyw
        real(rp) :: chord, cx, cy, xi, best, tx, ty, ds, nhx, nhy, d1, d2
        integer  :: n, i, ib, i_nxt, i_prv, i_le, i_te, dir

        allocate(xw(nx), yw(nx), nxw(nx), nyw(nx))
        call gather_wall(xw, yw, nxw, nyw)
        call check_body_wall('wall probes', xw, yw)
        call wall_chord(xw, yw, i_le, i_te, chord)
        cx = (xw(i_te)-xw(i_le))/chord
        cy = (yw(i_te)-yw(i_le))/chord

        wloc = 0.0_rp
        do n = 1,prb_nwall
           dir  = side_dir(xw, yw, i_le, i_te, prb_side(n))
           best = huge(1.0_rp); ib = 0
           i = mod(i_le-1+dir+nx,nx) + 1
           do while(i /= i_te)
              xi = ((xw(i)-xw(i_le))*cx + (yw(i)-yw(i_le))*cy)/chord
              if(abs(xi-prb_xc(n)) < best) then
                best = abs(xi-prb_xc(n)); ib = i; wloc(4,n) = xi
              endif
              i = mod(i-1+dir+nx,nx) + 1
           enddo
           if(ib == 0) then
             if(rank == root) print*, ' ERROR: wall probe', n, ' not found on side', prb_side(n)
             call secure_stop
           endif
           wloc(1:3,n) = [real(ib,rp), xw(ib), yw(ib)]

           if(sy == 1 .and. ib >= sx .and. ib <= ex) then
             i_nxt = mod(ib,nx) + 1
             i_prv = mod(ib-2+nx,nx) + 1
             tx  = 0.5_rp*(xw(i_nxt) - xw(i_prv))
             ty  = 0.5_rp*(yw(i_nxt) - yw(i_prv))
             ds  = sqrt(tx**2 + ty**2)
             tx  = tx/ds
             ty  = ty/ds
             nhx =  ty
             nhy = -tx
             if(nhx*(xx(ib,1)-xw(ib)) + nhy*(yy(ib,1)-yw(ib)) < 0.0_rp) then
               nhx = -nhx; nhy = -nhy
             endif
             if(tx*cx + ty*cy < 0.0_rp) then
               tx = -tx; ty = -ty
             endif
             d1 = nhx*(xx(ib,1)-xw(ib)) + nhy*(yy(ib,1)-yw(ib))
             d2 = nhx*(xx(ib,2)-xw(ib)) + nhy*(yy(ib,2)-yw(ib))
             if(d1 <= 0.0_rp .or. d2 <= d1) then
               print*, ' ERROR: wall probe', n, ': wall-normal distances not increasing at i =', ib, d1, d2
               call secure_stop
             endif
             nwl = nwl + 1
             wl_n(nwl)   = n
             wl_i(nwl)   = ib
             wl_g(:,nwl) = [d2/(d1*(d2-d1)), -d1/(d2*(d2-d1)), tx, ty]
           endif
        enddo
        deallocate(xw, yw, nxw, nyw)

        return
end subroutine locate_wall_probes

subroutine locate_field_probes(floc)
        implicit none
        real(rp), dimension(4,prb_max), intent(out) :: floc
        real(rp) :: d, d_lcl, d_min
        integer  :: n, i, j, id_lcl, id

        floc = 0.0_rp
        do n = 1,prb_nfield
           d_lcl = huge(1.0_rp); id_lcl = huge(1)
           do j = sy,ey
              do i = sx,ex
                 d = (xx(i,j)-prb_x(n))**2 + (yy(i,j)-prb_y(n))**2
                 if(d < d_lcl) then
                   d_lcl = d; id_lcl = i + (j-1)*nx
                 endif
              enddo
           enddo
           call MPI_allreduce(d_lcl, d_min, 1, MPI_RP, MPI_MIN, mpi_comm_cart, mpi_err)
           if(d_lcl /= d_min) id_lcl = huge(1)
           call MPI_allreduce(id_lcl, id, 1, MPI_INTEGER, MPI_MIN, mpi_comm_cart, mpi_err)
           i = mod(id-1,nx) + 1
           j = (id-1)/nx + 1

           floc(:,n) = 0.0_rp
           if(i >= sx .and. i <= ex .and. j >= sy .and. j <= ey) then
             if(cart_coord(3) == 0) floc(3:4,n) = [xx(i,j), yy(i,j)]
             nfl = nfl + 1
             fl_n(nfl) = n
             fl_i(nfl) = i
             fl_j(nfl) = j
           endif
           call MPI_allreduce(MPI_IN_PLACE, floc(3:4,n), 2, MPI_RP, MPI_SUM, mpi_comm_cart, mpi_err)
           floc(1:2,n) = [real(i,rp), real(j,rp)]
        enddo

        return
end subroutine locate_field_probes

subroutine sample_probes
        implicit none
        real(rp) :: s1, s2, s3, s4, s5, s6, mu_w, dudn, dvdn, sc
        integer  :: m, n, i, j, k, c

        call StartProfRange('sample_probes')

        if(nwl > 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector
          do m = 1,nwl
             i  = wl_i(m)
             s1 = 0.0_rp; s2 = 0.0_rp
             !$acc loop seq
             do k = sz,ez
                mu_w = mu_inf*laminar_viscosity(0.5_rp*(T(i,0,k) + T(i,1,k)), Tref, vis_flag)
                dudn = wl_g(1,m)*U(i,1,k) + wl_g(2,m)*U(i,2,k)
                dvdn = wl_g(1,m)*V(i,1,k) + wl_g(2,m)*V(i,2,k)
                s1   = s1 + 0.5_rp*(P(i,0,k) + P(i,1,k))
                s2   = s2 + mu_w*(dudn*wl_g(3,m) + dvdn*wl_g(4,m))
             enddo
             wl_s(1,m) = s1
             wl_s(2,m) = s2
          enddo
          !$acc end parallel
          !$acc update host(wl_s)
        endif

        if(nfl > 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector
          do m = 1,nfl
             i  = fl_i(m)
             j  = fl_j(m)
             s1 = 0.0_rp; s2 = 0.0_rp; s3 = 0.0_rp
             s4 = 0.0_rp; s5 = 0.0_rp; s6 = 0.0_rp
             !$acc loop seq
             do k = sz,ez
                s1 = s1 + phi(i,j,k,1)
                s2 = s2 + U(i,j,k)
                s3 = s3 + V(i,j,k)
                s4 = s4 + W(i,j,k)
                s5 = s5 + P(i,j,k)
                s6 = s6 + T(i,j,k)
             enddo
             fl_s(1,m) = s1
             fl_s(2,m) = s2
             fl_s(3,m) = s3
             fl_s(4,m) = s4
             fl_s(5,m) = s5
             fl_s(6,m) = s6
          enddo
          !$acc end parallel
          !$acc update host(fl_s)
        endif

        s_lcl = 0.0_rp
        do m = 1,nwl
           n = wl_n(m)
           s_lcl(2*n-1:2*n) = wl_s(:,m)
        enddo
        do m = 1,nfl
           c = 2*prb_nwall + 6*(fl_n(m)-1)
           s_lcl(c+1:c+6) = fl_s(:,m)
        enddo
        call MPI_reduce(s_lcl, s_gbl, ncol-1, MPI_RP, MPI_SUM, root, mpi_comm_cart, mpi_err)

        if(rank == root) then
          sc = 1.0_rp/real(nz,rp)
          n_bfr = n_bfr + 1
          bfr_it(n_bfr)  = it
          bfr(1,n_bfr)   = time*u_inf
          do n = 1,prb_nwall
             bfr(2*n  ,n_bfr) = (s_gbl(2*n-1)*sc - 1.0_rp)/q_inf
             bfr(2*n+1,n_bfr) =  s_gbl(2*n  )*sc/q_inf
          enddo
          do c = 2*prb_nwall+1,ncol-1
             bfr(c+1,n_bfr) = s_gbl(c)*sc
          enddo
          if(n_bfr == pb_bfr_sz) call flush_probes
        endif

        call EndProfRange

        return
end subroutine sample_probes

subroutine flush_probes
        implicit none
        integer :: unit, ios, n

        call StartProfRange('flush_probes')

        if(.not. probes_on .or. rank /= root .or. n_bfr == 0) then
          call EndProfRange
          return
        endif

        open(newunit=unit, file=trim(pb_file), status='old', position='append', iostat=ios)
        if(ios /= 0) then
          print*, ' ERROR: probes: cannot append to ', trim(pb_file)
          call MPI_Abort(mpi_comm_cart, 1, mpi_err)
        endif
        do n = 1,n_bfr
           write(unit,row_fmt) bfr_it(n), bfr(:,n)
        enddo
        close(unit)
        n_bfr = 0

        call EndProfRange

        return
end subroutine flush_probes

subroutine end_probes
        implicit none

        if(probes_on) then
          !$acc exit data delete(wl_i, wl_g, fl_i, fl_j, wl_s, fl_s)
        endif
        if(wm_comm /= MPI_COMM_NULL) then
          call MPI_File_close(wm_fh, mpi_err)
          call MPI_Type_free(wm_type, mpi_err)
          call MPI_Comm_free(wm_comm, mpi_err)
          !$acc exit data delete(wm_p)
        endif

        return
end subroutine end_probes

subroutine init_wall_map
        implicit none
        character(len=dl) :: fname
        integer :: color, wrank

        if(itWallMap <= 0) return

        select case(trim(bc(3)))
          case('adiabatic_wall', 'isothermal_wall', 'diabatic_wall')
          case default
            if(rank == root) print*, ' ERROR: itWallMap needs a no-slip wall on bc(3), found ', trim(bc(3))
            call secure_stop
        endselect

        wm_rec = 16 + 4*int(nx,8)*nz
        if(rank == root) write(*,'(A,I0,A,F0.1,A,F0.1,A)') ' Wall pressure map every ', itWallMap, ' it: ', &
                real(wm_rec,rp)/1.0e6_rp, ' MB per sample, ', &
                real(wm_rec,rp)*real(max(itmax-it,0)/itWallMap,rp)/1.0e9_rp, ' GB to itmax'

        color = MPI_UNDEFINED
        if(sy == 1) color = 1
        call MPI_Comm_split(mpi_comm_cart, color, rank, wm_comm, mpi_err)
        if(sy /= 1) return

        call MPI_Type_create_subarray(2, [nx, nz], [ex-sx+1, ez-sz+1], [sx-1, sz-1], &
                                      MPI_ORDER_FORTRAN, MPI_REAL, wm_type, mpi_err)
        call MPI_Type_commit(wm_type, mpi_err)
        allocate(wm_p(sx:ex,sz:ez))
        !$acc enter data create(wm_p)

        fname = 'DATA/'//trim(data_dir)//'/WALL_PRESSURE.bin'
        call MPI_Comm_rank(wm_comm, wrank, mpi_err)
        if(wrank == 0) call prepare_record_file(fname, [nx, nz], wm_rec, it, wm_off)
        call MPI_Bcast(wm_off, 1, MPI_OFFSET, 0, wm_comm, mpi_err)
        call MPI_File_open(wm_comm, fname, MPI_MODE_WRONLY, MPI_INFO_NULL, wm_fh, mpi_err)
        if(mpi_err /= MPI_SUCCESS) then
          print*, ' ERROR: cannot open ', trim(fname)
          call MPI_Abort(mpi_comm_cart, 1, mpi_err)
        endif
        call MPI_File_set_errhandler(wm_fh, MPI_ERRORS_ARE_FATAL, mpi_err)

        return
end subroutine init_wall_map

subroutine write_wall_map
        implicit none
        integer(8) :: it8
        real(8)    :: t8
        integer    :: i, k, wrank

        call StartProfRange('write_wall_map')

        if(wm_comm == MPI_COMM_NULL) then
          call EndProfRange
          return
        endif

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do k = sz,ez
           do i = sx,ex
              wm_p(i,k) = real(0.5_rp*(P(i,0,k) + P(i,1,k)),4)
           enddo
        enddo
        !$acc end parallel
        !$acc update host(wm_p)

        call MPI_Comm_rank(wm_comm, wrank, mpi_err)
        if(wrank == 0) then
          it8 = it
          t8  = time
          call MPI_File_write_at(wm_fh, wm_off    , it8, 1, MPI_INTEGER8        , MPI_STATUS_IGNORE, mpi_err)
          call MPI_File_write_at(wm_fh, wm_off + 8, t8 , 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpi_err)
        endif
        call MPI_File_set_view(wm_fh, wm_off + 16, MPI_REAL, wm_type, 'native', MPI_INFO_NULL, mpi_err)
        call MPI_File_write_all(wm_fh, wm_p, size(wm_p), MPI_REAL, MPI_STATUS_IGNORE, mpi_err)
        call MPI_File_set_view(wm_fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, mpi_err)
        wm_off = wm_off + wm_rec

        call EndProfRange

        return
end subroutine write_wall_map

subroutine gather_wall(xw, yw, nxw, nyw)
        implicit none
        real(rp), dimension(nx), intent(out) :: xw, yw, nxw, nyw
        integer :: i

        xw = 0.0_rp; yw = 0.0_rp; nxw = 0.0_rp; nyw = 0.0_rp
        if(sy == 1 .and. cart_coord(3) == 0) then
          do i = sx,ex
             xw(i)  = 0.5_rp*(xx(i,0) + xx(i,1))
             yw(i)  = 0.5_rp*(yy(i,0) + yy(i,1))
             nxw(i) = nx_eta(i,1)
             nyw(i) = ny_eta(i,1)
          enddo
        endif
        call MPI_allreduce(MPI_IN_PLACE,xw ,nx,MPI_RP,MPI_SUM,mpi_comm_cart,mpi_err)
        call MPI_allreduce(MPI_IN_PLACE,yw ,nx,MPI_RP,MPI_SUM,mpi_comm_cart,mpi_err)
        call MPI_allreduce(MPI_IN_PLACE,nxw,nx,MPI_RP,MPI_SUM,mpi_comm_cart,mpi_err)
        call MPI_allreduce(MPI_IN_PLACE,nyw,nx,MPI_RP,MPI_SUM,mpi_comm_cart,mpi_err)

        return
end subroutine gather_wall

subroutine check_body_wall(who, xw, yw)
        implicit none
        character(*)           , intent(in) :: who
        real(rp), dimension(nx), intent(in) :: xw, yw
        real(rp) :: ds
        integer  :: i

        select case(trim(bc(3)))
          case('adiabatic_wall', 'isothermal_wall', 'diabatic_wall')
          case default
            if(rank == root) print*, ' ERROR: ', who, ' needs a no-slip wall on bc(3), found ', trim(bc(3))
            call secure_stop
        endselect
        if(trim(bc(1)) /= 'periodic' .or. trim(bc(2)) /= 'periodic') then
          if(rank == root) print*, ' ERROR: ', who, ' needs a closed wall loop (bc(1) = bc(2) = periodic).'
          call secure_stop
        endif
        ds = 0.0_rp
        do i = 1,nx-1
           ds = max(ds, sqrt((xw(i+1)-xw(i))**2 + (yw(i+1)-yw(i))**2))
        enddo
        if(sqrt((xw(1)-xw(nx))**2 + (yw(1)-yw(nx))**2) > ds) then
          if(rank == root) print*, ' ERROR: ', who, ' needs a wall that closes on itself (O-grid body).'
          call secure_stop
        endif

        return
end subroutine check_body_wall

subroutine wall_chord(xw, yw, i_le, i_te, chord)
        implicit none
        real(rp), dimension(nx), intent(in)  :: xw, yw
        integer                , intent(out) :: i_le, i_te
        real(rp)               , intent(out) :: chord
        real(rp) :: d2, d2_max
        integer  :: i, j

        d2_max = -1.0_rp
        i_le = 1; i_te = 1
        do i = 1,nx
           do j = i+1,nx
              d2 = (xw(j)-xw(i))**2 + (yw(j)-yw(i))**2
              if(d2 > d2_max) then
                d2_max = d2; i_le = i; i_te = j
              endif
           enddo
        enddo
        if(xw(i_te) < xw(i_le)) then
          i = i_le; i_le = i_te; i_te = i
        endif
        chord = sqrt(d2_max)

        return
end subroutine wall_chord

function side_dir(xw, yw, i_le, i_te, side) result(dir)
        implicit none
        real(rp), dimension(nx), intent(in) :: xw, yw
        integer                , intent(in) :: i_le, i_te, side
        integer  :: dir, i
        real(rp) :: cx, cy, sd

        cx = xw(i_te) - xw(i_le)
        cy = yw(i_te) - yw(i_le)
        sd = 0.0_rp
        i  = mod(i_le,nx) + 1
        do while(i /= i_te)
           sd = sd + cx*(yw(i)-yw(i_le)) - cy*(xw(i)-xw(i_le))
           i  = mod(i,nx) + 1
        enddo
        dir = 1
        if(side*sd < 0.0_rp) dir = -1

        return
end function side_dir

end module output_probes
