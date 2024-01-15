module storage_module
use parameters_module
use mpi_module
use mesh_module
use allocate_module
use iso_c_binding
implicit none
integer            :: i,j,k,iv

! integers
integer :: err = 0
integer :: istop = 0
integer, parameter :: nvAve1D = 20
integer, parameter :: nVAve2D = 35
integer, parameter :: nVAve2D_aux = 20
integer, parameter :: nVWMLESData = 5

! time variables
real(rp)  :: time, t_restart

! conservative variables
real(rp), allocatable, dimension(:,:,:,:) :: phi 
real(rp), allocatable, dimension(:,:,:,:) :: RHS, phi_n

! Auxiliary variables
real(rp), allocatable, dimension(:,:,:) :: U          !< u-velocity
real(rp), allocatable, dimension(:,:,:) :: V          !< v-velocity
real(rp), allocatable, dimension(:,:,:) :: W          !< w-velocity
real(rp), allocatable, dimension(:,:,:) :: P          !< thermodynamical pressure
real(rp), allocatable, dimension(:,:,:) :: T          !< temperature
real(rp), allocatable, dimension(:,:,:) :: VIS        !< molecular viscosity
real(rp), allocatable, dimension(:,:,:) :: LMD        !< molecular diffusivity
real(rp), allocatable, dimension(:,:,:) :: DIV        !< velocity divergency
real(rp), allocatable, dimension(:,:,:) :: SSensor    !< shock sensor

real(rp), allocatable, dimension(:,:,:) :: tilde_op_x
real(rp), allocatable, dimension(:,:,:) :: tilde_op_y
real(rp), allocatable, dimension(:,:,:) :: tilde_op_z

real(rp), allocatable, dimension(:,:) :: pri_1D_x
real(rp), allocatable, dimension(:,:) :: pri_1D_y
real(rp), allocatable, dimension(:,:) :: pri_1D_z

real(rp), allocatable, dimension(:,:) :: phi_arr_x
real(rp), allocatable, dimension(:,:) :: phi_arr_y
real(rp), allocatable, dimension(:,:) :: phi_arr_z

real(rp), allocatable, dimension(:,:) :: flx_arr_x
real(rp), allocatable, dimension(:,:) :: flx_arr_y
real(rp), allocatable, dimension(:,:) :: flx_arr_z

real(rp), allocatable, dimension(:,:) :: flx_x
real(rp), allocatable, dimension(:,:) :: flx_y
real(rp), allocatable, dimension(:,:) :: flx_z

integer(1), allocatable, dimension(:) :: ishock_x
integer(1), allocatable, dimension(:) :: ishock_y
integer(1), allocatable, dimension(:) :: ishock_z

real(rp), allocatable, dimension(:) :: phirand, psirand
#ifdef AMD
real(c_double), dimension(:), allocatable :: psirandnum, phirandnum
type(c_ptr)                               :: psirandnumptr, phirandnumptr
integer(c_size_t)                         :: psirandsize, phirandsize
#endif

real(rp), allocatable, dimension(:) :: sponge_x


! statistics
real(rp), dimension(:,:)  , allocatable :: vmean1D
real(rp), dimension(:,:,:), allocatable :: vmean2D
real(rp), dimension(:,:,:), allocatable :: vmean2D_aux

! wall modelled stats
real(rp), dimension(:,:,:), allocatable :: WMLES_DATA_LW,WMLES_DATA_UW,WMLES_DATA
real(rp), dimension(:,:)  , allocatable :: vmean1D_wmles
real(rp), dimension(:)    , allocatable :: vmean0D_wmles

integer(1), allocatable, dimension(:,:,:) :: weno_flag
integer(1), allocatable, dimension(:,:,:) :: weno_flag_xyz

integer(1), parameter :: weno_smooth = 1
integer(1), parameter :: weno_shock  = 0

contains
subroutine init_GL_variables
        !=========================================================================================
        !
        !   ... ---o---o---o---o---o---      ...     ---o---o---o---o---o--- ...
        !   ...   -3  -2  -1   0   1                    n  n+1 n+2 n+3 n+4   ...
        !          |       |   ^   |----nodi interni----|   ^   |       |
        !          |-- GN--|   BC                           BC  |--GN --|
        !
        !=========================================================================================
        implicit none
        integer :: ltot

        ! azzero il tempo
        time = 0._rp
        t_restart = 0._rp
        
        ! evolutionary variables
        call AllocateReal(phi  ,lbx,ubx,lby,uby,lbz,ubz,1,5)
        call AllocateReal(RHS  ,lbx,ubx,lby,uby,lbz,ubz,1,5)
        call AllocateReal(phi_n,lbx,ubx,lby,uby,lbz,ubz,1,5)

        ! alloco le variabili di campo ausiliarie
        call AllocateReal(U  ,lbx,ubx,lby,uby,lbz,ubz)
        call AllocateReal(V  ,lbx,ubx,lby,uby,lbz,ubz)
        call AllocateReal(W  ,lbx,ubx,lby,uby,lbz,ubz)
        call AllocateReal(P  ,lbx,ubx,lby,uby,lbz,ubz)
        call AllocateReal(T  ,lbx,ubx,lby,uby,lbz,ubz)
        call AllocateReal(DIV,lbx,ubx,lby,uby,lbz,ubz)

        ltot = fd_order/2
        call AllocateReal(tilde_op_x,lbx,ubx,1,ltot,1,5)
        call AllocateReal(tilde_op_y,lby,uby,1,ltot,1,5)
        call AllocateReal(tilde_op_z,lbz,ubz,1,ltot,1,5)

        call AllocateReal(pri_1D_x,lbx,ubx,1,6)
        call AllocateReal(pri_1D_y,lby,uby,1,6)
        call AllocateReal(pri_1D_z,lbz,ubz,1,6)

        call AllocateReal(phi_arr_x,lbx,ubx,1,5)
        call AllocateReal(phi_arr_y,lby,uby,1,5)
        call AllocateReal(phi_arr_z,lbz,ubz,1,5)

        call AllocateReal(flx_arr_x,lbx,ubx,1,5)
        call AllocateReal(flx_arr_y,lby,uby,1,5)
        call AllocateReal(flx_arr_z,lbz,ubz,1,5)

        call AllocateReal(flx_x,lbx,ubx,1,5)
        call AllocateReal(flx_y,lby,uby,1,5)
        call AllocateReal(flx_z,lbz,ubz,1,5)

        call AllocateInteger(ishock_x,sx-1,ex)
        call AllocateInteger(ishock_y,sy-1,ey)
        call AllocateInteger(ishock_z,sz-1,ez)

        call AllocateReal(VIS,lbx,ubx,lby,uby,lbz,ubz)
        call AllocateReal(LMD,lbx,ubx,lby,uby,lbz,ubz)

        if(hybrid_weno) then
          call AllocateInteger(weno_flag,lbx,ubx,lby,uby,lbz,ubz)
          weno_flag = weno_smooth

          call AllocateInteger(weno_flag_xyz,lbx,ubx,lby,uby,lbz,ubz)
                
          call AllocateReal(SSensor,lbx,ubx,lby,uby,lbz,ubz)
        endif

        call AllocateReal(phirand,1,16)
        call AllocateReal(psirand,1,20)
#ifdef AMD
        allocate(phirandnum(16))
        allocate(psirandnum(20))
        phirandsize = 16
        psirandsize = 20
#endif

        if(wmles) then
          call AllocateReal(WMLES_DATA_LW,lbx,ubx,lbz,ubz,1,nvWmlesData)
          call AllocateReal(WMLES_DATA_UW,lbx,ubx,lbz,ubz,1,nvWmlesData)
          call AllocateReal(WMLES_DATA   ,lbx,ubx,lbz,ubz,1,nvWmlesData)
        endif

        return        
end subroutine init_GL_variables


subroutine init_Statistics_fields
! --------------------------------------------------------------------------------------------
!       This subroutine allocate reduced variables to exploit statistical periodic directions 
! --------------------------------------------------------------------------------------------

        implicit none
        character(1), dimension(3) :: pchar
        character(3)               :: pcharAll
        
        if(StFlg) then

          pchar = 'F'
          if(periods(1)) pchar(1) = 'T'
          if(periods(2)) pchar(2) = 'T'
          if(periods(3)) pchar(3) = 'T'

          pcharAll = pchar(1)//pchar(2)//pchar(3)

          selectcase(pcharAll)

          case('TFF') ! x   is periodic

          case('FTF') ! y   is periodic

          case('FFT') ! z   is periodic
            call AllocateReal(vmean2D    ,lbx,ubx,lby,uby,1,nvAve2D)
            call AllocateReal(vmean2D_aux,lbx,ubx,lby,uby,1,nvAve2D_aux)
            
            if(wmles) then
              call AllocateReal(vmean1D_wmles,lbx,ubx,1,nvWmlesData)
            endif

          case('TTF') ! xy are periodic

          case('FTT') ! yz are periodic
            call AllocateReal(vmean2D    ,lbx,ubx,lby,uby,1,nvAve2D)
            call AllocateReal(vmean2D_aux,lbx,ubx,lby,uby,1,nvAve2D_aux)

          case('TFT') ! xz are periodic
            call AllocateReal(vmean1D,lby,uby,1,nvAve1D)
            if(wmles) then
              call AllocateReal(vmean0D_wmles,1,nvWmlesData)
            endif

          endselect

        endif
        return
end subroutine init_statistics_fields




subroutine init_FD_coefficients
! ------------------------------------------------------------------------
!
!       This subroutine initializes finite difference coefficients for all
!       the solver. 
!       
! -----------------------------------------------------------------------
        implicit none

        central_fd_order = fd_order
          bward_fd_order = fd_order
          fward_fd_order = fd_order

        if    (fd_order == 2) then
                fd_L = 0
                fd_R = 1

        elseif(fd_order == 4) then
                fd_L = -1
                fd_R =  2

        elseif(fd_order == 6) then
                fd_L = -2
                fd_R =  3

        endif

        call AllocateReal(bward_1,-bward_fd_order,0)
        call AllocateReal(bward_2,-bward_fd_order-1,0)
        call AllocateReal(fward_1,0,fward_fd_order)
        call AllocateReal(fward_2,0,fward_fd_order+1)
        call AllocateReal(mid_point_lele,fd_L,fd_R)

        central_1 = 0.0_rp
        central_2 = 0.0_rp

        central_1_one_half = 0.0_rp
        central_2_one_half = 0.0_rp

        if    (fd_order == 2) then
                
                central_1(-1) = - 0.5_rp
                central_1( 0) =   0.0_rp
                central_1( 1) = - central_1(-1)

                central_2(-1) =   1.0_rp
                central_2( 0) = - 2.0_rp
                central_2( 1) =   central_2(-1)

                central_1_one_half(0) = 0.5_rp
                central_1_one_half(1) = central_1_one_half(0)

                mid_point_lele(0) = 0.5_rp
                mid_point_lele(1) = mid_point_lele(0)

                bward_1 = (/0.5_rp, -2._rp, 1.5_rp/)
                bward_2 = (/-1._rp, 4._rp, -5._rp, 2._rp/)

                fward_1 = (/-1.5_rp, 2._rp, -0.5_rp/)
                fward_2 = (/2._rp,-5._rp,4._rp,-1._rp/)

                central_2_one_half(0) = - 1.0_rp
                central_2_one_half(1) = - central_2_one_half(0)

        elseif(fd_order == 4) then

                central_1(-2) =   1.0_rp/12.0_rp
                central_1(-1) = - 2.0_rp/ 3.0_rp
                central_1( 0) =   0.0_rp
                central_1( 1) = - central_1(-1)
                central_1( 2) = - central_1(-2)

                central_2(-2) = - 1.0_rp/12.0_rp
                central_2(-1) =   4.0_rp/ 3.0_rp
                central_2( 0) = - 5.0_rp/ 2.0_rp
                central_2( 1) =   central_2(-1)
                central_2( 2) =   central_2(-2)

                central_1_one_half(-1) = - 1.0_rp/12.0_rp
                central_1_one_half( 0) =   7.0_rp/12.0_rp
                central_1_one_half( 1) = central_1_one_half( 0)
                central_1_one_half( 2) = central_1_one_half(-1)

                mid_point_lele(-1) = -1.0_rp/16.0_rp
                mid_point_lele( 0) = +9.0_rp/16.0_rp
                mid_point_lele( 1) = mid_point_lele( 0)
                mid_point_lele( 2) = mid_point_lele(-1)

                central_2_one_half(-1) =  1.0_rp/12.0_rp
                central_2_one_half( 0) = -5.0_rp/ 4.0_rp
                central_2_one_half( 1) = -central_2_one_half( 0)
                central_2_one_half( 2) = -central_2_one_half(-1)
                
                bward_1 = (/3._rp,-16._rp,+36._rp,-48._rp,+25._rp/)/(12._rp)
                bward_2 = (/-10._rp,+61._rp,-156._rp,+214._rp,-154._rp,+45._rp/)/(12._rp)

                fward_1 = (/-25._rp,+48._rp,-36._rp,+16._rp,-3._rp/)/(12._rp)
                fward_2 = (/45._rp,-154._rp,+214._rp,-156._rp,+61._rp,-10._rp/)/(12._rp)

        elseif(fd_order == 6) then

                central_1(-3) = - 1.0_rp/60.0_rp
                central_1(-2) =   3.0_rp/20.0_rp
                central_1(-1) = - 3.0_rp/ 4.0_rp
                central_1( 0) =   0.0_rp
                central_1( 1) = - central_1(-1)
                central_1( 2) = - central_1(-2)
                central_1( 3) = - central_1(-3)

                central_2(-3) =    1.0_rp/90.0_rp
                central_2(-2) = -  3.0_rp/20.0_rp
                central_2(-1) =    3.0_rp/ 2.0_rp
                central_2( 0) = - 49.0_rp/18.0_rp
                central_2( 1) = central_2(-1)
                central_2( 2) = central_2(-2)
                central_2( 3) = central_2(-3)

                central_1_one_half(-2) =   1.0_rp/60.0_rp
                central_1_one_half(-1) = - 2.0_rp/15.0_rp
                central_1_one_half( 0) =  37.0_rp/60.0_rp
                central_1_one_half( 1) = central_1_one_half( 0)
                central_1_one_half( 2) = central_1_one_half(-1)
                central_1_one_half( 3) = central_1_one_half(-2)

                mid_point_lele(-2) =  3.0_rp /256.0_rp
                mid_point_lele(-1) = -25.0_rp/256.0_rp
                mid_point_lele( 0) =  75.0_rp/128.0_rp
                mid_point_lele( 1) = mid_point_lele( 0)
                mid_point_lele( 2) = mid_point_lele(-1)
                mid_point_lele( 3) = mid_point_lele(-2)

                central_2_one_half(-2) = - 1.0_rp/90.0_rp
                central_2_one_half(-1) =   5.0_rp/36.0_rp
                central_2_one_half( 0) = -49.0_rp/36.0_rp
                central_2_one_half( 1) = - central_2_one_half( 0)
                central_2_one_half( 2) = - central_2_one_half(-1)
                central_2_one_half( 3) = - central_2_one_half(-2)

                bward_1 = (/1._rp/6._rp, -6._rp/5._rp, 15._rp/4._rp, -20._rp/3._rp, 15._rp/2._rp, -6._rp, 49._rp/20._rp/)
                bward_2 = (/-7._rp/10._rp, 1019._rp/180._rp, -201._rp/10._rp, 41._rp, -949._rp/18._rp, 879._rp/20._rp, &
                            -223._rp/10._rp, 469._rp/90._rp/)

                fward_1 = (/-49._rp/20._rp, 6._rp, -15._rp/2._rp, 20._rp/3._rp, -15._rp/4._rp,6._rp/5._rp, -1._rp/6._rp/)
                fward_2 = (/469._rp/90._rp, -223._rp/10._rp, 879._rp/20._rp, -949._rp/18._rp, 41._rp, -201._rp/10._rp, &
                            1019._rp/180._rp, -7._rp/10._rp/)


        else
                print*, ' Finite difference order ', fd_order, ' is not implemented.'
        endif

        if(abs(sum(central_1)) > toll_equality) stop 'central_1 FD coefficients are unbalanced'
        if(abs(sum(central_2)) > toll_equality) stop 'central_2 FD coefficients are unbalanced'

        if(abs(sum(central_1_one_half)-1.0_rp) > toll_equality) stop 'central_1_one_half are unbalanced'
        if(abs(sum(central_2_one_half)    ) > toll_equality) stop 'central_2_one_half are unbalanced'

        if(abs(sum(bward_1)) > toll_equality) stop 'bward_1 FD coefficients are unbalanced'
        if(abs(sum(bward_2)) > toll_equality) stop 'bward_2 FD coefficients are unbalanced'
        if(abs(sum(fward_1)) > toll_equality) stop 'fward_1 FD coefficients are unbalanced'
        if(abs(sum(fward_2)) > toll_equality) stop 'fward_2 FD coefficients are unbalanced'

        if(abs(sum(mid_point_lele)-1.0_rp) > toll_equality) stop 'mid_point_lele are unbalanced'

        return
end subroutine init_FD_coefficients


subroutine init_RK_coefficients
        implicit none
        
        allocate(a_rk(n_step), b_rk(n_step), c_rk(n_step), stat = err)

        if(err .ne. 0) then
          if(rank == root) print*, ' Allocation error in Runge-Kutta coefficients.'
          call secure_stop
        endif
        
        if(n_step == 3) then

                a_rk = (/0.0_rp, 0.75_rp, 1.0_rp/3.0_rp/)
                b_rk = (/1.0_rp, 0.25_rp, 2.0_rp/3.0_rp/)
                c_rk = (/0.0_rp, 0.50_rp, 0.50_rp   /)

        else
                if(rank == root) print*, ' Runge Kutta coefficients for ', n_step, ' steps are not implemented'
                call secure_stop
        endif
        return
end subroutine init_RK_coefficients



subroutine init_sponge()

  implicit none

  real(rp), parameter :: Ls = 1.0_rp
  real(rp), parameter :: As = 20.0_rp
  real(rp)            :: xini, xend, i_sizex, sponge
  integer             :: i

  call AllocateReal(sponge_x,lbx,ubx)

  xini = xmax - Ls
  xend = xmax
  i_sizex = 1.0_rp/(xend - xini)

  do i = lbx, ubx
    sponge_x(i) = 0.0_rp
    if(x(i) > xini) then

      sponge = (x(i) - xini)*i_sizex
      sponge = sponge * sponge
      sponge_x(i) = As * sponge * sponge

    end if
  enddo

end subroutine init_sponge




subroutine end_all_variables
        implicit none

        deallocate(bward_1  , bward_2          , &
                   fward_1  , fward_2          , &

                   a_rk, b_rk, c_rk            , &

                   phi, RHS, phi_n             , &
                   U, V, W, P, T, stat = err)

        if(err .ne. 0) then
          if(rank == root) print*, " WARNING: Deallocation error in end_all_variables"
        endif

        if(allocated(x)) deallocate(x)
        if(allocated(y)) deallocate(y)
        if(allocated(z)) deallocate(z)

        return
end subroutine end_all_variables












end module storage_module
