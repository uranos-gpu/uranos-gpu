module df_module
use parameters_module, only: rp
use mpi_module, only: sy,ey,sz,ez
use random_module
use profiling_module

implicit none

integer, parameter    , public :: DF_N      = 60
real(rp), dimension(3), public :: DF_dimMin = (/0.05_rp, 0.05_rp, 0.05_rp/)
real(rp), dimension(3), public :: DF_dimMax = (/0.25_rp, 0.20_rp, 0.25_rp/)

real(rp), dimension(:,:)  , allocatable, public :: DF_ylen
real(rp), dimension(:,:)  , allocatable, public :: DF_zlen
real(rp), dimension(:)    , allocatable, public :: DF_uconv
real(rp), dimension(:,:,:), allocatable, public :: DF_By
real(rp), dimension(:,:,:), allocatable, public :: DF_Bz

real(rp), dimension(:,:,:)  , allocatable, public :: DF_RnD2D
integer(8)                               , public :: DF_ndraw = 0
real(rp), dimension(:,:,:)  , allocatable, public :: DF_LundMatrix
real(rp), dimension(:,:,:)  , allocatable, public :: DF_fy

real(rp) :: DF_xlen1 = 0.50_rp
real(rp) :: DF_xlen2 = 0.15_rp
real(rp) :: DF_xlen3 = 0.15_rp

contains

subroutine DFInputData_TBL(ReTau,Mach,Prandtl,y_gbl,ny)

        use parameters_module     , only: data_dir, Trat
        use mpi_module            , only: rank
        use real_to_integer_module, only: locate
        use interpolation_module  , only: polint
        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection
#ifdef DEBUG
        use FileModule
#endif

        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: y_gbl
        real(rp)                           , intent(in)    :: Mach, Prandtl, ReTau
        integer                            , intent(in)    :: ny

        character(50), dimension(:)    , allocatable :: dtFile
        real(rp)     , dimension(:)    , allocatable :: DtReth
        real(rp)     , dimension(:,:,:), allocatable :: DtData, DtDataInt
        real(rp)     , dimension(:,:,:), allocatable :: ReyStress
        real(rp)     , dimension(:)    , allocatable :: tmp0, tmp1, tmp2

        real(rp)     , dimension(0:ny) :: UPlusInc,UplusCmp,RhoY, TmpY
        real(rp)     , dimension(3,3)  :: stress, A
        character(50)                  :: database, lFile
        integer                        :: dtUnit = 12, err = 0

        integer                        :: l, lines, m, j, jj, jjj, ltau, lltau
        real(rp)                       :: yy, uui,vvi,wwi,uvi, duu, dvv, dww, duv
        real(rp)                       :: rWall, tWall, rRatio, uEdge
#ifdef DEBUG
        type(FileType) :: DF_inputs, LundFile
#endif

        selectcase(nint(Retau))

          case(0:1272)
            database   = trim('libs/database_schlatter/')

            lines = 512
            allocate(dtFile(1:10), dtReTh(1:10))
            allocate(DtData(10,7,lines), stat=err)
            if(err .ne. 0) stop ' Allocation error in ComputeLundMatrixFromDataBase'

            dtFile( 1) = trim('vel_0670_dns.prof')
            dtFile( 2) = trim('vel_1000_dns.prof')
            dtFile( 3) = trim('vel_1410_dns.prof')
            dtFile( 4) = trim('vel_2000_dns.prof')
            dtFile( 5) = trim('vel_2540_dns.prof')
            dtFile( 6) = trim('vel_3030_dns.prof')
            dtFile( 7) = trim('vel_3270_dns.prof')
            dtFile( 8) = trim('vel_3630_dns.prof')
            dtFile( 9) = trim('vel_3970_dns.prof')
            dtFile(10) = trim('vel_4060_dns.prof')

          case(1273:)
            stop ' Empy database for Reynolds stresses at Retau > 1273.0_rp'

        endselect
        DtData = 0.0_rp
        do l = 1, size(dtFile)

           lFile = trim(database)//trim(dtFile(l))
           open(unit = dtUnit, file = lFile, status = 'old', iostat = err)
           if(err .ne. 0) stop ' Error opening Reynolds stress file'

           read(dtUnit,*) dtReTh(l)
           do j = 1,lines
              read(dtUnit,*,iostat=err) dtData(l,:,j)
              if(err .ne. 0) stop ' Error in DFInputData'
           enddo

           do j = 1, lines
              dtData(l,4,j) = dtData(l,4,j)**2/dtData(l,3,lines)**2
              dtData(l,5,j) = dtData(l,5,j)**2/dtData(l,3,lines)**2
              dtData(l,6,j) = dtData(l,6,j)**2/dtData(l,3,lines)**2
              dtData(l,7,j) = dtData(l,7,j)   /dtData(l,3,lines)**2
           enddo

           close(dtUnit)

        enddo

        allocate(tmp0(1:lines), dtDataInt(10,7,0:ny), stat=err)
        if(err .ne. 0) stop ' Allocatation error in DFInputData'

        m = 2
        do    l = 1,size(dtFile)
           do j = 0,ny

              yy = y_gbl(j)
              dtDataInt(l,1,j) = yy

              tmp0 = dtData(l,1,:)
              call locate(tmp0,1,lines,yy,jj)
              jjj = min(max(jj-(m-1)/2,1),lines+1-m)

              allocate(tmp1(jjj:jjj+m-1), tmp2(jjj:jjj+m-1))
              tmp1 = dtData(l,1,jjj:jjj+m-1)

              tmp2 = dtData(l,4,jjj:jjj+m-1)
              call polint(tmp1,tmp2,m,yy,uui,duu)

              tmp2 = dtData(l,5,jjj:jjj+m-1)
              call polint(tmp1,tmp2,m,yy,vvi,dvv)

              tmp2 = dtData(l,6,jjj:jjj+m-1)
              call polint(tmp1,tmp2,m,yy,wwi,dww)

              tmp2 = dtData(l,7,jjj:jjj+m-1)
              call polint(tmp1,tmp2,m,yy,uvi,duv)

              dtDataInt(l,4,j) = uui
              dtDataInt(l,5,j) = vvi
              dtDataInt(l,6,j) = wwi
              dtDataInt(l,7,j) = uvi

              deallocate(tmp1, tmp2)

           enddo
        enddo
        deallocate(tmp0)

        call MuskerProfile(ny,y_gbl,ReTau,UplusInc)
        call CompressibleCorrection(ny,UPlusInc,UPlusCmp,RhoY,TmpY,tWall,rWall)

        uEdge = maxval(UPlusCmp)
        do j = 0,ny
           DF_uconv(j) = min(max(UPlusCmp(j),10.0_rp)/uEdge,1.0_rp)
        enddo

        m = 1
        call locate(dtReTh,1,size(dtFile),ReTau,ltau)
        lltau = min(max(ltau-(m-1)/2,1),size(dtFile)+1-m)

        allocate(ReyStress(3,3,-3:ny+4),stat=err)
        if(err.ne.0) stop ' Allocation error in Reynolds Stress'

        ReyStress = 0.0_rp
        do j = 0,ny

           call polint(dtReTh(lltau),dtDataInt(lltau:lltau+m-1,4,j),m,retau,uui,duu)
           call polint(dtReTh(lltau),dtDataInt(lltau:lltau+m-1,5,j),m,retau,vvi,dvv)
           call polint(dtReTh(lltau),dtDataInt(lltau:lltau+m-1,6,j),m,retau,wwi,dww)
           call polint(dtReTh(lltau),dtDataInt(lltau:lltau+m-1,7,j),m,retau,uvi,duv)

           rRatio      = rWall/RhoY(j)

           Stress(1,:) = [uui, uvi, 0.0_rp]*rRatio
           Stress(2,:) = [uvi, vvi, 0.0_rp]*rRatio
           Stress(3,:) = [0.0_rp, 0.0_rp, wwi]*rRatio

           ReyStress(:,:,j) = Stress(:,:)

        enddo
        allocate(DF_LundMatrix(3,3,0:ny),stat=err)
        if(err .ne. 0) stop ' Allocation error of Lunds Matrix'

        DF_LundMatrix = 0.0_rp
        do j = 0,ny

          if(ReyStress(1,1,j) < 1.0E-08_rp) then
            a(:,:) = 0.0_rp
          else

            a(1,1) = sqrt(abs(ReyStress(1,1,j)))
            a(1,2) = 0.0_rp
            a(1,3) = 0.0_rp
            a(2,1) = ReyStress(2,1,j)/a(1,1)
            a(2,2) = sqrt(abs(ReyStress(2,2,j) - a(2,1)**2))
            a(2,3) = 0.0_rp
            a(3,1) = ReyStress(3,1,j)/a(1,1)
            a(3,2) = (ReyStress(3,2,j) - a(2,1)*a(3,1))/a(2,2)
            a(3,3) = sqrt(abs(ReyStress(3,3,j) - a(3,1)**2 - a(3,2)**2))

          endif

          DF_LundMatrix(:,:,j) = a(:,:)

        enddo

        allocate(DF_fy(3, sy-DF_N:ey+DF_N, sz-DF_N:ez+DF_N), stat = err)
        if(err .ne. 0) stop ' Allocation error of DF_fy'

#ifdef DEBUG
        DF_inputs%name ='ReynoldsStress'
        DF_inputs%dir  = trim(data_dir)//'/DF_INPUTS'
        call OpenNewFile(DF_inputs,0)
        do j = 0,ny
           write(df_inputs%unit,*) y_gbl(j), ReyStress(1,1,j), &
           ReyStress(2,2,j), ReyStress(3,3,j), ReyStress(1,2,j)
        enddo
        call CloseFile(DF_inputs)

        LundFile%name = 'lundMatrix_'//trim(str(rank))
        LundFile%dir  = trim(data_dir)//'/DF_INPUTS'
        call OpenNewFile(lundFile,0)
        do j = 0,ny
           write(lundFile%unit,'(7e18.6)') y_gbl(j), DF_LundMatrix(1,1,j), &
                                        DF_LundMatrix(2,2,j), &
                                        DF_LundMatrix(3,3,j), &
                                        DF_LundMatrix(2,1,j), &
                                        DF_LundMatrix(3,1,j), &
                                        DF_LundMatrix(3,2,j)
        enddo
        call CloseFile(LundFile)
#endif

        deallocate(dtFile, dtReth, dtData, dtDataInt, ReyStress)

        return
end subroutine DFInputData_TBL

subroutine DFRandomField2D(ny,nz)
        use random_module, only: philox_normal4

        implicit none
        integer     , intent(in)    :: ny,nz

        real(rp)  :: mean_x, mean_y, mean_z
        real(rp)  :: rmsq_x, rmsq_y, rmsq_z
        real(rp)  :: irms_x, irms_y, irms_z
        integer   :: j,k,l
        real(rp)  :: irNz, rnd1, rnd2, rnd3, rnd4
        integer(8):: ndrw_lo, ndrw_hi

        call StartProfRange('DFRandomField2D')

        ndrw_lo = iand(DF_ndraw, 4294967295_8)
        ndrw_hi = iand(ishft(DF_ndraw,-32), 4294967295_8)

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do    k = 1,nz
           do j = 1-DF_N,ny+DF_N

              call philox_normal4(int(j+DF_N,8), int(k,8), ndrw_lo, ndrw_hi, &
                                  rnd1, rnd2, rnd3, rnd4)

              DF_Rnd2D(1,j,k) = rnd1
              DF_Rnd2D(2,j,k) = rnd2
              DF_Rnd2D(3,j,k) = rnd3

           enddo
        enddo
        !$acc end parallel

        DF_ndraw = DF_ndraw + 1_8

        irNz = 1.0_rp/real(nz,rp)

        !$acc parallel default(present)
        !$acc loop gang, vector
        do j = 1-DF_N,ny+DF_N

           mean_x = 0.0_rp
           mean_y = 0.0_rp
           mean_z = 0.0_rp
           rmsq_x = 0.0_rp
           rmsq_y = 0.0_rp
           rmsq_z = 0.0_rp

           !$acc loop seq
           do k = 1,nz

              rnd1 = DF_Rnd2D(1,j,k)
              rnd2 = DF_Rnd2D(2,j,k)
              rnd3 = DF_Rnd2D(3,j,k)

              mean_x = mean_x + rnd1
              mean_y = mean_y + rnd2
              mean_z = mean_z + rnd3
              rmsq_x = rmsq_x + rnd1*rnd1
              rmsq_y = rmsq_y + rnd2*rnd2
              rmsq_z = rmsq_z + rnd3*rnd3
           enddo
           mean_x = mean_x*irNz
           mean_y = mean_y*irNz
           mean_z = mean_z*irNz
           rmsq_x = rmsq_x*irNz
           rmsq_y = rmsq_y*irNz
           rmsq_z = rmsq_z*irNz
           irms_x  = 1.0_rp/(sqrt(rmsq_x - mean_x*mean_x))
           irms_y  = 1.0_rp/(sqrt(rmsq_y - mean_y*mean_y))
           irms_z  = 1.0_rp/(sqrt(rmsq_z - mean_z*mean_z))

           !$acc loop seq
           do k = 1,nz
              DF_Rnd2D(1,j,k) = (DF_Rnd2D(1,j,k) - mean_x)*irms_x
              DF_Rnd2D(2,j,k) = (DF_Rnd2D(2,j,k) - mean_y)*irms_y
              DF_Rnd2D(3,j,k) = (DF_Rnd2D(3,j,k) - mean_z)*irms_z
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = 1,DF_N
           do    j = 1-DF_N,ny+DF_N
              do l = 1,3
                 DF_Rnd2D(l,j,nz+k) = DF_Rnd2D(l,j,k)
                 DF_Rnd2D(l,j,1 -k) = DF_Rnd2D(l,j,nz+1-k)
              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine DFRandomField2D

subroutine DFIntegralLenght_TBL(y_gbl,ny)
        use parameters_module, only: data_dir, ReTau, ymin
        use mpi_module       , only: rank
#ifdef DEBUG
        use FileModule
#endif

        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: y_gbl
        integer                            , intent(in)    :: ny

        real(rp) :: zlenou_u, zlenou_v, zlenou_w
        real(rp) :: zlenin_u, zlenin_v, zlenin_w
        real(rp) :: ftany, ylim
        integer  :: j

#ifdef DEBUG
        type(FileType) :: IntLenFile
#endif

        zlenou_u = 0.40_rp
        zlenou_v = 0.30_rp
        zlenou_w = 0.40_rp
        zlenin_u = min(150.0_rp/ReTau,zlenou_u)
        zlenin_v = min( 75.0_rp/ReTau,zlenou_v)
        zlenin_w = min(150.0_rp/ReTau,zlenou_w)

        do j=0,ny
           ftany = 0.5_rp*(1.0_rp+tanh((y_gbl(j)-0.2_rp)/0.03_rp))
           DF_zlen(1,j)  = zlenin_u+ftany*(zlenou_u-zlenin_u)
           DF_zlen(2,j)  = zlenin_v+ftany*(zlenou_v-zlenin_v)
           DF_zlen(3,j)  = zlenin_w+ftany*(zlenou_w-zlenin_w)

        enddo
        do j = 0,ny
           ylim = max(0.41_rp*max(y_gbl(j)-ymin,0.0_rp), 15.0_rp/ReTau)
           DF_ylen(:,j) = min(0.7_rp*DF_zlen(:,j), ylim)
        enddo

        DF_xlen1 = 2.40_rp
        DF_xlen2 = 0.30_rp
        DF_xlen3 = 0.30_rp

#ifdef DEBUG
        intLenFile%name = 'integralLen'//trim(str(rank))
        intLenFile%dir  = trim(data_dir)//'/DF_INTEGRAL_LENGTH'
        call OpenNewFile(intLenFile,0)
        do j = 0,ny
           write(intLenFile%unit,*) y_gbl(j), DF_ylen(:,j), DF_zlen(:,j)
        enddo
        call CloseFile(intLenFile)
#endif
        return
end subroutine DFIntegralLenght_TBL

subroutine DFCoefficients(y_gbl,ny,gbl_min_step)
        use parameters_module, only: rp, pi, data_dir

#ifdef DEBUG
        use FileModule
        use mpi_module, only: rank
#endif

        implicit none
        real(rp)    , dimension(:), allocatable, intent(in)    :: y_gbl
        real(rp)    , dimension(:)             , intent(in)    :: gbl_min_step
        integer                                , intent(in)    :: ny

        real(rp), dimension(3) :: sumby, sumbz
        real(rp)               :: dy, dz
        integer  :: m, j,jj
#ifdef DEBUG
        type(FileType) :: DFCoeffFIle
#endif

        do j = 1,ny

           dy = y_gbl(j) - y_gbl(j-1)
           dz = gbl_min_step(3)

           sumby = 0.0_rp
           do jj = -DF_N,DF_N
              DF_By(:,jj,j) = exp(-pi*abs(jj)/(DF_yLen(:,j)/dy))

              sumby(:) = sumby(:) + DF_By(:,jj,j)**2
           enddo
           do m = 1,3
              DF_By(m,:,j) = DF_By(m,:,j)/sqrt(sumby(m))
           enddo

           sumbz = 0.0_rp
           do jj = -DF_N,DF_N
              DF_Bz(:,jj,j) = exp(-pi*abs(jj)/(DF_zLen(:,j)/dz))

              sumbz(:) = sumbz(:) + DF_Bz(:,jj,j)**2
           enddo
           do m = 1,3
              DF_Bz(m,:,j) = DF_Bz(m,:,j)/sqrt(sumbz(m))
           enddo

        enddo

#ifdef DEBUG
        DFCoeffFile%name = 'Bcoefficient_'//trim(str(rank))
        DFCoeffFIle%dir  = trim(data_dir)//'/DF_COEFFICIENT'
        call OpenNewFile(DFCoeffFile,0)
        do j = 0, ny
           do jj = -DF_N,DF_N
              write(DFCoeffFile%unit,*) y_gbl(j), jj, DF_By(:,jj,j), DF_Bz(:,jj,j)
           enddo
        enddo
        call CloseFile(DFCoeffFile)
#endif

        return
end subroutine DFcoefficients

subroutine DFConvolution2D(sy,ey,sz,ez,vf)
        implicit none

        integer                                , intent(in)    :: sy, ey, sz, ez
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: vf

        integer :: l, j,k,jj,kk

        call StartProfRange('DFConvolution2D')

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz-DF_N,ez+DF_N
           do    j = sy-DF_N,ey+DF_N
              do l = 1,3
                 DF_fy(l,j,k) = 0.0_rp
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do          k = sz-DF_N,ez+DF_N
           do       j = sy,ey
              do    l = 1,3
                 do jj = -DF_N,DF_N
                    DF_fy(l,j,k) = DF_fy(l,j,k) + DF_By(l,jj,j)*DF_RnD2D(l,j+jj,k)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do    k = sz,ez
           do j = sy,ey
              do l = 1,3
                 vf(l,j,k) = 0.0_rp
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do          k = sz,ez
           do       j = sy,ey
              do    l = 1,3
                 do kk = -DF_N,DF_N
                    vf(l,j,k) = vf(l,j,k) + DF_Bz(l,kk,j)*DF_fy(l,j,k+kk)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine DFConvolution2D

subroutine DFCastroTimeCorrelation(ik,c_rk,u_inf,dt,sy,ey,sz,ez,vf_old,vf_new)

        use parameters_module, only: rp, pi

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: vf_new, vf_old
        real(rp), dimension(3)                 , intent(in)    :: c_rk
        real(rp)                               , intent(in)    :: u_inf, dt
        integer                                , intent(in)    :: ik
        integer                                , intent(in)    :: sy,ey, sz, ez

        real(rp) :: cx, cy, cz, uc
        real(rp) :: e1x, e1y, e1z, sqx, sqy, sqz
        real(rp) :: dtOld
        integer  :: j,k,ikk

        call StartProfRange('turbulent_inflow')

        dtold = 1.0_rp/3.0_rp*dt

        cx = -0.5_rp*pi*dtOld/DF_xlen1
        cy = -0.5_rp*pi*dtOld/DF_xlen2
        cz = -0.5_rp*pi*dtOld/DF_xlen3

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do k    = sz,ez
           do j = sy,ey
              uc  = u_inf*DF_uconv(j)
              e1x = exp(cx*uc)
              e1y = exp(cy*uc)
              e1z = exp(cz*uc)
              sqx = sqrt(1._rp-e1x*e1x)
              sqy = sqrt(1._rp-e1y*e1y)
              sqz = sqrt(1._rp-e1z*e1z)

              vf_new(1,j,k) = vf_old(1,j,k)*e1x + vf_new(1,j,k)*sqx
              vf_new(2,j,k) = vf_old(2,j,k)*e1y + vf_new(2,j,k)*sqy
              vf_new(3,j,k) = vf_old(3,j,k)*e1z + vf_new(3,j,k)*sqz

              vf_old(1,j,k) = vf_new(1,j,k)
              vf_old(2,j,k) = vf_new(2,j,k)
              vf_old(3,j,k) = vf_new(3,j,k)
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine DFCastroTimeCorrelation

subroutine DFEnforceReynoldsStresses2D(sy,ey,sz,ez,vf,uf)

        implicit none
        integer                                , intent(in)    :: sy, ey, sz, ez
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: vf
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: uf

        integer :: j, k

        call StartProfRange('DFEnforceReynoldsStresses2D')

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do    k = sz,ez
           do j = sy,ey
              uf(1,j,k) = 0.0_rp
              uf(2,j,k) = 0.0_rp
              uf(3,j,k) = 0.0_rp
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do    k = sz,ez
           do j = sy,ey
              uf(1,j,k) = uf(1,j,k) + DF_LundMatrix(1,1,j)*vf(1,j,k)
              uf(2,j,k) = uf(2,j,k) + DF_LundMatrix(2,1,j)*vf(1,j,k)
              uf(3,j,k) = uf(3,j,k) + DF_LundMatrix(3,1,j)*vf(1,j,k)

              uf(1,j,k) = uf(1,j,k) + DF_LundMatrix(1,2,j)*vf(2,j,k)
              uf(2,j,k) = uf(2,j,k) + DF_LundMatrix(2,2,j)*vf(2,j,k)
              uf(3,j,k) = uf(3,j,k) + DF_LundMatrix(3,2,j)*vf(2,j,k)

              uf(1,j,k) = uf(1,j,k) + DF_LundMatrix(1,3,j)*vf(3,j,k)
              uf(2,j,k) = uf(2,j,k) + DF_LundMatrix(2,3,j)*vf(3,j,k)
              uf(3,j,k) = uf(3,j,k) + DF_LundMatrix(3,3,j)*vf(3,j,k)
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine DFEnforceReynoldsStresses2D

end module df_module
