module df_module
use parameters_module, only: rp
use mpi_module, only: sy,ey,sz,ez
use random_module
use profiling_module
#ifdef _OPENACC
#ifdef NVIDIA
use curand
use curand_device
use openacc_curand
#endif
#ifdef AMD
use iso_c_binding
use hipfort_rocrand
#endif
#endif

implicit none

integer, parameter    , public :: DF_N      = 60 !30                !< filter width
real(rp), dimension(3), public :: DF_dimMin = (/0.05_rp, 0.05_rp, 0.05_rp/)  !< filter min Length
real(rp), dimension(3), public :: DF_dimMax = (/0.25_rp, 0.20_rp, 0.25_rp/)  !< filter max Length

real(rp), dimension(:,:)  , allocatable, public :: DF_ylen
real(rp), dimension(:,:)  , allocatable, public :: DF_zlen
real(rp), dimension(:,:,:), allocatable, public :: DF_By
real(rp), dimension(:,:,:), allocatable, public :: DF_Bz

real(rp), dimension(:,:,:)  , allocatable, public :: DF_RnD2D
real(rp), dimension(:,:,:,:), allocatable, public :: DF_RnD3D
real(rp), dimension(:,:,:)  , allocatable, public :: DF_LundMatrix
real(rp), dimension(:,:,:)  , allocatable, public :: DF_fy

real(rp) :: DF_xlen1 = 0.50_rp
real(rp) :: DF_xlen2 = 0.15_rp
real(rp) :: DF_xlen3 = 0.15_rp

#ifdef AMD
real(c_double), dimension(:), allocatable, public :: rnd2Dnum
type(c_ptr)                              , public :: rnd2Dnumptr
integer(c_size_t)                        , public :: rnd2Dnumsize
integer       , dimension(3)             , public :: rnd2Dnumshape
#endif

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

        ! local declarations
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
        real(rp)                       :: rWall, tWall, rRatio
#ifdef DEBUG
        type(FileType) :: DF_inputs, LundFile
#endif

        ! 
        ! === select the database
        !
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
        !
        ! === read the data  in the database file
        !
        DtData = 0.0_rp
        do l = 1, size(dtFile)
        
           ! open the file
           lFile = trim(database)//trim(dtFile(l))
           open(unit = dtUnit, file = lFile, status = 'old', iostat = err)
           if(err .ne. 0) stop ' Error opening Reynolds stress file'
              
           ! read data
           read(dtUnit,*) dtReTh(l)
           do j = 1,lines
              read(dtUnit,*,iostat=err) dtData(l,:,j)
              if(err .ne. 0) stop ' Error in DFInputData'
           enddo
                
           ! rescale some data
           do j = 1, lines
              dtData(l,4,j) = dtData(l,4,j)**2/dtData(l,3,lines)**2
              dtData(l,5,j) = dtData(l,5,j)**2/dtData(l,3,lines)**2
              dtData(l,6,j) = dtData(l,6,j)**2/dtData(l,3,lines)**2
              dtData(l,7,j) = dtData(l,7,j)   /dtData(l,3,lines)**2
           enddo

           close(dtUnit)

        enddo

        !
        ! === interpolation of the database on the global grid
        !
        allocate(tmp0(1:lines), dtDataInt(10,7,0:ny), stat=err)
        if(err .ne. 0) stop ' Allocatation error in DFInputData'

        m = 2 ! order of interpolation
        do    l = 1,size(dtFile)
           do j = 0,ny
        
              ! find the location
              yy = y_gbl(j)
              dtDataInt(l,1,j) = yy
        
              tmp0 = dtData(l,1,:)
              call locate(tmp0,1,lines,yy,jj)
              jjj = min(max(jj-(m-1)/2,1),lines+1-m)
        
              ! data interpolation
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

              ! save to DataInt
              dtDataInt(l,4,j) = uui
              dtDataInt(l,5,j) = vvi
              dtDataInt(l,6,j) = wwi
              dtDataInt(l,7,j) = uvi

              deallocate(tmp1, tmp2)

           enddo
        enddo
        deallocate(tmp0)

        !
        ! === interpolation of the Reynolds stress for the actual ReTau
        !
        ! getting the inflow nominal profile
        call MuskerProfile(ny,y_gbl,ReTau,UplusInc)
        call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)
        
        ! get the position in the database
        m = 1 ! order of interpolation (ERA posto a 3)
        call locate(dtReTh,1,size(dtFile),ReTau,ltau) 
        lltau = min(max(ltau-(m-1)/2,1),size(dtFile)+1-m)
        
        ! interpolate the database on the actual ReTau
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
        !
        ! === filling the lund Matrix inside every proc
        !
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
            !
            a(2,1) = ReyStress(2,1,j)/a(1,1)
            a(2,2) = sqrt(abs(ReyStress(2,2,j) - a(2,1)**2))
            a(2,3) = 0.0_rp
            !
            a(3,1) = ReyStress(3,1,j)/a(1,1)
            a(3,2) = (ReyStress(3,2,j) - a(2,1)*a(3,1))/a(2,2)
            a(3,3) = sqrt(abs(ReyStress(3,3,j) - a(3,1)**2 - a(3,2)**2))

          endif
        
          DF_LundMatrix(:,:,j) = a(:,:)

        enddo

        allocate(DF_fy(3, sy-DF_N:ey+DF_N, sz-DF_N:ez+DF_N), stat = err)
        if(err .ne. 0) stop ' Allocation error of DF_fy'

#ifdef DEBUG
        ! write input Reynolds stresses
        DF_inputs%name ='ReynoldsStress'
        DF_inputs%dir  = trim(data_dir)//'/DF_INPUTS'
        call OpenNewFile(DF_inputs,0)
        do j = 0,ny
           write(df_inputs%unit,*) y_gbl(j), ReyStress(1,1,j), &
           ReyStress(2,2,j), ReyStress(3,3,j), ReyStress(1,2,j)
        enddo
        call CloseFile(DF_inputs)

        ! write Lund's Matrix components
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


subroutine DFInputData_HTURB(y_gbl,ny)

#ifdef DEBUG
        use parameters_module, only: data_dir
        use mpi_module       , only: rank
        use FileModule
#endif
        use parameters_module, only: u_inf, turbulent_intensity

        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: y_gbl
        integer                            , intent(in)    :: ny

        ! local declarations
        real(rp), dimension(:,:,:), allocatable :: ReyStress
        real(rp), dimension(3,3)                :: a

        real(rp) :: trb_int
        integer  :: j, err = 0
#ifdef DEBUG
        type(FileType) :: DF_inputs, LundFile
#endif

        trb_int = turbulent_intensity * u_inf**2

        !
        ! ==== compute Reynolds stresses
        !
        allocate(ReyStress(3,3,-3:ny+4), stat = err)
        if(err.ne.0) stop ' Allocation error in DFInputData_HTURB'

        ReyStress = 0.0_rp
        do j = 0,ny

           ReyStress(1,1,j) = trb_int
           ReyStress(2,2,j) = trb_int
           ReyStress(3,3,j) = trb_int

        enddo
        !
        ! ==== filling lund matrix
        !
        allocate(DF_LundMatrix(3,3,0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in DFInputData_HTURB'

        DF_LundMatrix = 0.0_rp
        do j = 0,ny
                
           a(1,1) = sqrt(ReyStress(1,1,j))
           a(1,2) = 0.0_rp
           a(1,3) = 0.0_rp
           !
           a(2,1) = ReyStress(2,1,j)/a(1,1)
           a(2,2) = sqrt(ReyStress(2,2,j) - a(2,1)**2)
           a(2,3) = 0.0_rp
           !
           a(3,1) = ReyStress(3,1,j)/a(1,1)
           a(3,2) = (ReyStress(3,2,j) - a(2,1)*a(3,1))/a(2,2)
           a(3,3) = sqrt(ReyStress(3,3,j) - a(3,1)**2 - a(3,2)**2)

           DF_LundMAtrix(:,:,j) = a(:,:)

        enddo

#ifdef DEBUG
        ! write input Reynolds stresses
        DF_inputs%name ='ReynoldsStress'
        DF_inputs%dir  = trim(data_dir)//'/DF_INPUTS'
        call OpenNewFile(DF_inputs,0)
        do j = 0,ny
           write(df_inputs%unit,*) y_gbl(j), ReyStress(1,1,j), &
           ReyStress(2,2,j), ReyStress(3,3,j), ReyStress(1,2,j)
        enddo
        call CloseFile(DF_inputs)

        ! write Lund's Matrix components
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
        
        deallocate(ReyStress)

        allocate(DF_fy(3, sy-DF_N:ey+DF_N, sz-DF_N:ez+DF_N), stat = err)
        if(err .ne. 0) stop ' Allocation error of DF_fy'

        return
end subroutine DFInputData_HTURB


subroutine DFRandomField2D(ny,nz)

        use random_module
#ifdef AMD
        use iso_c_binding
#endif

        implicit none
        integer     , intent(in)    :: ny,nz
        
        ! local 
        real(rp)  :: mean_x, mean_y, mean_z
        real(rp)  :: rmsq_x, rmsq_y, rmsq_z
        real(rp)  :: rms_x, rms_y, rms_z
        integer   :: m,i,j,k,l
        real(rp)  :: irNz

        real(rp), parameter :: s  = 0.449871_rp
        real(rp), parameter :: t  = -0.386595_rp
        real(rp), parameter :: a  = 0.19600_rp
        real(rp), parameter :: b = 0.25472_rp
        real(rp), parameter :: half = 0.5_rp
        real(rp), parameter :: r1 = 0.27597_rp
        real(rp), parameter :: r2 = 0.27846_rp
        real(rp)            :: u, v, x, y, q

#ifndef AMD
        real(rp), parameter       :: mean = 0.0_rp, sdev = 1.0_rp
#else
        real(c_double), parameter :: mean = 0.0, sdev = 1.0
#endif
        integer, dimension(3) :: dim_
        integer(8)            :: err, sz

        call StartProfRange('DFRandomField2D')

#ifdef _OPENACC
#ifdef NVIDIA
        dim_ = shape(DF_Rnd2D)
        sz = dim_(1)*dim_(2)*dim_(3)

        !$acc host_data use_device(DF_Rnd2D)
        err = curandGenerateNormalDouble(generator, DF_Rnd2D, sz, mean, sdev)
        if (err.ne.0) print*,"Error in curandGenerateNormalDouble: ", err
        !$acc end host_data
#endif
#ifdef AMD
        !$acc host_data use_device(rnd2Dnum, rnd2Dnumptr)
        rnd2Dnumptr = c_loc(rnd2Dnum(1))
        
        err = rocrand_generate_normal_double(generator, rnd2Dnumptr, rnd2Dnumsize, mean, sdev)
        if (err.ne.0) print*,"Error in rocrand_generate_normal_double: ", err
        !$acc end host_data

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k=1,rnd2Dnumshape(3)
          do j=1,rnd2Dnumshape(2)
            do i=1,rnd2Dnumshape(1)
               DF_Rnd2D(i,j-DF_N,k-DF_N) = rnd2Dnum(i+(j-1)*rnd2Dnumshape(1)+(k-1)*rnd2Dnumshape(1)*rnd2Dnumshape(2))
            enddo
          enddo
        enddo
        !$acc end parallel
#endif
#else

        ! === compute random field
        do       k = 1,nz
           do    j = 1-DF_N,ny+DF_N
              do m = 1,3

                 do
                   call random_number(u)
                   call random_number(v)
                   v = 1.7156_rp * (v - half)
                 
                   ! Evaluate the quadratic form
                   x = u - s
                   y = ABS(v) - t
                   q = x**2 + y*(a*y - b*x)
                 
                   ! Accept P if inside inner ellipse
                   if (q < r1) exit
                   ! Reject P if outside outer ellipse
                   if (q > r2) cycle
                   ! Reject P if outside acceptance region
                   if (v**2 < -4.0_rp*LOG(u)*u**2) exit
                 enddo

                 DF_Rnd2D(m,j,k) = v/u

              enddo
           enddo
        enddo

#endif
                
        irNz = 1.0_rp/real(nz,rp)

        ! === remove the mean computed spanwise
        !$acc parallel default(present)
        !$acc loop gang, vector
        do j = 1-DF_N,ny+DF_N

           mean_x = 0.0_rp
           mean_y = 0.0_rp
           mean_z = 0.0_rp
           !
           rmsq_x = 0.0_rp
           rmsq_y = 0.0_rp
           rmsq_z = 0.0_rp

           !$acc loop seq
           do k = 1,nz
              mean_x = mean_x + DF_Rnd2D(1,j,k)
              mean_y = mean_y + DF_Rnd2D(2,j,k)
              mean_z = mean_z + DF_Rnd2D(3,j,k)
              !
              rmsq_x = rmsq_x + DF_Rnd2D(1,j,k)**2
              rmsq_y = rmsq_y + DF_Rnd2D(2,j,k)**2
              rmsq_z = rmsq_z + DF_Rnd2D(3,j,k)**2
           enddo
           mean_x = mean_x*irNz
           mean_y = mean_y*irNz
           mean_z = mean_z*irNz
           !
           rmsq_x = rmsq_x*irNz
           rmsq_y = rmsq_y*irNz
           rmsq_z = rmsq_z*irNz
           !
           rms_x  = sqrt(rmsq_x - mean_x**2)
           rms_y  = sqrt(rmsq_y - mean_y**2)
           rms_z  = sqrt(rmsq_z - mean_z**2)
        
           !$acc loop seq
           do k = 1,nz
              DF_Rnd2D(1,j,k) = (DF_Rnd2D(1,j,k) - mean_x)/rms_x
              DF_Rnd2D(2,j,k) = (DF_Rnd2D(2,j,k) - mean_y)/rms_y
              DF_Rnd2D(3,j,k) = (DF_Rnd2D(3,j,k) - mean_z)/rms_z
           enddo
        enddo
        !$acc end parallel

        ! === apply periodicity spanwise

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








subroutine DFRandomField3D(nx,ny,nz,gn)

        use random_module, only: random_normal
#ifdef DEBUG
        use FileModule
        use mpi_module       , only: rank
        use parameters_module, only: data_dir
#endif

        implicit none
        integer     , intent(in)    :: nx, ny, nz,gn
        
        ! local 
        real(rp), dimension(3) :: mean, rmsq, rms
        integer                :: m, i,j,k
        real(rp)               :: r
#ifdef DEBUG
        type(FileType) :: randomField3D
#endif
        !
        ! === compute random field
        !
        do k = 1,nz
           do j = 1-DF_N,ny+DF_N
              do i = 1-gn,nx+gn+1
                 do m = 1,3

                    call random_normal(r)
                    DF_Rnd3D(m,i,j,k) = r

                 enddo
              enddo
           enddo
        enddo
        !
        ! === remove the mean computed spanwise
        !
        do i = 1-gn,nx+gn+1
           do j = 1-DF_N,ny+DF_N

              mean = 0.0_rp
              rmsq = 0.0_rp
              do k = 1,nz
                 do m = 1,3
                    mean(m) = mean(m) + DF_Rnd3D(m,i,j,k)
                    rmsq(m) = rmsq(m) + DF_Rnd3D(m,i,j,k)**2
                 enddo
              enddo
              mean = mean/real(nz,rp)
              rmsq = rmsq/real(nz,rp)
              rms  = sqrt(rmsq - mean**2)

              do k = 1,nz
                 do m = 1,3
                    DF_Rnd3D(m,i,j,k) = (DF_Rnd3D(m,i,j,k) - mean(m))/rms(m)
                 enddo
              enddo
           enddo
        enddo
        !
        ! === apply periodicity spanwise
        !
        do k = 1,DF_N
           do j = 1-DF_N,ny+DF_N
              do i = 1-GN,nx+gn+1
                 do m = 1,3
                    DF_Rnd3D(m,i,j,nz+k) = DF_Rnd3D(m,i,j,k)
                    DF_Rnd3D(m,i,j,1 -k) = DF_Rnd3D(m,i,j,nz+1-k)
                 enddo
              enddo
           enddo
        enddo

#ifdef DEBUG
        randomField3D%name = 'randomField_'//trim(str(rank))
        randomField3D%dir  = trim(data_dir)//'/INIT_RANDOM_FIELD'
        call OpenNewFile(randomField3D,0)
        do k =1-DF_N,nz+DF_N
           write(randomField3D%unit,*) k, DF_Rnd3D(1,1,1,k)
        enddo
        call CloseFile(randomField3D)
#endif

        return
end subroutine DFRandomField3D











subroutine DFIntegralLenght_TBL(y_gbl,ny)
! ------------------------------------------------------------------------------------------
!       This subroutine compute the integral lenght scales for the digital filtering process
!
!       INPUT : y                !< mpi local y grid coordinate
!               sy, ey           !< local grid start/end indices
!               gbl_min_step     !< minimum grid step of the whole domain
!
! ----------------------------------------------------------------------------
        use parameters_module, only: data_dir, ReTau
        use mpi_module       , only: rank
#ifdef DEBUG
        use FileModule
#endif

        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: y_gbl
        integer                            , intent(in)    :: ny

        ! local declarations
        real(rp) :: zlenou_u, zlenou_v, zlenou_w
        real(rp) :: zlenin_u, zlenin_v, zlenin_w
        real(rp) :: ftany
        integer  :: j

#ifdef DEBUG
        type(FileType) :: IntLenFile
#endif

        !Outer z scale for u,v,w
        zlenou_u = 0.40_rp
        zlenou_v = 0.30_rp
        zlenou_w = 0.40_rp
        !Inner z scale for u,v,w
        zlenin_u = min(150.0_rp/ReTau,zlenou_u)
        zlenin_v = min( 75.0_rp/ReTau,zlenou_v)
        zlenin_w = min(150.0_rp/ReTau,zlenou_w)

        do j=0,ny
           ftany = 0.5_rp*(1.0_rp+tanh((y_gbl(j)-0.2_rp)/0.03_rp)) ! blending function
           DF_zlen(1,j)  = zlenin_u+ftany*(zlenou_u-zlenin_u)
           DF_zlen(2,j)  = zlenin_v+ftany*(zlenou_v-zlenin_v)
           DF_zlen(3,j)  = zlenin_w+ftany*(zlenou_w-zlenin_w)

           !DF_zlen(1,j) = min(DF_dimMin(1) + (y_gbl(j)/0.3_rp)*0.20_rp,DF_dimMax(1))
           !DF_zlen(2,j) = min(DF_dimMin(2) + (y_gbl(j)/0.3_rp)*0.15_rp,DF_dimMax(2))
           !DF_zlen(3,j) = min(DF_dimMin(3) + (y_gbl(j)/0.3_rp)*0.20_rp,DF_dimMax(3))
        enddo
        DF_ylen = 0.7_rp*DF_zlen

        DF_xlen1 = 0.80_rp
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




subroutine DFIntegralLenght_HTURB(Lz)
        implicit none
        real(rp)    , intent(in)    :: Lz

        DF_ylen = 0.25_rp*Lz
        DF_zlen = 0.25_rp*Lz
        
        return
end subroutine DFIntegralLenght_HTURB



subroutine DFCoefficients(y_gbl,ny,gbl_min_step)
! -----------------------------------------------------------------------
!
!       Computation of the Digital filters exponetial coefficients
!
! ---------------------------------------------------------------
        use parameters_module, only: rp, pi, data_dir

#ifdef DEBUG
        use FileModule 
        use mpi_module, only: rank
#endif

        implicit none
        real(rp)    , dimension(:), allocatable, intent(in)    :: y_gbl
        real(rp)    , dimension(:)             , intent(in)    :: gbl_min_step
        integer                                , intent(in)    :: ny

        ! local declarations
        real(rp), dimension(3) :: sumby, sumbz
        real(rp)               :: dy, dz
        integer  :: m, j,jj
#ifdef DEBUG
        type(FileType) :: DFCoeffFIle
#endif
        
        
        do j = 1,ny
                
           dy = y_gbl(j) - y_gbl(j-1)
           dz = gbl_min_step(3)

           ! === y filter coefficients
           sumby = 0.0_rp
           do jj = -DF_N,DF_N
              DF_By(:,jj,j) = exp(-pi*abs(jj)/(DF_yLen(:,j)/dy))

              sumby(:) = sumby(:) + DF_By(:,jj,j)**2
           enddo
           do m = 1,3
              DF_By(m,:,j) = DF_By(m,:,j)/sqrt(sumby(m))
           enddo

           ! === z filter coefficients
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

        ! y convolution
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

        ! z convolution

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
        

subroutine DFConvolution3D(sx,ex,sy,ey,sz,ez,vf)

#ifdef DEBUG
        use FileModule 
        use parameters_module, only: data_dir
        use mpi_module       , only: rank
#endif
        implicit none
        integer                                  , intent(in)    :: sx, ex, sy, ey, sz, ez
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: vf

        ! local declarations
        real(rp), allocatable, dimension(:,:,:,:) :: fy
        integer :: m,j,k,jj,kk,err = 0
#ifdef DEBUG
        type(FileType) :: velFluc
#endif

        allocate(fy(3, sx:ex, sy-DF_N:ey+DF_N, sz-DF_N:ez+DF_N), stat = err)
        if(err .ne. 0) stop ' Allocation error in DFConvolution2D'

        ! y convolution
        fy = 0.0_rp
        do k = sz-DF_N,ez+DF_N
           do j = sy,ey
              do jj = -DF_N,DF_N
                 do m = 1,3
                    fy(m,sx:ex,j,k) = fy(m,sx:ex,j,k) + DF_By(m,jj,j)*DF_RnD3D(m,sx:ex,j+jj,k)
                 enddo
              enddo
           enddo
        enddo

        ! z convolution
        vf = 0.0_rp
        do k = sz,ez
           do j = sy,ey
              do kk = -DF_N,DF_N
                 do m = 1,3
                    vf(m,sx:ex,j,k) = vf(m,sx:ex,j,k) + DF_Bz(m,kk,j)*fy(m,sx:ex,j,k+kk)
                 enddo
              enddo
           enddo
        enddo

        deallocate(fy)

#ifdef DEBUG
        velfluc%name = 'velocity_fluctuations'//trim(str(rank))
        velfluc%dir  = trim(data_dir)//'/DF_VELOCITY_FLUCTUATIONS'
        call OpenNewFile(velfluc,0)
        do j = sy,ey
           write(velfluc%unit,*) j, vf(:,1,j,1)
        enddo
        call CloseFile(velFluc)
#endif

        return
end subroutine DFConvolution3D


subroutine DFCastroTimeCorrelation(ik,c_rk,u_inf,dt,sy,ey,sz,ez,vf_old,vf_new)

        use parameters_module, only: rp, pi

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: vf_new, vf_old
        real(rp), dimension(3)                 , intent(in)    :: c_rk
        real(rp)                               , intent(in)    :: u_inf, dt
        integer                                , intent(in)    :: ik
        integer                                , intent(in)    :: sy,ey, sz, ez

        ! local declarations
        real(rp) :: tlen_x, tlen_y, tlen_z
        real(rp) :: exparg1_x, exparg1_y, exparg1_z
        real(rp) :: exparg2_x, exparg2_y, exparg2_z
        real(rp) :: exp1_x, exp1_y, exp1_z
        real(rp) :: sqrtexp_x, sqrtexp_y, sqrtexp_z
        real(rp) :: dtOld
        integer  :: j,k,ikk

        call StartProfRange('turbulent_inflow')

        ikk = mod(ik,3)-1
        ikk = mod(ikk+3,3)
        if(ikk == 0) ikk = 3
        dtold = c_rk(ikk)*dt

        dtold = 1.0_rp/3.0_rp*dt

        tlen_x    = DF_xlen1/u_inf
        tlen_y    = DF_xlen2/u_inf
        tlen_z    = DF_xlen3/u_inf

        exparg1_x = -0.5_rp*pi*dtOld/tlen_x
        exparg1_y = -0.5_rp*pi*dtOld/tlen_y
        exparg1_z = -0.5_rp*pi*dtOld/tlen_z

        exparg2_x = -pi*dtOld/tlen_x
        exparg2_y = -pi*dtOld/tlen_y
        exparg2_z = -pi*dtOld/tlen_z

        exp1_x    = exp(exparg1_x) 
        exp1_y    = exp(exparg1_y) 
        exp1_z    = exp(exparg1_z) 

        sqrtexp_x = sqrt(1._rp-exp(exparg2_x))
        sqrtexp_y = sqrt(1._rp-exp(exparg2_y))
        sqrtexp_z = sqrt(1._rp-exp(exparg2_z))
        
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do k    = sz,ez
           do j = sy,ey
              vf_new(1,j,k) = vf_old(1,j,k)*exp1_x + vf_new(1,j,k)*sqrtexp_x
              vf_new(2,j,k) = vf_old(2,j,k)*exp1_y + vf_new(2,j,k)*sqrtexp_y
              vf_new(3,j,k) = vf_old(3,j,k)*exp1_z + vf_new(3,j,k)*sqrtexp_z

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

        ! local declarations
        integer :: j, k

        call StartProfRange('DFEnforceReynoldsStresses2D')

        ! perform the multiplication with Lund's Matrix

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






subroutine DFEnforceReynoldsStresses3D(sx,ex,sy,ey,sz,ez,vf,uf)

        implicit none
        integer                                  , intent(in)    :: sx, ex, sy, ey, sz, ez
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: vf
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: uf

        ! local declarations
        integer :: i,j, k, m, l 

        ! perform the multiplication with Lund's Matrix
        uf = 0.0_rp
        do k = sz,ez
           do j = sy,ey
              do i = sx,ex
                 do m = 1,3
                    do l = 1,3
                       uf(m,i,j,k) = uf(m,i,j,k) + DF_LundMatrix(m,l,j)*vf(l,i,j,k)
                    enddo
                 enddo
              enddo
           enddo
        enddo

        return
end subroutine DFEnforceReynoldsStresses3D

























end module df_module
