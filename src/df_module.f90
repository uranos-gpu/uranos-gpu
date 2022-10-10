module df_module
use parameters_module, only: rp
use mpi_module, only: sy,ey,sz,ez
use random_module
#ifdef _OPENACC
use openacc_curand
#endif

implicit none

type dfData
  integer                :: N      = 60 !30                !< filter width
  real(rp), dimension(3) :: dimMin = (/0.05_rp, 0.05_rp, 0.05_rp/)  !< filter min Length
  real(rp), dimension(3) :: dimMax = (/0.25_rp, 0.20_rp, 0.25_rp/)  !< filter max Length
        
  real(rp), dimension(3)                  :: xlen = (/0.5_rp, 0.15_rp, 0.15_rp/)
  real(rp), dimension(:,:)  , allocatable :: ylen
  real(rp), dimension(:,:)  , allocatable :: zlen
  real(rp), dimension(:,:,:), allocatable :: By
  real(rp), dimension(:,:,:), allocatable :: Bz

  real(rp), dimension(:,:,:)  , allocatable :: RnD2D
  real(rp), dimension(:,:,:,:), allocatable :: RnD3D
  real(rp), dimension(:,:,:)  , allocatable :: LundMatrix
  real(rp), dimension(:,:,:)  , allocatable :: fy
endtype DfData
  
type(DfData), public :: DF
contains

subroutine DFInputData_TBL(ReTau,Mach,Prandtl,y_gbl,ny,DF)
        
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
        type(DfData)                       , intent(inout) :: DF
        
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
        allocate(DF%LundMatrix(3,3,0:ny),stat=err)
        if(err .ne. 0) stop ' Allocation error of Lunds Matrix'

        DF%LundMatrix = 0.0_rp
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
        
          DF%LundMatrix(:,:,j) = a(:,:)

        enddo

        allocate(DF%fy(3, sy-DF%N:ey+DF%N, sz-DF%N:ez+DF%N), stat = err)
        if(err .ne. 0) stop ' Allocation error of DF%fy'

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
           write(lundFile%unit,'(7e18.6)') y_gbl(j), DF%LundMatrix(1,1,j), &
                                        DF%LundMatrix(2,2,j), &
                                        DF%LundMatrix(3,3,j), &
                                        DF%LundMatrix(2,1,j), &
                                        DF%LundMatrix(3,1,j), &
                                        DF%LundMatrix(3,2,j)
        enddo
        call CloseFile(LundFile)
#endif

        deallocate(dtFile, dtReth, dtData, dtDataInt, ReyStress)

        return
end subroutine DFInputData_TBL


subroutine DFInputData_HTURB(y_gbl,ny,DF)

#ifdef DEBUG
        use parameters_module, only: data_dir
        use mpi_module       , only: rank
        use FileModule
#endif

        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: y_gbl
        integer                            , intent(in)    :: ny
        type(DfData)                       , intent(inout) :: DF

        ! local declarations
        real(rp), dimension(:,:,:), allocatable :: ReyStress
        real(rp), dimension(3,3)                :: a

        real(rp), parameter :: trb_int = 1.0E-4_rp
        integer             :: j, err = 0
#ifdef DEBUG
        type(FileType) :: DF_inputs, LundFile
#endif

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
        allocate(DF%LundMatrix(3,3,0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in DFInputData_HTURB'

        DF%LundMatrix = 0.0_rp
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

           DF%LundMAtrix(:,:,j) = a(:,:)

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
           write(lundFile%unit,'(7e18.6)') y_gbl(j), DF%LundMatrix(1,1,j), &
                                        DF%LundMatrix(2,2,j), &
                                        DF%LundMatrix(3,3,j), &
                                        DF%LundMatrix(2,1,j), &
                                        DF%LundMatrix(3,1,j), &
                                        DF%LundMatrix(3,2,j)
        enddo
        call CloseFile(LundFile)
#endif
        
        deallocate(ReyStress)

        return
end subroutine DFInputData_HTURB


subroutine DFRandomField2D(ny,nz,DF)

        use random_module, only: random_normal

        implicit none
        type(DfData), intent(inout) :: DF
        integer     , intent(in)    :: ny,nz
        
        ! local 
        real(rp)  :: mean_x, mean_y, mean_z
        real(rp)  :: rmsq_x, rmsq_y, rmsq_z
        real(rp)  :: rms_x, rms_y, rms_z
        integer   :: m,j,k, l
        real(rp)  :: irNz

        real(rp), parameter :: s  = 0.449871_rp
        real(rp), parameter :: t  = -0.386595_rp
        real(rp), parameter :: a  = 0.19600_rp
        real(rp), parameter :: b = 0.25472_rp
        real(rp), parameter :: half = 0.5_rp
        real(rp), parameter :: r1 = 0.27597_rp
        real(rp), parameter :: r2 = 0.27846_rp
        real(rp)            :: u, v, x, y, q


#ifdef _OPENACC

        !$acc parallel num_gangs(1) vector_length(1) private(cudaState) &
        !$acc default(present)
        !$acc loop collapse(3) 
        do       k = 1,nz
           do    j = 1-DF%N,ny+DF%N
              do m = 1,3
                      
                 DF%Rnd2D(m,j,k) = curand_normal(cudaState)

              enddo
           enddo
        enddo
        !$acc end parallel
        
#else


        ! === compute random field
        do       k = 1,nz
           do    j = 1-DF%N,ny+DF%N
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

                 DF%Rnd2D(m,j,k) = v/u

              enddo
           enddo
        enddo

#endif
                
        irNz = 1.0_rp/real(nz,rp)

        ! === remove the mean computed spanwise
        !$acc parallel default(present)
        !$acc loop gang, vector
        do j = 1-DF%N,ny+DF%N

           mean_x = 0.0_rp
           mean_y = 0.0_rp
           mean_z = 0.0_rp
           !
           rmsq_x = 0.0_rp
           rmsq_y = 0.0_rp
           rmsq_z = 0.0_rp

           !$acc loop seq
           do k = 1,nz
              mean_x = mean_x + DF%Rnd2D(1,j,k)
              mean_y = mean_y + DF%Rnd2D(2,j,k)
              mean_z = mean_z + DF%Rnd2D(3,j,k)
              !
              rmsq_x = rmsq_x + DF%Rnd2D(1,j,k)**2
              rmsq_y = rmsq_y + DF%Rnd2D(2,j,k)**2
              rmsq_z = rmsq_z + DF%Rnd2D(3,j,k)**2
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
              DF%Rnd2D(1,j,k) = (DF%Rnd2D(1,j,k) - mean_x)/rms_x
              DF%Rnd2D(2,j,k) = (DF%Rnd2D(2,j,k) - mean_y)/rms_y
              DF%Rnd2D(3,j,k) = (DF%Rnd2D(3,j,k) - mean_z)/rms_z
           enddo
        enddo
        !$acc end parallel

        ! === apply periodicity spanwise

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = 1,DF%N
           do    j = 1-DF%N,ny+DF%N
              do l = 1,3
                 DF%Rnd2D(l,j,nz+k) = DF%Rnd2D(l,j,k)
                 DF%Rnd2D(l,j,1 -k) = DF%Rnd2D(l,j,nz+1-k)
              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine DFRandomField2D








subroutine DFRandomField3D(nx,ny,nz,gn,DF)

        use random_module, only: random_normal
#ifdef DEBUG
        use FileModule
        use mpi_module       , only: rank
        use parameters_module, only: data_dir
#endif

        implicit none
        type(DfData), intent(inout) :: DF
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
           do j = 1-DF%N,ny+DF%N
              do i = 1-gn,nx+gn+1
                 do m = 1,3

                    call random_normal(r)
                    DF%Rnd3D(m,i,j,k) = r

                 enddo
              enddo
           enddo
        enddo
        !
        ! === remove the mean computed spanwise
        !
        do i = 1-gn,nx+gn+1
           do j = 1-DF%N,ny+DF%N

              mean = 0.0_rp
              rmsq = 0.0_rp
              do k = 1,nz
                 do m = 1,3
                    mean(m) = mean(m) + DF%Rnd3D(m,i,j,k)
                    rmsq(m) = rmsq(m) + DF%Rnd3D(m,i,j,k)**2
                 enddo
              enddo
              mean = mean/real(nz,rp)
              rmsq = rmsq/real(nz,rp)
              rms  = sqrt(rmsq - mean**2)

              do k = 1,nz
                 do m = 1,3
                    DF%Rnd3D(m,i,j,k) = (DF%Rnd3D(m,i,j,k) - mean(m))/rms(m)
                 enddo
              enddo
           enddo
        enddo
        !
        ! === apply periodicity spanwise
        !
        do k = 1,DF%N
           do j = 1-DF%N,ny+DF%N
              do i = 1-GN,nx+gn+1
                 do m = 1,3
                    DF%Rnd3D(m,i,j,nz+k) = DF%Rnd3D(m,i,j,k)
                    DF%Rnd3D(m,i,j,1 -k) = DF%Rnd3D(m,i,j,nz+1-k)
                 enddo
              enddo
           enddo
        enddo

#ifdef DEBUG
        randomField3D%name = 'randomField_'//trim(str(rank))
        randomField3D%dir  = trim(data_dir)//'/INIT_RANDOM_FIELD'
        call OpenNewFile(randomField3D,0)
        do k =1-DF%N,nz+DF%N
           write(randomField3D%unit,*) k, DF%Rnd3D(1,1,1,k)
        enddo
        call CloseFile(randomField3D)
#endif

        return
end subroutine DFRandomField3D











subroutine DFIntegralLenght_TBL(y_gbl,ny,DF)
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
        type(dfData)                       , intent(inout) :: DF
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
           DF%zlen(1,j)  = zlenin_u+ftany*(zlenou_u-zlenin_u)
           DF%zlen(2,j)  = zlenin_v+ftany*(zlenou_v-zlenin_v)
           DF%zlen(3,j)  = zlenin_w+ftany*(zlenou_w-zlenin_w)

           !DF%zlen(1,j) = min(DF%dimMin(1) + (y_gbl(j)/0.3_rp)*0.20_rp,DF%dimMax(1))
           !DF%zlen(2,j) = min(DF%dimMin(2) + (y_gbl(j)/0.3_rp)*0.15_rp,DF%dimMax(2))
           !DF%zlen(3,j) = min(DF%dimMin(3) + (y_gbl(j)/0.3_rp)*0.20_rp,DF%dimMax(3))
        enddo
        DF%ylen = 0.7_rp*DF%zlen

        DF%xlen(1) = 0.80_rp
        DF%xlen(2) = 0.30_rp
        DF%xlen(3) = 0.30_rp

#ifdef DEBUG
        intLenFile%name = 'integralLen'//trim(str(rank))
        intLenFile%dir  = trim(data_dir)//'/DF_INTEGRAL_LENGTH'
        call OpenNewFile(intLenFile,0)
        do j = 0,ny
           write(intLenFile%unit,*) y_gbl(j), Df%ylen(:,j), Df%zlen(:,j)
        enddo
        call CloseFile(intLenFile)
#endif
        return
end subroutine DFIntegralLenght_TBL




subroutine DFIntegralLenght_HTURB(Lz,DF)
        implicit none
        type(dfData), intent(inout) :: DF
        real(rp)    , intent(in)    :: Lz

        DF%ylen = 0.25_rp*Lz
        DF%zlen = 0.25_rp*Lz
        
        return
end subroutine DFIntegralLenght_HTURB



subroutine DFCoefficients(y_gbl,ny,gbl_min_step,DF)
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
        type(dfData)                           , intent(inout) :: DF

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
           do jj = -DF%N,DF%N
              DF%By(:,jj,j) = exp(-pi*abs(jj)/(DF%yLen(:,j)/dy))

              sumby(:) = sumby(:) + DF%By(:,jj,j)**2
           enddo
           do m = 1,3
              DF%By(m,:,j) = DF%By(m,:,j)/sqrt(sumby(m))
           enddo

           ! === z filter coefficients
           sumbz = 0.0_rp
           do jj = -DF%N,DF%N
              DF%Bz(:,jj,j) = exp(-pi*abs(jj)/(DF%zLen(:,j)/dz))

              sumbz(:) = sumbz(:) + DF%Bz(:,jj,j)**2
           enddo
           do m = 1,3
              DF%Bz(m,:,j) = DF%Bz(m,:,j)/sqrt(sumbz(m))
           enddo

        enddo

#ifdef DEBUG
        DFCoeffFile%name = 'Bcoefficient_'//trim(str(rank))
        DFCoeffFIle%dir  = trim(data_dir)//'/DF_COEFFICIENT'
        call OpenNewFile(DFCoeffFile,0)
        do j = 0, ny
           do jj = -DF%N,DF%N
              write(DFCoeffFile%unit,*) y_gbl(j), jj, DF%By(:,jj,j), DF%Bz(:,jj,j)
           enddo
        enddo
        call CloseFile(DFCoeffFile)
#endif
        
        return
end subroutine DFcoefficients


subroutine DFConvolution2D(sy,ey,sz,ez,DF,vf)
        implicit none
        
        integer                                , intent(in)    :: sy, ey, sz, ez
        type(DfData)                           , intent(inout) :: DF
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: vf

        integer :: l, j,k,jj,kk

        ! y convolution
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz-DF%N,ez+DF%N
           do    j = sy-DF%N,ey+DF%N
              do l = 1,3
                 DF%fy(l,j,k) = 0.0_rp
              enddo
           enddo
        enddo
        !$acc end parallel
        
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          k = sz-DF%N,ez+DF%N
           do       j = sy,ey
              do    jj = -DF%N,DF%N
                 do l = 1,3
                    DF%fy(l,j,k) = DF%fy(l,j,k) + DF%By(l,jj,j)*DF%RnD2D(l,j+jj,k)
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
        !$acc loop gang, vector collapse(4)
        do          k = sz,ez
           do       j = sy,ey
              do    kk = -DF%N,DF%N
                 do l = 1,3
                    vf(l,j,k) = vf(l,j,k) + DF%Bz(l,kk,j)*DF%fy(l,j,k+kk)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine DFConvolution2D
        

subroutine DFConvolution3D(sx,ex,sy,ey,sz,ez,DF,vf)

#ifdef DEBUG
        use FileModule 
        use parameters_module, only: data_dir
        use mpi_module       , only: rank
#endif
        implicit none
        integer                                  , intent(in)    :: sx, ex, sy, ey, sz, ez
        type(DfData)                             , intent(inout) :: DF
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: vf

        ! local declarations
        real(rp), allocatable, dimension(:,:,:,:) :: fy
        integer :: m,j,k,jj,kk,err = 0
#ifdef DEBUG
        type(FileType) :: velFluc
#endif

        allocate(fy(3, sx:ex, sy-DF%N:ey+DF%N, sz-DF%N:ez+DF%N), stat = err)
        if(err .ne. 0) stop ' Allocation error in DFConvolution2D'

        ! y convolution
        fy = 0.0_rp
        do k = sz-DF%N,ez+DF%N
           do j = sy,ey
              do jj = -DF%N,DF%N
                 do m = 1,3
                    fy(m,sx:ex,j,k) = fy(m,sx:ex,j,k) + DF%By(m,jj,j)*DF%RnD3D(m,sx:ex,j+jj,k)
                 enddo
              enddo
           enddo
        enddo

        ! z convolution
        vf = 0.0_rp
        do k = sz,ez
           do j = sy,ey
              do kk = -DF%N,DF%N
                 do m = 1,3
                    vf(m,sx:ex,j,k) = vf(m,sx:ex,j,k) + DF%Bz(m,kk,j)*fy(m,sx:ex,j,k+kk)
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

        ikk = mod(ik,3)-1
        ikk = mod(ikk+3,3)
        if(ikk == 0) ikk = 3
        dtold = c_rk(ikk)*dt

        dtold = 1.0_rp/3.0_rp*dt

        tlen_x    = DF%xlen(1)/u_inf
        tlen_y    = DF%xlen(2)/u_inf
        tlen_z    = DF%xlen(3)/u_inf

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

        return
end subroutine DFCastroTimeCorrelation

        

subroutine DFEnforceReynoldsStresses2D(sy,ey,sz,ez,vf,DF,uf)

        implicit none
        integer                                , intent(in)    :: sy, ey, sz, ez
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: vf
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: uf
        type(DfData)                           , intent(inout) :: DF

        ! local declarations
        integer :: j, k, m, l 

        ! perform the multiplication with Lund's Matrix

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz,ez
           do    j = sy,ey
              do m = 1,3
                 uf(m,j,k) = 0.0_rp
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          k = sz,ez
           do       j = sy,ey
              do    m = 1,3
                 do l = 1,3
                    uf(m,j,k) = uf(m,j,k) + DF%LundMatrix(m,l,j)*vf(l,j,k)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine DFEnforceReynoldsStresses2D






subroutine DFEnforceReynoldsStresses3D(sx,ex,sy,ey,sz,ez,vf,DF,uf)

        implicit none
        integer                                  , intent(in)    :: sx, ex, sy, ey, sz, ez
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: vf
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: uf
        type(DfData)                             , intent(inout) :: DF

        ! local declarations
        integer :: i,j, k, m, l 

        ! perform the multiplication with Lund's Matrix
        uf = 0.0_rp
        do k = sz,ez
           do j = sy,ey
              do i = sx,ex
                 do m = 1,3
                    do l = 1,3
                       uf(m,i,j,k) = uf(m,i,j,k) + DF%LundMatrix(m,l,j)*vf(l,i,j,k)
                    enddo
                 enddo
              enddo
           enddo
        enddo

        return
end subroutine DFEnforceReynoldsStresses3D

























end module df_module
