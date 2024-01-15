module post_stat_storage_module
use parameters_module
use FileModule
use mpi_module
use storage_module


type(DirType) :: MainDir
type(DirType) :: WmlesDir

integer       :: f         !< file index
integer       :: file_ini  !< index first file  
integer       :: file_end  !< index last  file
character(dl) :: ifile
character(dl) :: ifile_wmles
character(3)  :: pcharAll

real(rp), allocatable, dimension(:,:) :: FavreVelU2D
real(rp), allocatable, dimension(:,:) :: FavreVelV2D
real(rp), allocatable, dimension(:,:) :: FavreVelW2D

real(rp), allocatable, dimension(:,:) :: FavreMach2D

real(rp), allocatable, dimension(:,:) :: MeanMuTot2D

real(rp), allocatable, dimension(:,:) :: FavreRUU_2D
real(rp), allocatable, dimension(:,:) :: FavreRVV_2D
real(rp), allocatable, dimension(:,:) :: FavreRWW_2D
real(rp), allocatable, dimension(:,:) :: FavreRUV_2D

real(rp), allocatable, dimension(:,:) :: FavreUrms2D
real(rp), allocatable, dimension(:,:) :: FavreVrms2D
real(rp), allocatable, dimension(:,:) :: FavreWrms2D

real(rp), allocatable, dimension(:,:) :: FavreEkT_2D

real(rp), allocatable, dimension(:,:) :: FavreOMX_2D
real(rp), allocatable, dimension(:,:) :: FavreOMY_2D
real(rp), allocatable, dimension(:,:) :: FavreOMZ_2D
real(rp), allocatable, dimension(:,:) :: Favre_absOM

real(rp), allocatable, dimension(:,:) :: Favre_Ux_2D
real(rp), allocatable, dimension(:,:) :: Favre_Vx_2D
real(rp), allocatable, dimension(:,:) :: Favre_Uy_2D
real(rp), allocatable, dimension(:,:) :: Favre_Vy_2D

real(rp), allocatable, dimension(:,:) :: RMS_Sens_2D
real(rp), allocatable, dimension(:,:) :: RMS_Wflg_2D

! BOUNDARY LAYER stats
real(rp), allocatable, dimension(:,:) :: StreamWiseStats


contains
subroutine GetPeriodicDirs
        implicit none
        character(1), dimension(3) :: pchar

        pchar = 'F'
        if(bc(1) == 'periodic') pchar(1) = 'T'
        if(bc(3) == 'periodic') pchar(2) = 'T'
        if(bc(5) == 'periodic') pchar(3) = 'T'
        pcharAll = pchar(1)//pchar(2)//pchar(3)

        selectcase(pcharAll)
          case('FFT')
          write(*,'(A)') ' The problem is z-periodic'

          case('TFT')
          write(*,'(A)') ' The problem is xz-periodic'

          case('FTT')
          write(*,'(A)') ' The problem is yz-periodic'

        endselect

        return
end subroutine GetPeriodicDirs

subroutine InitStatsFields
        implicit none
        integer :: nsw = 21
        integer :: err = 0

        ! init indices
        sx = 1; ex = nx
        sy = 1; ey = ny
        sz = 1; ez = nz

        lbx = sx - GN ; ubx = ex + GN
        lby = sy - GN ; uby = ey + GN
        lbz = sz - GN ; ubz = ez + GN
        
        ! init statistical field
        selectcase(pcharAll)
          case('FFT','FTT') 
          call AllocateReal(vmean2D,lbx,ubx,lby,uby,1,nvAve2D)

          call AllocateReal(FavreVelU2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreVelV2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreVelW2D,lbx,ubx,lby,uby)

          call AllocateReal(FavreMach2D,lbx,ubx,lby,uby)

          call AllocateReal(MeanMuTot2D,lbx,ubx,lby,uby)

          call AllocateReal(FavreRUU_2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreRVV_2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreRWW_2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreRUV_2D,lbx,ubx,lby,uby)


          call AllocateReal(FavreUrms2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreVrms2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreWrms2D,lbx,ubx,lby,uby)

          call AllocateReal(FavreEkT_2D,lbx,ubx,lby,uby)

          call AllocateReal(FavreOMX_2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreOMY_2D,lbx,ubx,lby,uby)
          call AllocateReal(FavreOMZ_2D,lbx,ubx,lby,uby)
          call AllocateReal(Favre_absOM,lbx,ubx,lby,uby)

          call AllocateReal(Favre_Ux_2D,lbx,ubx,lby,uby)
          call AllocateReal(Favre_Uy_2D,lbx,ubx,lby,uby)
          call AllocateReal(Favre_Vx_2D,lbx,ubx,lby,uby)
          call AllocateReal(Favre_Vy_2D,lbx,ubx,lby,uby)

          call AllocateReal(RMS_Sens_2D,lbx,ubx,lby,uby)
          call AllocateReal(RMS_Wflg_2D,lbx,ubx,lby,uby)

          case('TFT')
          call AllocateReal(vmean1D,lby,uby,1,nvAve1D)

        endselect

        selectcase(ic)
        case('turbulent_BL','swbli')
          if(wmles) then
            nsw = nsw + nvWmlesData
            call AllocateReal(vmean1D_wmles,lbx,ubx,1,nvWmlesData)
          endif
          call AllocateReal(StreamWiseStats,lbx,ubx,1,nsw)
        endselect

        return
end subroutine InitStatsFields


subroutine DestroyStatsFields
        implicit none

        call DeallocateReal(vmean1D)
        call DeallocateReal(vmean2D)

        return
end subroutine DestroyStatsFields





end module post_stat_storage_module
