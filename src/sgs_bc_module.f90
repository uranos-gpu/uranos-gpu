module sgs_bc_module
use parameters_module, only: rp
use mpi_module
use profiling_module

implicit none
private
real(rp), parameter :: Pr_turb = 0.9_rp

public TrbVisBC, &
       ZeroGradientTrbVis, &
       WallResolvedTrbVis

contains
subroutine TrbVisBC(VIS,LMD)

        use bc_module, only: face_type

        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: VIS, LMD

        type(face_type), dimension(6) :: all_bound
        type(face_type)               :: bound
        real(rp)                      :: tWall
        integer                       :: f

        call StartProfRange("TrbVisBC")

        if(.not.mpi_flag) my_neighbour = MPI_PROC_NULL

        all_bound%node = (/ sx , ex , sy , ey , sz ,  ez/)
        all_bound%norm = (/ -1 ,  1 , -1 ,  1 , -1 ,  1 /)
        all_bound%face = (/ 'W', 'E', 'S', 'N', 'B', 'F'/)

        do f = 1, 6
           if(my_neighbour(f) == MPI_PROC_NULL) then

             bound%node = all_bound(f)%node
             bound%norm = all_bound(f)%norm
             bound%face = all_bound(f)%face

             selectcase(trim(bc(f)))
               case('adiabatic_wall', 'diabatic_wall')
                 tWall = Trat*T_recovery
                 call WallResolvedTrbVis(bound,tWall,VIS,LMD)

               case('isothermal_wall')
                 tWall = 1.0_rp
                 call WallResolvedTrbVis(bound,tWall,VIS,LMD)

               case default
                 call ZeroGradientTrbVis(bound,VIS,LMD)

             endselect

           endif
        enddo

        call EndProfRange

        return
end subroutine TrbVisBC

subroutine ZeroGradientTrbVis(b,VIS,LMD)

        use bc_module, only: face_type

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS,LMD
        type(face_type)                        , intent(in)    :: b

        integer :: i,j,k,bnode,bnorm

        call StartProfRange("ZeroGradientTrbVis")

        bnode = b%node
        bnorm = b%norm

        selectcase(b%face)

          case('W','E')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(3)
            do          k = sz,ez
               do       j = sy,ey
                  do    i = 1,GN
                     VIS(bnode+bnorm*i,j,k) = VIS(bnode-bnorm*(i-1),j,k)
                     LMD(bnode+bnorm*i,j,k) = LMD(bnode-bnorm*(i-1),j,k)
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('S','N')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(3)
            do          k = sz,ez
               do       j = 1,GN
                  do    i = sx,ex
                     VIS(i,bnode+bnorm*j,k) = VIS(i,bnode-bnorm*(j-1),k)
                     LMD(i,bnode+bnorm*j,k) = LMD(i,bnode-bnorm*(j-1),k)
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('B','F')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(3)
            do          k = 1,GN
               do       j = sy,ey
                  do    i = sx,ex
                     VIS(i,j,bnode+bnorm*k) = VIS(i,j,bnode-bnorm*(k-1))
                     LMD(i,j,bnode+bnorm*k) = LMD(i,j,bnode-bnorm*(k-1))
                  enddo
               enddo
            enddo
            !$acc end parallel

        endselect

        call EndProfRange

        return
end subroutine ZeroGradientTrbVis

subroutine WallResolvedTrbVis(b,tWall,VIS,LMD)

        use fluid_functions_module, only: laminar_viscosity
        use bc_module             , only: face_type
        use parameters_module     , only: k_inf, iPrandtl, thermally_perfect_gas
        use thermodynamics_module , only: cp_local

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS,LMD
        real(rp)                               , intent(in)    :: tWall
        type(face_type)                        , intent(in)    :: b

        real(rp) :: mWall, kWall, twoM, twoK
        integer  :: i,j,k,bnode,bnorm

        call StartProfRange("WallResolvedTrbVis")

        mWall = laminar_viscosity(tWall,Tref,vis_flag)
        if(thermally_perfect_gas > 0) then
           kWall = cp_local(tWall) * iPrandtl * mWall
        else
           kWall = k_inf * mWall
        endif
        twoM  = 2.0_rp*mWall
        twoK  = 2.0_rp*kWall

        bnode = b%node
        bnorm = b%norm

        selectcase(b%face)

          case('W','E')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(3)
            do          k = sz,ez
               do       j = sy,ey
                  do    i = 1,GN
                     VIS(bnode+bnorm*i,j,k) = twoM - VIS(bnode-bnorm*(i-1),j,k)
                     LMD(bnode+bnorm*i,j,k) = twoK - LMD(bnode-bnorm*(i-1),j,k)
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('S','N')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(3)
            do          k = sz,ez
               do       j = 1,GN
                  do    i = sx,ex
                     VIS(i,bnode+bnorm*j,k) = twoM - VIS(i,bnode-bnorm*(j-1),k)
                     LMD(i,bnode+bnorm*j,k) = twoK - LMD(i,bnode-bnorm*(j-1),k)
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('B','F')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(3)
            do          k = 1,GN
               do       j = sy,ey
                  do    i = sx,ex
                     VIS(i,j,bnode+bnorm*k) = twoM - VIS(i,j,bnode-bnorm*(k-1))
                     LMD(i,j,bnode+bnorm*k) = twoK - LMD(i,j,bnode-bnorm*(k-1))
                  enddo
               enddo
            enddo
            !$acc end parallel

          case default
            print*, ' WallResolvedTrbVis is not implemented for face ', &
                      trim(b%face)
            stop

        endselect

        call EndProfRange

        return
end subroutine WallResolvedTrbVis

end module sgs_bc_module
