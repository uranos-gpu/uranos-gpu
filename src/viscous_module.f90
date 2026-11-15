module viscous_module
use profiling_module
use parameters_module
use mpi_module
use mesh_module
use storage_module
use viscous_shared_module
use viscous_laplacian_module
use viscous_staggered_module
use viscous_conservative_module

implicit none
private
public viscous_fluxes

contains

subroutine viscous_fluxes(dims)

        implicit none
        integer, intent(in) :: dims

        call StartProfRange("viscous_fluxes")

        if    (diffusion_scheme == 'laplacian') then
            call viscous_flux_3D_laplacian
            call viscous_heatflux_3D_laplacian
            call DivergencyGradient3D(U,V,W,VIS,RHS)
        elseif(diffusion_scheme == 'staggered') then
            call viscous_flux_3D_staggered
            if(non_orthogonal) call cross_flux_3D_staggered
            call heatflux_3D_staggered
            call DivergencyGradient3D(U,V,W,VIS,RHS)
        elseif(diffusion_scheme == 'conservative') then
            call viscous_faces_conservative
            if(inplane_axis_aligned) then
                call viscous_flux_3D_conservative_xi_cart
                call viscous_flux_3D_conservative_eta_cart
                if(dims == 3) call viscous_flux_3D_conservative_zta_cart
            else
                call viscous_flux_3D_conservative_xi_curv
                call viscous_flux_3D_conservative_eta_curv
                if(dims == 3) call viscous_flux_3D_conservative_zta_curv
            endif
            if(LES .and. yoshizawa) call divergence_3D_conservative
        else
            if(rank == root) print*, "Diffusion scheme ", trim(diffusion_scheme), " is not implemented"
            call secure_stop
        endif

        if(LES) then
          if(yoshizawa) call ComputeSGSIsotropicTerm()
        end if

        call EndProfRange

        return
end subroutine viscous_fluxes

end module viscous_module
