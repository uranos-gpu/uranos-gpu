module nscbc_boundary_conditions
use parameters_module
use mpi_module
use mesh_module
use storage_module
implicit none

private
public refl_wall

contains

subroutine refl_wall(phi,bnode,bnorm,bface,Twall)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        integer                                  , intent(in)    :: bnode
        integer                                  , intent(in)    :: bnorm
        character(1)                             , intent(in)    :: bface
        real(rp)                                 , intent(in)    :: Twall

        real(rp) :: rho, uu, vv, ww, tt, h, qq, cc, c, ci
        real(rp) :: b1, b2, b3
        real(rp) :: p_rho, p_e, etot, df
        real(rp) :: ttW,uuW,vvW,wwW,rrW,rEW,yW,y1,y2,idy
        integer  :: i,j,k,j0,j1,j2,m,mm
        real(rp), dimension(5) :: dw_dn, dwc_dn, ev
        real(rp), dimension(5,5) :: el, er
        integer  :: node, norm

        node = bnode
        norm = bnorm
        j0 = bnode + bnorm
        j1 = bnode
        j2 = bnode - bnorm

        ttW = TWall
        uuW = 0.0_rp
        vvW = 0.0_rp
        wwW = 0.0_rp
        yW = 0.5_rp*(y(j1) + y(j0))
        y1 = abs(y(j1)-yW)
        y2 = abs(y(j2)-yW)
        idy = 1._rp/(abs(y(j1) - y(j0)))

        selectcase(bface)
        case('S','N')

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2) private(dw_dn,dwc_dn,ev,el,er)
        do        k = sz,ez
           do     i = sx,ex


               rrW = 0.5_rp*(phi(i,j0,k,1)+phi(i,j1,k,1))
               rEW = rrW*ttW/(gamma0-1._rp)
        
               ! wall derivative of conservative variables
               dw_dn(1) = phi(i,j1,k,1)-phi(i,j0,k,1) 
               dw_dn(2) = phi(i,j1,k,2)-phi(i,j0,k,2) 
               dw_dn(3) = phi(i,j1,k,3)-phi(i,j0,k,3) 
               dw_dn(4) = phi(i,j1,k,4)-phi(i,j0,k,4) 
               dw_dn(5) = phi(i,j1,k,5)-phi(i,j0,k,5) 

               !Compute eigenvectors
               rho = rrW
               uu  = uuW
               vv  = vvW
               ww  = wwW
               tt  = ttw
               h   = cp*tt
               qq  = 0.5_rp*(uu*uu  +vv*vv + ww*ww)
               cc  = gamma0*tt
               c   = sqrt(cc)
               ci  =  1._rp/c
               p_rho     = tt 
               p_e = rho*(gamma0-1._rp)
               etot = h - tt
               
               b3 = etot - rho*p_rho/p_e
               b2 = p_e/(rho*cc)
               b1 = p_rho/cc - b2*(etot - 2._rp*qq)
               
               el(1,1) =   0.5_rp * (b1     + vv * ci)
               el(2,1) =  -0.5_rp * (b2 * uu         )
               el(3,1) =  -0.5_rp * (b2 * vv +     ci)
               el(4,1) =  -0.5_rp * (b2 * ww         )
               el(5,1) =   0.5_rp * b2
               el(1,2) =   1._rp - b1
               el(2,2) =   b2 * uu
               el(3,2) =   b2 * vv
               el(4,2) =   b2 * ww
               el(5,2) =  -b2
               el(1,3) =   0.5_rp * (b1     - vv * ci)
               el(2,3) =  -0.5_rp * (b2 * uu         )
               el(3,3) =  -0.5_rp * (b2 * vv -     ci)
               el(4,3) =  -0.5_rp * (b2 * ww         )
               el(5,3) =   0.5_rp * b2
               el(1,4) =  -ww 
               el(2,4) =   0._rp 
               el(3,4) =   0._rp
               el(4,4) =   1._rp 
               el(5,4) =   0._rp
               el(1,5) =   uu 
               el(2,5) =  -1._rp 
               el(3,5) =   0._rp
               el(4,5) =   0._rp 
               el(5,5) =   0._rp

               er(1,1) =  1._rp
               er(2,1) =  1._rp
               er(3,1) =  1._rp
               er(4,1) =  0._rp
               er(5,1) =  0._rp
               er(1,2) =  uu
               er(2,2) =  uu
               er(3,2) =  uu
               er(4,2) =  0._rp 
               er(5,2) = -1._rp 
               er(1,3) =  vv - c
               er(2,3) =  vv
               er(3,3) =  vv + c
               er(4,3) =  0._rp
               er(5,3) =  0._rp
               er(1,4) =  ww
               er(2,4) =  ww
               er(3,4) =  ww
               er(4,4) =  1._rp  
               er(5,4) =  0._rp 
               er(1,5) =  h  - vv * c
               er(2,5) =  b3 
               er(3,5) =  h  + vv * c
               er(4,5) =  ww     
               er(5,5) =  -uu    

               !$acc loop seq
               do m=1,5
                dwc_dn(m) = 0._rp
                !$acc loop seq
                 do mm=1,5 
                  dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
                 enddo
               enddo

               ! Compute eigenvalues
               ev(1) = vv-c
               ev(2) = vv
               ev(3) = vv+c
               ev(4) = ev(2)
               ev(5) = ev(2)

               !$acc loop seq
               do m=1,5
                dwc_dn(m) = ev(m) * dwc_dn(m)
               enddo

               dwc_dn(2) = 0._rp
               if(bface.eq.'S') then
                 dwc_dn(3) = dwc_dn(1)
               elseif(bface.eq.'N') then
                 dwc_dn(1) = dwc_dn(3)
               endif
               dwc_dn(4) = 0._rp
               dwc_dn(5) = 0._rp

               !$acc loop seq
               do m=1,5
                df = 0._rp
                !$acc loop seq
                do mm=1,5
                 df = df + er(mm,m) * dwc_dn(mm)
                enddo
                RHS(i,j1,k,m) = RHS(i,j1,k,m) - df*idy !- RHS(i,j0,k,1) 
               enddo

           enddo
        enddo
        !$acc end parallel
        
        case default 
          print*, 'refl_wall is not implemented for face ', trim(bface)
          stop

        endselect

        return
end subroutine refl_wall
        


end module nscbc_boundary_conditions
