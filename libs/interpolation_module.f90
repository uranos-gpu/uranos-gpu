module interpolation_module
! -----------------------------------------------------------------------------------------
!
!  This module produces an interface to interpolate the value of a generic field in a point.
!  USAGE:
!  - interpolate a field in 2 dimensions:
!
!       field_in_pt = interpl_field_2D(pt,x,y,field)     !< for 3D field
!       field_in_pt = interpl_field_2D(pt,x,y,field,id)  !< for 3D field with co-dimension (id)
!
! -----------------------------------------------------------------------------------------

private
public interpl_field_1D, interpl_field_2D, interpl_field_3D, polint

interface interpl_field_1D
  module procedure interpl_1Dfield_1D
end interface

interface interpl_field_2D
  module procedure interpl_3Dfield_2D, interpl_4Dfield_2D
endinterface

interface interpl_field_3D
  module procedure interpl_4Dfield_3D
endinterface

contains
function interpl_1Dfield_1D(pt,x,field1D) result(field_pt)

        use parameters_module       , only: rp
        use mpi_module              , only: lbx,ubx
        use real_to_integer_module  , only: optFloor

        implicit none
        real(rp)                           , intent(in) :: pt
        real(rp), allocatable, dimension(:), intent(in) :: x
        real(rp), allocatable, dimension(:), intent(in) :: field1D
        real(rp)                                        :: field_pt

        ! local declarations
        integer  :: il
        real(rp) :: x0, x1, phi0, phi1

        ! ===== find the integer position of the cell
        il = optFloor(x,lbx,ubx,pt)
       
        ! ===== compute the line passing through two points
        x0   = x(il)
        x1   = x(il+1)

        phi0 = field1D(il)
        phi1 = field1D(il+1)

        ! ===== interpolate the field linearly
        field_pt = phi0 + (phi1 - phi0)/(x1 - x0) * (pt - x0)
        
        return
end function interpl_1Dfield_1D



function interpl_3Dfield_2D(pt,x,y,k,field3D) result(field_pt)
! -------------------------------------------------------------------
!
!       This function interpolate the value of a 3D field in a point.
!
!       IN : pt          !< point in which we want to interpolate the field
!            k           !< plane index
!            field3D     !< field we want to interpolate
!       OUT: field_pt    !> interpolated value
!
! -------------------------------------------------------------------

        use parameters_module       , only: rp
        use mpi_module              , only: lbx, ubx, lby, uby
        use real_to_integer_module  , only: optFloor
        use matrix_inversion_module , only: matinv4

        implicit none
        real(rp), dimension(3)                 , intent(in) :: pt
        real(rp), dimension(:)    , allocatable, intent(in) :: x, y
        real(rp), dimension(:,:,:), allocatable, intent(in) :: field3D
        integer                                , intent(in) :: k
        real(rp)                                            :: field_pt

        ! local declarations
        integer, dimension(2)   :: il
        integer, dimension(4)   :: ic_, jc_
        real(rp),dimension(4)   :: phi_i,c
        real(rp),dimension(4,4) :: V
        real(rp)                :: xi, yj
        integer                 :: i0, i1, j0, j1, i,j,ib
        
        ! ===== find the integer position of the cell
        il(1) = optFloor(x,lbx,ubx,pt(1))
        il(2) = optFloor(y,lby,uby,pt(2))

        ! ===== store the nodes in a array
        i0 = il(1); i1 = il(1)+1
        j0 = il(2); j1 = il(2)+1

        ic_ = (/i0,i1,i0,i1/)
        jc_ = (/j0,j0,j1,j1/)
        
        ! ===== interpolating the value of the cell
        do ib = 1,4

           i = ic_(ib)
           j = jc_(ib)

           xi = x(i)
           yj = y(j)

           V(ib,:)   = [xi*yj, xi, yj, 1.0_rp]
           phi_i(ib) = field3D(i,j,k)

        enddo

        c = matmul(matinv4(V),phi_i)

        field_pt = c(1)*pt(1)*pt(2) + c(2)*pt(1) + c(3)*pt(2) + c(4)

        return
end function interpl_3Dfield_2D



function interpl_4Dfield_2D(pt,x,y,k,field4D,id) result(field_pt)
! -------------------------------------------------------------------
!
!       This function interpolate the value of a 4D field in a point.
!
!       IN : id          !< index of the field co-dimension
!            pt          !< point in which we want to interpolate the field
!            k           !< plane index
!            field4D     !< field we want to interpolate
!       OUT: field_pt    !> interpolated value
!
! -------------------------------------------------------------------

        use parameters_module       , only: rp
        use mpi_module              , only: lbx,ubx, lby,uby
        use real_to_integer_module  , only: optFloor
        use matrix_inversion_module , only: matinv4

        implicit none
        real(rp), dimension(3)                   , intent(in) :: pt
        real(rp), dimension(:)      , allocatable, intent(in) :: x,y
        real(rp), dimension(:,:,:,:), allocatable, intent(in) :: field4D
        integer                                  , intent(in) :: id, k
        real(rp)                                              :: field_pt

        ! local declarations
        integer, dimension(2)   :: il
        integer, dimension(4)   :: ic_, jc_
        real(rp),dimension(4)   :: phi_i,c
        real(rp),dimension(4,4) :: V
        real(rp)                :: xi, yj
        integer                 :: i0, i1, j0, j1, i, j, ib
        
        ! ===== find the integer position of the cell
        il(1) = optFloor(x,lbx,ubx,pt(1))
        il(2) = optFloor(y,lby,uby,pt(2))

        ! ===== store the nodes in a array
        i0 = il(1); i1 = il(1)+1
        j0 = il(2); j1 = il(2)+1

        ic_ = (/i0,i1,i0,i1/)
        jc_ = (/j0,j0,j1,j1/)

        ! ===== interpolating the value of the cell
        do ib = 1,4

           i = ic_(ib)
           j = jc_(ib)

           xi = x(i)
           yj = y(j)

           V(ib,:)   = [xi*yj, xi, yj, 1.0_rp]
           phi_i(ib) = field4D(i,j,k,id)

        enddo

        c = matmul(matinv4(V),phi_i)

        field_pt = c(1)*pt(1)*pt(2) + c(2)*pt(1) + c(3)*pt(2) + c(4)

        return
end function interpl_4Dfield_2D




function interpl_4Dfield_3D(pt,x,sx,ex,y,sy,ey,z,sz,ez,field4D,id) result(field_pt)
! -------------------------------------------------------------------
!
!       This function interpolate the value of a 4D field in a point (xp,yp,zp).
!
!       IN : pt          !< point in which we want to interpolate the field
!            field4D     !< field we want to interpolate
!            id          !< index of the field
!       OUT: field_pt    !> interpolated value
!
! -------------------------------------------------------------------

        use parameters_module       , only: rp
        use real_to_integer_module  , only: optFloor
        use matrix_inversion_module , only: gauss_elimination

        implicit none
        real(rp), dimension(3)                   , intent(in) :: pt
        real(rp), dimension(:)      , allocatable, intent(in) :: x, y, z
        real(rp), dimension(:,:,:,:), allocatable, intent(in) :: field4D
        integer                                  , intent(in) :: id
        integer                                  , intent(in) :: sx,ex,sy,ey,sz,ez
        real(rp)                                              :: field_pt

        ! local declarations
        integer, dimension(3)   :: il
        integer, dimension(8)   :: ic_, jc_, kc_
        real(rp),dimension(8)   :: phi_i,c, pt_array
        real(rp),dimension(8,8) :: V
        real(rp)                :: xi, yj, zk
        integer                 :: i0, i1,j0,j1,k0,k1,i,j,k,ib
        
        ! ===== find the integer position of the cell
        il(1) = optFloor(x,sx,ex,pt(1))
        il(2) = optFloor(y,sy,ey,pt(2))
        il(3) = optFloor(z,sz,ez,pt(3))

        ! ===== store the nodes in a array
        i0 = il(1); i1 = il(1)+1
        j0 = il(2); j1 = il(2)+1
        k0 = il(3); k1 = il(3)+1

        ic_ = (/i0,i1,i1,i0,i0,i1,i1,i0/)
        jc_ = (/j0,j0,j1,j1,j0,j0,j1,j1/)
        kc_ = (/k0,k0,k0,k0,k1,k1,k1,k1/)
        
        ! ===== interpolating the value of the cell
        do ib = 1,8

           i = ic_(ib)
           j = jc_(ib)
           k = kc_(ib)

           xi = x(i)
           yj = y(j)
           zk = z(k)

           V(ib,:)   = [xi*yj*zk, xi*yj, xi*zk, yj*zk, xi, yj, zk, 1.0_rp]
           phi_i(ib) = field4D(i,j,k,id)

        enddo

        call gauss_elimination(V,phi_i,c,8)

        pt_array = (/pt(1)*pt(2)*pt(3), pt(1)*pt(2), pt(1)*pt(3), pt(2)*pt(3), pt(1), pt(2), pt(3), 1.0_rp/)

        field_pt = dot_product(c,pt_array)

        return
end function interpl_4Dfield_3D



      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      use parameters_module       , only: rp

      implicit none

      INTEGER , PARAMETER :: NMAX=10
      INTEGER :: N
      REAL(rp) :: X,Y,DY
      REAL(rp) , DIMENSION(N) :: XA,YA
      REAL(rp) , DIMENSION(NMAX) ::C,D

      real(rp) :: DEN,HO,DIF,DIFT,HP,W
      integer :: I,NS,M

      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
!         IF(DEN.EQ.0._rp) PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
!







end module interpolation_module
