module dirac_delta_module
use parameters_module, only: rp, pi

integer, parameter :: nbox = 2

interface interpl_diracdelta
  module procedure interpl_diracdelta_3D, interpl_diracdelta_2D
end interface

contains
function interpl_diracdelta_3D(pt,x,y,z,field) result(field_in_pt)

        use real_to_integer_module, only: nearest_integer_opt

        implicit none
        real(rp), dimension(3)                 , intent(in) :: pt
        real(rp), dimension(:)    , allocatable, intent(in) :: x,y,z
        real(rp), dimension(:,:,:), allocatable, intent(in) :: field
        real(rp)                                            :: field_in_pt

        ! local declaration
        real(rp), dimension(3) :: step, dlta
        integer , dimension(3) :: sb,eb
        integer , dimension(3) :: il
        integer                :: i,j,k
        real(rp)               :: delta_dv, xi, yj, zk

        sb(1) = lbound(x,1) + 4
        sb(2) = lbound(y,1) + 4
        sb(3) = lbound(z,1) + 4

        eb(1) = ubound(x,1) - 4
        eb(2) = ubound(y,1) - 4
        eb(3) = ubound(z,1) - 4

        il(1) = nearest_integer_opt(x,sb(1),eb(1),pt(1))
        il(2) = nearest_integer_opt(y,sb(2),eb(2),pt(2))
        il(3) = nearest_integer_opt(z,sb(3),eb(3),pt(3))
        
        ! init field_in_pt
        field_in_pt = 0.0_rp
        do k = il(3)-nbox, il(3)+nbox
        
           zk      = z(k)
           step(3) = zk - z(k-1)
           dlta(3) = delta(zk - pt(3), step(3))

           do j = il(2)-nbox, il(2)+nbox

              yj      = y(j)
              step(2) = yj - y(j-1)
              dlta(2) = delta(yj - pt(2), step(2))

              do i = il(1)-nbox, il(1)+nbox

                 xi      = x(i)
                 step(1) = xi - x(i-1)
                 dlta(1) = delta(xi - pt(1), step(1))
                
                 delta_dv = dlta(1) * step(1) * &
                            dlta(2) * step(2) * &
                            dlta(3) * step(3)

                 field_in_pt = field_in_pt + field(i,j,k) * delta_dv

              enddo
           enddo
        enddo


        return
end function interpl_diracdelta_3D



function interpl_diracdelta_2D(pt,x,y,k,field) result(field_in_pt)

        use real_to_integer_module, only: nearest_integer_opt

        implicit none
        real(rp), dimension(3)                 , intent(in) :: pt
        real(rp), dimension(:)    , allocatable, intent(in) :: x,y
        real(rp), dimension(:,:,:), allocatable, intent(in) :: field
        integer                                , intent(in) :: k
        real(rp)                                            :: field_in_pt

        ! local declaration
        real(rp), dimension(2) :: step, dlta
        integer , dimension(2) :: sb,eb
        integer , dimension(2) :: il
        integer                :: i,j
        real(rp)               :: delta_dv, xi, yj

        sb(1) = lbound(x,1) + 4
        sb(2) = lbound(y,1) + 4

        eb(1) = ubound(x,1) - 4
        eb(2) = ubound(y,1) - 4

        il(1) = nearest_integer_opt(x,sb(1),eb(1),pt(1))
        il(2) = nearest_integer_opt(y,sb(2),eb(2),pt(2))
        
        ! init field_in_pt
        field_in_pt = 0.0_rp
        do j = il(2)-nbox, il(2)+nbox

           yj      = y(j)
           step(2) = yj - y(j-1)
           dlta(2) = delta(yj - pt(2), step(2))

           do i = il(1)-nbox, il(1)+nbox

              xi      = x(i)
              step(1) = xi - x(i-1)
              dlta(1) = delta(xi - pt(1), step(1))
             
              delta_dv = dlta(1) * step(1) * &
                         dlta(2) * step(2)

              field_in_pt = field_in_pt + field(i,j,k) * delta_dv

           enddo
        enddo
        return
end function interpl_diracdelta_2D






function delta_Roma(x,h) result(DiracDelta)
!-----------------------------------------------------------------------
! function that evaluates Dirac delta in a point r (1 point needed)
!
! REFERENCE : Roma et. al, 
!             "An adaptive version of immersed boundary ..."
!             Journal of Computational Physics, pag. 519
!-----------------------------------------------------------------------
        implicit none
        real(rp), intent(in) :: x              !< dirac delta argument
        real(rp), intent(in) :: h              !< grid space
        real(rp)             :: abs_r          !< phi argument
        real(rp)             :: un_h           !< inverse of grid step
        real(rp)             :: DiracDelta     !< Dirac delta
        
        ! calculate arguments
        un_h = 1._rp/h
        abs_r = abs(x*un_h)
        
        ! calculate Dirac delta function
        DiracDelta = 0._rp
        if((abs_r.ge.0.5_rp).and.(abs_r.le.1.5_rp)) then
                
                DiracDelta = un_h * 1._rp/6._rp * (5._rp - 3._rp*abs_r - sqrt(-3._rp*(1._rp-abs_r)**2 + 1._rp))
        
        elseif((abs_r.le.0.5_rp)) then

                DiracDelta = un_h * 1._rp/3._rp * (1._rp + sqrt(-3._rp*abs_r**2 + 1))

        endif
        
return
end function delta_Roma


function delta_LiWang(x,h) result(DiracDelta)
!-----------------------------------------------------------------------
! function that evaluates Dirac delta in a point r (2 points needed)
!
! REFERENCE : A. Li Wang et. al, 
!             "An immersed boundary method for fluid structure interection ..."
!             Journal of Computational Physics, pag. 137
!-----------------------------------------------------------------------
        implicit none
        real(rp), intent(in) :: x              !< dirac delta argument
        real(rp), intent(in) :: h              !< grid space
        real(rp)             :: abs_r          !< lambda argument
        real(rp)             :: un_h           !< inverse of grid step
        real(rp)             :: DiracDelta     !< Dirac delta
        
        ! calculate lambda argument
        un_h = 1._rp/h
        abs_r = abs(un_h*x)
        
        ! calculate lambda function
        DiracDelta = 0._rp
        if((abs_r.ge.0._rp).and.(abs_r.lt.1._rp))then
                
                DiracDelta = un_h * 1._rp/8._rp * (3._rp - 2._rp*abs_r + sqrt(1._rp + 4._rp*abs_r - 4._rp*abs_r**2))
        
        elseif((abs_r.ge.1._rp).and.(abs_r.lt.2._rp))then
        
                DiracDelta = un_h * 1._rp/8._rp * (5._rp - 2._rp*abs_r - sqrt(-7._rp + 12._rp*abs_r - 4._rp*abs_r**2)) 

        endif
return
end function delta_LiWang


function delta_phi3star(x,h) result(DiracDelta)
!-----------------------------------------------------------------------
! function that evaluates Dirac delta in a point r (2 points needed)
!
! REFERENCE : Xiaolei Yang et. al, 
!             "A smoothing technique for discrete delta functions ..."
!             Journal of Computational Physics, pag. 7821
!-----------------------------------------------------------------------
        implicit none
        real(rp), intent(in) :: x              !< dirac delta argument
        real(rp), intent(in) :: h              !< grid space
        real(rp)             :: abs_r          !< lambda argument
        real(rp)             :: un_h           !< inverse of grid step
        real(rp)             :: DiracDelta     !< Dirac delta

        ! calculate lambda argument
        un_h = 1._rp/h
        abs_r = abs(un_h*x)

        DiracDelta = 0._rp
        if((abs_r.ge.0._rp).and.(abs_r.lt.1._rp))then

                DiracDelta = 17._rp/48._rp + sqrt(3._rp)*pi/(108._rp) + abs_r/4._rp - abs_r**2/4._rp    &
                           + (1._rp-2*abs_r)/(16._rp) * sqrt(-12._rp*abs_r**2 + 12._rp*abs_r+1._rp)  &
                           - sqrt(3._rp)/12._rp*asin(sqrt(3._rp)/2._rp*(2*abs_r-1._rp))

        elseif((abs_r.ge.1._rp).and.(abs_r.lt.2._rp))then

                DiracDelta = 55._rp/48._rp - sqrt(3._rp)*pi/108._rp - 13._rp*abs_r/12._rp + abs_r**2/4._rp &
                           + (2*abs_r-3)/48._rp * sqrt(-12._rp*abs_r**2 + 36._rp*abs_r - 23._rp)  &
                           + sqrt(3._rp)/36._rp*asin(sqrt(3._rp)/2._rp*(2*abs_r-3._rp))

        endif

        DiracDelta = un_h * DiracDelta

        return
end function delta_phi3star


function delta(x,h) result(DiracDelta)
!-----------------------------------------------------------------------
! function that evaluates Dirac delta in a point r (2 points needed)
!
! REFERENCE : Qiu et al.
!             'A boundary condition-enforced immersed boundary method 
!              for compressible viscous flows'
!             Computers and fluids
!-----------------------------------------------------------------------
        implicit none
        real(rp), intent(in) :: x              !< dirac delta argument
        real(rp), intent(in) :: h              !< grid space
        real(rp)             :: abs_r          !< lambda argument
        real(rp)             :: un_h           !< inverse of grid step
        real(rp)             :: DiracDelta     !< Dirac delta
        
        ! calculate lambda argument
        un_h = 1._rp/h
        abs_r = abs(un_h*x)

        DiracDelta = 0.0_rp
        if    (abs_r .le. 0.5_rp) then

                DiracDelta = 3._rp/8._rp + pi/32._rp - abs_r**2/4._rp

        elseif((abs_r .gt. 0.5_rp) .and. (abs_r.le.1.5_rp)) then

                DiracDelta = 1._rp/4._rp+1._rp/8._rp*(1._rp-abs_r)*sqrt(-2._rp+8._rp*abs_r-4._rp*abs_r**2) - &
                             1._rp/8._rp*asin(sqrt(2._rp)*(abs_r -1._rp))

        elseif((abs_r.gt.1.5_rp) .and. (abs_r.le.2.5_rp)) then

                DiracDelta = 17._rp/16._rp - pi/64._rp + ((abs_r-2._rp)/16._rp)*(sqrt(-14._rp+16._rp*abs_r-4._rp*abs_r**2)) + &
                           & 1._rp/16._rp*asin(sqrt(2._rp)*(abs_r -2._rp)) + abs_r**2/8._rp - 3._rp*abs_r/4._rp

        endif

        DiracDelta = DiracDelta*un_h

        return
endfunction delta

function my_delta(x,h) result(DiracDelta)
        implicit none
        real(rp), intent(in) :: x              !< dirac delta argument
        real(rp), intent(in) :: h              !< grid space
        real(rp)             :: abs_r          !< lambda argument
        real(rp)             :: i_h            !< inverse of grid step
        real(rp)             :: DiracDelta     !< Dirac delta
        real(rp) :: int_
        real(rp) :: n = 2.5_rp

        i_h = 1.0_rp/h
        abs_r = abs(x*i_h)
        int_ = sqrt(pi/n)

        DiracDelta = i_h * exp(-n*abs_r**2)/int_

        return
end function my_delta









end module dirac_delta_module

      
