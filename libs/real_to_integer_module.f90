module real_to_integer_module

implicit none
contains

integer function RealToInt(x,lbx,ubx,xi) result(minIndex)
        !$acc routine seq
        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:), intent(in) :: x
        real(rp)                           , intent(in) :: xi
        integer                            , intent(in) :: lbx,ubx

        real(rp) :: minDistance, distance
        integer  :: left, right, mid

        left        = lbx
        right       = ubx
        minIndex    = left
        minDistance = abs(xi - x(left))

        do while (left <= right)

           mid = (left + right)/2

           distance = abs(xi - x(mid))

           if(distance < minDistance) then
             minDistance = distance
             minIndex = mid
           endif

           if(xi < x(mid)) then
             right = mid - 1
           else
             left  = mid + 1
           endif

        enddo

        return
end function RealToInt

function optFloor(x,lb,ub,xl) result(il)

        use parameters_module, only: rp

        implicit none
        real(rp), dimension(:), allocatable, intent(in) :: x
        integer                            , intent(in) :: lb, ub
        real(rp)                           , intent(in) :: xl
        integer                                         :: il

        integer            :: lo, hi, mid

        lo = lb
        hi = ub
        il = lo

        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x(mid) <= xl) then
                lo = mid
            else
                hi = mid
            end if
        end do

        if (x(il) <= xl) then
            il = lo
        else
            il = lo - 1
        end if

        return
end function optFloor

subroutine locate(x,sx,ex,xl,j)
        use parameters_module, only: rp

        implicit none
        integer                   , intent(in)  :: sx,ex
        real(rp), dimension(sx:ex), intent(in)  :: x
        real(rp)                  , intent(in)  :: xl
        integer                   , intent(out) :: j

        integer :: ju, jl, jm

        jl = sx-1
        ju = ex+1

     10 if(ju-jl .gt. 1) then
          jm = (ju+jl)/2

          if((x(ex).gt.x(sx)) .eqv. (xl.gt.x(jm))) then
            jl = jm
          else
            ju = jm
          endif

          goto 10

        endif

        j = jl

        return
end subroutine locate

end module real_to_integer_module
