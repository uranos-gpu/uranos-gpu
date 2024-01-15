module real_to_integer_module
#ifdef TIME
use performance_module
#endif

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

        ! local declarations
        integer            :: lo, hi, mid

        lo = lb
        hi = ub
        il = lo

        ! Binary search
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x(mid) <= xl) then
                lo = mid
            else
                hi = mid
            end if
        end do

        ! Determine the floor index
        if (x(il) <= xl) then
            il = lo
        else
            il = lo - 1
        end if

        return
end function optFloor




function locFloor(x,lb,ub,xl) result(il)
        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:), intent(in) :: x
        integer                            , intent(in) :: lb, ub
        real(rp)                           , intent(in) :: xl
        integer                                         :: il

        real(rp) :: d, d_old
        integer  :: i

        il    = lb
        d_old = abs(xl - x(lb))

        do i = lb+1,ub

           d = abs(xl - x(i))
           if(d < d_old) then
              il = i
              d_old = d
           endif

        enddo
        
        ! get the floor
        if(x(il) > xl) il = il-1

        return
end function locFloor



function locNearest(x,lb,ub,xl) result(il)

        use parameters_module, only: rp

        implicit none
        real(rp), dimension(:), allocatable, intent(in) :: x
        integer                            , intent(in) :: lb, ub
        real(rp)                           , intent(in) :: xl
        integer                                         :: il

        real(rp) :: d, d_old
        integer  :: i

        il    = lb
        d_old = abs(xl - x(lb))

        do i = lb+1,ub

           d = abs(xl - x(i))
           if(d < d_old) then
              il = i
              d_old = d
           endif

        enddo

        return
end function locNearest









subroutine locate(x,sx,ex,xl,j)
! ---------------------------------------------------------------------------
!
!       Computation of the nearest integer position in respect of a point
!       
!       INPUT:  x        !< coordinate along wich we wan to find the position
!               sx       !< starting index of the cooerdinate
!               ex       !< ending   index of the cooerdinate
!               xl       !< real point we wanto to locate
!       OUTPUT: j        !< nearest integer position of xl
!       
! ---------------------------------------------------------------------------
        use parameters_module, only: rp

        implicit none
        integer                   , intent(in)  :: sx,ex
        real(rp), dimension(sx:ex), intent(in)  :: x
        real(rp)                  , intent(in)  :: xl
        integer                   , intent(out) :: j

        ! local 
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

function nearest_integer_opt(x,sx,ex,pt) result(ix)

        use parameters_module, only: rp

        implicit none
        real(rp), dimension(:), allocatable, intent(in) :: x
        integer                            , intent(in) :: sx, ex
        real(rp)                           , intent(in) :: pt
        integer                                         :: ix
#ifdef TIME
        call mpi_stime(s_nei_time)
#endif
        
        ix = minloc( abs(x(sx:ex) - pt), 1) + sx - 1 

#ifdef TIME
        call mpi_etime(s_nei_time,t_nei_calls,t_nei_time)
#endif
        return
end function nearest_integer_opt



end module real_to_integer_module
