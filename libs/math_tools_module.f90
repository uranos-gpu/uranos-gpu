module math_tools_module
use parameters_module, only: rp

implicit none

interface newton_raphson
  module procedure newton_raphson_0par, newton_raphson_2par, newton_raphson_3par
end interface

interface det
  module procedure det3
end interface


interface QuickSort
  module procedure QuickSort2D
end interface



contains
function MatVec(A,b,n) result(c)
        implicit none
        integer , intent(in) :: n 
        real(rp), intent(in) :: A(n,n), b(n)
        real(rp)             :: c(n)
        integer              :: i,j

        do i = 1,n
           c(i) = 0.0_rp
           do j = 1,n
              c(i) = c(i) + A(i,j)*b(j)
           enddo
        enddo

        return
end function MatVec





subroutine newton_raphson_0par(f,itmax, toll, x0, xnew)
! ------------------------------------------------------------------------------------------
!       
!       This subroutine implements the Newton-Rapson algorithm for a scalar function f(x).
!       INPUT:  f        !< funtion
!               itmax    !< number of iterations
!               toll     !< tollerance to convergency
!               x0       !< first guess point
!       OUTPUT: xnew     !> solution
!
! ------------------------------------------------------------------------------------------
        implicit none

        interface 
          function f(x) result(fx)
            integer , parameter  :: rp = 8
            real(rp), intent(in) :: x
            real(rp)             :: fx

          end function
        end interface

        integer , intent(in)  :: itmax 
        real(rp), intent(in)  :: x0   
        real(rp), intent(in)  :: toll    
        real(rp), intent(out) :: xnew

        ! local declarations
        integer  :: iter                     !< iter variable
        real(rp) :: ctoll                    !< calculated tolerance
        real(rp) :: xold1, xold2, hk
        
        ! initializin the loop
        xold1 = x0
        xold2 = x0 + 1.E-14_rp
        ctoll = 2._rp*toll
        iter  = 0

        do while ((ctoll.ge.toll).and.(iter.le.itmax))
           iter = iter + 1
           
           ! calculting finite difference
           hk = (f(xold1) - f(xold2))/(xold1 - xold2)
           
           ! update x with Newton-Raphson formula
           xnew = xold1 - f(xold1)/hk

           ! update tollerance
           ctoll = abs(xnew-xold1)

           ! update loop points
           xold2 = xold1
           xold1 = xnew
        enddo

        return
end subroutine newton_raphson_0par



subroutine newton_raphson_2par(f,par1,par2,itmax, toll, x0, xnew)
! ------------------------------------------------------------------------------------------
!       
!       This subroutine implements the Newton-Rapson algorithm for a scalar function f(x).
!       INPUT:  f        !< funtion
!               itmax    !< number of iterations
!               toll     !< tollerance to convergency
!               x0       !< first guess point
!       OUTPUT: xnew     !> solution
!
! ------------------------------------------------------------------------------------------
        implicit none

        interface 
          function f(x,par1,par2) result(fx)
            integer , parameter  :: rp = 8
            real(rp), intent(in) :: x
            real(rp), intent(in) :: par1
            real(rp), intent(in) :: par2
            real(rp)             :: fx

          end function
        end interface

        integer , intent(in)  :: itmax 
        real(rp), intent(in)  :: x0   
        real(rp), intent(in)  :: par1
        real(rp), intent(in)  :: par2
        real(rp), intent(in)  :: toll    
        real(rp), intent(out) :: xnew

        ! local declarations
        integer  :: iter                     !< iter variable
        real(rp) :: ctoll                    !< calculated tolerance
        real(rp) :: xold1, xold2, hk
        
        ! initializin the loop
        xold1 = x0
        xold2 = x0 + 1.E-14_rp
        ctoll = 2._rp*toll
        iter  = 0

        do while ((ctoll.ge.toll).and.(iter.le.itmax))
           iter = iter + 1
           
           ! calculting finite difference
           hk = (f(xold1,par1,par2) - f(xold2,par1,par2))/(xold1 - xold2)
           
           ! update x with Newton-Raphson formula
           xnew = xold1 - f(xold1,par1,par2)/hk

           ! update tollerance
           ctoll = abs(xnew-xold1)

           ! update loop points
           xold2 = xold1
           xold1 = xnew
        enddo

        return
end subroutine newton_raphson_2par





subroutine newton_raphson_3par(f,par1,par2,par3,itmax, toll, x0, xnew)
! ------------------------------------------------------------------------------------------
!       
!       This subroutine implements the Newton-Rapson algorithm for a scalar function f(x).
!       INPUT:  f        !< funtion
!               itmax    !< number of iterations
!               toll     !< tollerance to convergency
!               x0       !< first guess point
!       OUTPUT: xnew     !> solution
!
! ------------------------------------------------------------------------------------------
        implicit none

        interface 
          function f(x,par1,par2,par3) result(fx)
            integer , parameter  :: rp = 8
            real(rp), intent(in) :: x
            real(rp), intent(in) :: par1
            real(rp), intent(in) :: par2
            real(rp), intent(in) :: par3
            real(rp)             :: fx

          end function
        end interface

        integer , intent(in)  :: itmax 
        real(rp), intent(in)  :: x0   
        real(rp), intent(in)  :: par1
        real(rp), intent(in)  :: par2
        real(rp), intent(in)  :: par3
        real(rp), intent(in)  :: toll    
        real(rp), intent(out) :: xnew

        ! local declarations
        integer  :: iter                     !< iter variable
        real(rp) :: ctoll                    !< calculated tolerance
        real(rp) :: xold1, xold2, hk
        
        ! initializin the loop
        xold1 = x0
        xold2 = x0 + 1.E-14_rp
        ctoll = 2._rp*toll
        iter  = 0

        do while ((ctoll.ge.toll).and.(iter.le.itmax))
           iter = iter + 1
           
           ! calculting finite difference
           hk = (f(xold1,par1,par2,par3) - f(xold2,par1,par2,par3))/(xold1 - xold2)
           
           ! update x with Newton-Raphson formula
           xnew = xold1 - f(xold1,par1,par2,par3)/hk

           ! update tollerance
           ctoll = abs(xnew-xold1)

           ! update loop points
           xold2 = xold1
           xold1 = xnew
        enddo

        return
end subroutine newton_raphson_3par



function dist(a,b) result(d)
! ---------------------------------------------------
!       Computation of the distance between to points
! ---------------------------------------------------
        implicit none
        real(rp), dimension(:), intent(in) :: a, b
        real(rp)                           :: d

        d = sqrt(dot_product(a-b,a-b))

        return
end function dist


function cross_product(a,b) result(c)
! -------------------------------------------------
!       Computation of the cross product c = a ^ b
! -------------------------------------------------
        implicit none
        real(rp), dimension(3), intent(in) :: a,b
        real(rp), dimension(3)             :: c

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)

        return
end function cross_product

function norm(a) result(norm_a)
        implicit none
        real(rp), dimension(:), intent(in) :: a
        real(rp)                           :: norm_a

        norm_a = sqrt(dot_product(a,a))

        return
end function norm

function arctan2(y,x) result(theta)

        implicit none
        real(rp), intent(in) :: y
        real(rp), intent(in) :: x
        real(rp)             :: theta
        
        real(rp), parameter  :: pi   = 2.0_rp*asin(1.0_rp)

        theta = atan2(y,x)
        if(theta < 0.0_rp) theta = theta + 2*pi

        return
end function arctan2

function cotanh(x) result(fx)
! ------------------------------------
!       Hyperbolic co-tangent
! ------------------------------------
        implicit none
        real(rp), intent(in) :: x
        real(rp)             :: fx

        fx = 1.0_rp/(tanh(x))

        return
end function cotanh

function sech(x) result(fx)
! ------------------------------------
!       Hyperbolic secant
! ------------------------------------
        implicit none
        real(rp), intent(in) :: x
        real(rp)             :: fx

        fx = 1.0_rp/(cosh(x))

        return
end function sech

function csch(x) result(fx)
! ------------------------------------
!       Hyperbolic co-secant
! ------------------------------------
        implicit none
        real(rp), intent(in) :: x
        real(rp)             :: fx

        fx = 1.0_rp/(sinh(x))

        return
end function csch


elemental function DegToRad(deg) result(rad)
        implicit none
        real(rp), intent(in) :: deg
        real(rp), parameter  :: pi = 2.0_rp*asin(1.0_rp)
        real(rp)             :: rad
        
        rad = deg*pi/180.0_rp

        return
end function DegToRad 


elemental function RadToDeg(rad) result(deg)
        implicit none
        real(rp), intent(in) :: rad
        real(rp), parameter  :: pi = 2.0_rp*asin(1.0_rp)
        real(rp)             :: deg
        
        deg = rad*180.0_rp/pi

        return
end function RadToDeg


pure function twoptsline_coeff(p1,p2) result(c)
! -------------------------------------------------
!       a1*x + a2*y + a3 = 0
! -------------------------------------------------
        implicit none
        real(rp), dimension(:), intent(in) :: p1
        real(rp), dimension(:), intent(in) :: p2
        real(rp), dimension(3)             :: c

        c(1) = p1(2) - p2(2)
        c(2) = p2(1) - p1(1)
        c(3) = p1(1)*p2(2) - p2(1)*p1(2)

        return
end function twoptsline_coeff

pure function closest_point_on_line(p0,c) result(x)
! ------------------------------------------------
!       Computing the nearest point x in respect of a point p0
!       liing on a line with c1*x+c2y+c3 = 0
! ------------------------------------------------
        implicit none
        real(rp), dimension(:), intent(in) :: p0
        real(rp), dimension(:), intent(in) :: c
        real(rp), dimension(3)             :: x

        real(rp) :: denom

        denom = 1.0_rp/(c(1)**2 + c(2)**2)

        x(1) = (c(2)*( c(2)*p0(1) - c(1)*p0(2)) - c(1)*c(3))*denom
        x(2) = (c(1)*(-c(2)*p0(1) + c(1)*p0(2)) - c(2)*c(3))*denom
        x(3) = p0(3)

        return
end function closest_point_on_line


function trace(A) result(trac)
        implicit none
        real(rp), dimension(:,:), intent(in) :: A
        real(rp)                             :: trac
        integer                              :: i

        trac = 0.0_rp
        do i = 1, size(A,1)
           trac = trac + A(i,i)
        enddo

        return
end function trace

function det3(A) result(det)
        implicit none
        real(rp), dimension(3,3), intent(in) :: A
        real(rp)                             :: det

        det =   A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))  &
              - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))  &
              + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

        if(abs(det) < 10.0E-14_rp) stop ' detA_3 fails! the metrix is singular'

        return
end function det3

!function eigenvalues(A,n) result(lambda)
!! -------------------------------------------------------------------
!!       This subroutine compute the eigenvalues of a Symmetric matrix 
!! -------------------------------------------------------------------
!        implicit none
!        integer                 , intent(in) :: n
!        real(rp), dimension(n,n), intent(in) :: A
!        real(rp), dimension(n)               :: lambda
!        
!        ! local declarations
!        real(rp), dimension(:), allocatable :: tmp
!        real(rp), dimension(n)              :: swap
!        integer                             :: err = 0, lwork, i
!        external                            :: DSYEV
!
!        lwork = max(1,3*n-1)
!        allocate(tmp(lwork))
!
!        call dsyev('N','U',n,A,n,lambda,TMP,lwork,err)
!
!        ! === swap lambda in ascendent order
!        swap(:) = lambda(:)
!        do i = 1,n
!           lambda(i) = swap(n+1-i)
!        enddo
!
!        if(err /= 0) stop ' FATAL ERROR! LAPACK eigenvalues failed!'
!        
!        deallocate(tmp)
!        return
!end function eigenvalues


subroutine Eigenvalues_Jacobi(M,lambda,toll,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
        implicit none
        integer , intent(in)  :: n
        real(rp), intent(in)  :: M(n,n)
        real(rp), intent(in)  :: toll
        real(rp), intent(out) :: lambda(n)

        ! local declarations
        real(rp) :: i_nr, aji, aik, ajk, aki, akj
        logical  :: mask(n)
        real(rp) :: A(n,n)
        real(rp) :: tmp(n)
        real(rp) :: b2, bar
        real(rp) :: beta, coeff, c, s, cs, sc
        integer  :: i, j, k, l, it, itmax = 10
        
        A = M
        
        ! find the sum of all off-diagonal elements (squared)
        b2 = 0.0_rp
        do j=1,n
          do i=j+1,n
             b2 = b2 + 2*a(i,j)*a(i,j)
          end do
        end do
        
        if (b2 < toll) return

        ! average for off-diagonal elements /2
        i_nr = 1.0_rp/real(n*n,rp)
        bar = 0.5_rp*b2*i_nr
        
        it = 1
        iter:do while (b2.gt.toll .and. it < itmax)

          do i=1,n-1
            do j=i+1,n

              aji = a(j,i)

              if (aji*aji < bar) cycle  ! do not touch small elements
              b2 = b2 - 2.0_rp*aji*aji
              bar = 0.5_rp*b2*i_nr

              ! calculate coefficient c and s for Givens matrix
              beta = (a(j,j)-a(i,i))/(2.0_rp*aji)
              coeff = 0.5_rp*beta/sqrt(1.0_rp+beta*beta)
              s = sqrt(max(0.5_rp+coeff,0.0_rp))
              c = sqrt(max(0.5_rp-coeff,0.0_rp))

              ! recalculate rows i and j
              do k=1,n
                aik = a(i,k)
                ajk = a(j,k)
                cs =  c*aik+s*ajk
                sc = -s*aik+c*ajk
                a(i,k) = cs
                a(j,k) = sc
              end do

              ! new matrix a_{k+1} from a_{k}
              do k=1,n
                aki = a(k,i)
                akj = a(k,j)
                cs =  c*aki+s*akj
                sc = -s*aki+c*akj
                a(k,i) = cs
                a(k,j) = sc
              end do

            end do
          end do

          it = it + 1

        end do iter

        ! sort eigenvalues
        do i = 1,n
           tmp(i) = A(i,i)
        enddo

        mask = .true.
        do i = 1,n
           l = maxloc(tmp,1,mask)
           lambda(i) = tmp(l)
           mask(l)   = .false.
        enddo

        return
end subroutine Eigenvalues_Jacobi




function smooth_step(x0,x1,x) result(Sx)
! --------------------------------------------------
!
!       This function provides a sigmoid smooth step
!       INPUT: x0       !< start position of the step
!              x1       !< end   position of the step
!       OUTPUT:Sx       !> step x
!
! --------------------------------------------------
        implicit none
        real(rp), intent(in) :: x0, x1, x
        real(rp)             :: Sx, x_tmp

        ! rescale x
        x_tmp = (x-x0)/(x1-x0)

        ! computing smooth step
        if(x_tmp < 1.0_rp) then
          Sx = 6*x_tmp**5 - 15*x_tmp**4 + 10*x_tmp**3
        else
          Sx = 1.0_rp
        endif

        return
end function smooth_step


subroutine QuickSort2D(v,vsort)
        implicit none
        real(rp), dimension(:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:)  , allocatable, intent(inout) :: vsort

        ! local declaration
        logical, dimension(:,:), allocatable :: mask
        integer, dimension(2)                :: id
        integer                              :: s1, e1, s2, e2, i, vsize
        
        s1 = lbound(v,1)
        e1 = ubound(v,1)
        s2 = lbound(v,2)
        e2 = ubound(v,2)
        vsize = (e1 - s1 + 1) * (e2 - s2 + 1)

        allocate(mask(s1:e1, s2:e2))
        mask = .true.

        do i = 1,vsize
           vsort(i) = minval(v,mask)
           id       = minloc(abs(v - vsort(i)))
           mask(id(1),id(2)) = .false.
        enddo


        deallocate(mask)
        return
end subroutine QuickSort2D


end module math_tools_module


