module matrix_inversion_module

implicit none
private
public matinv3, matinv4, gauss_elimination, tdma, inverse

contains

pure function matinv3(A) result(B)
! ----------------------------------------------------------------
! Performs a direct calculation of the inverse of a 3×3 matrix.
! ----------------------------------------------------------------
    integer , parameter  :: rp = 8
    real(rp), intent(in) :: A(3,3)   !! Matrix
    real(rp)             :: B(3,3)   !! Inverse matrix
    real(rp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0_rp/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
        
    return
end function matinv3






pure function matinv4(A) result(B)
! ----------------------------------------------------------------
! Performs a direct calculation of the inverse of a 4×4 matrix.
! ----------------------------------------------------------------
    implicit none
    integer , parameter  :: rp = 8 
    real(rp), intent(in) :: A(4,4)   !! Matrix
    real(rp)             :: B(4,4)   !! Inverse matrix
    real(rp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
    1.0_rp/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))

    return
end function matinv4


subroutine gauss_elimination(a,b,x,n)
!---------------------------------------------------------------
!
!       Optimised solutions to a system of linear equations A*x=b
!
!       INPUT:     
!               A(n,n) !< system matrix
!               b(n)   !< system right hand side
!               n      !< system rank
!       OUTPUT:
!               x(n)   !< system solution
!
!---------------------------------------------------------------
        implicit none
        
        integer , parameter   :: rp = 8
        integer , intent(in ) :: n
        real(rp), intent(in ) :: a(n,n)
        real(rp), intent(in ) :: b(n)
        real(rp), intent(out) :: x(n)
        
        ! local declarations
        real(rp) :: a2(n,n)
        integer  :: i, j, k, p
        real(rp) :: row(n)
        real(rp) :: t
        
        a2(1:n,1:n) = a(1:n,1:n)
        x(1:n) = b(1:n)
        
        do k = 1, n
          !
          !  Find the maximum element in column I.
          !
          p = k
        
          do i = k + 1, n
            if ( abs ( a2(p,k) ) < abs ( a2(i,k) ) ) then
              p = i
            end if
          end do
        
          if (abs(a2(p,k)) < 1.0E-14_rp ) stop ' In Gauss elimination the matrix is singular'
          !
          !  Switch rows K and P.
          !
          if ( k /= p ) then
        
            row(1:n) = a2(k,1:n)
            a2(k,1:n) = a2(p,1:n)
            a2(p,1:n) = row(1:n)
        
            t    = x(k)
            x(k) = x(p)
            x(p) = t
        
          end if
          !
          !  Scale the pivot row.
          !
          a2(k,k+1:n) = a2(k,k+1:n) / a2(k,k)
          x(k) = x(k) / a2(k,k)
          a2(k,k) = 1.0_rp
          !
          !  Use the pivot row to eliminate lower entries in that column.
          !
          do i = k + 1, n
            if (abs(a2(i,k))> 1.0E-14_rp) then
              t = - a2(i,k)
              a2(i,k) = 0.0_rp
              a2(i,k+1:n) = a2(i,k+1:n) + t * a2(k,k+1:n)
              x(i) = x(i) + t * x(k)
            end if
          end do
        
        end do
        !
        !  Back solve.
        !
        do j = n, 2, -1
          x(1:j-1) = x(1:j-1) - a2(1:j-1,j) * x(j)
        end do
        
        return
end subroutine gauss_elimination

subroutine tdma(a,b,c,r,x,s,e)
        !$acc routine seq

        implicit none
        integer , parameter      :: rp = 8
        integer , intent(in)     :: s,e
        real(rp), dimension(s:e) :: a,b,c,r,x    

        integer :: i

        ! forward elimination phase
        do i=s+1,e
                b(i) = b(i) - a(i)/b(i-1)*c(i-1)
                r(i) = r(i) - a(i)/b(i-1)*r(i-1)
        end do
        ! backward substitution phase 
        x(e) = r(e)/b(e)
        do i=e-1,s,-1
                x(i) = (r(i)-c(i)*x(i+1))/b(i)
        end do

        return
end subroutine tdma



 subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer , parameter :: rp = 8
integer  :: n
real(rp) :: a(n,n), c(n,n)
real(rp) :: L(n,n), U(n,n), b(n), d(n), x(n)
real(rp) :: coeff
integer  :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0_rp
U=0.0_rp
b=0.0_rp

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0_rp
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0_rp
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0_rp
end do
end subroutine inverse




!subroutine lapack_gauss_elimination(A,b,x,n)
!
!        implicit none
!        integer, parameter                    :: rp = 8
!        integer                 , intent(in)  :: n
!        real(rp), dimension(n,n), intent(in)  :: A
!        real(rp), dimension(n)  , intent(in)  :: b
!        real(rp), dimension(n)  , intent(out) :: x
!
!        ! local declaration
!        real(rp), dimension(n,n) :: Ainv
!        real(rp), dimension(n)   :: work  ! work array for LAPACK
!        integer , dimension(n)   :: ipiv  ! pivot indices
!        integer , dimension(2)   :: err = 0
!        
!        ! External procedures defined in LAPACK
!        external DGETRF !< LU decomposition
!        external DGETRI !< A^(-1) from LU decomp
!        
!        ! prevent lapack to overwrite
!        Ainv = A
!        
!        call DGETRF(n, n, Ainv, n, ipiv, err(1))
!        call DGETRI(n, Ainv, n, ipiv, work, n, err(2))
!
!        x = matmul(Ainv,b)
!
!        if(sum(err) /= 0) stop ' Error in gauss_elimination'
!        
!        return
!end subroutine lapack_gauss_elimination





end module matrix_inversion_module
