module random_module
use parameters_module, only: rp
#ifdef _OPENACC
#ifdef NVIDIA
use curand
use curand_device
use openacc_curand
#endif
#ifdef AMD
use iso_c_binding
use hipfort_rocrand
#endif
#endif

implicit none
private
integer, public, allocatable, dimension(:,:,:)   :: seedField2DYZ
integer, public, allocatable, dimension(:,:,:,:) :: seedField3D

#ifdef _OPENACC
#ifdef NVIDIA
type(curandGenerator), public :: generator
#endif
#ifdef AMD
type(c_ptr),       public :: generator
#endif
#endif

public init_DefaultSeed, rnd,  mpi_random_field2D, random_normal

contains
subroutine init_DefaultSeed
#ifdef _OPENACC
      implicit none
      integer(8), parameter :: seed = 1234
      integer(8)            :: err

#ifdef NVIDIA
      err = curandCreateGenerator(generator, CURAND_RNG_PSEUDO_XORWOW)
      err = curandSetPseudoRandomGeneratorSeed (generator, seed)

      if (err.ne.0) print*,"Error in curandCreateGenerator: ",err
#endif
#ifdef AMD
      err = rocrand_Create_Generator(generator, ROCRAND_RNG_PSEUDO_XORWOW)
      if (err.ne.0) print*,"Error in rocrand_Create_Generator: ",err

      err = rocrand_Set_Seed(generator, seed)
      if (err.ne.0) print*,"Error in rocrand_Set_Seed: ",err
#endif

#else

        implicit none
        integer, allocatable, dimension(:) :: DefaultSeed
        integer                            :: seedSz, err = 0

        
        ! get seed size
        call random_seed(size=seedSz)

        allocate(DefaultSeed(1:seedSz),stat = err)
        if(err .ne. 0) stop ' Allocation error in init_DefaultSeed'

        DefaultSeed(:) = 1

        call random_seed(put=DefaultSeed)

        deallocate(DefaultSeed)
        return

#endif
end subroutine init_DefaultSeed



subroutine mpi_random_field2D(time,ny,nz,periodic,GN,Rf,comm) 
! -------------------------------------------------------------------
!       Computation of a three-dimensional random normal field
! -------------------------------------------------------------------
        use mpi
        
        implicit none
        integer, parameter :: rp = 8

        integer                                , intent(in)    :: ny, nz
        integer                                , intent(in)    :: GN
        real(rp)                               , intent(in)    :: time
        integer                                , intent(in)    :: comm
        logical , dimension(3)                 , intent(in)    :: periodic
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: Rf

        ! local declaration
        integer , dimension(:), allocatable :: seed
        real(rp), dimension(:), allocatable :: lcl_mean, gbl_mean
        real(rp), dimension(:), allocatable :: lcl_rmsq, gbl_rmsq
        integer                             :: lcl_npts, gbl_npts
        integer , dimension(3) :: lb,ub
        integer                :: sy, sz, ey, ez
        real(rp)               :: r
        integer                :: j,k,m,l,iTime, err=0
        integer                :: jj, kk
        !
        ! === Get field dimensions
        !
        lb(:) = (/lbound(Rf,1), lbound(Rf,2), lbound(Rf,3)/)
        ub(:) = (/ubound(Rf,1), ubound(Rf,2), ubound(Rf,3)/)

        sy = lb(2) + GN; ey = ub(2) - GN
        sz = lb(3) + GN; ez = ub(3) - GN

        lcl_npts = (ey-sy+1)*(ez-sz+1)
        
        call random_seed(size=m)
        allocate(seed(1:m))
        
        ! variable to obtain time-decorrelation
        iTime = nint(time)

        !
        ! === Generate a MPI safe random field
        !
        do k    = lb(3), ub(3)
           do j = lb(2), ub(2)
                
              ! exploit y periodicity
              jj = j
              if(periodic(2)) then
                jj = mod(j,ny)
                jj = mod(jj+ny,ny)
              endif

              ! exploit z periodicity
              kk = k
              if(periodic(3)) then
                kk = mod(k,nz)
                kk = mod(kk+nz,nz)
              endif
              
              ! put the seed as a funcion of the glb grid location
              seed(:) = jj + (kk-1)*ny + iTime
              
              call random_seed(put=seed)
              
              do l = lb(1), ub(1)

                 call random_normal(r)
                 Rf(l,j,k) = r

              enddo
          enddo
        enddo

        deallocate(seed)
        !
        ! === Compute the mean
        !
        allocate(lcl_mean(ub(1)-lb(1)+1),gbl_mean(ub(1)-lb(1)+1))

        do l = lb(1), ub(1)
           lcl_mean(l) = sum(Rf(l,sy:ey, sz:ez))
        enddo

        call MPI_allreduce(lcl_mean, gbl_mean, size(lcl_mean), MPI_DOUBLE_PRECISION, MPI_SUM, comm, err)
        call MPI_allreduce(lcl_npts, gbl_npts, 1             , MPI_INTEGER         , MPI_SUM, comm, err)
        
        gbl_mean = gbl_mean/real(gbl_npts,rp)

        !
        ! === Compute the RMSQ
        !
        allocate(lcl_rmsq(ub(1)-lb(1)+1),gbl_rmsq(ub(1)-lb(1)+1))
        
        do l = lb(1),ub(1)
           lcl_rmsq(l) = sum((Rf(l,sy:ey, sz:ez) - gbl_mean(l))**2)
        enddo

        call MPI_allreduce(lcl_rmsq, gbl_rmsq, size(lcl_rmsq), MPI_DOUBLE_PRECISION, MPI_SUM, comm, err)

        gbl_rmsq = sqrt(gbl_rmsq/real(gbl_npts,rp))
        
        !
        ! === subtract the mean and rescale by the RMSQ
        !
        do l = lb(1), ub(1)
           Rf(l,:,:) = (Rf(l,:,:) - gbl_mean(l))/gbl_rmsq(l)
        enddo

        deallocate(lcl_mean, gbl_mean) 
        deallocate(lcl_rmsq, gbl_rmsq) 

        return
end subroutine mpi_random_field2D







subroutine random_normal(ran_norm)
! -----------------------------------------------------------------
!       This subroutine generate random number Gaussian distributed
! -----------------------------------------------------------------

        implicit none
        integer, parameter    :: rp = 8
        real(rp), intent(out) :: ran_norm
        
        ! Local variables
        real(rp), parameter :: s    = 0.449871_rp, t  = -0.386595_rp, a  = 0.19600_rp, b = 0.25472_rp
        real(rp), parameter :: half = 0.5_rp     , r1 = 0.27597_rp  , r2 = 0.27846_rp
        real(rp)            :: u, v, x, y, q
        
        do
          call random_number(u)
          call random_number(v)
          v = 1.7156_rp * (v - half)
        
          ! Evaluate the quadratic form
          x = u - s
          y = ABS(v) - t
          q = x**2 + y*(a*y - b*x)
        
          ! Accept P if inside inner ellipse
          if (q < r1) exit
          ! Reject P if outside outer ellipse
          if (q > r2) cycle
          ! Reject P if outside acceptance region
          if (v**2 < -4.0_rp*LOG(u)*u**2) exit
        enddo
        
        ! Return ratio of P's coordinates as the normal deviate
        ran_norm = v/u

        return
end subroutine random_normal



function rnd() result(r)
! -------------------------------------------------------------------------
!       this function provides a random generator with mean = 0
! -------------------------------------------------------------------------
        use parameters_module, only: rp

        implicit none
        real(rp) :: r

        call random_number(r)
        r = r - 0.5_rp

        return
end function rnd




end module random_module





































