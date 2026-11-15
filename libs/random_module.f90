module random_module
use parameters_module, only: rp, pi
#ifdef _OPENACC
#ifdef NVIDIA
use curand
use curand_device
use openacc_curand
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
#endif

public init_DefaultSeed, rnd,  mpi_random_field2D, random_normal
public philox_normal4

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

#else

        implicit none
        integer, allocatable, dimension(:) :: DefaultSeed
        integer                            :: seedSz, err = 0

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
        use mpi

        implicit none
        integer, parameter :: rp = 8

        integer                                , intent(in)    :: ny, nz
        integer                                , intent(in)    :: GN
        real(rp)                               , intent(in)    :: time
        integer                                , intent(in)    :: comm
        logical , dimension(3)                 , intent(in)    :: periodic
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: Rf

        integer , dimension(:), allocatable :: seed
        real(rp), dimension(:), allocatable :: lcl_mean, gbl_mean
        real(rp), dimension(:), allocatable :: lcl_rmsq, gbl_rmsq
        integer                             :: lcl_npts, gbl_npts
        integer , dimension(3) :: lb,ub
        integer                :: sy, sz, ey, ez
        real(rp)               :: r
        integer                :: j,k,m,l,iTime, err=0
        integer                :: jj, kk
        lb(:) = (/lbound(Rf,1), lbound(Rf,2), lbound(Rf,3)/)
        ub(:) = (/ubound(Rf,1), ubound(Rf,2), ubound(Rf,3)/)

        sy = lb(2) + GN; ey = ub(2) - GN
        sz = lb(3) + GN; ez = ub(3) - GN

        lcl_npts = (ey-sy+1)*(ez-sz+1)

        call random_seed(size=m)
        allocate(seed(1:m))

        iTime = nint(time)

        do k    = lb(3), ub(3)
           do j = lb(2), ub(2)

              jj = j
              if(periodic(2)) then
                jj = mod(j,ny)
                jj = mod(jj+ny,ny)
              endif

              kk = k
              if(periodic(3)) then
                kk = mod(k,nz)
                kk = mod(kk+nz,nz)
              endif

              seed(:) = jj + (kk-1)*ny + iTime

              call random_seed(put=seed)

              do l = lb(1), ub(1)

                 call random_normal(r)
                 Rf(l,j,k) = r

              enddo
          enddo
        enddo

        deallocate(seed)
        allocate(lcl_mean(ub(1)-lb(1)+1),gbl_mean(ub(1)-lb(1)+1))

        do l = lb(1), ub(1)
           lcl_mean(l) = sum(Rf(l,sy:ey, sz:ez))
        enddo

        call MPI_allreduce(lcl_mean, gbl_mean, size(lcl_mean), MPI_DOUBLE_PRECISION, MPI_SUM, comm, err)
        call MPI_allreduce(lcl_npts, gbl_npts, 1             , MPI_INTEGER         , MPI_SUM, comm, err)

        gbl_mean = gbl_mean/real(gbl_npts,rp)

        allocate(lcl_rmsq(ub(1)-lb(1)+1),gbl_rmsq(ub(1)-lb(1)+1))

        do l = lb(1),ub(1)
           lcl_rmsq(l) = sum((Rf(l,sy:ey, sz:ez) - gbl_mean(l))**2)
        enddo

        call MPI_allreduce(lcl_rmsq, gbl_rmsq, size(lcl_rmsq), MPI_DOUBLE_PRECISION, MPI_SUM, comm, err)

        gbl_rmsq = sqrt(gbl_rmsq/real(gbl_npts,rp))

        do l = lb(1), ub(1)
           Rf(l,:,:) = (Rf(l,:,:) - gbl_mean(l))/gbl_rmsq(l)
        enddo

        deallocate(lcl_mean, gbl_mean)
        deallocate(lcl_rmsq, gbl_rmsq)

        return
end subroutine mpi_random_field2D

subroutine random_normal(ran_norm)

        implicit none
        integer, parameter    :: rp = 8
        real(rp), intent(out) :: ran_norm

        real(rp), parameter :: s    = 0.449871_rp, t  = -0.386595_rp, a  = 0.19600_rp, b = 0.25472_rp
        real(rp), parameter :: half = 0.5_rp     , r1 = 0.27597_rp  , r2 = 0.27846_rp
        real(rp)            :: u, v, x, y, q

        do
          call random_number(u)
          call random_number(v)
          v = 1.7156_rp * (v - half)

          x = u - s
          y = ABS(v) - t
          q = x**2 + y*(a*y - b*x)

          if (q < r1) exit
          if (q > r2) cycle
          if (v**2 < -4.0_rp*LOG(u)*u**2) exit
        enddo

        ran_norm = v/u

        return
end subroutine random_normal

function rnd() result(r)
        use parameters_module, only: rp

        implicit none
        real(rp) :: r

        call random_number(r)
        r = r - 0.5_rp

        return
end function rnd

subroutine philox_normal4(c1,c2,c3,c4, z1,z2,z3,z4)
!$acc routine seq
        implicit none
        integer(8), intent(in)  :: c1,c2,c3,c4
        real(rp)  , intent(out) :: z1,z2,z3,z4

        integer(8), parameter :: MASK = 4294967295_8
        integer(8), parameter :: M0   = 3528531795_8
        integer(8), parameter :: M1   = 3449720151_8
        integer(8), parameter :: W0   = 2654435769_8
        integer(8), parameter :: W1   = 3144134277_8
        integer(8), parameter :: KY0  = 1234_8
        integer(8), parameter :: KY1  = 0_8

        real(rp), parameter :: I32  = 2.3283064365386963e-10_rp
        real(rp), parameter :: HI32 = 1.1641532182693481e-10_rp

        integer(8) :: b0,b1,b2,b3, k0,k1
        integer(8) :: alo,ahi,blo,bhi,t0,t1,t2,t3,mid
        integer(8) :: hi0,lo0,hi1,lo1
        real(rp)   :: u1,u2,u3,u4, r12, r34, a12, a34
        integer    :: n

        b0 = iand(c1,MASK); b1 = iand(c2,MASK)
        b2 = iand(c3,MASK); b3 = iand(c4,MASK)
        k0 = KY0          ; k1 = KY1

        do n = 1,10

           alo = iand(M0,65535_8); ahi = ishft(M0,-16)
           blo = iand(b0,65535_8); bhi = ishft(b0,-16)
           t0  = alo*blo; t1 = ahi*blo; t2 = alo*bhi; t3 = ahi*bhi
           mid = t1 + t2 + ishft(t0,-16)
           lo0 = ior(ishft(iand(mid,65535_8),16), iand(t0,65535_8))
           hi0 = iand(t3 + ishft(mid,-16), MASK)

           alo = iand(M1,65535_8); ahi = ishft(M1,-16)
           blo = iand(b2,65535_8); bhi = ishft(b2,-16)
           t0  = alo*blo; t1 = ahi*blo; t2 = alo*bhi; t3 = ahi*bhi
           mid = t1 + t2 + ishft(t0,-16)
           lo1 = ior(ishft(iand(mid,65535_8),16), iand(t0,65535_8))
           hi1 = iand(t3 + ishft(mid,-16), MASK)

           b0 = ieor(ieor(hi1,b1),k0)
           b1 = lo1
           b2 = ieor(ieor(hi0,b3),k1)
           b3 = lo0

           k0 = iand(k0 + W0, MASK)
           k1 = iand(k1 + W1, MASK)

        enddo

        u1 = real(b0,rp)*I32 + HI32
        u2 = real(b1,rp)*I32 + HI32
        u3 = real(b2,rp)*I32 + HI32
        u4 = real(b3,rp)*I32 + HI32

        r12 = sqrt(-2.0_rp*log(u1)); a12 = 2.0_rp*pi*u2
        r34 = sqrt(-2.0_rp*log(u3)); a34 = 2.0_rp*pi*u4

        z1 = r12*cos(a12)
        z2 = r12*sin(a12)
        z3 = r34*cos(a34)
        z4 = r34*sin(a34)

        return
end subroutine philox_normal4

end module random_module
