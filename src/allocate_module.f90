module allocate_module
use parameters_module, only: rp

implicit none
private
public :: AllocateReal,   DeallocateReal
public :: AllocateInteger, DeallocateInteger

interface AllocateReal
  module procedure AllocateReal1D
  module procedure AllocateReal2D
  module procedure AllocateReal3D
  module procedure AllocateReal4D
end interface AllocateReal

interface DeallocateReal
  module procedure DeallocateReal1D
  module procedure DeallocateReal2D
  module procedure DeallocateReal3D
  module procedure DeallocateReal4D
end interface DeallocateReal

interface AllocateInteger
  module procedure AllocateInteger_1_1D
  module procedure AllocateInteger_1_2D
  module procedure AllocateInteger_1_3D
  module procedure AllocateInteger_1_4D
  module procedure AllocateInteger_4_1D
  module procedure AllocateInteger_4_2D
  module procedure AllocateInteger_4_3D
  module procedure AllocateInteger_4_4D
end interface AllocateInteger

interface DeallocateInteger
  module procedure DeallocateInteger_1_1D
  module procedure DeallocateInteger_1_2D
  module procedure DeallocateInteger_1_3D
  module procedure DeallocateInteger_1_4D
  module procedure DeallocateInteger_4_1D
  module procedure DeallocateInteger_4_2D
  module procedure DeallocateInteger_4_3D
  module procedure DeallocateInteger_4_4D
end interface DeallocateInteger

contains

subroutine AllocateReal1D(var, s1,e1)
        implicit none
        real(rp), allocatable, dimension(:), intent(inout) :: var
        integer                            , intent(in)    :: s1,e1

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0.0_rp

        return
end subroutine AllocateReal1D

subroutine AllocateReal2D(var, s1,e1, s2,e2)
        implicit none
        real(rp), allocatable, dimension(:,:), intent(inout) :: var
        integer                              , intent(in)    :: s1,e1
        integer                              , intent(in)    :: s2,e2

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0.0_rp

        return
end subroutine AllocateReal2D

subroutine AllocateReal3D(var, s1,e1, s2,e2, s3,e3)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: var
        integer                                , intent(in)    :: s1,e1
        integer                                , intent(in)    :: s2,e2
        integer                                , intent(in)    :: s3,e3

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2, s3:e3), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0.0_rp

        return
end subroutine AllocateReal3D

subroutine AllocateReal4D(var, s1,e1, s2,e2, s3,e3, s4,e4)
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: var
        integer                                  , intent(in)    :: s1,e1
        integer                                  , intent(in)    :: s2,e2
        integer                                  , intent(in)    :: s3,e3
        integer                                  , intent(in)    :: s4,e4

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2, s3:e3, s4:e4), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0.0_rp

        return
end subroutine AllocateReal4D

subroutine DeallocateReal1D(var)
        implicit none
        real(rp), allocatable, dimension(:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateReal1D

subroutine DeallocateReal2D(var)
        implicit none
        real(rp), allocatable, dimension(:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateReal2D

subroutine DeallocateReal3D(var)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateReal3D

subroutine DeallocateReal4D(var)
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateReal4D

subroutine AllocateInteger_1_1D(var, s1,e1)
        implicit none
        integer(1), allocatable, dimension(:), intent(inout) :: var
        integer                              , intent(in)    :: s1,e1

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_1

        return
end subroutine AllocateInteger_1_1D

subroutine AllocateInteger_1_2D(var, s1,e1, s2,e2)
        implicit none
        integer(1), allocatable, dimension(:,:), intent(inout) :: var
        integer                                , intent(in)    :: s1,e1
        integer                                , intent(in)    :: s2,e2

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_1

        return
end subroutine AllocateInteger_1_2D

subroutine AllocateInteger_1_3D(var, s1,e1, s2,e2, s3,e3)
        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: var
        integer                                  , intent(in)    :: s1,e1
        integer                                  , intent(in)    :: s2,e2
        integer                                  , intent(in)    :: s3,e3

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2, s3:e3), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_1

        return
end subroutine AllocateInteger_1_3D

subroutine AllocateInteger_1_4D(var, s1,e1, s2,e2, s3,e3, s4,e4)
        implicit none
        integer(1), allocatable, dimension(:,:,:,:), intent(inout) :: var
        integer                                    , intent(in)    :: s1,e1
        integer                                    , intent(in)    :: s2,e2
        integer                                    , intent(in)    :: s3,e3
        integer                                    , intent(in)    :: s4,e4

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2, s3:e3, s4:e4), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_1

        return
end subroutine AllocateInteger_1_4D

subroutine AllocateInteger_4_1D(var, s1,e1)
        implicit none
        integer(4), allocatable, dimension(:), intent(inout) :: var
        integer                              , intent(in)    :: s1,e1

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_4

        return
end subroutine AllocateInteger_4_1D

subroutine AllocateInteger_4_2D(var, s1,e1, s2,e2)
        implicit none
        integer(4), allocatable, dimension(:,:), intent(inout) :: var
        integer                                , intent(in)    :: s1,e1
        integer                                , intent(in)    :: s2,e2

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_4

        return
end subroutine AllocateInteger_4_2D

subroutine AllocateInteger_4_3D(var, s1,e1, s2,e2, s3,e3)
        implicit none
        integer(4), allocatable, dimension(:,:,:), intent(inout) :: var
        integer                                  , intent(in)    :: s1,e1
        integer                                  , intent(in)    :: s2,e2
        integer                                  , intent(in)    :: s3,e3

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2, s3:e3), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_4

        return
end subroutine AllocateInteger_4_3D

subroutine AllocateInteger_4_4D(var, s1,e1, s2,e2, s3,e3, s4,e4)
        implicit none
        integer(4), allocatable, dimension(:,:,:,:), intent(inout) :: var
        integer                                    , intent(in)    :: s1,e1
        integer                                    , intent(in)    :: s2,e2
        integer                                    , intent(in)    :: s3,e3
        integer                                    , intent(in)    :: s4,e4

        integer :: err = 0

        if(.not.allocated(var)) allocate(var(s1:e1, s2:e2, s3:e3, s4:e4), stat=err)
        if(err /= 0) stop ' Allocation error'

        var = 0_4

        return
end subroutine AllocateInteger_4_4D

subroutine DeallocateInteger_1_1D(var)
        implicit none
        integer(1), allocatable, dimension(:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_1_1D

subroutine DeallocateInteger_1_2D(var)
        implicit none
        integer(1), allocatable, dimension(:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_1_2D

subroutine DeallocateInteger_1_3D(var)
        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_1_3D

subroutine DeallocateInteger_1_4D(var)
        implicit none
        integer(1), allocatable, dimension(:,:,:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_1_4D

subroutine DeallocateInteger_4_1D(var)
        implicit none
        integer(4), allocatable, dimension(:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_4_1D

subroutine DeallocateInteger_4_2D(var)
        implicit none
        integer(4), allocatable, dimension(:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_4_2D

subroutine DeallocateInteger_4_3D(var)
        implicit none
        integer(4), allocatable, dimension(:,:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_4_3D

subroutine DeallocateInteger_4_4D(var)
        implicit none
        integer(4), allocatable, dimension(:,:,:,:), intent(inout) :: var

        integer :: err = 0

        if(allocated(var)) deallocate(var, stat=err)
        if(err /= 0) stop ' Deallocation error'

        return
end subroutine DeallocateInteger_4_4D

end module allocate_module
