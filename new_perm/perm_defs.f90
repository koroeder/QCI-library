module perm_defs
   use prec
   implicit none
   !> #permutational checks
   !> perm file
   character(len=80) :: permfile = "perm.allow"
   !> energy tolerance for alignment
   real(kind=real64) :: edifftol = 1.0D-6
   !> orbital cutoff for atoms
   real(kind=real64) :: orbittol = 0.3D0
   !> LPERMDIST alignment threshold
   real(kind=real64) :: localpermcut = 0.5D0
   !> alignment cutoff
   real(kind=real64) :: localpermcut2 = 5.0D0
   !> ???
   integer :: localpermneigh = 11
   !> max number of atoms in secondary sets
   integer :: maxnsets = 49
   !> maximum number of random rotations to try in minpermdist
   integer :: nranrot = 0

   !> number of permutational groups 
   integer :: npermgroup = 0

   integer, allocatable :: npermsize(:)
   integer, allocatable :: permgroup(:)
   integer, allocatable :: bestperm(:)
   integer, allocatable :: nsets(:)
   integer, allocatable :: sets(:,:,:)

   logical, allocatable :: permutable(:), permutable2(:)
   integer, allocatable :: ingroup(:)

end module perm_defs