! ============================================================================
! test_lpermdist.f90
!
! Minimal functional test for the lpermdist module (subroutine lopermdist).
!
! Setup: a synthetic 5-atom system with one permutable pair (atoms 1,2 --
! think two equivalent hydrogens) and three fixed/distinguishable "anchor"
! atoms (3,4,5). Structure A is built from structure B by:
!   1. swapping the labels of the permutable pair
!   2. applying a known rigid rotation (25 deg about z) + translation
! With no noise added, an exact solver should recover distance ~ 0 and
! nmove == 2 (only the permuted pair should differ from the identity).
!
! This is a black-box numerical check on lopermdist's *output* (distance,
! nmove) -- it deliberately does not assert on the exact contents of
! newperm, since the routine's own header comment ("COORDSA becomes ...
! without the permutations") appears to contradict what the final COORDSA
! assignment actually does (see accompanying documentation, Known Issue).
! distance/nmove are unambiguous regardless of that question.
!
! ============================================================================
program test_lpermdist
   use qcikeys,      only: natoms, debug
   use perm_defs
   use lpermdist, only: lopermdist
   use qciprec
   implicit none


   real(kind=real64), allocatable :: coordsb(:), coordsa(:)
   real(kind=real64) :: distance, dist2
   real(kind=real64) :: rmatbest(3,3)
   integer :: nmove
   integer, allocatable :: newperm(:)

   real(kind=real64) :: rawb(3,5)
   real(kind=real64) :: theta, ct, st, tx, ty, tz
   integer :: i
   logical :: pass

   natoms = 5
   debug  = .true.

   allocate(coordsb(3*natoms), coordsa(3*natoms), newperm(natoms))

   ! --- Reference structure B --------------------------------------------
   rawb(:,1) = (/  1.0d0,  0.6d0,  0.0d0 /)   ! permutable atom "a"
   rawb(:,2) = (/  1.0d0, -0.6d0,  0.0d0 /)   ! permutable atom "b"
   rawb(:,3) = (/  0.0d0,  0.0d0,  0.0d0 /)   ! fixed anchor
   rawb(:,4) = (/ -1.0d0,  0.2d0,  0.4d0 /)   ! fixed anchor
   rawb(:,5) = (/ -1.2d0, -0.8d0, -0.5d0 /)   ! fixed anchor

   do i = 1, natoms
      coordsb(3*(i-1)+1:3*(i-1)+3) = rawb(:,i)
   end do

   ! --- Mobile structure A: swap the permutable pair, then apply a known
   !     rigid rotation (about z) and translation ------------------------
   theta = 25.0d0 * (3.14159265358979d0 / 180.0d0)
   !theta = 0.0d0
   ct = cos(theta); st = sin(theta)
   tx = 0.3d0; ty = -0.2d0; tz = 0.5d0
   !tx = 0.0d0; ty = 0.0d0; tz = 0.0d0

   call place(1, rawb(:,2))   ! slot 1 gets B's atom-2 coordinates (swapped)
   call place(2, rawb(:,1))   ! slot 2 gets B's atom-1 coordinates (swapped)
   call place(3, rawb(:,3))
   call place(4, rawb(:,4))
   call place(5, rawb(:,5))

   ! --- Permutation-group setup: one group {1,2}, no obligatory sets -----
   npermgroup = 1
   allocate(npermsize(npermgroup)); npermsize(1) = 2
   allocate(permgroup(2));          permgroup    = (/ 1, 2 /)
   allocate(nsets(npermgroup));     nsets(1)     = 0
   allocate(sets(1,1,1));           sets(1,1,1)  = 0
   allocate(permutable(natoms));    permutable   = .false.
   allocate(permutable2(natoms));   permutable2  = .false.
   permutable(1) = .true.
   permutable(2) = .true.
   allocate(ingroup(natoms));       ingroup      = 0
   ingroup(1) = 1
   ingroup(2) = 1

   localpermneigh = 3        ! force nats < natoms -> exercises standard_orient
   localpermcut   = 1.0d0    ! generous accept threshold for this noiseless test
   localpermcut2  = 3.0d0    ! generous neighbour-inclusion distance cutoff

   write(*,*) "Input coords"
   write(*,*) '--- coordsb  ---'
   do i = 1, natoms
      write(*,'(A,I0,A,3F12.6)') 'B(', i, ') = ', coordsb(3*(i-1)+1:3*(i-1)+3)
   end do
 
   write(*,*) '--- coordsa  ---'
   do i = 1, natoms
      write(*,'(A,I0,A,3F12.6)') 'A(', i, ') = ', coordsa(3*(i-1)+1:3*(i-1)+3)
   end do


   ! --- Run --------------------------------------------------------------
   call lopermdist(coordsb, coordsa, distance, dist2, rmatbest, nmove, newperm)

   write(*,*) '============================================================'
   write(*,*) 'newperm  = ', newperm
   write(*,*) 'nmove    = ', nmove
   write(*,*) 'distance = ', distance
   write(*,*) 'dist2    = ', dist2

   pass = (distance < 1.0d-6) .and. (nmove == 2)

   if (pass) then
      write(*,*) 'RESULT: PASS -- permutation + rigid alignment recovered exactly'
   else
      write(*,*) 'RESULT: FAIL -- expected distance ~ 0 and nmove == 2, see values above'
   end if
   write(*,*) '============================================================'
   
   write(*,*) '--- coordsb (reference, unchanged) ---'
   do i = 1, natoms
      write(*,'(A,I0,A,3F12.6)') 'B(', i, ') = ', coordsb(3*(i-1)+1:3*(i-1)+3)
   end do
 
   write(*,*) '--- coordsa (returned, post-alignment) ---'
   do i = 1, natoms
      write(*,'(A,I0,A,3F12.6)') 'A(', i, ') = ', coordsa(3*(i-1)+1:3*(i-1)+3)
   end do

   
contains

   subroutine place(islot, bvec)
      integer, intent(in) :: islot
      real(kind=real64), intent(in) :: bvec(3)
      real(kind=real64) :: v(3)
      v(1) = ct*bvec(1) - st*bvec(2) + tx
      v(2) = st*bvec(1) + ct*bvec(2) + ty
      v(3) = bvec(3) + tz
      coordsa(3*(islot-1)+1 : 3*(islot-1)+3) = v
   end subroutine place

end program test_lpermdist