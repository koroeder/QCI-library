!note: keep nsize instead of natoms for helper functions (so we can roatate small groups)

module align
   use qciprec
   implicit none

   contains

      !> Align final image to the start image 
      subroutine align_endpoints()
         use qcikeys, only: natoms, e2e_dist, qciambert, qcihiret
         use INTERPOLATION_KEYS, only: xstart, xfinal
         use lpermdist, only: lopermdist, mindist

         implicit none
         real(kind=real64) :: newxf(3*natoms)
         real(kind=real64) :: dist2, dworst, rmat(3,3)
         integer :: nmove, newperm(natoms)


         ! make sure endpoints are both centred at the origin
         call move_to_origin(natoms, xstart)
         call move_to_origin(natoms, xfinal)

         if (qcihiret) then
            call mindist(xstart,xfinal,rmat,e2e_dist,dworst,newxf)
            xfinal = newxf
         else if (qciambert) then
            call lopermdist(xstart,xfinal,e2e_dist,dist2,rmat,nmove,newperm, .false.)
         end if
      end subroutine align_endpoints

      !> Align first nsize xb coords to xa
      subroutine alignxbtoa(xa, xb, nsize)
         use qcikeys, only: natoms
         use qciprec
         use lpermdist, only: get_rot_matrix, apply_rot
         implicit none 
         real(kind = real64), intent(in) :: xa(3*natoms)    !< reference coords
         real(kind = real64), intent(inout) :: xb(3*natoms) !< coords to be aligned
         integer, intent(in) :: nsize                       !< number of atoms in xa and xb
         real(kind = real64) :: ra(3*nsize), rb(3*nsize), rb_rot(3*nsize)
         real(kind = real64) :: cxa(3), cxb(3), rmat(3,3)
         real(kind = real64) :: dist, r0(3), r1
         integer :: i, j, k
  
         ra(1:3*nsize) = xa(1:3*nsize)
         rb(1:3*nsize) = xb(1:3*nsize)
         
         !centre ra
         call move_to_origin(nsize,ra,cxa)
         !centre rb around origin
         call move_to_origin(nsize,rb,cxb)
         
         !align coordinates
         call get_rot_matrix(nsize, ra, rb, rmat, dist)
         call apply_rot(nsize, rb, rmat, rb_rot )
             
         ! move coordinates to ra frame of reference for maximum alignment
         call move_coords(nsize,rb_rot,-cxa)
         ! write coordinates back to xb to be returned
         xb(1:3*nsize) = rb_rot(1:3*nsize)
      
      end subroutine alignxbtoa


      subroutine move_to_origin(nsize, x, cox_out)
         
         implicit none
         integer, intent(in) :: nsize
         real(kind=real64),intent(inout) :: x(3*nsize)
         real(kind=real64),intent(out), optional :: cox_out(3)
         real(kind=real64) :: cox(3)
         integer :: idx

         cox(1:3) = 0.0d0
         do idx = 1,nsize
            cox(1) = cox(1) + x(3*idx-2)
            cox(2) = cox(2) + x(3*idx-1)            
            cox(3) = cox(3) + x(3*idx)
         end do

         cox(1:3) = cox(1:3)/nsize

         do idx = 1,nsize
            x(3*idx-2) = x(3*idx-2) - cox(1)
            x(3*idx-1) = x(3*idx-1) - cox(2)
            x(3*idx)   = x(3*idx)   - cox(3)
         end do

         if (present(cox_out)) cox_out(1:3) = cox(1:3)

      end subroutine move_to_origin


      !> Translate image by vector -CX
      subroutine move_coords(nsize,x,cx)
         implicit none
         integer, intent(in) :: nsize
         real(kind = real64), intent(inout) :: x(3*nsize)
         real(kind = real64), intent(in) :: cx(3) !< translational vector
         integer :: i, j, idx
         
         do i=1,nsize
            do j=1,3
               idx = 3*(i-1)+j
               x(idx) = x(idx) - cx(j)
            end do
         end do
      end subroutine move_coords 
      
end module align