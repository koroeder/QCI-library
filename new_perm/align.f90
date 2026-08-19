module align
   use qciprec
   implicit none

   contains

      subroutine align_endpoints()
         use qcikeys, only: natoms
         use INTERPOLATION_KEYS, only: xstart, xfinal, dist
         implicit none
         real(kind=real64) :: newxf(3*natoms)
         real(kind=real64) :: dist2, dworst, rmat(3,3)
         integer :: nmove, nperm(natoms)

         ! make sure endpoints are both centred at the origin
         call move_to_origin(xstart)
         call move_to_origin(xfinal)

         if (hiret) then
            call mindist(xstart,xfinal,rmat,dist,dworst,newxf)
            xf = newxf
         else if (ambert) then
            call lpermdist(xstart,xfinal,dist,dist2,rmat,nmove,newperm)
         end if
      end subroutine align_endpoints

      subroutine move_to_origin(x)
         use qcikeys, only: natoms
         implicit none
         real(kind=real64) :: x(3*natoms)
         real(kind=real64) :: cox(3)
         integer :: idx

         cox(1:3) = 0.0d0
         do idx = 1,natoms
            cox(1) = cox(1) + x(3*idx-2)
            cox(2) = cox(2) + x(3*idx-1)            
            cox(3) = cox(3) + x(3*idx)
         end do

         cox(1:3) = cox(1:3)/natoms

         do idx = 1,natoms
            x(3*idx-2) = x(3*idx-2) - cox(1)
            x(3*idx-1) = x(3*idx-1) - cox(2)
            x(3*idx)   = x(3*idx)   - cox(3)
         end do
      end subroutine move_to_origin
end module align