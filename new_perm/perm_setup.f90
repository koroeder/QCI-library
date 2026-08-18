module perm_setup
   use prec
   use perm_defs
   implicit none

   
   contains
      subroutine alloc_perm()
         call dealloc_perm()
         allocate(npermsize(natoms))
         allocate(permgroup(natoms))
         allocate(nsets(natoms))
         allocate(sets(natoms,npermgroup,maxnsets))
         allocate(bestperm(natoms))
         allocate(permutable(natoms),permutable2(natoms),ingroup(natoms))
      end subroutine alloc_perm

      subroutine dealloc_perm()
         if (allocated(npermsize)) deallocate(npermsize)
         if (allocated(permgroup)) deallocate(permgroup)
         if (allocated(nsets)) deallocate(nsets)
         if (allocated(sets)) deallocate(sets)   
         if (allocated(bestperm)) deallocate(bestperm)   
         if (allocated(permutable)) deallocate(permutable) 
         if (allocated(permutable2)) deallocate(permutable2) 
         if (allocated(ingroup)) deallocate(ingroup) 
      end subroutine dealloc_perm

      subroutine setup_perm()
         use defs
         use utils_io, only: file_exist, getunit
         use logger, only: log_message
         implicit none
         integer :: funit
         integer :: ndummy, j1, j2, j3

         if (.not.file_exist("perm.allow")) then
            call log_message(3, "In setup_perm, perm.allow file not found - terminating")
            stop
         end if
         
         funit = getunit()
         open(unit=funit,file=permfile,status='old')
         !read first line with number of groups   
         read(funit,*) npermgroup

         call alloc_perm()
         npermsize(:) = 0
         permgroup(:) = 0
         nsets(:) = 0
         sets(:,:,:) = 0
         ! The above dimensions are fixed at NATOMS because:
         ! (a) Atoms were not allowed to appear in more than one group.
         ! (b) The maximum number of pair exchanges associated with a group is three.


         ndummy = 1
         do j1=1,npermgroup
            read(funit,*) npermsize(j1),nsets(j1)
            ! Sanity checks!
            if (nsets(j1).gt.maxnsets) then
               call log_message(3, "perm_Setup> number of secondary sets exceeds limit")
               stop
            end if
            if (ndummy+npermsize(j1)-1.gt.natoms) then              
               call log_message(3, " perm_Setup> number of atoms to be permuted in all groups is > 3*number of atoms")
               stop
            endif
            read(funit,*) permgroup(ndummy:ndummy+npermsize(j1)-1), ((sets(permgroup(j3),j1,j2), &
                              j3=ndummy,ndummy+npermsize(j1)-1), j2=1,nsets(j1))
            ndummy=ndummy+npermsize(j1)
         enddo
         close(funit)

         ! And another sanity check!
         do j1=1,ndummy
            do j2=j1+1,ndummy
               if (permgroup(j2).eq.permgroup(j1)) then
                  call log_message(3, " perm_setup> atom appears in more than one group")
                  stop
               endif
            enddo
         enddo

         permutable(1:natoms) = .false.
         permutable2(1:natoms) = .false.
         ingroup(1:natom) = 0 !is this the correct default initialisation?
         ! Set up some auxiliary data that is needed in lpermdist
         ndummy=1
         do j1=1,npermgroup
            do j2=1,npermsize(j1)
               permutable(permgroup(ndummy+j2-1))=.true.
               ingroup(permgroup(ndummy+j2-1))=j1
            enddo
            ! This block is needed to allow associated permutational atoms to ve recognised as permutable later.
            ! This is not necessary of the associated permutable atons appear in perm.allow in their own list,
            ! as for methy and amine hydrogens, but it is necessary for associated H atoms in phenylalanine,
            ! as these cannot permute without the associated C atoms, etc. The choice of primary and associated
            ! sets for phenylalanine is arbitrary.
            if (nsets(j1).gt.0) then
               do j2=1,npermsize(j1)
                  do j3=1,nsets(j1)
                     permutable2(sets(permgroup(ndummy+j2-1),j1,j3))=.true.  
                  enddo
               enddo
            endif
            ndummy=ndummy+npermsize(j1)
         enddo
      end subroutine setup_perm
end module perm_setup