!
! Set up start and finish atom indices for each group once in intlbfgs
! STARTGROUP(DOGROUP) and ENDGROUP(DOGROUP)
!

SUBROUTINE QCIOPTPERM(COORDSB,COORDSA,NATOMS,DOGROUP,PERM,STARTGROUP,IM)
USE KEY,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS
IMPLICIT NONE
INTEGER NATOMS
INTEGER PERM(NATOMS), DOGROUP, J1, STARTGROUP(NPERMGROUP), J3, I1, I2, I3, IM
DOUBLE PRECISION COORDSA(3*NATOMS), COORDSB(3*NATOMS), D123, D132, D231, D213, D312, D321, DMIN, D12, D21
LOGICAL IDENTITY

IDENTITY=.TRUE.
DO J1=1,NATOMS
   PERM(J1)=J1
ENDDO

IF (NPERMSIZE(DOGROUP).EQ.3) THEN
   I1=PERMGROUP(STARTGROUP(DOGROUP))
   I2=PERMGROUP(STARTGROUP(DOGROUP)+1)
   I3=PERMGROUP(STARTGROUP(DOGROUP)+2)
! permutation 1 2 3
   D123=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &      COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &      COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I3-1)+3)-&
 &      COORDSA(3*(I3-1)+3))**2
   DMIN=D123
! permutation 1 3 2
   D132=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &      COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &      COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I3-1)+3)-&
 &      COORDSA(3*(I2-1)+3))**2
   IF (D132.LT.DMIN) THEN 
      DMIN=D132
      PERM(I1)=I1
      PERM(I2)=I3
      PERM(I3)=I2
      IDENTITY=.FALSE.
   ENDIF
! permutation 2 1 3
   D213=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &      COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &      COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I3-1)+3)-&
 &      COORDSA(3*(I3-1)+3))**2
   IF (D213.LT.DMIN) THEN 
      DMIN=D213
      PERM(I1)=I2
      PERM(I2)=I1
      PERM(I3)=I3
      IDENTITY=.FALSE.
   ENDIF
! permutation 2 3 1
   D231=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &      COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &      COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I3-1)+3)-&
 &      COORDSA(3*(I1-1)+3))**2
   IF (D231.LT.DMIN) THEN 
      DMIN=D231
      PERM(I1)=I2
      PERM(I2)=I3
      PERM(I3)=I1
      IDENTITY=.FALSE.
   ENDIF
! permutation 3 2 1
   D321=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &      COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &      COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I3-1)+3)-&
 &      COORDSA(3*(I1-1)+3))**2
   IF (D321.LT.DMIN) THEN 
      DMIN=D321
      PERM(I1)=I3
      PERM(I2)=I2
      PERM(I3)=I1
      IDENTITY=.FALSE.
   ENDIF
! permutation 3 1 2
   D312=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &      COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &      COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I3-1)+3)-&
 &      COORDSA(3*(I2-1)+3))**2
   IF (D312.LT.DMIN) THEN 
      DMIN=D312
      PERM(I1)=I3
      PERM(I2)=I1
      PERM(I3)=I2
      IDENTITY=.FALSE.
   ENDIF
!  IF (.NOT.IDENTITY) WRITE(*,'(A,2I6,6G15.5,3I6)') ' qcioptperm> images,D123,D132,D213,D231,D312,D321,permutation: ', &
!&               IM,IM+1,D123,D132,D213,D231,D312,D321,PERM(I1),PERM(I2),PERM(I3)
   WRITE(*,'(A,2I6,6G15.5,3I6)') ' qcioptperm> images,D123,D132,D213,D231,D312,D321,permutation: ', &
 &               IM,IM+1,D123,D132,D213,D231,D312,D321,PERM(I1),PERM(I2),PERM(I3)
ELSE IF (NPERMSIZE(DOGROUP).EQ.2) THEN
! permutation 1 2
   I1=PERMGROUP(STARTGROUP(DOGROUP))
   I2=PERMGROUP(STARTGROUP(DOGROUP)+1)
   D12=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &     COORDSA(3*(I1-1)+3))**2 &
 &    +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &     COORDSA(3*(I2-1)+3))**2 
   DMIN=D12
! permutation 2 1
   D21=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I1-1)+3)-&
 &     COORDSA(3*(I2-1)+3))**2 &
 &    +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I2-1)+3)-&
 &     COORDSA(3*(I1-1)+3))**2 
   IF (D21.LT.DMIN) THEN 
      DMIN=D21
      PERM(I1)=I2
      PERM(I2)=I1
      IF (NSETS(DOGROUP).GT.0) THEN
         DO J3=1,NSETS(DOGROUP) ! can be 2 or 3
!           WRITE(*,'(A,4I6)') 'I1,I2,sets(I1,j3),sets(I2,j3)=',I1,I2,sets(I1,j3),sets(I2,j3)
            PERM(SETS(I1,DOGROUP,J3))=SETS(I2,DOGROUP,J3)
            PERM(SETS(I2,DOGROUP,J3))=SETS(I1,DOGROUP,J3)
         ENDDO
      ENDIF

   ENDIF
   IF (.NOT.IDENTITY) THEN
      WRITE(*,'(A,2I6,2G15.5,2I6)') ' qcioptperm> images,D12,D21,permutation: ',IM,IM+1,D12,D21,PERM(I1),PERM(I2)
      IF (NSETS(DOGROUP).GT.0) WRITE(*,'(A,6I6)') ' qcioptperm> associated permutation: ', &
  &                         (PERM(SETS(I1,DOGROUP,J3)),PERM(SETS(I2,DOGROUP,J3)),J3=1,NSETS(DOGROUP))  
   ENDIF
ELSE
   WRITE(*,'(A,I6,A,I6)') ' qcioptperm> Unknown group size ',NPERMSIZE(DOGROUP),' for group ',DOGROUP
   STOP
ENDIF

END SUBROUTINE QCIOPTPERM