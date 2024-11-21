MODULE ADDREMOVE_IMAGES
   USE QCIPREC
   IMPLICIT NONE

   CONTAINS
      !TODO are the energy checks needed here?
      ! I am not sure we actually use them before we test the energy anyway - it's not a lot of overhead but a little

      SUBROUTINE ADD_IMAGE(IDXMAX, ETOTAL, RMS)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES, MUPDATE
         USE MOD_INTCOORDS, ONLY: ALLOC_INTCOORDS, ALLOC_PREVCOORDS, XYZ, EEE, GGG
         USE INTERPOLATION_KEYS, ONLY: ALLOC_STEPTAKING, SEARCHSTEP, GDIF
         USE REPULSION, ONLY: CHECKREP
         USE CONSTR_E_GRAD, ONLY: CONGRAD
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IDXMAX
         REAL(KIND=REAL64), INTENT(OUT) :: ETOTAL, RMS
         REAL(KIND=REAL64), ALLOCATABLE :: TEMPXYZ(:), TEMPSTEP(:,:), TEMPGDIF(:,:)
         INTEGER :: OFFSET1, OFFSET2, OFFSET3, J1, NIMAGEDUMMY, NIMGPREV

         !allocate temporary arrays to save current values
         ALLOCATE(TEMPXYZ(3*NATOMS*(NIMAGES+2)))
         ALLOCATE(TEMPSTEP(0:MUPDATE,1:3*NATOMS*NIMAGES))
         ALLOCATE(TEMPGDIF(0:MUPDATE,1:3*NATOMS*NIMAGES))
         !save XYZ coordinates
         TEMPXYZ = XYZ
         !increase number of images by 1
         NIMGPREV = NIMAGES
         NIMAGES = NIMAGES+1
         !reallocate the interpolation coordinates
         !allocate calls deallocate first
         CALL ALLOC_INTCOORDS()
         CALL ALLOC_PREVCOORDS()
         !write xyz back into correct array
         OFFSET1 = 3*NATOMS*(IDXMAX-1)
         OFFSET2 = 3*NATOMS*IDXMAX
         OFFSET3 = 3*NATOMS*(IDXMAX+1)
         XYZ(1:OFFSET2) = TEMPXYZ(1:OFFSET2)
         XYZ(OFFSET2+1:OFFSET3) = (TEMPXYZ(OFFSET1+1:OFFSET2)+TEMPXYZ(OFFSET2+1:OFFSET3))/2.0D0
         XYZ(OFFSET3+1:3*NATOMS*(NIMAGES+2)) = TEMPXYZ(OFFSET2+1:3*NATOMS*(NIMAGES+1))
         !save searchstep and gdif
         TEMPSTEP = SEARCHSTEP
         TEMPGDIF = GDIF
         !reallocate all associated variables now
         CALL ALLOC_STEPTAKING()
         !write search step and gdif back into the correct variables

         DO J1=1,MUPDATE
            !new version of setting up SEARCHSTEP and GDIF
            ! if IDMAX is 1, then we have to add the new image to the start
            IF (IDXMAX.EQ.1) THEN
               SEARCHSTEP(J1,1:3*NATOMS) = TEMPSTEP(J1,1:3*NATOMS)
               SEARCHSTEP(J1,3*NATOMS+1:3*NATOMS*NIMAGES) = TEMPSTEP(J1,1:3*NATOMS*NIMGPREV)
               GDIF(J1,1:3*NATOMS) = TEMPGDIF(J1,1:3*NATOMS)
               GDIF(J1,3*NATOMS+1:3*NATOMS*NIMAGES) = TEMPGDIF(J1,1:3*NATOMS*NIMGPREV)
            ! if IDMAX is bigger equal to NIMAGES, the new image is the final image
            ELSE IF (IDXMAX.EQ.NIMAGES) THEN
               SEARCHSTEP(J1,1:3*NATOMS*NIMGPREV) = TEMPSTEP(J1,1:3*NATOMS*NIMGPREV)
               SEARCHSTEP(J1,3*NATOMS*NIMGPREV+1:3*NATOMS*NIMAGES) = TEMPSTEP(J1,3*NATOMS*(NIMGPREV-1)+1:3*NATOMS*NIMGPREV)
               GDIF(J1,1:3*NATOMS*NIMGPREV) = TEMPGDIF(J1,1:3*NATOMS*NIMGPREV)
               GDIF(J1,3*NATOMS*NIMGPREV+1:3*NATOMS*NIMAGES) = TEMPGDIF(J1,3*NATOMS*(NIMGPREV-1)+1:3*NATOMS*NIMGPREV)
            ELSE
               SEARCHSTEP(J1,1:OFFSET2) = TEMPSTEP(J1,1:OFFSET2)
               SEARCHSTEP(J1,OFFSET2+1:OFFSET3) = TEMPSTEP(J1,OFFSET1+1:OFFSET2)
               SEARCHSTEP(J1,OFFSET3+1:3*NATOMS*NIMAGES) = TEMPSTEP(J1,OFFSET2+1:3*NATOMS*NIMGPREV)
               GDIF(J1,1:OFFSET2) = TEMPGDIF(J1,1:OFFSET2)
               GDIF(J1,OFFSET2+1:OFFSET3) = TEMPGDIF(J1,OFFSET1+1:OFFSET2)
               GDIF(J1,OFFSET3+1:3*NATOMS*NIMAGES) = TEMPGDIF(J1,OFFSET2+1:3*NATOMS*NIMGPREV)
            END IF
         END DO

         ! before we continue check repulsion neighbour list
         CALL CHECKREP(XYZ,0,1)
         ! call congrad routine
         CALL CONGRAD(ETOTAL, XYZ, GGG, EEE, RMS)

         DEALLOCATE(TEMPXYZ,TEMPSTEP,TEMPGDIF)
      END SUBROUTINE ADD_IMAGE


      SUBROUTINE REMOVE_IMAGE(IDXMIN,ETOTAL,RMS)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES, MUPDATE
         USE MOD_INTCOORDS, ONLY: ALLOC_INTCOORDS, ALLOC_PREVCOORDS, XYZ, EEE, GGG
         USE INTERPOLATION_KEYS, ONLY: ALLOC_STEPTAKING, SEARCHSTEP, GDIF
         USE REPULSION, ONLY: CHECKREP
         USE CONSTR_E_GRAD, ONLY: CONGRAD
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IDXMIN
         REAL(KIND=REAL64), INTENT(OUT) :: ETOTAL, RMS
         REAL(KIND=REAL64), ALLOCATABLE :: TEMPXYZ(:), TEMPSTEP(:,:), TEMPGDIF(:,:)
         INTEGER :: OFFSET1, OFFSET2, OFFSET3, J1

         !allocate temporary arrays to save current values
         ALLOCATE(TEMPXYZ(3*NATOMS*(NIMAGES+2)))
         ALLOCATE(TEMPSTEP(0:MUPDATE,1:3*NATOMS*NIMAGES))
         ALLOCATE(TEMPGDIF(0:MUPDATE,1:3*NATOMS*NIMAGES))
         !save XYZ coordinates
         TEMPXYZ = XYZ
         !decrease number of images by 1
         NIMAGES = NIMAGES-1
         !reallocate the interpolation coordinates
         !allocate calls deallocate first
         CALL ALLOC_INTCOORDS()
         CALL ALLOC_PREVCOORDS()
         !write xyz back into correct array
         OFFSET1 = 3*NATOMS*(IDXMIN-1)
         OFFSET2 = 3*NATOMS*IDXMIN
         OFFSET3 = 3*NATOMS*(IDXMIN-2)
         XYZ(1:OFFSET1) = TEMPXYZ(1:OFFSET1)
         XYZ(OFFSET1+1:3*NATOMS*(NIMAGES+2)) = TEMPXYZ(OFFSET2+1:3*NATOMS*(NIMAGES+3))
         !save searchstep and gdif
         TEMPSTEP = SEARCHSTEP
         TEMPGDIF = GDIF
         !reallocate all associated variables now
         CALL ALLOC_STEPTAKING()        
         !write search step and gdif back into the correct variables
         DO J1=1,MUPDATE
            SEARCHSTEP(J1,1:OFFSET3)=TEMPSTEP(J1,1:OFFSET3)
            SEARCHSTEP(J1,OFFSET3+1:3*NATOMS*NIMAGES) = TEMPSTEP(J1,OFFSET1+1:3*NATOMS*(NIMAGES-1))
            
            GDIF(J1,1:OFFSET1)=TEMPGDIF(J1,1:OFFSET1)
            GDIF(J1,OFFSET3+1:3*NATOMS*NIMAGES) = TEMPGDIF(J1,OFFSET1+1:3*NATOMS*(NIMAGES-1))
         END DO

         ! before we continue check repulsion neighbour list
         CALL CHECKREP(XYZ,0,1)
         ! call congrad routine
         CALL CONGRAD(ETOTAL, XYZ, GGG, EEE, RMS)

         DEALLOCATE(TEMPXYZ,TEMPSTEP,TEMPGDIF)

      END SUBROUTINE REMOVE_IMAGE

END MODULE ADDREMOVE_IMAGES