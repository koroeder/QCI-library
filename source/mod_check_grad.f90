! Module containing the adjusted congrad routines to check the gradient

MODULE MOD_CHECK_GRAD
    USE QCIPREC
    USE CONSTR_E_GRAD
    CONTAINS
 
       SUBROUTINE CHECK_GRAD(XYZ)
          USE QCIKEYS, ONLY: NIMAGES, NATOMS, CHECKCONINT
          REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates

          REAL(KIND = REAL64) :: GGG(3*NATOMS*(NIMAGES+2)), GMINUS(3*NATOMS*(NIMAGES+2)), GPLUS(3*NATOMS*(NIMAGES+2)) ! gradient for each atom in each image
          REAL(KIND = REAL64) :: EEE(NIMAGES+2), EEEMINUS(NIMAGES+2), EEEPLUS(NIMAGES+2)             ! energy for each image
          REAL(KIND = REAL64) :: ETOTAL, EMINUS, EPLUS ! overall energy
          REAL(KIND = REAL64) :: RMS, RMSMINUS, RMSPLUS ! total force
 
          REAL(KIND = REAL64) :: XMINUS(3*NATOMS*(NIMAGES+2)), XPLUS(3*NATOMS*(NIMAGES+2))
          REAL(KIND = REAL64) :: DX=1.0D-2
          REAL(KIND = REAL64) :: DFA, DFN
          INTEGER :: J, K

          XMINUS = XYZ - DX
          XPLUS = XMINUS + 2*DX

          !XMINUS= XYZ
          !XPLUS = XYZ
          !XMINUS(3*NATOMS+1:3*NATOMS*(NIMAGES+1)) = XMINUS(3*NATOMS+1:3*NATOMS*(NIMAGES+1)) - DX
          !XPLUS(3*NATOMS+1:3*NATOMS*(NIMAGES+1)) = XPLUS(3*NATOMS+1:3*NATOMS*(NIMAGES+1)) + DX

          ! call correct congrad routine
          IF (CHECKCONINT) THEN
             CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
             CALL CONGRAD2(EMINUS, XMINUS, GMINUS, EEEMINUS, RMSMINUS)
             CALL CONGRAD2(EPLUS, XPLUS, GPLUS, EEEPLUS, RMSPLUS)
          ELSE
             CALL CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
             CALL CONGRAD1(EMINUS, XMINUS, GMINUS, EEEMINUS, RMSMINUS)
             CALL CONGRAD1(EPLUS, XPLUS, GPLUS, EEEPLUS, RMSPLUS)
          END IF
        
          DFN = (EPLUS-EMINUS)/(2.0D0*DX)
          DFA = RMS
          WRITE(*,*) 'Gradient numerical  = ', DFN
          WRITE(*,*) 'Gradient analytical = ', DFA
          WRITE(*,*) 'Difference in gradients = ', ABS(DFN-DFA)
          WRITE(*,*) 'Total energies: ', ETOTAL, EPLUS, EMINUS
          DO J=1,NIMAGES+2
             DFN = (EEEPLUS(J)-EEEMINUS(J))/(2.0D0*DX)
             DFA = 0.0D0
             DO K=1,3*NATOMS
                DFA = DFA + GGG((3*NATOMS)*(J-1)+K)**2
             END DO
             DFA = SQRT(DFA/(3*NATOMS))
             WRITE(*,*) "Image ", J, " numerical: ", DFN, " analytical: ", DFA, " difference: ", ABS(DFN-DFA)
          END DO
       END SUBROUTINE CHECK_GRAD
END MODULE