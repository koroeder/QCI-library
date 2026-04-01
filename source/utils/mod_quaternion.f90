!> Module containing quaternion operations
MODULE QUATERNIONS
    USE QCIPREC
    
    IMPLICIT NONE
    
    
    CONTAINS



   SUBROUTINE SLERP_INTERPOLATION(Q0, Q1, ALPHA, Q_INTERP)
    IMPLICIT NONE
    
    REAL(KIND=8), INTENT(IN) :: Q0(4), Q1(4)      ! Quaternions (w, x, y, z)
    REAL(KIND=8), INTENT(IN) :: ALPHA             ! Interpolation parameter [0,1]
    REAL(KIND=8), INTENT(OUT) :: Q_INTERP(4)      ! Interpolated quaternion
    
    REAL(KIND=8) :: DOT_PRODUCT, THETA, SIN_THETA
    REAL(KIND=8) :: Q1_ADJUSTED(4)
    REAL(KIND=8) :: WEIGHT0, WEIGHT1
    
    ! Compute dot product
    DOT_PRODUCT = Q0(1)*Q1(1) + Q0(2)*Q1(2) + Q0(3)*Q1(3) + Q0(4)*Q1(4)
    
    ! Adjust Q1 if dot product is negative (take shorter path)
    Q1_ADJUSTED = Q1
    IF(DOT_PRODUCT < 0.0D0) THEN
        Q1_ADJUSTED = -Q1_ADJUSTED
        DOT_PRODUCT = -DOT_PRODUCT
    END IF
    
    ! Clamp dot product to avoid numerical issues
    DOT_PRODUCT = MAX(-1.0D0, MIN(1.0D0, DOT_PRODUCT))
    
    ! Calculate angle between quaternions
    THETA = ACOS(DOT_PRODUCT)
    SIN_THETA = SIN(THETA)
    
    ! Avoid division by zero
    IF(ABS(SIN_THETA) < 1.0D-6) THEN
        ! Quaternions are very close, use linear interpolation
        Q_INTERP = (1.0D0 - ALPHA) * Q0 + ALPHA * Q1_ADJUSTED
    ELSE
        ! Standard SLERP formula
        WEIGHT0 = SIN((1.0D0 - ALPHA) * THETA) / SIN_THETA
        WEIGHT1 = SIN(ALPHA * THETA) / SIN_THETA
        Q_INTERP = WEIGHT0 * Q0 + WEIGHT1 * Q1_ADJUSTED
    END IF
    
    ! Normalize result
    CALL NORMALIZE_QUATERNION(Q_INTERP)

    END SUBROUTINE SLERP_INTERPOLATION

    SUBROUTINE MATRIX_TO_QUATERNION(ROTATION_MATRIX, Q)
        USE QCIPREC
        IMPLICIT NONE
        REAL(KIND=REAL64), INTENT(IN) :: ROTATION_MATRIX(3, 3)
        REAL(KIND=REAL64), INTENT(OUT) :: Q(4)  ! (w, x, y, z)
        
        REAL(KIND=REAL64) :: TRACE, S
        INTEGER :: I
        
        TRACE = ROTATION_MATRIX(1,1) + ROTATION_MATRIX(2,2) + ROTATION_MATRIX(3,3)
        
        IF(TRACE > 0.0D0) THEN
            ! w is largest
            S = SQRT(TRACE + 1.0D0) * 2.0D0
            Q(1) = 0.25D0 * S
            Q(2) = (ROTATION_MATRIX(3,2) - ROTATION_MATRIX(2,3)) / S
            Q(3) = (ROTATION_MATRIX(1,3) - ROTATION_MATRIX(3,1)) / S
            Q(4) = (ROTATION_MATRIX(2,1) - ROTATION_MATRIX(1,2)) / S
        ELSE IF((ROTATION_MATRIX(1,1) > ROTATION_MATRIX(2,2)) .AND. &
                (ROTATION_MATRIX(1,1) > ROTATION_MATRIX(3,3))) THEN
            ! x is largest
            S = SQRT(1.0D0 + ROTATION_MATRIX(1,1) - ROTATION_MATRIX(2,2) - &
                     ROTATION_MATRIX(3,3)) * 2.0D0
            Q(1) = (ROTATION_MATRIX(3,2) - ROTATION_MATRIX(2,3)) / S
            Q(2) = 0.25D0 * S
            Q(3) = (ROTATION_MATRIX(1,2) + ROTATION_MATRIX(2,1)) / S
            Q(4) = (ROTATION_MATRIX(1,3) + ROTATION_MATRIX(3,1)) / S
        ELSE IF(ROTATION_MATRIX(2,2) > ROTATION_MATRIX(3,3)) THEN
            ! y is largest
            S = SQRT(1.0D0 + ROTATION_MATRIX(2,2) - ROTATION_MATRIX(1,1) - &
                     ROTATION_MATRIX(3,3)) * 2.0D0
            Q(1) = (ROTATION_MATRIX(1,3) - ROTATION_MATRIX(3,1)) / S
            Q(2) = (ROTATION_MATRIX(1,2) + ROTATION_MATRIX(2,1)) / S
            Q(3) = 0.25D0 * S
            Q(4) = (ROTATION_MATRIX(2,3) + ROTATION_MATRIX(3,2)) / S
        ELSE
            ! z is largest
            S = SQRT(1.0D0 + ROTATION_MATRIX(3,3) - ROTATION_MATRIX(1,1) - &
                     ROTATION_MATRIX(2,2)) * 2.0D0
            Q(1) = (ROTATION_MATRIX(2,1) - ROTATION_MATRIX(1,2)) / S
            Q(2) = (ROTATION_MATRIX(1,3) + ROTATION_MATRIX(3,1)) / S
            Q(3) = (ROTATION_MATRIX(2,3) + ROTATION_MATRIX(3,2)) / S
            Q(4) = 0.25D0 * S
        END IF
        
        CALL NORMALIZE_QUATERNION(Q)
    END SUBROUTINE MATRIX_TO_QUATERNION


    SUBROUTINE QUATERNION_TO_MATRIX(Q, ROTATION_MATRIX)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: Q(4)              ! (w, x, y, z)
        REAL(KIND=8), INTENT(OUT) :: ROTATION_MATRIX(3, 3)
        
        REAL(KIND=8) :: W, X, Y, Z
        
        W = Q(1)
        X = Q(2)
        Y = Q(3)
        Z = Q(4)
        
        ROTATION_MATRIX(1,1) = 1.0D0 - 2.0D0*(Y**2 + Z**2)
        ROTATION_MATRIX(1,2) = 2.0D0*(X*Y - W*Z)
        ROTATION_MATRIX(1,3) = 2.0D0*(X*Z + W*Y)
        
        ROTATION_MATRIX(2,1) = 2.0D0*(X*Y + W*Z)
        ROTATION_MATRIX(2,2) = 1.0D0 - 2.0D0*(X**2 + Z**2)
        ROTATION_MATRIX(2,3) = 2.0D0*(Y*Z - W*X)
        
        ROTATION_MATRIX(3,1) = 2.0D0*(X*Z - W*Y)
        ROTATION_MATRIX(3,2) = 2.0D0*(Y*Z + W*X)
        ROTATION_MATRIX(3,3) = 1.0D0 - 2.0D0*(X**2 + Y**2)

    END SUBROUTINE QUATERNION_TO_MATRIX

    SUBROUTINE NORMALIZE_QUATERNION(Q)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(INOUT) :: Q(4)
        REAL(KIND=8) :: NORM
            
        NORM = SQRT(Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2)
        IF(NORM > 1.0D-12) THEN
            Q = Q / NORM
        ELSE
            Q = [1.0D0, 0.0D0, 0.0D0, 0.0D0]  ! Default to identity
        END IF
            
    END SUBROUTINE NORMALIZE_QUATERNION

END MODULE QUATERNIONS