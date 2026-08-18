!     Interface to spjv.f for calculating minimum distance
!     of two atomic configurations with respect to
!     particle permutations.
!     The function permdist determines the distance or weight function,
!     minperm is the main routine.
!
!         Tomas Oppelstrup, Jul 10, 2003
!         tomaso@nada.kth.se
!

!     This is the main routine for minimum distance calculation.
!     Given two coordinate vectors p,q of particles each, return
!     the minimum distance in dist, and the permutation in perm.
!     perm is an integer vector such that
!       p(i) <--> q(perm(i))
!     i.e.
!       sum(i=1,n) permdist(p(i), q(perm(i))) == dist

module minperm_mod
   use prec, only: int64, real64
   implicit none
   ! Save the largest arrays between iterations to reduce allocations.
   !     cc, kk:  Sparse matrix of distances
   integer(kind=int64), allocatable, dimension(:) :: cc, kk

   private
   public :: minperm, minperm_deallocate

   contains
      subroutine minperm(n, p, q, perm, dist, worstdist, worstradius)
         implicit none

         !     Input
         !       n  : System size
         !       p,q: Coordinate vectors (n particles)
         !       pbc: Periodic boundary conditions?
         integer, intent(in) :: n
         real(kind=real64) :: p(3*n), q(3*n), worstdist, worstradius
         !remove--> sx, sy,sz, s
         !remove--> logical pbc

         !     Output
         !       perm: Permutation so that p(i) <--> q(perm(i))
         !       dist: Minimum attainable distance
         integer, intent(out) :: perm(n)
         real(kind=real64),intent(out) :: dist
         real(kind=real64) :: dummy
         !     Parameters
         !       scale : Precision
         !       maxnei: Maximum number of closest neighbours
         real(kind=real64), parameter :: scale = 1.0d6
         integer, parameter :: maxnei = 60

         !     Internal variables
         !     first:
         !       Sparse matrix of distances
         !     first(i):
         !       Beginning of row i in data,index vectors
         !     kk(first(i)..first(i+1)-1):
         !       Column indexes of existing elements in row i
         !     cc(first(i)..first(i+1)-1):
         !       Matrix elements of row i
         integer(kind=int64) :: first(n+1), x(n), y(n)
         integer(kind=int64) :: u(n), v(n), h
         integer(kind=int64) :: m, i, j, k, l, l2, t, a, i3, j3
         integer(kind=int64) :: n8, sz8, d
         integer(kind=int64) :: ndone, j1, j2

         if (.not. allocated(kk)) allocate(kk(n*maxnei))
         if (size(kk) .ne. n*maxnei) then
            deallocate(kk)
            allocate(kk(n*maxnei))
         end if
         if (.not. allocated(cc)) allocate(cc(n*maxnei))
         if (size(cc) .ne. n*maxnei) then
            deallocate(cc)
            allocate(cc(n*maxnei))
         end if

         m = maxnei
         if(n .le. maxnei) m = n
         sz8 = m*n
         n8 = n

         do i=0,n
            first(i+1) = i*m + 1
         enddo

         if(m .eq. n) then
   !     Compute the full matrix...

            do i=1,n
               k = first(i)-1
               do j=1,n
                  cc(k+j) = permdist(p(3*i-2), q(3*j-2), s, pbc)*scale
                  kk(k+j) = j
               enddo
            enddo
         else
         !     We need to store the distances of the maxnei closeest neighbors
         !     of each particle. The following builds a heap to keep track of
         !     the maxnei closest neighbours seen so far. It might be more
         !     efficient to use quick-select instead... (This is definately
         !     true in the limit of infinite systems.)
            do i=1,n
               k = first(i)-1
               do j=1,m
                  cc(k+j) = permdist(p(3*i-2), q(3*j-2), s, pbc)*scale
                  kk(k+j) = j
                  l = j
10                if(l .le. 1) goto 11
                  l2 = l/2
                  if(cc(k+l2) .lt. cc(k+l)) then
                     h = cc(k+l2)
                     cc(k+l2) = cc(k+l)
                     cc(k+l) = h
                     t = kk(k+l2)
                     kk(k+l2) = kk(k+l)
                     kk(k+l) = t
                     l = l2
                     goto 10
                  endif
11             enddo
            
               do j=m+1,n
                  d = permdist(p(3*i-2), q(3*j-2), s, pbc)*scale
                  if(d .lt. cc(k+1)) then
                     cc(k+1) = d
                     kk(k+1) = j
                     l = 1
20                   l2 = 2*l
                     if(l2+1 .gt. m) goto 21
                     if(cc(k+l2+1) .gt. cc(k+l2)) then
                        a = k+l2+1
                     else
                        a = k+l2
                     endif
                     if(cc(a) .gt. cc(k+l)) then
                        h = cc(a)
                        cc(a) = cc(k+l)
                        cc(k+l) = h
                        t = kk(a)
                        kk(a) = kk(k+l)
                        kk(k+l) = t
                        l = a-k
                        goto 20
                     endif
21                   if (l2 .le. m) THEN ! split IF statements to avoid a segmentation fault
                        if (cc(k+l2) .gt. cc(k+l)) then
                           h = cc(k+l2)
                           cc(k+l2) = cc(k+l)
                           cc(k+l) = h
                           t = kk(k+l2)
                           kk(k+l2) = kk(k+l)
                           kk(k+l) = t
                        endif
                     endif
                  endif
               enddo
            enddo    
         end if 

         !     Call bipartite matching routine
         call jovosap(n8, sz8, cc, kk, first, x, y, u, v, h)

         if(h .lt. 0) then
            !     If initial guess correct, deduce solution distance
            !     which is not done in jovosap
            h = 0
            do i=1,n
               j = first(i)
30             if (j.gt.n*maxnei) then
   !              PRINT '(A,I6,A)','minperm> WARNING A - matching failed'
                  do J1=1,n
                     perm(J1)=J1
                  enddo
                  return
               endif
               if(kk(j) .ne. x(i)) then
                  j = j + 1
                  goto 30
               endif
               h = h + cc(j)
            enddo
         endif

         do i=1,n
            perm(i) = x(i)
            if (perm(i).gt.n) perm(i)=n
            if (perm(i).lt.1) perm(i)=1
         enddo

         dist = dble(h) / scale

         worstdist=-1.0d0
         do i=1,n
         ! DUMMY=(p(3*(i-1)+1)-q(3*(perm(i)-1)+1))**2+(p(3*(i-1)+2)-q(3*(perm(i)-1)+2))**2+(p(3*(i-1)+3)-q(3*(perm(i)-1)+3))**2
            dummy=permdist(p(3*i-2),q(3*(perm(i)-1)+1),s,pbc)
            if (dummy.gt.worstdist) then
               worstdist=dummy 
               worstradius=p(3*(i-1)+1)**2+p(3*(i-1)+2)**2+p(3*(i-1)+3)**2
            endif
         enddo
         worstdist=sqrt(worstdist)
         worstradius=max(sqrt(worstradius),1.0d0)
      end subroutine minperm 
      
      !     permdist is the distance or weight function. It is coded
      !     separately for clarity. Just hope that the compiler
      !     knows how to to do proper inlining!
      !     Input
      !       p,q: Coordinates
      real(kind=real64) function permdist(p, q)
         implicit none
         real(kind=real64) , intent(in) :: p(3), q(3)
         real(kind=real64)  :: d

         d = (q(1) - p(1))**2+(q(2) - p(2))**2+(q(3) - p(3))**2
         permdist = d
      end function permdist

!     The following routine performs weighted bipartite matching for
!     for a sparse non-negative integer weight matrix.
!     The original source is
!         http://www.magiclogic.com/assignment.html
!     A publication reference can be found on the above homepage and
!     in a comment below
!     

      SUBROUTINE JOVOSAP(N,SZ,CC,KK,FIRST,X,Y,U,V,H)
      IMPLICIT NONE
      INTEGER*8 N, SZ
      INTEGER*8 CC(SZ),KK(SZ),FIRST(N+1),X(N),Y(N),U(N),V(N)
      INTEGER*8 H,CNT,L0,T,T0,TD,V0,VJ,DJ
      INTEGER*8 LAB(N),D(N),FREE(N),TODO(N)
      LOGICAL OK(N)
      INTEGER*8 J, I, J0, L, J1, MIN, K, I0
      INTEGER*8 BIGINT

!     I don't know how to make g77 read integer*8 constants/parameters.
!       PARAMETER (BIGINT = 10**12) does not work(!)
!     nor does
!       PARAMETER (BIGINT = 1000000000000)
!     but this seems to be ok:
      BIGINT = 10**9
      BIGINT = BIGINT * 1000

!
! THIS SUBROUTINE SOLVES THE SPARSE LINEAR ASSIGNMENT PROBLEM
! ACCORDING 
!
!   "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear   
!    Assignment Problems," Computing 38, 325-340, 1987
!   
!   by
!   
!   R. Jonker and A. Volgenant, University of Amsterdam.
!
!
! INPUT PARAMETERS :
! N = NUMBER OF ROWS AND COLUMNS
! C = WEIGHT MATRIX
!
! OUTPUT PARAMETERS
! X = COL ASSIGNED TO ROW
! Y = ROW ASSIGNED TO COL
! U = DUAL ROW VARIABLE
! V = DUAL COLUMN VARIABLE
! H = VALUE OF OPTIMAL SOLUTION
!
! INITIALIZATION

!     Next line added by tomaso@nada.kth.se, to enable detection
!     of solutions being equivalent to the initial guess

!
!  If Y(:) is initialised to zero then we see segmentation faults if 
!  a Y element is unset, etc.
!

      Y(1:N) = 0
      X(1:N) = 0
      TODO(1:N)=0
      h = -1
      J1 = 0
      DO 10 J=1,N
         V(J)=BIGINT
   10 CONTINUE
      DO 20 I=1,N
         X(I)=0
         DO 15 T=FIRST(I),FIRST(I+1)-1
            J=KK(T)
            IF (CC(T).LT.V(J)) THEN
              V(J)=CC(T)
              Y(J)=I
            END IF
   15    CONTINUE
   20 CONTINUE
      DO 30 J=1,N
         J0=N-J+1
         I=Y(J0)
         IF (I.EQ.0) THEN
!           PRINT '(A,I6,A)','minperm> WARNING B - matching failed'
            RETURN
         ENDIF
         IF (X(I).NE.0) THEN
           X(I)=-ABS(X(I))
           Y(J0)=0
         ELSE
           X(I)=J0
         END IF
   30 CONTINUE
      L=0
      DO 40 I=1,N
         IF (X(I).EQ.0) THEN
           L=L+1
           FREE(L)=I
           GOTO 40
         END IF
         IF (X(I).LT.0) THEN
           X(I)=-X(I)
         ELSE
           J1=X(I)
           MIN=BIGINT
           DO 31 T=FIRST(I),FIRST(I+1)-1
              J=KK(T)
              IF (J.EQ.J1) GOTO 31
              IF (CC(T)-V(J).LT.MIN) MIN=CC(T)-V(J)
   31      CONTINUE
           V(J1)=V(J1)-MIN
         END IF
   40 CONTINUE
! IMPROVE THE INITIAL SOLUTION
      CNT=0
      IF (L.EQ.0) GOTO 1000
   41 L0=L
      K=1
      L=0
   50 I=FREE(K)
      K=K+1
      V0=BIGINT
      VJ=BIGINT
      DO 42 T=FIRST(I),FIRST(I+1)-1
         J=KK(T)
         H=CC(T)-V(J)
         IF (H.LT.VJ) THEN
           IF (H.GE.V0) THEN
             VJ=H
             J1=J
           ELSE
             VJ=V0
             V0=H
             J1=J0
             J0=J
           END IF
         END IF
   42 CONTINUE
      I0=Y(J0)
      IF (V0.LT.VJ) THEN
        V(J0)=V(J0)-VJ+V0
      ELSE
        IF (I0.EQ.0) GOTO 43
        J0=J1
        I0=Y(J1)
      END IF
      IF (I0.EQ.0) GOTO 43
      IF (V0.LT.VJ) THEN
        K=K-1
        FREE(K)=I0
      ELSE
        L=L+1
        FREE(L)=I0
      END IF
   43 X(I)=J0
      Y(J0)=I
      IF (K.LE.L0) GOTO 50
      CNT=CNT+1
      IF ((L.GT.0).AND.(CNT.LT.2)) GOTO 41
! AUGMENTATION PART
      L0=L
      DO 90 L=1,L0
         DO 51 J=1,N
            OK(J)=.FALSE.
            D(J)=BIGINT
   51    CONTINUE
         MIN=BIGINT
         I0=FREE(L)
         TD=N
         DO 52 T=FIRST(I0),FIRST(I0+1)-1
            J=KK(T)
            DJ=CC(T)-V(J)
            D(J)=DJ
            LAB(J)=I0
            IF (DJ.LE.MIN) THEN
              IF (DJ.LT.MIN) THEN
                MIN=DJ
                K=1
                TODO(1)=J
              ELSE
                K=K+1
                TODO(K)=J
              END IF
            END IF
   52    CONTINUE
         DO 53 H=1,K
            J=TODO(H)
            IF (J.EQ.0) THEN
!              PRINT '(A,I6,A)','minperm> WARNING C - matching failed'
               RETURN
            ENDIF
            IF (Y(J).EQ.0) GOTO 80
            OK(J)=.TRUE.
   53    CONTINUE
! REPEAT UNTIL A FREE ROW HAS BEEN FOUND
   60    IF (K.EQ.0) THEN
!           PRINT '(A,I6,A)','minperm> WARNING D - matching failed'
            RETURN
         ENDIF
         J0=TODO(K)
         K=K-1
         I=Y(J0)
         TODO(TD)=J0
         TD=TD-1
         T0=FIRST(I)
         T=T0-1
   61    T=T+1
         IF (KK(T).NE.J0) GOTO 61
         H=CC(T)-V(J0)-MIN
         DO 62 T=T0,FIRST(I+1)-1
            J=KK(T)
            IF (.NOT. OK(J)) THEN
              VJ=CC(T)-H-V(J)
              IF (VJ.LT.D(J)) THEN
                D(J)=VJ
                LAB(J)=I
                IF (VJ.EQ.MIN) THEN
                  IF (Y(J).EQ.0) GOTO 70
                  K=K+1
                  TODO(K)=J
                  OK(J)=.TRUE.
                END IF
              END IF
            END IF
   62    CONTINUE
         IF (K.NE.0) GOTO 60
         MIN=BIGINT-1
         DO 63 J=1,N
            IF (D(J).LE.MIN) THEN
              IF (.NOT. OK(J)) THEN
                IF (D(J).LT.MIN) THEN
                  MIN=D(J)
                  K=1
                  TODO(1)=J
                ELSE
                  K=K+1
                  TODO(K)=J
                END IF
              END IF
            END IF
   63    CONTINUE
         DO 64 J0=1,K
            J=TODO(J0)
            IF (Y(J).EQ.0) GOTO 70
            OK(J)=.TRUE.
   64    CONTINUE
         GOTO 60
   70    IF (MIN.EQ.0) GOTO 80
         DO 71 K=TD+1,N
            J0=TODO(K)
            V(J0)=V(J0)+D(J0)-MIN
   71    CONTINUE
   80    I=LAB(J)
         Y(J)=I
         K=J
         J=X(I)
         X(I)=K
         IF (I0.NE.I) GOTO 80
   90 CONTINUE
      H=0
      DO 100 I=1,N
         J=X(I)
         T=FIRST(I)-1
  101    T=T+1
         IF (T.GT.SZ) THEN
            PRINT '(A,I6,A)','minperm> WARNING D - atom ',I,' not matched - maximum number of neighbours too small?'
            RETURN
         ENDIF
         IF (KK(T).NE.J) GOTO 101
         DJ=CC(T)
         U(I)=DJ-V(J)
         H=H+DJ
  100 CONTINUE

 1000 END SUBROUTINE JOVOSAP

      subroutine minperm_deallocate()

          implicit none

          if (allocated(kk)) deallocate(kk)
          if (allocated(cc)) deallocate(cc)

      end subroutine minperm_deallocate
end module minperm_mod
