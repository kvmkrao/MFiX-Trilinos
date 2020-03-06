      SUBROUTINE MFIXinterfacemu(VNAME, VNO, VAR, A_M, B_M, ITMAX, TOL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
!      USE param1
!      USE geometry
!      USE indices
      USE compar
!      USE sendrecv
!      USE leqsol
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
       CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
       INTEGER :: VNO
! variable
        DOUBLE PRECISION :: VAR(DIMENSION_3)
! Septadiagonal matrix A_m
        DOUBLE PRECISION:: A_m(DIMENSION_3, -3:3)

! Vector b_m
        DOUBLE PRECISION :: B_m(DIMENSION_3)
! maximum number of iterations (generally leq_it)
        INTEGER :: ITMAX, IER
! error indicator
        INTEGER :: ntrows
	real*8 :: TOL
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! OVERRELAXATION FACTOR
!        DOUBLE PRECISION, PARAMETER :: OMEGA = 1.0  !1.2
!        integer :: iidebug
!        parameter( iidebug = 0 )
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Variable
! Indices
!        INTEGER :: IJK
        INTEGER :: ITER

        DOUBLE PRECISION oAm
!       DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Anew
        DOUBLE PRECISION, allocatable, dimension(:,:) :: Anew,pos
        integer, allocatable, dimension(:) :: locgl


!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
	INTEGER, PARAMETER :: cols = 7
	DOUBLE PRECISION   :: aijmax
!        DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: Bnew

!       added 
!        DOUBLE PRECISION :: AVar(DIMENSION_3)
        INTEGER :: I, J, K, IJK, IJK_GL, ie, dimred, nstart,nend
!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
!        integer :: class, interval
!        integer :: j_start(2), j_end(2)
!       added 
!        allocate(Bnew(DIMENSION_3))
	allocate(Anew(DIMENSION_3,cols))
	allocate(pos(DIMENSION_3,cols-1))
!-----------------------------------------------

        DO IJK = ijkstart3, ijkend3

!         IF(.NOT.IS_ON_myPE_owns(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
	  aijmax = maxval(abs(A_M(ijk,:)) )
          OAM = 1.0d0/aijmax
	  if(aijmax.gt.1e-05) then
!	   A_M(IJK,:)*OAM
!	  OAM = 1.0d0/A_M(IJK,0)
!          A_M(IJK, 0) = 1.0d0
	  A_M(IJK, 0) = A_M(IJK, 0)*OAM
	  A_M(IJK,-3) = A_M(IJK,-3)*OAM
          A_M(IJK,-2) = A_M(IJK,-2)*OAM
          A_M(IJK,-1) = A_M(IJK,-1)*OAM
          A_M(IJK, 1) = A_M(IJK, 1)*OAM
          A_M(IJK, 2) = A_M(IJK, 2)*OAM
	  A_M(IJK, 3) = A_M(IJK, 3)*OAM
          B_M(IJK)    = B_M(IJK)*OAM

!         if(abs(A_M(IJK,-2)).lt.1) 
	  Anew(IJK,1) =  A_M(IJK,-3)
	  Anew(IJK,2) =  A_M(IJK,-1)
!         if(abs(A_M(IJK,-1)).lt.1) 
	  Anew(IJK,3) =  A_M(IJK,-2)
          Anew(IJK,4) =  A_M(IJK,0) ! -1.0d0
!         if(abs(A_M(IJK, 1)).lt.1) 
	  Anew(IJK,5) =  A_M(IJK,2)
!         if(abs(A_M(IJK, 2)).lt.1) 
	  Anew(IJK,6) =  A_M(IJK,1)
	  Anew(IJK,7) =  A_M(IJK,3)
	  locgl(IJK)  = IJK_GL

!	 if(ier.lt.2) then
	         pos(IJK,1) = KM_OF(IJK) - IJK
	         pos(IJK,2) = IM_OF(IJK) - IJK
                 pos(IJK,3) = JM_OF(IJK) - IJK 
                 pos(IJK,4) = JP_OF(IJK) - IJK
                 pos(IJK,5) = IP_OF(IJK) - IJK 
                 pos(IJK,6) = KP_OF(IJK) - IJK 
!	 end if 

!        write(80+myPE,*)IJK,pos(IJK,1),pos(IJK,2),pos(IJK,3),pos(IJK,4), & 
!        pos(IJK,5),pos(IJK,6)
!        write(90+myPE,*)IJK,Anew(IJK,1),Anew(IJK,2),Anew(IJK,3),Anew(IJK,4), & 
!        Anew(IJK,5),Anew(IJK,6),Anew(IJK,7), B_M(IJK)

!	 end if
	 end if
        ENDDO

!       CALL FWrapper(Anew,Bnew,VAR,pos,cols,DIMENSION_3,IER)
!       CALL FWrapper(Anew,B_M,VAR,pos,cols,DIMENSION_3,IER)
        CALL FWrapper_ext(B_M,VAR,pos,Anew,locgl,cols,DIMENSION_3,ijkmax3,ITMAX,TOL)
        END 
