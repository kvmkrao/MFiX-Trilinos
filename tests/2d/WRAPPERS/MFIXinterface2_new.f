      SUBROUTINE MFIXinterfacemu(VNAME, VNO, VAR, A_M, B_M, ITMAX, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
!     USE param1
!     USE geometry
!     USE indices
      USE compar
!     USE sendrecv
!     USE leqsol
      USE functions
!     USE mpi_utility

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
!       CHARACTER(LEN=*), INTENT(IN) :: Vname
       CHARACTER(LEN=*) :: Vname
! variable number (not really used here; see calling subroutine)
       INTEGER :: VNO
! variable
!        DOUBLE PRECISION, INTENT(INOUT) :: VAR(DIMENSION_3)
        DOUBLE PRECISION :: VAR(DIMENSION_3)
! Septadiagonal matrix A_m
!        DOUBLE PRECISION, INTENT(INOUT):: A_m(DIMENSION_3, -3:3)
        DOUBLE PRECISION :: A_m(DIMENSION_3, -3:3)

! Vector b_m
!        DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
        DOUBLE PRECISION :: B_m(DIMENSION_3)
! maximum number of iterations (generally leq_it)
        INTEGER :: ITMAX
! error indicator
        INTEGER :: IER,ntrows
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

!        DOUBLE PRECISION, allocatable, dimension(:,:) :: Anew
!        DOUBLE PRECISION, allocatable, dimension(:,:) :: pos
!        DOUBLE PRECISION, allocatable, dimension(:) :: Bn,xn
!	integer, allocatable, dimension(:) :: locgl,locijk

!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
	INTEGER, PARAMETER :: cols = 7
	DOUBLE PRECISION   :: aijmax
        INTEGER :: I, J, K, IJK, IJK_GL, ie, dimred,penn, nxstart,nxend


	DOUBLE PRECISION :: Anew(DIMENSION_3,cols) 
	DOUBLE PRECISION :: pos(DIMENSION_3,cols-1) 
	integer :: locgl(DIMENSION_3)
	integer :: locijk(DIMENSION_3)
	DOUBLE PRECISION :: Bn(DIMENSION_3)
	DOUBLE PRECISION :: xn(DIMENSION_3)

	character(len=7):: string

	dimred = funijk(iend2,jend2,kend2) - funijk(istart2,jstart2,kend2) - JMAX2 + 1

!       allocate(Anew(dimred,cols))
!       allocate(pos(dimred,cols-1))
!       allocate(locgl(dimred))
!       allocate(locijk(dimred))
!       allocate(Bn(dimred))
!       allocate(xn(dimred))

	nxstart = istart2
	nxend   = iend2 
	if(myPE.gt.0) nxstart = istart2 + 1
	if(mype.ne.(numPEs-1)) nxend = iend2 - 1 

	ie = 0 
            do k = kstart2,kend2
               do i = nxstart,nxend
                  do j = jstart2,jend2
!		IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
                     IJK = funijk(i,j,k)
                     IJK_GL= funijk_gl(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     OAM = one/aijmax
                     A_M(IJK,:) = A_M(IJK,:)*OAM
                     B_M(IJK) = B_M(IJK)*OAM
!                     IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
	  ie = ie + 1 
          Anew(ie,1) =  A_M(IJK,-3)
          Anew(ie,2) =  A_M(IJK,-1)
          Anew(ie,3) =  A_M(IJK,-2)
          Anew(ie,4) =  A_M(IJK,0) ! -1.0d0
          Anew(ie,5) =  A_M(IJK,2)
          Anew(ie,6) =  A_M(IJK,1)
          Anew(ie,7) =  A_M(IJK,3)
	  Bn(ie)     =  B_M(IJK)
	  locgl(ie)  =  IJK_GL- funijk(istart2,jstart2,kstart2) + 1
	  locijk(ie) =  IJK
!	  write(800+myPE,*) i,j,k,ie,locgl(ie),ijk,ijk_gl,(A_M(ijk,:)), B_M(ijk)
           pos(ie,1) = KM_OF(IJK) - IJK
           pos(ie,2) = IM_OF(IJK) - IJK
           pos(ie,3) = JM_OF(IJK) - IJK
           pos(ie,4) = JP_OF(IJK) - IJK
           pos(ie,5) = IP_OF(IJK) - IJK
           pos(ie,6) = KP_OF(IJK) - IJK

        write(80+myPE,*)ie,IJK,i,j,k,locgl(ie),IJK_GL,pos(ie,1),pos(ie,2),pos(ie,3),pos(ie,4), & 
        pos(ie,5),pos(ie,6)

        write(90+myPE,*)ie,IJK, locgl(ie),Anew(ie,1),Anew(ie,2),Anew(ie,3),Anew(ie,4), & 
        Anew(ie,5),Anew(ie,6),Anew(ie,7), Bn(ie)
                  enddo
               enddo
            enddo
	penn = ie

       if(myPE.eq.(numPEs-1)) ntrows  =  FUNIJK_GL(iend2,jend2,kend2) - &  !maxval(locgl)
                                         FUNIJK(istart2,jstart2,kstart2) + 1 

        call MPI_BCAST(ntrows,1,MPI_INTEGER,(numPEs-1), & 
        MPI_COMM_WORLD,mpierr) 

	do i=1,ie
	write(800+myPE,*) Anew(i,1),Anew(i,2),Anew(i,3),Anew(i,4), & 
        Anew(i,5),Anew(i,6),Anew(i,7), Bn(i)
	end do 

!	if(ie.ne.dimred) then 
!	write(*,*) "MFIXinterface:array sizes does match", ie, dimred
!	stop 
!	end if 

!	write(*,*) " max_dimension", myPE, ijkstart3, ijkend3, ntrows !,maxval(lglob)
!	CALL FWrapper(Anew,Bn,  xn,pos,lglob,cols,ie,ntrows)

	write(*,*) "MFIX_wrapper", myPE, ie, DIMENSION_3, ntrows,ijkmax2,ijkmax3,penn
!	stop
!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ntrows )
	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ijkmax2)

!	call cwrapper(Bn, xn, pos, Anew,locgl, cols, ie, ijkmax2, string)
!	ie = 0
!$omp parallel do default(shared) private(ijk,i,j,k,ie)
!	do k = kstart2,kend2
!               do i = istart2,iend2
!                  do j = jstart2,jend2
!                     IJK = funijk(i,j,k)
!		     ie = ie +1 
!                     VAR(IJK)  = xn(ie)
!                  enddo
!               enddo
!            enddo

	do k= 1, ie
       	        IJK = locijk(k)
	        VAR(IJK) = xn(k)
	end do 

        END


