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
        DOUBLE PRECISION, allocatable, dimension(:,:) :: Anew
        DOUBLE PRECISION, allocatable, dimension(:,:) :: pos
        DOUBLE PRECISION, allocatable, dimension(:) :: Bn,xn
	integer, allocatable, dimension(:) :: locgl,locijk

!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
	INTEGER, PARAMETER :: cols = 7
	DOUBLE PRECISION   :: aijmax

!       added 
!       DOUBLE PRECISION :: AVar(DIMENSION_3)
        INTEGER :: I, J, K, IJK, IJK_GL, ie, dimred, nxstart,nxend
	integer :: nistart,niend
!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
!        integer :: class, interval
!        integer :: j_start(2), j_end(2)
!       added 
!        allocate(Bnew(DIMENSION_3))

!	real*8 :: Anew(DIMENSION_3,cols) 
!	real*8 :: pos(DIMENSION_3,cols) 
!	integer:: lglob(DIMENSION_3)
!	real*8 :: Bn(DIMENSION_3)
!	real*8 :: xn(DIMENSION_3)

!	allocate(Anew(DIMENSION_3,cols))
!	allocate(pos(DIMENSION_3,cols-1))
!	allocate(lglob(DIMENSION_3))
!	allocate(Bn(DIMENSION_3))
!	allocate(xn(DIMENSION_3))
!-----------------------------------------------

!	dimred = funijk(iend2,jend2,kend2) - funijk(istart2,jstart2,kend2) - JMAX2 + 1
!	dimred =  ijkmax2 !DIMENSION_3  !
	nistart = istart1
	niend   = iend1 
	if(myPE.eq.0) nistart = nistart - 2
	if(mype.eq.(numPEs-1)) niend = niend + 2
	dimred = funijk(niend,jend1,kend1) - funijk(nistart,1,kstart1)  + 2  !1

       allocate(Anew(dimred,cols))
       allocate(pos(dimred,cols-1))
       allocate(locgl(dimred))
       allocate(locijk(dimred))
       allocate(Bn(dimred))
       allocate(xn(dimred))
!	nxstart = istart2
!	nxend   = iend2
!	if(myPE.eq.0) nxstart = istart - 1   
!	if(myPE.gt.0) nxstart = istart2 + 1
!	if(mype.ne.(numPEs-1)) nxend = iend2 - 1 

	ie = 0
!$omp parallel do default(shared) private(ijk,ijk_gl,i,j,k,oam,aijmax, ie)
            do k = kstart1,kend1
               do i = nistart,niend
                  do j = 1, jmax2 !jend1
!		IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
                     IJK = funijk(i,j,k)
                     IJK_GL= funijk_gl(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     OAM = one/aijmax
                     A_M(IJK,:) = A_M(IJK,:)*OAM
                     B_M(IJK) = B_M(IJK)*OAM
!	  IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
!	  if((i.ge.nxstart).and.(i.le.nxend)) then
	  ie = ie + 1 
!         if(abs(A_M(IJK,-2)).lt.1) 
          Anew(ie,1) =  A_M(IJK,-3)
          Anew(ie,2) =  A_M(IJK,-1)
          Anew(ie,3) =  A_M(IJK,-2)
          Anew(ie,4) =  A_M(IJK,0) ! -1.0d0
          Anew(ie,5) =  A_M(IJK,2)
          Anew(ie,6) =  A_M(IJK,1)
          Anew(ie,7) =  A_M(IJK,3)
	  Bn(ie)     =  B_M(IJK)
	  if(i.lt.1)        Bn(ie) = - B_M(IJK)
	  if(i.gt.iend+1)   Bn(ie) = - B_M(IJK)
	  locgl(ie)  =  i * JMAX2 + j ! (i-2) * JMAX2 + j  !IJK_GL - 2* JMAX2 !funijk(istart1,jstart1,kstart1) + 1
	  locijk(ie) =  IJK
          pos(ie,1) =   KM_OF(IJK) - IJK
          pos(ie,2) =   IM_OF(IJK) - IJK
          pos(ie,3) =   JM_OF(IJK) - IJK
          pos(ie,4) =   JP_OF(IJK) - IJK
          pos(ie,5) =   IP_OF(IJK) - IJK
          pos(ie,6) =   KP_OF(IJK) - IJK

        write(80+myPE,*)ie,IJK,i,j,k,locgl(ie),IJK_GL,pos(ie,1),pos(ie,2),pos(ie,3),pos(ie,4), & 
        pos(ie,5),pos(ie,6)

        write(90+myPE,*)ie,IJK,locgl(ie),Anew(ie,1),Anew(ie,2),Anew(ie,3),Anew(ie,4), & 
        Anew(ie,5),Anew(ie,6),Anew(ie,7), Bn(ie)
                  enddo
               enddo
            enddo

!       if(myPE.eq.(numPEs-1)) ntrows  = maxval(locgl) ! FUNIJK_GL(iend2,jend2,kend2) - & 
!                                        !FUNIJK(istart2,jstart2,kstart2) + 1 
!        call MPI_BCAST(ntrows,1,MPI_INTEGER,(numPEs-1), MPI_COMM_WORLD,mpierr) 

!	if(myPE.gt.0) then 
!	do i = 1, JMAX2
!          Anew(i,1) = 0.0d0  
!          Anew(i,2) = 0.0d0  
!          Anew(i,3) = 0.0d0  
!          Anew(i,4) = 0.0d0  ! -1.0d0
!          Anew(i,5) = 0.0d0  
!          Anew(i,6) = 0.0d0  
!          Anew(i,7) = 0.0d0  
!          Bn(i)     = 0.0d0  
!        enddo
!	end if

	write(*,*) "istart", myPE, istart1,istart2,istart3,istart4, maxval(locgl),ijkmax3,ijkmax2
	write(*,*) "iend", myPE, iend1,  iend2,  iend3,  iend4, maxval(locgl),ijkmax3,ijkmax2
!	stop 
	if(ie.ne.dimred) then 
	write(*,*) "MFIXinterface:array sizes does match", ie, dimred
	stop 
	end if 

	call MPI_BARRIER(MPI_COMM_WORLD, mpierr) 
!	write(*,*) " max_dimension", myPE, ijkstart3, ijkend3, ntrows !,maxval(lglob)
!	CALL FWrapper(Anew,Bn,  xn,pos,lglob,cols,ie,ntrows)

!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,dimred, ntrows)
!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ijkmax3, jmax2)

!	if(myPE.eq.0) ie = ie - JMAX2*2 

!	stop
	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ijkmax3, jmax2)

!	stop

!	ie = 0
!$omp parallel do default(shared) private(ijk,i,j,k,ie)
!	do k = kstart1,kend1
!               do i = istart1,iend1
!                  do j = 1,jend1
!!		     IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
!		     ie = ie + 1
!                     IJK = funijk(i,j,k)
!                     VAR(IJK)  = xn(ie)
!                  enddo
!               enddo
!            enddo

!	stop
	do k= 1, ie
       	IJK = locijk(k)
	        VAR(IJK) = xn(k)
	end do 

        END


