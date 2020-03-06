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
        INTEGER :: I, J, K, IJK, IJK_GL, ie, dimred,penn
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

	dimred = funijk(iend2,jend2,kend2) - funijk(istart2,jstart2,kend2) + 1

       allocate(Anew(dimred,cols))
       allocate(pos(dimred,cols-1))
       allocate(locgl(dimred))
       allocate(locijk(dimred))
       allocate(Bn(dimred))
       allocate(xn(dimred))

!	ie = 0 
!       if(myPE.eq.0) then
!!$omp parallel do default(shared) private(ijk,ijk_gl,i,j,k,oam,aijmax, ie)
!	i = istart3
!        do k = kstart3,kend3
!          do j = jstart3,jend3
!          IJK = funijk(i,j,k)
!          IJK_GL= funijk_gl(i,j,k)
!          aijmax = maxval(abs(A_M(ijk,:)) )
!          OAM = one/aijmax
!          A_M(IJK,:) = A_M(IJK,:)*OAM
!          B_M(IJK) = B_M(IJK)*OAM
!          ie = ie + 1
!          Anew(ie,1) =  A_M(IJK,-3)
!          Anew(ie,2) =  A_M(IJK,-1)
!          Anew(ie,3) =  A_M(IJK,-2)
!          Anew(ie,4) =  A_M(IJK,0) ! -1.0d0
!          Anew(ie,5) =  A_M(IJK,2)
!          Anew(ie,6) =  A_M(IJK,1)
!          Anew(ie,7) =  A_M(IJK,3)
!          Bn(ie)     =  B_M(IJK)
!          locgl(ie)  =  IJK_GL !- funijk(istart2,jstart2,kstart2) + 1
!          locijk(ie) =  IJK
!!          write(800+myPE,*) i,j,k,ie,locgl(ie),locijk(ie),ijk_gl,(A_M(ijk,:)), B_M(ijk)
!          write(900,*) i,j,k,ie,locgl(ie),locijk(ie),ijk_gl,(A_M(ijk,:)), B_M(ijk)
!          pos(ie,1) = KM_OF(IJK) - IJK
!          pos(ie,2) = IM_OF(IJK) - IJK
!          pos(ie,3) = JM_OF(IJK) - IJK
!          pos(ie,4) = JP_OF(IJK) - IJK
!          pos(ie,5) = IP_OF(IJK) - IJK
!          pos(ie,6) = KP_OF(IJK) - IJK
!	  end do 
!	end do 
!	end if


!	write(*,*) JMAX, JMAX2, JMAX3
!	stop

	ie = 0 
!$omp parallel do default(shared) private(ijk,ijk_gl,i,j,k,oam,aijmax, ie)
            do k = kstart2,kend2
               do i = istart2,iend2
                  do j = jstart2,jend2
		IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
                     IJK = funijk(i,j,k)
                     IJK_GL= funijk_gl(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     OAM = one/aijmax
                     A_M(IJK,:) = A_M(IJK,:)*OAM
                     B_M(IJK) = B_M(IJK)*OAM
	  ie = ie + 1 
!         if(abs(A_M(IJK,-2)).lt.1) 
          Anew(ie,1) =  A_M(IJK,-3)
          Anew(ie,2) =  A_M(IJK,-1)
!         if(abs(A_M(IJK,-1)).lt.1) 
          Anew(ie,3) =  A_M(IJK,-2)
          Anew(ie,4) =  A_M(IJK,0) ! -1.0d0
!         if(abs(A_M(IJK, 1)).lt.1) 
          Anew(ie,5) =  A_M(IJK,2)
!         if(abs(A_M(IJK, 2)).lt.1) 
          Anew(ie,6) =  A_M(IJK,1)
          Anew(ie,7) =  A_M(IJK,3)
	  Bn(ie)     =  B_M(IJK)
	  locgl(ie)  =  IJK_GL - JMAX2 !funijk(istart2,jstart2,kstart2) + 1 
!	  if(myPE.gt.0) locgl(ie) = locgl(ie) - 2 
	  locijk(ie) =  IJK
	  write(800+myPE,*) i,j,k,ie,locgl(ie),locijk(ie),ijk_gl,(A_M(ijk,:)), B_M(ijk)
           pos(ie,1) = KM_OF(IJK) - IJK
           pos(ie,2) = IM_OF(IJK) - IJK
           pos(ie,3) = JM_OF(IJK) - IJK
           pos(ie,4) = JP_OF(IJK) - IJK
           pos(ie,5) = IP_OF(IJK) - IJK
           pos(ie,6) = KP_OF(IJK) - IJK

!        write(80+myPE,*)ie,IJK,i,j,k,locgl(ie),IJK_GL,pos(ie,1),pos(ie,2),pos(ie,3),pos(ie,4), & 
!        pos(ie,5),pos(ie,6)

!        write(90+myPE,*)ie,IJK,Anew(ie,1),Anew(ie,2),Anew(ie,3),Anew(ie,4), & 
!        Anew(ie,5),Anew(ie,6),Anew(ie,7), Bn(ie)
                  enddo
               enddo
            enddo

       if(myPE.eq.(numPEs-1)) ntrows  = maxval(locgl)
!FUNIJK_GL(iend2,jend2,kend2) - & 
                                        !FUNIJK(istart2,jstart2,kstart2) + 1 
        call MPI_BCAST(ntrows,1,MPI_INTEGER,(numPEs-1), & 
        MPI_COMM_WORLD,mpierr) 
	write(*,*) "ntrows" , ntrows-2

	call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
!	if(ie.ne.dimred) then 
!	write(*,*) "MFIXinterface:array sizes does match", ie, dimred
!	stop 
!	end if 

!	write(*,*) " max_dimension", myPE, ijkstart3, ijkend3, ntrows !,maxval(lglob)
!	CALL FWrapper(Anew,Bn,  xn,pos,lglob,cols,ie,ntrows)

!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,dimred, ntrows)

	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ntrows)

!	ie = 0
!!$omp parallel do default(shared) private(ijk,i,j,k,ie)
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
!	VAR(IJK) = xn(k)
	end do 

        END

