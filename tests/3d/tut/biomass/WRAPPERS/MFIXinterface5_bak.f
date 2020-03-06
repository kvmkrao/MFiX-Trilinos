!       SUBROUTINE MFIXinterfacemu(VNAME, VNO, VAR, A_M, B_M, pos, ireds, ITMAX, TOL)
       SUBROUTINE MFIXinterfacemu(VNAME, VNO, VAR, A_M, B_M, ITMAX, TOL)
! ******************* Modules  **********************
      USE param
!     USE param1
!     USE geometry
!     USE indices
      USE compar
!     USE sendrecv
!     USE leqsol
      USE functions
!     USE mpi_utility
      USE funits
      USE indices
      USE sendrecv
      USE mpi_utility
      USE functions
! *************************************************** 

      IMPLICIT NONE
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
        INTEGER :: ntrows,class,ireds
	real*8 :: TOL
!-----------------------------------------------
! OVERRELAXATION FACTOR
!        DOUBLE PRECISION, PARAMETER :: OMEGA = 1.0  !1.2
!        integer :: iidebug
!        parameter( iidebug = 0 )
!-----------------------------------------------
!        INTEGER :: IJK
        INTEGER :: ITER, II

        DOUBLE PRECISION oAm
        DOUBLE PRECISION, allocatable, dimension(:,:) :: Anew
        DOUBLE PRECISION, allocatable, dimension(:,:) :: pos
        DOUBLE PRECISION, allocatable, dimension(:) :: Bn,xn
	integer, allocatable, dimension(:) :: locgl,locijk

!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
	INTEGER, PARAMETER :: cols = 7
	DOUBLE PRECISION   :: aijmax

        INTEGER :: I, J, K, IJK, IJK_GL, ie, dimred, nstart,nend

	real*8 r1,r2,r3,r4,r5,r6
!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
!        integer :: class, interval
!        integer :: j_start(2), j_end(2)
!       added 

!	real*8 :: Anew(DIMENSION_3,cols) 
!	real*8 :: pos(ireds,6)
!	integer:: lglob(DIMENSION_3)
!	real*8 :: Bn(DIMENSION_3)
!	real*8 :: xn(DIMENSION_3)

!-----------------------------------------------
	nstart = istart
	nend   = iend
	if(myPE.eq.0)         nstart = istart - 1 
	if(myPE.eq.(numPEs-1)) nend  = iend + 1 
	dimred = (nend-nstart+1) * (jend-jstart+1) * (kend-kstart+1)

!	write(*,*) myPE, nstart,nend, istart,istart1,istart2,istart3,istart4
!	write(*,*) myPE, nstart,nend, iend,iend1,iend2,iend3,iend4
!	write(*,*) myPE, nstart,nend, jstart,jstart1,jstart2,jstart3,jstart4
!	dimred = funijk(nend,jend,kend) - funijk(nstart,jstart,kstart)  + 1  !1

       allocate(Anew(dimred,cols))
       allocate(pos(dimred,cols-1))
       allocate(locgl(dimred))
       allocate(locijk(dimred))
       allocate(Bn(dimred))
       allocate(xn(dimred))

	write(*,*) "ireds", dimred
!	stop

	ie = 0
!$omp parallel do default(shared) private(ijk,ijk_gl,i,j,k,oam,aijmax, ie)

!	do k = kstart3,kend3
!               do i = istart3,iend3
!!               do i = istart4,iend4
!                  do j = jstart3, jend3 !jend1
!                     IJK = funijk(i,j,k)
!                     IJK_GL= funijk_gl(i,j,k)
!                      aijmax = maxval(abs(A_M(ijk,:)) )
!                      OAM = one/aijmax
!                      A_M(IJK,:) = A_M(IJK,:)*OAM
!                      B_M(IJK) = B_M(IJK)*OAM
!          ie = ie + 1
!          Anew(ie,1) =  A_M(IJK,-3)
!          Anew(ie,2) =  A_M(IJK,-1)
!          Anew(ie,3) =  A_M(IJK,-2)
!          Anew(ie,4) =  A_M(IJK,0) ! -1.0d0
!          Anew(ie,5) =  A_M(IJK,2)
!          Anew(ie,6) =  A_M(IJK,1)
!          Anew(ie,7) =  A_M(IJK,3)
!          Bn(ie)     =  B_M(IJK)
!          locgl(ie)  =  IJK_GL !(i*imax2)+ j +(k-1)*(imax2+2)*jmax2
!          locijk(ie) =  IJK
!	  IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE

!!          pos(ie,1) =   KM_OF(IJK) - IJK !0 
!!          pos(ie,2) =   - (IJK - IM_OF(IJK)) !-JMAX2
!!!          pos(ie,2) =   INCREMENT_FOR_MP(2,class)
!!          pos(ie,3) =   JM_OF(IJK) - IJK !-1
!!          pos(ie,4) =   JP_OF(IJK) - IJK !1
!!          pos(ie,5) =   IP_OF(IJK) - IJK !JAMX2
!!          pos(ie,6) =   KP_OF(IJK) - IJK !0
!
!!	  II = (J + C0 + I*C1 + K*C2)
!!          r1 = INCREMENT_FOR_MP(5,class)
!!	 II = (J + C0 + I*C1 + K*C2)
!!         CLASS = CELL_CLASS(IJK)
!!	 r1 = FUNIJK_GL(I,J,K-1) - FUNIJK_GL(I,J,K) !INCREMENT_FOR_MP(5,class)
!
!         write(80+myPE,*)ie,IJK,i,j,k,locgl(ie),pos(ie,1),pos(ie,2),pos(ie,3),pos(ie,4), & 
!        pos(ie,5),pos(ie,6)
!
!!        write(80+myPE,*)ie,IJK,i,j,k,locgl(ie),r1,r2,r3,r4,r5,r6 
!
!        write(90+myPE,*)ie,IJK,locgl(ie),Anew(ie,1),Anew(ie,2),Anew(ie,3),Anew(ie,4), & 
!        Anew(ie,5),Anew(ie,6),Anew(ie,7), Bn(ie)
!
!!         if(pos(ie,1).eq.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,1)
!!         if(pos(ie,2).eq.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,2)
!!         if(pos(ie,3).eq.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,3)
!!         if(pos(ie,4).eq.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,5)
!!         if(pos(ie,5).eq.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,6)
!!         if(pos(ie,6).eq.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,7)
!
!                  enddo
!               enddo
!            enddo

	do ijk=1,dimred
		Anew(i,1) = 0.0d0
		Anew(i,2) = 0.0d0
		Anew(i,3) = 0.0d0
		Anew(i,4) = -1.0d0
		Anew(i,5) = 0.0d0
		Anew(i,6) = 0.0d0
		Anew(i,7) = 0.0d0
	        Bn(i)      = 0.0d0
                locgl(ijk)  = funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk))
	end do 


	ie = 0 
        do ijk = ijkstart3,ijkend3
!	locijk(ijk) = ijk
!        locgl(ijk)  = funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk))
	if(I_OF(IJK).eq.0.or.I_OF(IJK).eq.13) then
	locijk(ijk) = ijk
        locgl(ijk)  = funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk))
	write(*,*) ijk, locijk(ijk), locgl(ijk)
	end if 
!       I_OF(IJK),J_OF(IJK), K_OF(IJK)

        IF(.NOT.IS_ON_myPE_owns(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
!  909	continue
         ie = ie + 1
          Anew(ijk,1) =  A_M(IJK,-3)
          Anew(ijk,2) =  A_M(IJK,-1)
          Anew(ijk,3) =  A_M(IJK,-2)
          Anew(ijk,4) =  A_M(IJK,0) ! -1.0d0
          Anew(ijk,5) =  A_M(IJK,2)
          Anew(ijk,6) =  A_M(IJK,1)
          Anew(ijk,7) =  A_M(IJK,3)
          Bn(ijk)     =  B_M(IJK)
	  pos(ijk,1) =   KM_OF(IJK) - IJK !0 
          pos(ijk,2) =   IM_OF(IJK) - IJK !-JMAX2
          pos(ijk,3) =   JM_OF(IJK) - IJK !-1
          pos(ijk,4) =   JP_OF(IJK) - IJK !1
          pos(ijk,5) =   IP_OF(IJK) - IJK !JAMX2
          pos(ijk,6) =   KP_OF(IJK) - IJK !0

          locgl(ijk)  =  funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk)) !- IMAX2 !(i*imax2)+ j +(k-1)*(imax2+2)*jmax2
          locijk(ijk) =  IJK
        write(900+myPE,*)ie, IJK, locgl(ijk), &  !funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk)), &
        I_OF(IJK),J_OF(IJK),K_OF(IJK),KM_OF(IJK)-IJK, IM_OF(IJK) - IJK, & 
        JM_OF(IJK) - IJK, JP_OF(IJK) - IJK, IP_OF(IJK) - IJK, KP_OF(IJK) - IJK
!	write(*,*) ijk, locijk(ijk), locgl(ijk) 
       end do

	do i=1,dimred 
	write(*,*) i, locijk(i), locgl(i)
	end do 

!	stop

!  888	format(6(i12","),6(1f12.1","))
!  999	format(3(i12","),8(f12.1","))

  888	format(6i12,6f12.1)
  999	format(3i12,8f12.1)

       if(myPE.eq.(numPEs-1)) ntrows  = maxval(locgl) ! FUNIJK_GL(iend2,jend2,kend2) - & 
                                        !FUNIJK(istart2,jstart2,kstart2) + 1 
        call MPI_BCAST(ntrows,1,MPI_INTEGER,(numPEs-1), MPI_COMM_WORLD,mpierr) 
	write(*,*) "i j k ", myPE, i, j, k, ntrows, ijkmax2, ijkmax3
!	stop

	write(*,*) "istart", myPE, nstart,istart1,istart2,istart3,istart4 !, maxval(locgl)
	write(*,*) "iend",   myPE, nend,  iend1,  iend2,  iend3,  iend4       !, maxval(locgl)

	write(*,*) "jstart", myPE, jstart,jstart1,jstart2,jstart3,jstart4 !, maxval(locgl)
	write(*,*) "jend",   myPE, jend,  jend1,  jend2,  jend3,  jend4       !, maxval(locgl)

	write(*,*) "kstart", myPE, kstart,kstart1,kstart2,kstart3,kstart4 !, maxval(locgl)
	write(*,*) "kend",   myPE, kend,  kend1,  kend2,  kend3,  kend4       !, maxval(locgl)

	write(*,*) "ijkmax", myPE, ijkmax1,ijkmax2,ijkmax3,ijkmax4, ntrows
!	stop 

!	if(ie.ne.dimred) then 
!	write(*,*) "MFIXinterface:array sizes does match", ie, dimred
!	stop 
!	end if 

!	call MPI_BARRIER(MPI_COMM_WORLD, mpierr) 
	write(*,*) " max_dimension", myPE, ijkstart3, ijkend3, ijkmax3, ijkmax2, ntrows,ireds !,maxval(lglob)

	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,dimred, ijkmax3,ITMAX, TOL)
!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ntrows, jmax2)

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

	do k= 1, ie
       	IJK = locijk(k)
	        VAR(IJK) = xn(k)
	        write(800+myPE,*) k,locgl(k),VAR(IJK)
	end do 

	call send_recv(var,2)
!	call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
        END
