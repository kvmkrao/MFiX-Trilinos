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

        DOUBLE PRECISION oAm,res1,tmp 
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
	nstart = kstart
	nend   = kend
	if(myPE.eq.0)         nstart = kstart - 1 
	if(myPE.eq.(numPEs-1)) nend  = kend + 1 
	dimred = (nend-nstart+1) * (jend-jstart+1) * (iend-istart+1)

!	write(*,*) "istart", myPE, nstart,nend, istart,istart1,istart2,istart3,istart4
!	write(*,*) "iend",   myPE, nstart,nend, iend,iend1,iend2,iend3,iend4
!	write(*,*) "jstart", myPE, jstart,jstart1,jstart2,jstart3,jstart4
!	write(*,*) "jend",   myPE, jend,jend1,jend2,jend3,jend4
!	write(*,*) "kstart", myPE, kstart,kstart1,kstart2,kstart3,kstart4
!	write(*,*) "kend",   myPE, kend,kend1,kend2,kend3,kend4

!	stop
	dimred = funijk(iend,jend,nend) - funijk(istart,jstart,nstart)  + 1  !1

       allocate(Anew(dimred,cols))
       allocate(pos(dimred,cols-1))
       allocate(locgl(dimred))
       allocate(locijk(dimred))
       allocate(Bn(dimred))
       allocate(xn(dimred))

	ie = 0 
        do k = nstart,nend
               do i = istart,iend
                  do j = jstart, jend !jend1
                     IJK = funijk(i,j,k)
                     IJK_GL= funijk_gl(i,j,k)
!                     aijmax = maxval(abs(A_M(ijk,:)) )
!                     OAM = one/aijmax
!                     A_M(IJK,:) = A_M(IJK,:)*OAM
!                     B_M(IJK) = B_M(IJK)*OAM
	             ie = ie + 1 
		     locijk(ie) = ijk
                     locgl(ie)  = funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk))
		      Anew(ie,1) =  A_M(IJK,-3)
                      Anew(ie,2) =  A_M(IJK,-1)
                      Anew(ie,3) =  A_M(IJK,-2)
                      Anew(ie,4) =  A_M(IJK, 0) ! -1.0d0
                      Anew(ie,5) =  A_M(IJK, 2)
                      Anew(ie,6) =  A_M(IJK, 1)
                      Anew(ie,7) =  A_M(IJK, 3)
                      Bn(ie)     =  B_M(IJK)
	              pos(ie,1) =   KM_OF(IJK) - IJK !0 
                      pos(ie,2) =   IM_OF(IJK) - IJK !-JMAX2
                      pos(ie,3) =   JM_OF(IJK) - IJK !-1
                      pos(ie,4) =   JP_OF(IJK) - IJK !1
                      pos(ie,5) =   IP_OF(IJK) - IJK !JAMX2
                      pos(ie,6) =   KP_OF(IJK) - IJK !0
!		  write(*,*) ijk, ie, locijk(ie), locgl(ie)
!                   write(900+myPE,*)ie, IJK, locgl(ie), &  !funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk)), &
!             I_OF(IJK),J_OF(IJK),K_OF(IJK),KM_OF(IJK)-IJK, IM_OF(IJK) - IJK, & 
!            JM_OF(IJK) - IJK, JP_OF(IJK) - IJK, IP_OF(IJK) - IJK, KP_OF(IJK) - IJK

!        write(80+myPE,*)IJK, IJK,i,j,k,IJK_GL,KM_OF(IJK) - IJK, IM_OF(IJK) - IJK, & 
!        JM_OF(IJK) - IJK, JP_OF(IJK) - IJK, IP_OF(IJK) - IJK, KP_OF(IJK) - IJK

!        write(90+myPE,999)IJK, IJK,IJK_GL,A_M(IJK,-3), A_M(IJK,-1), A_M(IJK,-2), A_M(IJK,0), & 
!!        write(90+myPE,999)IJK_GL ,A_M(IJK,-3,M), A_M(IJK,-1,M), A_M(IJK,-2,M), A_M(IJK,0,M), & 
!        A_M(IJK,2), A_M(IJK,1), A_M(IJK,3), B_M(IJK),VAR(IJK)

         if(pos(ie,1).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,1)
         if(pos(ie,2).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,2)
         if(pos(ie,3).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,3)
         if(pos(ie,4).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,5)
         if(pos(ie,5).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,6)
         if(pos(ie,6).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,7)

		  end do 
		end do 
	end do 

  999   format(3i6,9f12.6)

  888	format(6i12,6f12.1)

!       if(myPE.eq.(numPEs-1)) ntrows  = maxval(locgl) ! FUNIJK_GL(iend2,jend2,kend2) - & 
                                        !FUNIJK(istart2,jstart2,kstart2) + 1 
!        call MPI_BCAST(ntrows,1,MPI_INTEGER,(numPEs-1), MPI_COMM_WORLD,mpierr) 
!	write(*,*) "i j k ", myPE, i, j, k, ntrows, ijkmax2, ijkmax3
!	stop

!	write(*,*) "istart", myPE, nstart,istart1,istart2,istart3,istart4 !, maxval(locgl)
!	write(*,*) "iend",   myPE, nend,  iend1,  iend2,  iend3,  iend4       !, maxval(locgl)
!
!	write(*,*) "jstart", myPE, jstart,jstart1,jstart2,jstart3,jstart4 !, maxval(locgl)
!	write(*,*) "jend",   myPE, jend,  jend1,  jend2,  jend3,  jend4       !, maxval(locgl)
!
!	write(*,*) "kstart", myPE, kstart,kstart1,kstart2,kstart3,kstart4 !, maxval(locgl)
!	write(*,*) "kend",   myPE, kend,  kend1,  kend2,  kend3,  kend4       !, maxval(locgl)
!
!	write(*,*) "ijkmax", myPE, ijkmax1,ijkmax2,ijkmax3,ijkmax4, ntrows

!!	if(ie.ne.dimred) then 
!!	write(*,*) "MFIXinterface:array sizes does match", ie, dimred
!!	stop 
!!	end if 
!
!!	call MPI_BARRIER(MPI_COMM_WORLD, mpierr) 
	write(*,*) " max_dimension", myPE, ijkstart3, ijkend3, ijkmax3, ijkmax2, ntrows,ireds !,maxval(lglob)

	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,dimred, ijkmax3,ITMAX, TOL)
!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ntrows, jmax2)

	ie = 0
!$omp parallel do default(shared) private(ijk,i,j,k,ie)
	do k = nstart,nend
               do i = istart,iend
                  do j = jstart,jend
!		     IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
		     ie = ie + 1
                     IJK = funijk(i,j,k)
                     VAR(IJK)  = xn(ie)
                  enddo
               enddo
            enddo

!	do k= 1, ie
!       	     IJK = locijk(k)
!	     VAR(IJK) = xn(k)
!!	     write(800+myPE,*) k,locgl(k),VAR(IJK)
!	end do 

	call send_recv(var,2)

!	res1 = 0 
!	do k=1, ie
!	    IJK = locijk(k) 
!	    tmp = B_M(IJK) -A_M(IJK,-3)*VAR(KM_OF(IJK)) -A_M(IJK,-1)*VAR(IM_OF(IJK)) & 
!                                -A_M(IJK,-2)*VAR(JM_OF(IJK)) -A_M(IJK, 0)*VAR(IJK)        &
!	                        -A_M(IJK, 2)*VAR(JP_OF(IJK)) -A_M(IJK, 1)*VAR(IP_OF(IJK)) &
!	                        -A_M(IJK, 3)*VAR(KP_OF(IJK)) 
!	    res1 = res1 + abs(tmp)
!	    write(*,*) IJK, locgl(k), tmp, res1
!	end do 

        END
