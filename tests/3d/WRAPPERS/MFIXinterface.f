      SUBROUTINE MFIXinterfacemu(VNAME, VNO, VAR, A_M, B_M, ITMAX, TOL,dimred)
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
       use machine, only: wall_time

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
	INTEGER, PARAMETER :: cols = 7
	DOUBLE PRECISION   :: aijmax

        INTEGER :: I, J, K, IJK, IJK_GL, ie, dimred, nstart,nend
	double precision :: time1, time2, time3, time4 

	real*8 r1,r2,r3,r4,r5,r6
	character(len=7) :: string
!        integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1

	real*8 :: Anew(dimred,cols) 
	real*8 :: pos(dimred,cols-1)
	integer:: locgl(dimred)
	integer:: locijk(dimred)
	real*8 :: Bn(dimred)
	real*8 :: xn(dimred)
!-----------------------------------------------

	string = "CWRAP"
	if(NODESK.gt.1) then
	time1 = WALL_TIME()

	nstart = kstart
	nend   = kend
	if(myPE.eq.0)          nstart = kstart - 1 
	if(myPE.eq.(numPEs-1)) nend   = kend   + 1 

	if(numPEs.eq.1) then
	nstart = kstart
	nend   = kend 
	end if 

! Murali commented feb 
!!	ie = 0 
!!        do k = nstart,nend
!!               do i = istart,iend
!!                  do j = jstart, jend !jend1
!!                     IJK = funijk(i,j,k)
!!                     IJK_GL= funijk_gl(i,j,k)
!!!                     aijmax = maxval(abs(A_M(ijk,:)) )
!!!                     OAM = one/aijmax
!!!                     A_M(IJK,:) = A_M(IJK,:)*OAM
!!!                     B_M(IJK) = B_M(IJK)*OAM
!!	             ie = ie + 1 
!!		     locijk(ie) = ijk
!!!                     locgl(ie)  = funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk))
!!                      locgl(ie)  = funijk_gl(i,j,k)
!!		      Anew(ie,1) =  A_M(IJK,-3)
!!                      Anew(ie,2) =  A_M(IJK,-1)
!!                      Anew(ie,3) =  A_M(IJK,-2)
!!                      Anew(ie,4) =  A_M(IJK, 0) ! -1.0d0
!!                      Anew(ie,5) =  A_M(IJK, 2)
!!                      Anew(ie,6) =  A_M(IJK, 1)
!!                      Anew(ie,7) =  A_M(IJK, 3)
!!         !if(myPE.eq.0) write(*,*) "Mwrapper",ie, Anew(ie,1),Anew(ie,2),Anew(ie,3),Anew(ie,4)
!!                      Bn(ie)     =  B_M(IJK)
!!	              pos(ie,1) =   KM_OF(IJK) - IJK !0 
!!                      pos(ie,2) =   IM_OF(IJK) - IJK !-JMAX2
!!                      pos(ie,3) =   JM_OF(IJK) - IJK !-1
!!                      pos(ie,4) =   JP_OF(IJK) - IJK !1
!!                      pos(ie,5) =   IP_OF(IJK) - IJK !JAMX2
!!                      pos(ie,6) =   KP_OF(IJK) - IJK !0
!!
!!!         if(myPE.eq.0) write(*,*) "Mwrapper",ijk, KM_OF(IJK), IM_OF(IJK), JM_OF(IJK)
!!!         if(myPE.eq.0) write(*,*) "Mwrapper",ie, pos(ie,1), pos(ie,2), pos(ie,3)
!!!         if(myPE.eq.0) write(*,*) "Mwrapper",ijk,funijk_gl(i_of(ijk),j_of(ijk),k_of(ijk)),funijk_gl(i,j,k)
!!
!!         if(pos(ie,1).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,1)
!!         if(pos(ie,2).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,2)
!!         if(pos(ie,3).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,3)
!!         if(pos(ie,4).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,5)
!!         if(pos(ie,5).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,6)
!!         if(pos(ie,6).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,7)
!!		  end do 
!!		end do 
!!	end do 


  999   format(3i6,9f12.6)
  888	format(6i12,6f12.1)

!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,dimred, ijkmax3,ITMAX, TOL)
	time2 = WALL_TIME()
!	call cwrapper(Bn, xn, pos, Anew,locgl,cols,dimred, ijkmax3,ITMAX, TOL,string,A_M,DIMENSION_3)
	call cwrapper(B_M, xn, cols,dimred, ijkmax3,ITMAX, TOL,string,A_M,DIMENSION_3)
!	call cwrapper(B_M, xn, pos, Anew,locgl,cols,dimred, ijkmax3,ITMAX, TOL,string,A_M,DIMENSION_3)
!	call FWrapper_ext(Bn, xn, pos, Anew,locgl,cols,ie, ntrows, jmax2)
	time3 = WALL_TIME()

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

	if(numPEs.gt.1) call send_recv(var,2)
	time4 = WALL_TIME()

	if(myPE.eq.0) then
	write(*,*) "time: MFIX   wrappper  ", time4 - time1 
	write(*,*) "time: MFIX wrapper, form arrays ", time2 - time1 
	write(*,*) "time: sendrecv x array ", time4 - time3 
	write(*,*) "time: C, C++ wrapppers ", time3 - time2
	end if 

	else 

        ie = 0
            do k = kstart3,kend3
               do i = istart3,iend3
                  do j = jstart3,jend3
                     IJK = funijk(i,j,k)
                     IJK_GL= funijk_gl(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     OAM = one/aijmax
                     A_M(IJK,:) = A_M(IJK,:)*OAM
                     B_M(IJK) = B_M(IJK)*OAM
                     ie =  IJK !ie + 1 
                     Anew(ie,1) =  A_M(IJK,-3)
                     Anew(ie,2) =  A_M(IJK,-1)
                     Anew(ie,3) =  A_M(IJK,-2)
                     Anew(ie,4) =  A_M(IJK,0) ! -1.0d0
                     Anew(ie,5) =  A_M(IJK,2)
                     Anew(ie,6) =  A_M(IJK,1)
                     Anew(ie,7) =  A_M(IJK,3)
                     Bn(ie)     =  B_M(IJK)
                     locgl(ie)  =  IJK_GL !- funijk(istart2,jstart2,kstart2) + 1
                     locijk(ie) =  IJK
                      pos(ie,1) = KM_OF(IJK) - IJK
                      pos(ie,2) = IM_OF(IJK) - IJK
                      pos(ie,3) = JM_OF(IJK) - IJK
                      pos(ie,4) = JP_OF(IJK) - IJK
                      pos(ie,5) = IP_OF(IJK) - IJK
                      pos(ie,6) = KP_OF(IJK) - IJK
                     if(pos(ie,1).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,1)
                     if(pos(ie,2).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,2)
                     if(pos(ie,3).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,3)
                     if(pos(ie,4).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,5)
                     if(pos(ie,5).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,6)
                     if(pos(ie,6).eq.0.0) Anew(ie,4) = Anew(ie,4) + Anew(ie,7)

                  enddo
               enddo
            enddo

        !	call FWrapper_ext(Bn, xn, pos, Anew, locgl,cols,    ie, ijkmax3, itmax,tol)
        !	call FWrapper_ext(Bn, VAR, pos, Anew,locgl,cols,dimred, ijkmax3,ITMAX, TOL)
!	call cwrapper(Bn, VAR, pos, Anew,locgl,cols,dimred, ijkmax3,ITMAX, TOL,string)
	end if

        END
