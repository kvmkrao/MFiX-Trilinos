!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_LIN_EQ                                            C
!  Purpose: Interface for linear equation solver                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_LIN_EQ(VNAME, Vno, VAR, A_M, B_M, M, ITMAX,&
                              METHOD, SWEEP, TOL1, PC, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE residual
      USE toleranc
      USE leqsol
      use machine, only: wall_time
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number
!     Note: not really used beyond this subroutine. here it is
!     used for potentially adjusting the tolerances but it is
!     currently disabled code.
!     1 = pressure correction equation
!     2 = solids correction equation or gas/solids continuity
!     3 = gas/solids u-momentum
!     4 = gas/solids v-momentum
!     5 = gas/solids w-momentum
!     6 = temperature
!     7 = species
!     8 = granular temperature
!     9 = scalar, E_Turb_G, k_Turb_G
      INTEGER, INTENT(IN) :: Vno
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! phase index
      INTEGER, INTENT(IN) :: M
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! linear equation solver method (generally leq_method)
!     1 = sor
!     2 = bicgstab (default)
!     3 = gmres
!     5 = cg
      INTEGER, INTENT(IN) :: METHOD
! sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=4), INTENT(IN) :: SWEEP
! convergence tolerance for leq solver (leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL1
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) :: PC
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Adjust LEQ tolerance flag
      LOGICAL, PARAMETER :: adjust_leq_tol = .FALSE.
      LOGICAL, PARAMETER :: leq_tol_scheme1 = .FALSE.
! currently only used for gmres routine
      INTEGER, PARAMETER :: MAX_IT = 1
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
      DOUBLE PRECISION :: max_resid_local, tol_resid_max
! convergence tolerance for leq solver
      DOUBLE PRECISION :: TOL
! transpose of septadiaganol matrix A_M
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A_MT
! indices
      INTEGER :: II, IJK
! for constructing local character strings
      CHARACTER(LEN=80) :: LINE0, LINE1
!-----------------------------------------------
        DOUBLE PRECISION :: time1, time2,time3 

        CHARACTER(LEN=20) :: key
        INTEGER :: KEY_EXIT
        character(len=30) :: filename

!	write(*,*) "solve_lin_eq.f", VNAME, METHOD
!	 call write_parallel_info()

! Adjusting the tolerances
!        open(810,file='solv-time.txt', Access = 'append',status='unknown')
!	open(820,file='solv-timen.txt', Access = 'append',status='unknown')	
! ----------------------------------------------------------------
      IF(adjust_leq_tol) THEN
         max_resid_local = maxval(resid(:,M),1)
         tol_resid_max   = max(TOL_RESID, TOL_RESID_T, TOL_RESID_TH, TOL_RESID_X)
         IF(leq_tol_scheme1.AND.resid(Vno,M).LT.1.0D-1) THEN
            if(Vno.le.5) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID)
            elseif (Vno.eq.6) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_T)
            elseif (Vno.eq.7) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_X)
            elseif (Vno.eq.8) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_Th)
            endif
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, resid(Vno,M)
         ELSEIF(max_resid_local.LT.1.0D-1) THEN
            TOL = MAX(TOL1,TOL1*max_resid_local/TOL_RESID_max)
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, max_resid_local
         ENDIF
      ELSE
        TOL = TOL1
      ENDIF
! ----------------------------------------------------------------<<<


! Solve the linear system of equations
! ---------------------------------------------------------------->>>
      SELECT CASE (METHOD)
      CASE (1)
! SOR: Successive Sver Relaxation method from Templates
        CALL LEQ_SOR (VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), &
                      ITMAX, IER)

      CASE (2)
! BICGSTAB: BIConjugate Gradients STabilized method
         IF(do_transpose) THEN  ! mfix.dat keyword default=false
            allocate( A_mt(-3:3, ijkstart3:ijkend3 ))
!!$omp parallel do private(ijk,ii)
            DO ijk=ijkstart3,ijkend3
               do ii=-3,3
                  A_mt(ii,ijk) = A_m(ijk,ii,M)
               enddo
            ENDDO
            call leq_bicgst(VNAME, VNO, VAR, A_Mt(:,:), B_M(:,M), &
                            SWEEP, TOL, PC, ITMAX, IER)
            deallocate( A_mt )
         ELSE
           time1 = WALL_TIME()
            call leq_bicgs(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                           SWEEP, TOL, PC, ITMAX, IER)
         ENDIF


      CASE (3)
! GMRES: A Generalized Minimal RESidual Algorithm
         time1 = WALL_TIME()
         call leq_gmres(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                        SWEEP, TOL, ITMAX, MAX_IT, IER)
         write(810,*) VNAME, method, WALL_TIME() - time1

      CASE (4)
! Mix:
         IER = 0
         call leq_bicgs(VNAME,VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP,&
                       TOL, PC, ITMAX, IER)
         IF (IER .eq. -2) THEN
            IER = 0
            print*,'calling leq_gmres', Vname
            call leq_gmres(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                           SWEEP, TOL, ITMAX, MAX_IT, IER)
         ENDIF


      CASE (5)
! CG: Conjugate Gradients
         time1 = WALL_TIME()
         call leq_cg(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP,&
                     TOL, PC, ITMAX, IER)
        
        write(810,*) VNAME, method, WALL_TIME() - time1

        CASE(6)
        time1 = WALL_TIME()
!        CALL MFIXinterfacemu(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), ITMAX,IER)
!        write(810,*) VNAME, method, WALL_TIME() - time1
	call leq_bicgs(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                           SWEEP, TOL, PC, ITMAX, IER)

!       CASE (6) !- Disabled
! LSOR: Line Successive Over Relaxation method
!       CALL LEQ_LSOR(VNAME, VAR, A_M(:,:,M), B_M(:,M), ITMAX, IER)

!        call leq_bicgs(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
!                           SWEEP, TOL, PC, ITMAX, IER)

        write(810,*) VNAME, method, WALL_TIME() - time1 

      CASE DEFAULT
         LINE0(1:14) = 'SOLVE_LIN_EQ: '
         LINE0(15:80)= VName
         WRITE(LINE1,'(A, I2, A)') &
             'Error: LEQ_METHOD = ', METHOD, ' is invalid'
         CALL WRITE_ERROR(LINE0, LINE1, 1)
         CALL mfix_exit(myPE)
      END SELECT
! ----------------------------------------------------------------<<<

    

      RETURN
      END SUBROUTINE SOLVE_LIN_EQ


       SUBROUTINE WRITE_ABVAR1(A_M, B_M, var, file1)

!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31
!12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix

      USE geometry
      USE compar
      USE mpi_utility
      USE indices
      USE functions

      use machine

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Local index
        INTEGER          L
!
!                      cell index
        INTEGER          IJK
        INTEGER :: ITER, ie,nrows
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3)
!
!                      Source vector
      DOUBLE PRECISION b_m(DIMENSION_3)
!                      Source vector
      DOUBLE PRECISION var(DIMENSION_3)

      double precision, allocatable :: array1(:) , array2(:)   !//
      double precision, allocatable :: am(:,:)
      character(len=30) :: file1

     integer i, j, k
     if (myPE == PE_IO) then
         allocate (array1(ijkmax3))
         allocate (array2(ijkmax3))
         allocate (am(ijkmax3,-3:3))
      else
         allocate (array1(1))
         allocate (array2(1))
         allocate (am(1,-3:3))
      end if

      if (myPE == PE_IO) then
         CALL START_LOG
         WRITE (*,*)  ' Note : write_am_m is VERY inefficient '
!         WRITE (100+ITER, '(A,A)') &
!           '  IJK  I  J  K   b         s         w         p         e
!           ', &
!           '  n         t        Source     Variable'
        open(12, file=file1, access ='append',status='unknown')
!       open(12, file='lineqn',status='unknown')
      end if


      do L = -3,3
      call gather(a_m(:,L),array1,root)
        DO K = Kmin2, Kmax2
          DO I = Imin2, Imax2
           DO J = Jmin2, Jmax2
             IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
             IJK = FUNIJK_GL(IMAP_C(I),JMAP_C(J),KMAP_C(K))
!            IJK = FUNIJK_GL(I,J,K)
             if (myPE == PE_IO) am(ijk,l) = array1(ijk)
           END DO
          END DO
        END DO
      end do

      call gather(var(:),array1,root)
      call gather(b_m(:),array2,root)

      ie = 0
      DO K = Kmin2, Kmax2
        DO I = Imin2, Imax2
          DO J = Jmin2, Jmax2
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
!           IJK = FUNIJK_GL(I,J,K)
            IJK = FUNIJK_GL(IMAP_C(I),JMAP_C(J),KMAP_C(K))
            ie = ie +  1
            if (myPE == PE_IO) then
            WRITE(12,'(I8, 9(1X,1E20.10))') FUNIJK_IO(I,J,K),&
                     (AM(ijk,L),L=-3,3), array2(IJK),array1(IJK)
!           WRITE(*,'(2I8, 8(1X,G17.10))') ie, FUNIJK_IO(I,J,K),&
!                                    (AM(ijk,L),L=-3,3), array2(IJK)
            end if

          END DO
        END DO
      END DO

        if (myPE == PE_IO) CALL END_LOG


        deallocate (array1)    !//
        deallocate (array2)    !//
        deallocate(AM)
        if (myPE == PE_IO) close(12)
        RETURN

        END SUBROUTINE WRITE_ABVAR1

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization: array1,array2,i,j,k
!// 400 Added mpi_utility module and other global reduction (gather)
!calls

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_LIN_EQ                                            C
!  Purpose: Interface for linear equation solver                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_LIN_EQM(VNAME, Vno, VAR, A_M, B_M, M, ITMAX,&
                              METHOD, SWEEP, TOL1, PC, IER,NIT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE residual
      USE toleranc
      USE leqsol
      use machine, only: wall_time
      USE sendrecv
      USE functions
	use mpi_utility

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number
!     Note: not really used beyond this subroutine. here it is
!     used for potentially adjusting the tolerances but it is
!     currently disabled code.
!     1 = pressure correction equation
!     2 = solids correction equation or gas/solids continuity
!     3 = gas/solids u-momentum
!     4 = gas/solids v-momentum
!     5 = gas/solids w-momentum
!     6 = temperature
!     7 = species
!     8 = granular temperature
!     9 = scalar, E_Turb_G, k_Turb_G
      INTEGER, INTENT(IN) :: Vno
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! phase index
      INTEGER, INTENT(IN) :: M
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! linear equation solver method (generally leq_method)
!     1 = sor
!     2 = bicgstab (default)
!     3 = gmres
!     5 = cg
      INTEGER, INTENT(IN) :: METHOD
! sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=4), INTENT(IN) :: SWEEP
! convergence tolerance for leq solver (leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL1
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) :: PC
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Adjust LEQ tolerance flag
      LOGICAL, PARAMETER :: adjust_leq_tol = .FALSE.
      LOGICAL, PARAMETER :: leq_tol_scheme1 = .FALSE.
! currently only used for gmres routine
      INTEGER, PARAMETER :: MAX_IT = 1
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
      DOUBLE PRECISION :: max_resid_local, tol_resid_max
! convergence tolerance for leq solver
      DOUBLE PRECISION :: TOL
! transpose of septadiaganol matrix A_M
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A_MT
! indices
      INTEGER :: II, IJK,NIT, J , K, I, IJK_GL , isk
! for constructing local character strings
      CHARACTER(LEN=80) :: LINE0, LINE1
!-----------------------------------------------
        DOUBLE PRECISION :: time1, time2,time3, OAM,aijmax 

        CHARACTER(LEN=20) :: key
        INTEGER :: KEY_EXIT,ie,irarr
        character(len=30) :: filename
	real*8 tmp, res
!	integer inm1, inm2,inm3,in0,inp1,inp2,inp3,nstart,nend
	integer nstart,nend, dimred


!	write(*,*) "solve_lin_eqm.f", VNAME, METHOD

! Adjusting the tolerances
      open(810,file='solv-time.txt', Access = 'append',status='unknown')
      open(820,file='solv-timen.txt', Access = 'append',status='unknown')	
! ----------------------------------------------------------------
      IF(adjust_leq_tol) THEN
         max_resid_local = maxval(resid(:,M),1)
         tol_resid_max   = max(TOL_RESID, TOL_RESID_T, TOL_RESID_TH, TOL_RESID_X)
         IF(leq_tol_scheme1.AND.resid(Vno,M).LT.1.0D-1) THEN
            if(Vno.le.5) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID)
            elseif (Vno.eq.6) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_T)
            elseif (Vno.eq.7) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_X)
            elseif (Vno.eq.8) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_Th)
            endif
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, resid(Vno,M)
         ELSEIF(max_resid_local.LT.1.0D-1) THEN
            TOL = MAX(TOL1,TOL1*max_resid_local/TOL_RESID_max)
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, max_resid_local
         ENDIF
      ELSE
        TOL = TOL1
      ENDIF
! ----------------------------------------------------------------<<<


! Solve the linear system of equations
! ---------------------------------------------------------------->>>
       SELECT CASE (METHOD)
       CASE (1)
! SOR: Successive Sver Relaxation method from Templates
       CALL LEQ_SOR (VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), &
                      ITMAX, IER)

      CASE (2)
! BICGSTAB: BIConjugate Gradients STabilized method
       IF(do_transpose) THEN  ! mfix.dat keyword default=false
       allocate( A_mt(-3:3, ijkstart3:ijkend3 ))
!!$omp parallel do private(ijk,ii)
            DO ijk=ijkstart3,ijkend3
               do ii=-3,3
                  A_mt(ii,ijk) = A_m(ijk,ii,M)
               enddo
            ENDDO
            call leq_bicgst(VNAME, VNO, VAR, A_Mt(:,:), B_M(:,M), &
                            SWEEP, TOL, PC, ITMAX, IER)
            deallocate( A_mt )
         ELSE
           time1 = WALL_TIME()
!            call leq_bicgsmu(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
!                           SWEEP, TOL, PC, ITMAX, IER, NIT)

      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_owns(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         OAM = -ONE/A_M(IJK,0,M)
         A_M(IJK, 0,M) = -ONE
         A_M(IJK,-2,M) = A_M(IJK,-2,M)*OAM
         A_M(IJK,-1,M) = A_M(IJK,-1,M)*OAM
         A_M(IJK, 1,M) = A_M(IJK, 1,M)*OAM
         A_M(IJK, 2,M) = A_M(IJK, 2,M)*OAM
         A_M(IJK,-3,M) = A_M(IJK,-3,M)*OAM
         A_M(IJK, 3,M) = A_M(IJK, 3,M)*OAM
         B_M(IJK,M)    = B_M(IJK,M)*OAM
      ENDDO


	do k = kstart,kend
               do i = istart,iend
                  do j = jstart,jend
                     IJK = funijk(i,j,k)
                     IJK_GL= funijk_gl(i,j,k)
                     A_M(IJK,-3,M) =  0.166d0
                     A_M(IJK,-2,M) =  0.166d0  !0.1d0
                     A_M(IJK,-1,M) =  0.166d0
                     A_M(IJK,0,M)  = -1.0d0
                     A_M(IJK,1,M)  =  0.166d0
                     A_M(IJK,2,M)  =  0.166d0 ! 0.1d0
                     A_M(IJK,3,M)  =  0.166d0
		     B_M(IJK,M)    = A_M(IJK,-3,M)+A_M(IJK,-2,M)+A_M(IJK,-1,M)+ &    
                     A_M(IJK,0,M)+A_M(IJK,1,M)+A_M(IJK,2,M)+A_M(IJK,3,M)
                 end do 
               end do 
            end do 


	time1 = WALL_TIME()
	 call leq_bicgsmu(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
!	 call leq_bicgs(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                           SWEEP, TOL, PC, ITMAX, IER)

!	    call mpi_barrier(MPI_COMM_WORLD,mpierr)

	    do k = kstart3,kend3
               do i = istart3,iend3
                  do j = jstart3,jend3
                     IJK = funijk(i,j,k)
	             write(700+myPE,*) IJK,funijk_gl(I,J,K), VAR(IJK)
!	             write(700+myPE,*) IJK,VAR(IJK)
	          end do 
	       end do 
	    end do 

         ENDIF

      CASE (3)
! GMRES: A Generalized Minimal RESidual Algorithm
         time1 = WALL_TIME()
         call leq_gmres(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                        SWEEP, TOL, ITMAX, MAX_IT, IER)
         write(810,*) VNAME, method, WALL_TIME() - time1

      CASE (4)
! Mix:
         IER = 0
         call leq_bicgs(VNAME,VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP,&
                       TOL, PC, ITMAX, IER)
         IF (IER .eq. -2) THEN
            IER = 0
            print*,'calling leq_gmres', Vname
            call leq_gmres(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                           SWEEP, TOL, ITMAX, MAX_IT, IER)
         ENDIF


      CASE (5)
! CG: Conjugate Gradients
         time1 = WALL_TIME()
         call leq_cg(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP,&
                     TOL, PC, ITMAX, IER)
        
        write(810,*) VNAME, method, WALL_TIME() - time1

        CASE(6)

        if(NODESK.gt.1) then
        nstart = kstart
        nend   = kend
        if(myPE.eq.0)          nstart = kstart - 1
        if(myPE.eq.(numPEs-1)) nend   = kend   + 1

        if(numPEs.eq.1) then
        nstart = kstart
        nend   = kend
        end if

        dimred = funijk(iend,jend,nend) - funijk(istart,jstart,nstart)  + 1  !1

	else 
	dimred = DIMENSION_3
	end if 


	VAR(:) =  0.0d0 !B_M(:,M) 
          do k = kstart,kend
               do i = istart,iend
                  do j = jstart,jend
                     IJK = funijk(i,j,k)
                     IJK_GL= funijk_gl(i,j,k)

                     A_M(IJK,-3,M) =  0.166d0
                     A_M(IJK,-2,M) =  0.166d0  !0.1d0
                     A_M(IJK,-1,M) =  0.166d0
                     A_M(IJK,0,M)  = -1.0d0
                     A_M(IJK,1,M)  =  0.166d0
                     A_M(IJK,2,M)  =  0.166d0 ! 0.1d0
                     A_M(IJK,3,M)  =  0.166d0
		     B_M(IJK,M)    = A_M(IJK,-3,M)+A_M(IJK,-2,M)+A_M(IJK,-1,M)+ &    
                     A_M(IJK,0,M)+A_M(IJK,1,M)+A_M(IJK,2,M)+A_M(IJK,3,M)
                 end do 
               end do 
          end do 

	time1 = WALL_TIME()
	CALL MFIXinterfacemu(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), ITMAX,TOL,dimred)
	write(*,*) "Trilinos", VNAME, method, WALL_TIME() - time1
!	call mpi_barrier(MPI_COMM_WORLD,mpierr)

           do k = kstart3,kend3
               do i = istart3,iend3
                  do j = jstart3,jend3
                     IJK = funijk(i,j,k)
!                     write(700+myPE,*) IJK, VAR(IJK)
                     write(700+myPE,*) IJK, funijk_gl(i,j,k), VAR(IJK)
                  end do
               end do
            end do

!       CASE (6) !- Disabled
! LSOR: Line Successive Over Relaxation method


        write(810,*) VNAME, method, WALL_TIME() - time1 

      CASE DEFAULT
         LINE0(1:14) = 'SOLVE_LIN_EQ: '
         LINE0(15:80)= VName
         WRITE(LINE1,'(A, I2, A)') &
             'Error: LEQ_METHOD = ', METHOD, ' is invalid'
         CALL WRITE_ERROR(LINE0, LINE1, 1)
         CALL mfix_exit(myPE)
      END SELECT
! ----------------------------------------------------------------<<<

    
      RETURN
      END SUBROUTINE SOLVE_LIN_EQM

