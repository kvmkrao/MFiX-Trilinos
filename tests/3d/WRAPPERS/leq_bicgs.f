!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine LEQ_BICGS                                                C
!  Purpose: Solve system of linear system using BICGS method           C
!           Biconjugate gradients stabilized                           C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
! Handan Liu wrote below:               !Jan 22 2013
! The modification is as below:
!       Adding a loop of 2D RSRS sweep and parallelizing for OpenMP.
!       Splitting the existing 3D RSRS loop into two loops for OpenMP
!               due to data dependency
!       Adding openmp directives in all loops in leq_bicgs0.
!       Setting 'use_doloop = ture' to activate using do-loop.
!
      SUBROUTINE LEQ_BICGSmu(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                           TOL, PC, ITMAX, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE leqsol
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: Var
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3,-3:3), INTENT(INOUT) :: A_m
     DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3), INTENT(INOUT) :: A_m

! Vector b_m
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: B_m
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
! Note: this setting only seems to matter when leq_pc='line'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) ::  PC
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-------------------------------------------------
! Local Variables
!-------------------------------------------------

      if(PC.eq.'LINE') then   ! default
         call LEQ_BICGS0mu( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, IER )
      elseif(PC.eq.'DIAG') then
         call LEQ_BICGS0mu( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE1, IER )
      elseif(PC.eq.'NONE') then
         call LEQ_BICGS0mu( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE0, IER )
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) &
           'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_BICGSmu


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_BICGS0                                              C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_BICGS0mu(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                            TOL, ITMAX, MATVEC, MSOLVE, IER )

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE geometry
      USE compar
      USE mpi_utility
      USE sendrecv
      USE indices
      USE leqsol
      USE cutcell
      USE functions
      use machine, only: wall_time
	USE time_cpu
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments/procedure
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here-see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
!      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, INTENT(INOUT) :: A_m(ijkstart3:ijkend3,-3:3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3), INTENT(INOUT) :: A_m
! Vector b_m
!      DOUBLE PRECISION, INTENT(INOUT) :: B_m(ijkstart3:ijkend3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
! dummy arguments/procedures set as indicated
!     matvec->leq_matvec
! for preconditioner (leq_pc)
!    'line' msolve->leq_msolve  (default)
!    'diag' msolve->leq_msolve1
!    'none' msolve->leq_msolve0
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 0
      DOUBLE PRECISION, PARAMETER :: ratiotol = 0.2
      LOGICAL, PARAMETER :: do_unit_scaling = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION(:), allocatable :: R,Rtilde, P,Phat, Svec, Shat, Tvec,V

      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: &
                        alpha, beta, omega, rho
      DOUBLE PRECISION :: TxS, TxT, RtildexV, RtildexR, &
                          aijmax, oam
      DOUBLE PRECISION :: Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, j, k, ijk
      INTEGER :: iter
      DOUBLE PRECISION, DIMENSION(2) :: TxS_TxT
      DOUBLE PRECISION time1, time2, time21, time3, time4, time5,time6, time7
      DOUBLE PRECISION time0, time8, time9, time10, time11, time12,time01 
      LOGICAL            :: prt
      DOUBLE PRECISION timec1, timec2, timec01, timec3, timec4, timec5,timec6
      DOUBLE PRECISION timec0, timec8, timec9, timec10, timec7
      prt=.FALSE.
!-----------------------------------------------

	time0 = WALL_TIME()
	call CPU_TIME(timec0)

      allocate(R(DIMENSION_3))
      allocate(Rtilde(DIMENSION_3))
      allocate(P(DIMENSION_3))
      allocate(Phat(DIMENSION_3))
      allocate(Svec(DIMENSION_3))
      allocate(Shat(DIMENSION_3))
      allocate(Tvec(DIMENSION_3))
      allocate(V(DIMENSION_3))

      is_serial = numPEs.eq.1.and.is_serial

! these scalars should not be necessary to initialize but done as failsafe
      rnorm = ZERO
      rnorm0 = ZERO
      snorm = ZERO
      pnorm = ZERO

! initializing
      alpha(:)  = zero
      beta(:)   = zero
      omega(:)  = zero
      rho(:)    = zero

      use_doloop = .TRUE.   ! set true by Handan Liu for OpenMP do-loop

! Adding all by Handan Liu
! zero out R, Rtilde, P, Phat, Svec, Shat, Tvec, V
! --------------------------------
      if (use_doloop) then   ! mfix.dat keyword default=false
!$omp  parallel do default(shared) private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = zero
            Rtilde(ijk) = zero
            P(ijk) = zero
            Phat(ijk) = zero
            Svec(ijk) = zero
            Shat(ijk) = zero
            Tvec(ijk) = zero
            V(ijk) = zero
         enddo

      else
         R(:) = zero
         Rtilde(:) = zero
         P(:) = zero
         Phat(:) = zero
         Svec(:) = zero
         Shat(:) = zero
         Tvec(:) = zero
         V(:) = zero
      endif

      TOLMIN = EPSILON( one )

! Scale matrix to have unit diagonal
! ---------------------------------------------------------------->>>
      if (do_unit_scaling) then


         IF(RE_INDEXING) THEN  ! Loop only over active cells
!$omp parallel do default(shared) private(ijk,oam,aijmax)
            DO IJK = IJKSTART3,IJKEND3
               aijmax = maxval(abs(A_M(ijk,:)) )
               OAM = one/aijmax
               A_M(IJK,:) = A_M(IJK,:)*OAM
               B_M(IJK) = B_M(IJK)*OAM
            ENDDO

         ELSE

!$omp parallel do default(shared) private(ijk,i,j,k,oam,aijmax)
            do k = kstart2,kend2
               do i = istart2,iend2
                  do j = jstart2,jend2
                     IJK = funijk(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     OAM = one/aijmax
                     A_M(IJK,:) = A_M(IJK,:)*OAM
                     B_M(IJK) = B_M(IJK)*OAM
!	           write(90+myPE,999)IJK,funijk_gl(i,j,k),A_M(IJK,-3), A_M(IJK,-1), A_M(IJK,-2), & 
!                   A_M(IJK,0), A_M(IJK,2), A_M(IJK,1), A_M(IJK,3), B_M(IJK)
                  enddo
               enddo
            enddo
	 if(prt)write(*,*) "1. do_unit_scaling: no- reindexing"
         ENDIF

      endif

  999	format(2i6,8f10.5)

! ----------------------------------------------------------------<<<
! Compute initial residual (R = b-A*x) for Ax=b
!    assume initial guess in Var
!    rtilde = r
! ---------------------------------------------------------------->>>
      call MATVEC(Vname, Var, A_M, R)   ! returns R=A*Var

      if (use_doloop) then   ! mfix.dat keyword default=false
!$omp parallel do default(shared) private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = B_m(ijk) - R(ijk)
         enddo
      if(prt) write(*,*) "2. use_do loop"
      else
         R(:) = B_m(:) - R(:)
      endif

      call send_recv(R,nlayers_bicgs)

      if(is_serial) then
         Rnorm0 = zero
         if (use_doloop) then   ! mfix.dat keyword default=false
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm0)
            do ijk=ijkstart3,ijkend3
               Rnorm0 = Rnorm0 + R(ijk)*R(ijk)
            enddo
         else
            Rnorm0 = dot_product(R,R)
         endif
         Rnorm0 = sqrt( Rnorm0 )
      else
         Rnorm0 = sqrt( dot_product_par( R, R ) )
	if(prt) write(*,*) "3. parallel"
      endif

! determine an initial guess for the residual = residual + small random
! number (so it could be set to anything). note that since random_number
! is used to supply the guess, this line could potentially be the source
! of small differences between runs.  the random number is shifted below
! between -1 and 1 and then scaled by factor 1.0D-6*Rnorm0
      call random_number(Rtilde(:))

! Shift random number array to be consistent with case when RE_INDEXING is .FALSE.
       IF(RE_INDEXING) CALL SHIFT_DP_ARRAY(Rtilde)


      if (use_doloop) then   ! mfix.dat keyword default=false
!$omp parallel do default(shared) private(ijk)
         do ijk=ijkstart3,ijkend3
            Rtilde(ijk) = R(ijk) + (2.0d0*Rtilde(ijk)-1.0d0)*1.0d-6*Rnorm0
         enddo
	if(prt) write(*,*) "4. use_doloop"
      else
         Rtilde(:) = R(:) + (2.0d0*Rtilde(:)-1.0d0)*1.0d-6*Rnorm0
      endif

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
      endif
! ----------------------------------------------------------------<<<

	time01  = WALL_TIME()
	call CPU_TIME(timec01) 
	write(*,*) "Before main loop", time01-time0, timec01 - timec0


! Main loop
! ---------------------------------------------------------------->>>
      iter = 1
      do i=1,itmax

	time1  = WALL_TIME()
	call cpu_time(timec1) 
         if(is_serial) then
            if (use_doloop) then   ! mfix.dat keyword default=false
               RtildexR = zero
!$omp parallel do default(shared) private(ijk) reduction(+:RtildexR)
               do ijk=ijkstart3,ijkend3
                  RtildexR = RtildexR + Rtilde(ijk) * R(ijk)
               enddo
               rho(i-1) = RtildexR
            else
               rho(i-1) = dot_product( Rtilde, R )
            endif
         else
            rho(i-1) = dot_product_par( Rtilde, R )
           if(prt) write(*,*) "5  parallel "
         endif ! is_serial

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
               print*, 'leq_bicgs,',Vname,': rho(i-1) == 0 '
               ier = -2
	       stop
            else
! Method converged.  residual is already zero
! --------------------------------
               ier = 0
            endif
            call send_recv(var,2)
            return
         endif ! rho(i-1).eq.0

	time2  = WALL_TIME()
	call cpu_time(timec2) 
	if(myPE.eq.0) write(*,*) " MFIX: rho_i ", time2 - time1, timec2 - timec1 

         if (i .eq. 1) then
            if (use_doloop) then
!!$omp        parallel do private(ijk)
!$omp parallel do default(shared) private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = R(ijk)
               enddo
	    if(prt) write(*,*) "7. i = 1, use_doloop"
            else
               P(:) = R(:)
            endif
         else
            beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )

!           time21 = WALL_TIME()
!           write(*,*) " MFIX: beta ", time21 - time2

            if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = R(ijk) + beta(i-1)*( P(ijk) - omega(i-1)*V(ijk) )
               enddo
!	     write(*,*) "7 i /= 1 .no. use_doloop"
            else
               P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )
            endif
         endif ! i.eq.1

        time3 = WALL_TIME()
        call cpu_time(timec3) 
        if(myPE.eq.0) write(*,*) " MFIX: p_i ", time3 - time2, timec3 - timec2

! Solve A*Phat(:) = P(:)
! V(:) = A*Phat(:)
! --------------------------------
         call MSOLVE(Vname, P, A_m, Phat, CMETHOD) ! returns Phat
         call MATVEC(Vname, Phat, A_m, V)   ! returns V=A*Phat
	 if(prt) write(*,*) "8 msolve  matvec "

         if(is_serial) then
            if (use_doloop) then
               RtildexV = zero
!$omp  parallel do default(shared) private(ijk) reduction(+:RtildexV)
               do ijk=ijkstart3,ijkend3
                  RtildexV = RtildexV + Rtilde(ijk) * V(ijk)
               enddo
            else
               RtildexV = dot_product( Rtilde, V )
            endif
         else
            if(prt) write(*,*) "9 rtildexV"
            RtildexV = dot_product_par( Rtilde, V )
         endif ! is_serial

        time4 = WALL_TIME()
        call cpu_time(timec4) 
        if(myPE.eq.0) write(*,*) " MFIX: tilde v ", time4 - time3, timec4 - timec3

!         print*,'leq_bicgs, initial: ', Vname,' RtildexV ', RtildexV

! compute alpha
! --------------------------------
         alpha(i) = rho(i-1) / RtildexV
!        time5 = WALL_TIME()
!        write(*,*) " MFIX: alpha ", time5 - time4

! compute Svec
! --------------------------------
         if (use_doloop) then
!!$omp     parallel do private(ijk)
!$omp parallel do default(shared) private(ijk)
            do ijk=ijkstart3,ijkend3
               Svec(ijk) = R(ijk) - alpha(i) * V(ijk)
            enddo
	   if(prt) write(*,*) "10 svec "
!!$omp end parallel do
         else
            Svec(:) = R(:) - alpha(i) * V(:)
         endif ! use_doloop

	time6  = WALL_TIME()
	call cpu_time(timec6)
	if(myPE.eq.0) write(*,*) " MFIX: s ", time6 - time4, timec6 - timec4 

! Check norm of Svec(:); if small enough:
! set X(:) = X(:) + alpha(i)*Phat(:) and stop
! --------------------------------
         if(.not.minimize_dotproducts) then
            if(is_serial) then
               if (use_doloop) then
                  Snorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Snorm)
                  do ijk=ijkstart3,ijkend3
                     Snorm = Snorm + Svec(ijk) * Svec(ijk)
                  enddo
               else
                  Snorm = dot_product( Svec, Svec )
               endif
               Snorm = sqrt( Snorm )
            else
               if(prt) write(*,*) "11 snorm", Snorm
               Snorm = sqrt( dot_product_par( Svec, Svec ) )
            endif               ! is_serial
!            print*,'leq_bicgs, initial: ', Vname,' Snorm ', real(Snorm)


            if (Snorm <= TOLMIN) then
               if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
                  do ijk=ijkstart3,ijkend3
                     Var(ijk) = Var(ijk) + alpha(i)*Phat(ijk)
                  enddo
                  if(prt)   write(*,*) "12 Solution : Snorm <= tolmin"
               else
                  Var(:) = Var(:) + alpha(i)*Phat(:)
               endif            ! use_doloop

! Recompute residual norm
! --------------------------------
               if (idebugl >= 1) then
                  call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
!                  Rnorm = sqrt( dot_product_par( Var, Var ) )
!                  print*,'leq_bicgs, initial: ', Vname,' Vnorm ', Rnorm

                  if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
                     do ijk=ijkstart3,ijkend3
                        R(ijk) = B_m(ijk) - R(ijk)
                     enddo
		  if(prt) write(*,*) "13. idebugl >=1 & use_doloop"
                  else
                     R(:) = B_m(:) - R(:)
                  endif

                  if(is_serial) then
                     if (use_doloop) then
                        Rnorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm)
                        do ijk=ijkstart3,ijkend3
                           Rnorm = Rnorm + R(ijk)*R(ijk)
                        enddo
                     else
                        Rnorm =  dot_product( R, R )
                     endif
                     Rnorm = sqrt( Rnorm )
                  else
                     Rnorm = sqrt( dot_product_par( R, R ) )
		  if(prt)   write(*,*) "13_1. idebugl >=1 & parallel"
                  endif
!                  print*,'leq_bicgs, initial: ', Vname,' Rnorm ', Rnorm
               endif            ! idebugl >= 1
                           isConverged = .TRUE.
                           EXIT
            endif               ! end if (Snorm <= TOLMIN)
         endif                  ! end if (.not.minimize_dotproducts)



! Solve A*Shat(:) = Svec(:)
! Tvec(:) = A*Shat(:)
! --------------------------------
         call MSOLVE( Vname, Svec, A_m, Shat, CMETHOD)  ! returns Shat
         call MATVEC( Vname, Shat, A_m, Tvec )   ! returns Tvec=A*Shat
         if(prt) write(*,*) "14 msolve, matvec: svec, Shat"

         if(is_serial) then
            if (use_doloop) then
               TxS = zero
               TxT = zero
!$omp  parallel do default(shared) private(ijk) reduction(+:TxS,TxT)
               do ijk=ijkstart3,ijkend3
                  TxS = TxS + Tvec(ijk)  * Svec(ijk)
                  TxT = TxT + Tvec(ijk)  * Tvec(ijk)
               enddo
	      if(prt) write(*,*) "is_serial: use_doloop"
            else
               TxS = dot_product( Tvec, Svec )
               TxT = dot_product( Tvec, Tvec )
            endif
         else
            if(.not.minimize_dotproducts) then
               TxS = dot_product_par( Tvec, Svec )
               TxT = dot_product_par( Tvec, Tvec )
              if(prt)  write(*,*) "15 : Tvec Svec"
            else
               TxS_TxT = dot_product_par2(Tvec, Svec, Tvec, Tvec )
               TxS = TxS_TxT(1)
               TxT = TxS_TxT(2)
            endif
         endif

         IF(TxT.eq.Zero) TxT = SMALL_NUMBER

	time7 = WALL_TIME()
	call cpu_time(timec7) 
        if(myPE.eq.0) write(*,*) " MFIX: t ", time7 - time6, timec7 - timec6

! compute omega
! --------------------------------
         omega(i) = TxS / TxT
!	time8 = WALL_TIME()
!        write(*,*) " MFIX: omega_i ", time8 - time7

! compute new guess for Var
! --------------------------------
         if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
            do ijk=ijkstart3,ijkend3
               Var(ijk) = Var(ijk) +  &
                  alpha(i)*Phat(ijk) + omega(i)*Shat(ijk)
                  R(ijk) = Svec(ijk) - omega(i)*Tvec(ijk)
            enddo
	    if(prt) write(*,*) "16: Var"
!            if(prt) write(808,*) Var(:)
         else
            Var(:) = Var(:) +  &
               alpha(i)*Phat(:) + omega(i)*Shat(:)
               R(:) = Svec(:) - omega(i)*Tvec(:)
         endif

	time9 = WALL_TIME()
	call cpu_time(timec9) 
        if(myPE.eq.0) write(*,*) " MFIX: x_i, r_i",  time9 - time7, timec9 - timec7
        if(myPE.eq.0) write(*,*) " MFIX: total for", time9 - time1
        
! --------------------------------
         if(.not.minimize_dotproducts.or.(mod(iter,icheck_bicgs).eq.0)) then
            if(is_serial) then
               if (use_doloop) then
                  Rnorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm)
                  do ijk=ijkstart3,ijkend3
                     Rnorm = Rnorm + R(ijk) * R(ijk)
                  enddo
               else
                  Rnorm =  dot_product(R, R )
               endif
               Rnorm = sqrt( Rnorm )
            else
               Rnorm = sqrt( dot_product_par(R, R) )
               if(prt) write(*,*) "17 dot_product_par"
            endif               ! is_serial

            if (idebugl.ge.1) then
               if (myPE.eq.PE_IO) then
                  print*,'iter, Rnorm ', iter, Rnorm, Snorm
                  print*,'alpha(i), omega(i) ', alpha(i), omega(i)
                  print*,'TxS, TxT ', TxS, TxT
                  print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
               endif
            endif

!           call mfix_exit(myPE)

! Check convergence; continue if necessary
! for continuation, it is necessary that omega(i) .ne. 0
            isconverged = (Rnorm <= TOL*Rnorm0)

            if (isconverged) then
               iter_tot(vno) = iter_tot(vno) + iter + 1
               EXIT
            endif
         endif                  ! end if(.not.minimize_dotproducts)

! Advance the iteration count
         iter = iter + 1

      enddo   ! end do i=1,itmax
! end of linear solver loop
! ----------------------------------------------------------------<<<

      time10 = WALL_TIME()
      if(myPE.eq.0)write(*,*) " MFIX: end bicgstab",   time10 - time9, timec10 -timec9
      if(myPE.eq.0)write(*,*) " MFIX: total bicgstab", time10 - time01

	write(*,*) "MFIX: No of iterations" , iter,i

      if (idebugl >= 1) then
         call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
         if (use_doloop) then
!$omp  parallel do default(shared) private(ijk)
            do ijk=ijkstart3,ijkend3
               R(ijk) = R(ijk) - B_m(ijk)
            enddo
         else
            R(:) = R(:) - B_m(:)
         endif

         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm)
               do ijk=ijkstart3,ijkend3
                  Rnorm = Rnorm + R(ijk) * R(ijk)
               enddo
            else
               Rnorm = dot_product( R,R)
            endif
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R,R) )
         endif

         if(myPE.eq.0) print*,'leq_bicgs: final Rnorm ', Rnorm
         if(myPE.eq.0)  print*,'leq_bicgs ratio : ', Vname,' ',iter,  &
         ' L-2', Rnorm/Rnorm0
      endif   ! end if(idebugl >=1)

!      isconverged = (real(Rnorm) <= TOL*Rnorm0);
      if(.NOT.isConverged) isconverged = (real(Rnorm) <= TOL*Rnorm0);
!     write(*,*) '***',iter, isconverged, Rnorm, TOL, Rnorm0, myPE
      IER = 0
      if (.not.isconverged) then
         IER = -1
         iter_tot(vno) = iter_tot(vno) + iter
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif


!           do k = kstart2,kend2
!               do i = istart2,iend2
!                  do j = jstart2,jend2
!                    IJK = funijk(i,j,k)
!!                     aijmax = maxval(abs(A_M(ijk,:)) )
!!                     OAM = one/aijmax
!!                     A_M(IJK,:) = A_M(IJK,:)*OAM
!!                     B_M(IJK) = B_M(IJK)*OAM
!	           write(90+myPE,999)IJK,funijk_gl(i,j,k),A_M(IJK,-3), A_M(IJK,-1), A_M(IJK,-2), & 
!                   A_M(IJK,0), A_M(IJK,2), A_M(IJK,1), A_M(IJK,3), B_M(IJK)
!                  enddo
!               enddo
!            enddo


      call send_recv(var,2)

      deallocate(R)
      deallocate(Rtilde)
      deallocate(P)
      deallocate(Phat)
      deallocate(Svec)
      deallocate(Shat)
      deallocate(Tvec)
      deallocate(V)

        time11 = WALL_TIME()
        write(*,*) " call  bicgstab", time11 - time0

      return
      end subroutine LEQ_BICGS0mu

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine LEQ_BICGS                                                C
!  Purpose: Solve system of linear system using BICGS method           C
!           Biconjugate gradients stabilized                           C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
! Handan Liu wrote below:               !Jan 22 2013
! The modification is as below:
!       Adding a loop of 2D RSRS sweep and parallelizing for OpenMP.
!       Splitting the existing 3D RSRS loop into two loops for OpenMP
!               due to data dependency
!       Adding openmp directives in all loops in leq_bicgs0.
!       Setting 'use_doloop = ture' to activate using do-loop.
!
      SUBROUTINE LEQ_BICGS(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                           TOL, PC, ITMAX, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE leqsol
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: Var
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3,-3:3), INTENT(INOUT) :: A_m
     DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3), INTENT(INOUT) :: A_m

! Vector b_m
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: B_m
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
! Note: this setting only seems to matter when leq_pc='line'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) ::  PC
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!	integer :: NIT 
!-------------------------------------------------
! Local Variables
!-------------------------------------------------

      if(PC.eq.'LINE') then   ! default
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, IER)
      elseif(PC.eq.'DIAG') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE1, IER)
      elseif(PC.eq.'NONE') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE0, IER)
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) &
           'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_BICGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_BICGS0                                              C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_BICGS0(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                            TOL, ITMAX, MATVEC, MSOLVE, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE geometry
      USE compar
      USE mpi_utility
      USE sendrecv
      USE indices
      USE leqsol
      USE cutcell
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments/procedure
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here-see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
!      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, INTENT(INOUT) :: A_m(ijkstart3:ijkend3,-3:3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3), INTENT(INOUT) :: A_m
! Vector b_m
!      DOUBLE PRECISION, INTENT(INOUT) :: B_m(ijkstart3:ijkend3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
! dummy arguments/procedures set as indicated
!     matvec->leq_matvec
! for preconditioner (leq_pc)
!    'line' msolve->leq_msolve  (default)
!    'diag' msolve->leq_msolve1
!    'none' msolve->leq_msolve0
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 0
      DOUBLE PRECISION, PARAMETER :: ratiotol = 0.2
      LOGICAL, PARAMETER :: do_unit_scaling = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION(:), allocatable :: R,Rtilde, P,Phat, Svec, Shat, Tvec,V

      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: &
                        alpha, beta, omega, rho
      DOUBLE PRECISION :: TxS, TxT, RtildexV, RtildexR, &
                          aijmax, oam
      DOUBLE PRECISION :: Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, j, k, ijk, IJK_GL
      INTEGER :: iter
      DOUBLE PRECISION, DIMENSION(2) :: TxS_TxT
	INTEGER :: NITR
!-----------------------------------------------

      allocate(R(DIMENSION_3))
      allocate(Rtilde(DIMENSION_3))
      allocate(P(DIMENSION_3))
      allocate(Phat(DIMENSION_3))
      allocate(Svec(DIMENSION_3))
      allocate(Shat(DIMENSION_3))
      allocate(Tvec(DIMENSION_3))
      allocate(V(DIMENSION_3))

      is_serial = numPEs.eq.1.and.is_serial

! these scalars should not be necessary to initialize but done as failsafe
      rnorm = ZERO
      rnorm0 = ZERO
      snorm = ZERO
      pnorm = ZERO

! initializing
      alpha(:)  = zero
      beta(:)   = zero
      omega(:)  = zero
      rho(:)    = zero

          use_doloop = .TRUE.   ! set true by Handan Liu for OpenMP do-loop

! Adding all by Handan Liu
! zero out R, Rtilde, P, Phat, Svec, Shat, Tvec, V
! --------------------------------
      if (use_doloop) then   ! mfix.dat keyword default=false
!$omp  parallel do default(shared) private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = zero
            Rtilde(ijk) = zero
            P(ijk) = zero
            Phat(ijk) = zero
            Svec(ijk) = zero
            Shat(ijk) = zero
            Tvec(ijk) = zero
            V(ijk) = zero
         enddo

      else
         R(:) = zero
         Rtilde(:) = zero
         P(:) = zero
         Phat(:) = zero
         Svec(:) = zero
         Shat(:) = zero
         Tvec(:) = zero
         V(:) = zero
      endif

      TOLMIN = EPSILON( one )

! Scale matrix to have unit diagonal
! ---------------------------------------------------------------->>>
      if (do_unit_scaling) then


         IF(RE_INDEXING) THEN  ! Loop only over active cells
!$omp parallel do default(shared) private(ijk,oam,aijmax)
            DO IJK = IJKSTART3,IJKEND3
               aijmax = maxval(abs(A_M(ijk,:)) )
               OAM = one/aijmax
               A_M(IJK,:) = A_M(IJK,:)*OAM
               B_M(IJK) = B_M(IJK)*OAM
            ENDDO

         ELSE

!$omp parallel do default(shared) private(ijk,ijk_gl,i,j,k,oam,aijmax)
            do k = kstart2,kend2
               do i = istart2,iend2
                  do j = jstart2,jend2
                     IJK = funijk(i,j,k)
		     IJK_GL= funijk_gl(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     OAM = one/aijmax
                     A_M(IJK,:) = A_M(IJK,:)*OAM
                     B_M(IJK) = B_M(IJK)*OAM
!		if(NITR.eq.2) write(80+myPE,*) 'leq_bicgs.f', mypE, i,j,k,ijk,ijk_gl,(A_M(ijk,:)), B_M(ijk)
                  enddo
               enddo
            enddo

         ENDIF

      endif


! ----------------------------------------------------------------<<<


! Compute initial residual (R = b-A*x) for Ax=b
!    assume initial guess in Var
!    rtilde = r
! ---------------------------------------------------------------->>>
      call MATVEC(Vname, Var, A_M, R)   ! returns R=A*Var

      if (use_doloop) then   ! mfix.dat keyword default=false
!$omp parallel do default(shared) private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = B_m(ijk) - R(ijk)
         enddo
      else
         R(:) = B_m(:) - R(:)
      endif

      call send_recv(R,nlayers_bicgs)

      if(is_serial) then
         Rnorm0 = zero
         if (use_doloop) then   ! mfix.dat keyword default=false
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm0)
            do ijk=ijkstart3,ijkend3
               Rnorm0 = Rnorm0 + R(ijk)*R(ijk)
            enddo
         else
            Rnorm0 = dot_product(R,R)
         endif
         Rnorm0 = sqrt( Rnorm0 )
      else
         Rnorm0 = sqrt( dot_product_par( R, R ) )
      endif

! determine an initial guess for the residual = residual + small random
! number (so it could be set to anything). note that since random_number
! is used to supply the guess, this line could potentially be the source
! of small differences between runs.  the random number is shifted below
! between -1 and 1 and then scaled by factor 1.0D-6*Rnorm0
      call random_number(Rtilde(:))

! Shift random number array to be consistent with case when RE_INDEXING is .FALSE.
       IF(RE_INDEXING) CALL SHIFT_DP_ARRAY(Rtilde)


      if (use_doloop) then   ! mfix.dat keyword default=false
!$omp parallel do default(shared) private(ijk)
         do ijk=ijkstart3,ijkend3
            Rtilde(ijk) = R(ijk) + (2.0d0*Rtilde(ijk)-1.0d0)*1.0d-6*Rnorm0
         enddo
      else
         Rtilde(:) = R(:) + (2.0d0*Rtilde(:)-1.0d0)*1.0d-6*Rnorm0
      endif

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
      endif
! ----------------------------------------------------------------<<<


! Main loop
! ---------------------------------------------------------------->>>
      iter = 1
      do i=1,itmax

         if(is_serial) then
            if (use_doloop) then   ! mfix.dat keyword default=false
               RtildexR = zero
!$omp parallel do default(shared) private(ijk) reduction(+:RtildexR)
               do ijk=ijkstart3,ijkend3
                  RtildexR = RtildexR + Rtilde(ijk) * R(ijk)
               enddo
               rho(i-1) = RtildexR
            else
               rho(i-1) = dot_product( Rtilde, R )
            endif
         else
            rho(i-1) = dot_product_par( Rtilde, R )
         endif ! is_serial

!         print*,'leq_bicgs, initial: ', Vname,' rho(i-1) ', rho(i-1)

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
!               print*, 'leq_bicgs,',Vname,': rho(i-1) == 0 '
               ier = -2
            else
! Method converged.  residual is already zero
! --------------------------------
               ier = 0
            endif
            call send_recv(var,2)
            return
         endif ! rho(i-1).eq.0

         if (i .eq. 1) then
            if (use_doloop) then
!!$omp        parallel do private(ijk)
!$omp parallel do default(shared) private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = R(ijk)
               enddo
            else
               P(:) = R(:)
            endif
         else
            beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )
            if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = R(ijk) + beta(i-1)*( P(ijk) - omega(i-1)*V(ijk) )
               enddo
            else
               P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )
            endif
         endif ! i.eq.1


! Solve A*Phat(:) = P(:)
! V(:) = A*Phat(:)
! --------------------------------
         call MSOLVE(Vname, P, A_m, Phat, CMETHOD) ! returns Phat
         call MATVEC(Vname, Phat, A_m, V)   ! returns V=A*Phat

         if(is_serial) then
            if (use_doloop) then
               RtildexV = zero
!$omp  parallel do default(shared) private(ijk) reduction(+:RtildexV)
               do ijk=ijkstart3,ijkend3
                  RtildexV = RtildexV + Rtilde(ijk) * V(ijk)
               enddo
            else
               RtildexV = dot_product( Rtilde, V )
            endif
         else
            RtildexV = dot_product_par( Rtilde, V )
         endif ! is_serial

!         print*,'leq_bicgs, initial: ', Vname,' RtildexV ', RtildexV

! compute alpha
! --------------------------------
         alpha(i) = rho(i-1) / RtildexV

! compute Svec
! --------------------------------
         if (use_doloop) then
!!$omp     parallel do private(ijk)
!$omp parallel do default(shared) private(ijk)
            do ijk=ijkstart3,ijkend3
               Svec(ijk) = R(ijk) - alpha(i) * V(ijk)
            enddo
!!$omp end parallel do
         else
            Svec(:) = R(:) - alpha(i) * V(:)
         endif ! use_doloop


! Check norm of Svec(:); if small enough:
! set X(:) = X(:) + alpha(i)*Phat(:) and stop
! --------------------------------
         if(.not.minimize_dotproducts) then
            if(is_serial) then
               if (use_doloop) then
                  Snorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Snorm)
                  do ijk=ijkstart3,ijkend3
                     Snorm = Snorm + Svec(ijk) * Svec(ijk)
                  enddo
               else
                  Snorm = dot_product( Svec, Svec )
               endif
               Snorm = sqrt( Snorm )
            else
               Snorm = sqrt( dot_product_par( Svec, Svec ) )
            endif               ! is_serial
!            print*,'leq_bicgs, initial: ', Vname,' Snorm ', real(Snorm)


            if (Snorm <= TOLMIN) then
               if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
                  do ijk=ijkstart3,ijkend3
                     Var(ijk) = Var(ijk) + alpha(i)*Phat(ijk)
                  enddo
               else
                  Var(:) = Var(:) + alpha(i)*Phat(:)
               endif            ! use_doloop

! Recompute residual norm
! --------------------------------
               if (idebugl >= 1) then
                  call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
!                  Rnorm = sqrt( dot_product_par( Var, Var ) )
!                  print*,'leq_bicgs, initial: ', Vname,' Vnorm ', Rnorm

                  if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
                     do ijk=ijkstart3,ijkend3
                        R(ijk) = B_m(ijk) - R(ijk)
                     enddo
                  else
                     R(:) = B_m(:) - R(:)
                  endif

                  if(is_serial) then
                     if (use_doloop) then
                        Rnorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm)
                        do ijk=ijkstart3,ijkend3
                           Rnorm = Rnorm + R(ijk)*R(ijk)
                        enddo
                     else
                        Rnorm =  dot_product( R, R )
                     endif
                     Rnorm = sqrt( Rnorm )
                  else
                     Rnorm = sqrt( dot_product_par( R, R ) )
                  endif
!                  print*,'leq_bicgs, initial: ', Vname,' Rnorm ', Rnorm
               endif            ! idebugl >= 1
                           isConverged = .TRUE.
                           EXIT
            endif               ! end if (Snorm <= TOLMIN)
         endif                  ! end if (.not.minimize_dotproducts)



! Solve A*Shat(:) = Svec(:)
! Tvec(:) = A*Shat(:)
! --------------------------------
         call MSOLVE( Vname, Svec, A_m, Shat, CMETHOD)  ! returns Shat
         call MATVEC( Vname, Shat, A_m, Tvec )   ! returns Tvec=A*Shat

         if(is_serial) then
            if (use_doloop) then
               TxS = zero
               TxT = zero
!$omp  parallel do default(shared) private(ijk) reduction(+:TxS,TxT)
               do ijk=ijkstart3,ijkend3
                  TxS = TxS + Tvec(ijk)  * Svec(ijk)
                  TxT = TxT + Tvec(ijk)  * Tvec(ijk)
               enddo
            else
               TxS = dot_product( Tvec, Svec )
               TxT = dot_product( Tvec, Tvec )
            endif
         else
            if(.not.minimize_dotproducts) then
               TxS = dot_product_par( Tvec, Svec )
               TxT = dot_product_par( Tvec, Tvec )
            else
               TxS_TxT = dot_product_par2(Tvec, Svec, Tvec, Tvec )
               TxS = TxS_TxT(1)
               TxT = TxS_TxT(2)
            endif
         endif

         IF(TxT.eq.Zero) TxT = SMALL_NUMBER

! compute omega
! --------------------------------
         omega(i) = TxS / TxT

! compute new guess for Var
! --------------------------------
         if (use_doloop) then
!$omp parallel do default(shared) private(ijk)
            do ijk=ijkstart3,ijkend3
               Var(ijk) = Var(ijk) +  &
                  alpha(i)*Phat(ijk) + omega(i)*Shat(ijk)
                  R(ijk) = Svec(ijk) - omega(i)*Tvec(ijk)
            enddo
         else
            Var(:) = Var(:) +  &
               alpha(i)*Phat(:) + omega(i)*Shat(:)
               R(:) = Svec(:) - omega(i)*Tvec(:)
         endif


! --------------------------------
         if(.not.minimize_dotproducts.or.(mod(iter,icheck_bicgs).eq.0)) then
            if(is_serial) then
               if (use_doloop) then
                  Rnorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm)
                  do ijk=ijkstart3,ijkend3
                     Rnorm = Rnorm + R(ijk) * R(ijk)
                  enddo
               else
                  Rnorm =  dot_product(R, R )
               endif
               Rnorm = sqrt( Rnorm )
            else
               Rnorm = sqrt( dot_product_par(R, R) )
            endif               ! is_serial

            if (idebugl.ge.1) then
               if (myPE.eq.PE_IO) then
                  print*,'iter, Rnorm ', iter, Rnorm, Snorm
                  print*,'alpha(i), omega(i) ', alpha(i), omega(i)
                  print*,'TxS, TxT ', TxS, TxT
                  print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
               endif
            endif

!           call mfix_exit(myPE)

! Check convergence; continue if necessary
! for continuation, it is necessary that omega(i) .ne. 0
            isconverged = (Rnorm <= TOL*Rnorm0)

            if (isconverged) then
               iter_tot(vno) = iter_tot(vno) + iter + 1
               EXIT
            endif
         endif                  ! end if(.not.minimize_dotproducts)

! Advance the iteration count
         iter = iter + 1

      enddo   ! end do i=1,itmax
! end of linear solver loop
! ----------------------------------------------------------------<<<

!	if(myPE.eq.0) print*,'leq_bicgs: final Rnorm ', Rnorm
        if(myPE.eq.0)  print*,'leq_bicgs ratio : ', Vname,' ',iter,  &
         ' L-2', Rnorm/Rnorm0, tol


      if (idebugl >= 1) then
         call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
         if (use_doloop) then
!$omp  parallel do default(shared) private(ijk)
            do ijk=ijkstart3,ijkend3
               R(ijk) = R(ijk) - B_m(ijk)
            enddo
         else
            R(:) = R(:) - B_m(:)
         endif

         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!$omp parallel do default(shared) private(ijk) reduction(+:Rnorm)
               do ijk=ijkstart3,ijkend3
                  Rnorm = Rnorm + R(ijk) * R(ijk)
               enddo
            else
               Rnorm = dot_product( R,R)
            endif
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R,R) )
         endif

      endif   ! end if(idebugl >=1)

!      isconverged = (real(Rnorm) <= TOL*Rnorm0);
      if(.NOT.isConverged) isconverged = (real(Rnorm) <= TOL*Rnorm0);
!     write(*,*) '***',iter, isconverged, Rnorm, TOL, Rnorm0, myPE
      IER = 0
      if (.not.isconverged) then
         IER = -1
         iter_tot(vno) = iter_tot(vno) + iter
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif

      call send_recv(var,2)

  800	continue

      deallocate(R)
      deallocate(Rtilde)
      deallocate(P)
      deallocate(Phat)
      deallocate(Svec)
      deallocate(Shat)
      deallocate(Tvec)
      deallocate(V)

      return
      end subroutine LEQ_BICGS0

