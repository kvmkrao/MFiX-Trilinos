	subroutine FWrapper_ext(Bv, Var, aloc,   Am, lmap, nc, mr, nrows,maxit,tole)

!	subroutine FWrapper_ext(Bv, Var, Am,lglob,nc, mr, ntrows)
!       external cwrapper

	USE mpi_utility

!        include 'mpif.h'
!	integer myPE, numPEs,mpierr

	IMPLICIT NONE
!	integer myPE, numPEs,mpierr
        CHARACTER(len=7):: String
        integer :: nc,mr,i,nrows,j,maxit
!        real*8, INTENT(INOUT) :: Var(mr)
        real*8  :: Var(mr)
! Septadiagonal matrix A_m
!        real*8, INTENT(IN)  :: Am(mr,nc)
        real*8 :: Am(mr,nc)
! Vector b_m
!        real*8, INTENT(IN)  ::  Bv(mr)
        real*8 ::  Bv(mr)
!        real*8, INTENT(IN)  ::  aloc(mr,nc-1)
        real*8 ::  aloc(mr,nc-1)
!	integer, INTENT(IN) :: lmap(mr)
	integer :: lmap(mr)
	real*8  :: tole 

!       call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, mpierr )
!       call MPI_COMM_SIZE( MPI_COMM_WORLD, numPEs, mpierr )

        String = "HI"

!        Var(:)  = 0.0d0

!	write(*,*) "FWRAPPER_ext",  maxval(lmap), nrows, mr 
!	do i = 1, mr 
!	write(800+myPE,*) i,nc,lmap(i),(Am(i,j),j=1,7), Bv(i)
!	write(800+myPE,*) lmap(i), Bv(i)
!	end do 

!	write(*,*) "FWRAPPER", myPE, mr, nc, nrows,maxj
!	stop
        call cwrapper(Bv, Var, aloc, Am,lmap,nc, mr, nrows, maxit,tole, string)
	return
        end


      SUBROUTINE scaling(A_M, B_M)
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE leqsol
      USE cutcell
      USE functions
!      USE compar
      use machine, only: wall_time

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
        DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)

! Vector b_m
        DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
        DOUBLE PRECISION oAm
        DOUBLE PRECISION   :: aijmax

        INTEGER :: I, J, K, IJK, IJK_GL

!$omp parallel do default(shared) private(ijk,i,j,k,oam,aijmax)
            do k = kstart3,kend3
               do i = istart3,iend3
                  do j = jstart3,jend3
                     IJK = funijk(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     OAM = one/aijmax
                     A_M(IJK,:) = A_M(IJK,:)*OAM
                     B_M(IJK) = B_M(IJK)*OAM
                  end do
                end do
            end do

        return
        end


        subroutine xmlflt(check,string) ! , string)
        real*8 check
        integer length
        CHARACTER(LEN=20) :: string
!        string = "viscosity"
        length = LEN_TRIM(string)
!       write(*,*) "length", length
!       stop
        call xmlfltc(check,string,length)
        RETURN
        END


        subroutine xmlint(intvar,string) ! , string)
        integer intvar
        integer length
        CHARACTER(LEN=20) :: string
        length = LEN_TRIM(string)
        call xmlintc(intvar,string,length)
        RETURN
        END

        subroutine xmlstr(string1,string) ! , string)
        CHARACTER(LEN=50) :: string1
        integer length
        CHARACTER(LEN=20) :: string
        length = LEN_TRIM(string)
        call xmlstrc(string1,50, string,length)
        RETURN
        END

