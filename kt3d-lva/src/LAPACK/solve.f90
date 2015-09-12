!  solve.f90 
!
!  FUNCTIONS/SUBROUTINES exported from solve.dll:
!  solve      - AX = B
!
subroutine solve ( A, X, B, N, M )
    implicit none
  ! Variables
    integer, intent (in) :: N, M !Rank and Solutions

    real*8, dimension(N,N), intent (inout) :: A
    real*8, dimension(N,M), intent (inout) :: X
    real*8, dimension(N,M), intent (in   ) :: B
    integer, dimension(N) :: IPIV

    integer :: LDA, LDB !Leading dimension of A and B
    integer :: INFO  !return state of the program
    integer :: NRHS !number of righthandsides

    ! Body of solve
    LDA = N
    LDB = N
    NRHS = M
    X = B
    CALL DGESV( N, NRHS, A, LDA, IPIV, X, LDB, INFO ) 
end subroutine solve
