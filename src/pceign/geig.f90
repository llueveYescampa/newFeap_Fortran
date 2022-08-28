subroutine geig(g,h,d,p,t,nev,prt)
implicit none
integer  nev
double precision   g(*),h(*),d(*),p(nev,*),t(*)
logical  prt

!  Purpose: Solve general eigenproblem 'g*p = h*p*d'

!  Inputs:
!     g(*)     - Projected subspace array
!     h(*)     - Projected subspace array
!     t(*)     - Working vector
!     nev      - Size of subspace problem
!     prt      - Flag, output results of each iteration if true

!  Outputs:
!     d(*)     - Eigenvalues of subspace problem
!     p(nev,*) - Eigenvectors of subspace problem


    integer  ir

!  Compute choleski factors of 'h'

    if(prt) then 
      call wprojm(g,nev,1)
      call wprojm(h,nev,2)
    end if  
    call chlfac(h,nev)

!  Compute standard eigenvalue problem matrix 'c'

    call chlfwd(h,g,p,nev)

!  Perform eignfunction decomposition of 'c'

    call eisql(g,d,t,p,nev,ir)

!  Compute vectors of original problem

    call chlbac(h,p,nev)

end
