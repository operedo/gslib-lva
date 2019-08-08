
module error_report
    implicit none
    
    integer, public :: test1 !memory allocation test1
    
  contains
  
    !report a general error and terminate execution
    subroutine gen_error ( report, routine )
        character*(*), intent (in) :: report, routine
        write(*,10) trim(report), trim(routine)
    10  format(/,'ERROR: ',a,' IN ',a,/)
        STOP ; Return
    end subroutine gen_error
    
    !report an out of memory error
    subroutine out_of_mem ( routine )
        character*(*), intent (in) :: routine
        write(*,10) trim(routine)
    10  format(/,'ERROR OUT-OF-MEMORY IN ',a,/)
        STOP ; Return
    end subroutine out_of_mem
        
end module error_report


!same as johns, but ignores all the heassian stuff



module class_lvaf2
    use sort
    use grid_info
    use error_report
!    use group
    implicit none
    public :: dsqrd, &      !calculate anisotropic distance between 2 points
              fix_angles, & !correct angles
              get_rmatrix   !acquire the rotation matrix for a grid block
    
    
    !PRIVATE PARAMETERS
    real(kind=8), private, parameter :: PI      = 3.14159265358979D+00
    real(kind=8), private, parameter :: DEG2RAD = PI / 180.0D+00
    real(kind=8), private, parameter :: eps = 1.0e-4
    integer, private, parameter :: max_cp = 7 ! max_cp = 2^7
    integer, private, parameter :: dbg = 9 !debug file handle
    character*13, private, parameter :: debug_file = 'min_dsqrd.dbg' !debug file name
    
    
    !PRIVATE VARIABLES
    real(kind=8), private, allocatable, dimension (:,:) :: dsl !line segmented by a grid
    integer,      private, allocatable, dimension (:,:) :: bsl !block ids of segments in dsl
    integer,      private, allocatable, dimension ( : ) :: nsl !number of small segments per large segment
    logical, private :: debug_mode = .False. !debugging output, default if off
    
  contains
  
    subroutine debug_mode_on() !Turn on debug mode for this module
        logical :: isopen
        debug_mode = .True.
        inquire(dbg, opened = isopen)
        if(isopen) Return
        open(dbg, file = debug_file, status = 'unknown')
        write(dbg,10)
    10  format('Debug output for min_dsqrd in class_lvaf',//)
    end subroutine debug_mode_on
    subroutine debug_mode_off() !Turn off debug mode for this module
        logical :: isclosed
        debug_mode = .False.
        inquire(dbg, opened = isclosed)
        if(isclosed) return
        close(dbg)
    end subroutine debug_mode_off
  
    !Find the minimum anisotropic distance between two points.  This involves
    !starting with one internal point and refining the line after finding each
    !solution.
    real(kind=8) function min_dsqrd ( a, b ) result ( D )
        real(kind=8), dimension ( 3 ), intent (in) :: a, b !line endpoints
        real(kind=8), allocatable, dimension (:,:) :: L1, L2 !lines
        real(kind=8) :: Dn, D0
        integer :: i, j, k, n, test1
        
        !Determine a maximum number of line refinements to use
        n = block_count(a, b)
        if(n < 0) then
        call gen_error ( 'line endpoint beyond grid spec', 'min_dsqrd' )
        end if
        n = max( 1, min( int( alog(real(n)) / alog(2.0) ), max_cp ) )
        
        D = eps + 1.0 ; Dn = 0.0
        !Incrementally refine L and find the minimum anisotropic distance.
        j = 2
        allocate(L1(2,3), stat = test1)
        !if(test1) call out_of_mem ( 'min_dsqrd _1' )
        L1(1,:) = a ; L1(2,:) = b
        i = 0
        do while ( i <= n .and. D - Dn > eps )
            i = i + 1       !Increment counter
            j = j * 2 - 1   !Calculate # of control-points
            !Refine L2 from L1 (add mid-points)
            allocate(L2(j,3), stat = test1)
            !if(test1) call out_of_mem ( 'min_dsqrd _2' )
            do k = 1, j, 2 !Fill L2 with L1 and interpolated points
                L2(k,:) = L1((k + 1) / 2,:)
                if(k < j - 1)then
                    L2(k+1,:) = (L2(k,:) + L1((k + 3) / 2,:)) / 2.0
                endif
            enddo
            deallocate(L1)
            Dn = min_dsqrd_ ( L2, j, D ) !Minimize distance of L2
            if(i == 1) D0 = D
            !Replace L1 with L2 and continue refining
            allocate(L1(j,3), stat = test1)
            !if(test1) call out_of_mem ( 'min_dsqrd _3' )
            L1 = L2
            deallocate(L2)
        enddo
        !If debugging, write the line to file
        if(debug_mode) call dbg_print_path ( L1, D0, Dn )
        if(allocated(L1)) deallocate(L1)
    end function min_dsqrd
  
    !Minimize the distance along a line defined by n verices that exists
    !within an LVAF.  Return the covariance between the two points
    !Process: evaluate the gradient and hessian for the line, then
    !         use Newtons method to minimize the distance
    real(kind=8) function min_dsqrd_ ( L, n, D0 ) result ( D )
        use global
        integer, intent (in) :: n !# of vertices
        real(kind=8), dimension (n,3), intent (inout) :: L !line
        real(kind=8), intent (out) :: D0 !Initial distance
        !INTERNALS
        real(kind=8) :: t, Dt !Newton step parameter and distance for that step
        integer :: i, j
        logical :: outside,err
        
        !Initialization
        allocate(nsl( n ), stat = test1)
        !if(test1) call out_of_mem ( 'min_dsqrd  _1' )
        D = 0.0

        !Enter the Newton method process
        j = 1
            !segment the line

            call split_line ( L, n, outside )
            !Evaluate the distance, Jacobian and Hessian
            
            call fnc_eval ( n, L, D,err )
            
            n_seg=n_seg+1         
            if(max_seg>-1 .and. n_seg>max_seg  ) then
            
                D0 = D
            end if         
               
            if(err) then
                D0 = D
            end if
            
            
            if(j == 1) D0 = D
        deallocate(nsl)
    end function min_dsqrd_
    
    !Evaluate the distance, Jacobian and Hessian for a line
    !with n vertices.  The line must already have been segmented by
    !the grid and stored in dsl, bsl and nsl.
    subroutine fnc_eval ( n, Line, F, error )
        use global
        integer,               intent (in)  :: n !# of vertices
        real(kind=8),          intent (out) :: F !Function (distance)
        real(kind=8), dimension (n,3), intent (in) :: Line !line
        !INTERNALS
        real(kind=8), dimension (3,3) :: R !Rotation matrix
        real(kind=8), dimension ( 3 ) :: angles
        real(kind=8), dimension ( 2 ) :: ratios
        real(kind=8) :: t, dj
        integer :: i, j, k, l
        logical error
        error=.false.
        
        !Calculate the objective for each vertex in L using all segments
        !between L(i) and L(i-1)
        F  = 0.0  
        do i = 1, n - 1 !Loop over the large segments in Line
            l = (i - 1) * 3 + 1
            k = i * 3
            do j = nsl(i), nsl(i+1) - 1 !Loop over the small segments in dsl
                !set angles, anisotropic ratios and rotation matrix
                !for this small segment.
                angles = grid( bsl(j,1), bsl(j,2), bsl(j,3), 1:3 )
                ratios = grid( bsl(j,1), bsl(j,2), bsl(j,3), 4:5 )
                !call fix_angles  ( angles )
                call get_rmatrix ( angles, ratios, R )
                !Accumulate the distance
                dj = dsqrt( dsqrd_vect ( dsl(j,:), R ) )
                F = F + dj
            enddo
        enddo
        !Normalize the jacobian
        Return
    end subroutine fnc_eval
    
    !Create a new line that is defined by all intersections
    !of the current line with the underlying grid.  The new line
    !is stored in module variables dsl, bsl and nsl
    subroutine split_line ( L, n, outside )
        integer, intent (in) :: n
        real(kind=8), dimension (n,3), intent (in) :: L
        logical, intent (out) :: outside !line beyond grid indicator
        !INTERNALS
        integer :: i, j, k, nsa, nst
        integer, dimension ( 3 ) :: ida, idb
        real(kind=8), pointer, dimension (:,:) :: segs, tempverts
        integer,      pointer, dimension (:,:) :: segids, tempids
        
        !get the block indices of all vertices to count the
        !total number of intersections.  note that this is a max
        !on intersections.  If a vertex is outside the grid,
        !return with an error
        nsa = n
        do i = 2, n
            call get_index ( L( i ,:), ida, outside ) ; if(outside) Return
            call get_index ( L(i-1,:), idb, outside ) ; if(outside) Return
            nsa = nsa + sum(abs(ida - idb))
        enddo

        if(associated(tempverts)) deallocate(tempverts)
        if(associated(tempids))   deallocate(tempids)
        allocate(tempverts(nsa,3), tempids(nsa-1,3), stat = test1)
        !if(test1) call out_of_mem ( 'segment_lineset  _1' )
        
        !For each segment, calculate the intersection(s) and
        !add them to a new line segment
        tempverts(1,:) = L(1,:)
        nsl(1) = 1
        k = 1
        do i = 2, n
            !segment the line v(i)v(i-1)
            call split_segment ( L(i-1,:), L(i,:), nst, segs, segids, outside )
            
            !if there were any intersections, add them to tempverts
            if(nst > 2)then
                !store the indices of all small segments in the 
                !larger segment ( L(i-1) to L(i) )
                nsl(i) = nsl(i - 1) + nst - 1
                do j = 2, nst - 1
                    k = k + 1
                    tempverts(k,:) = segs(j,:)
                    tempids(k-1,:) = segids(j-1,:)
                enddo
            else
                nsl(i) = nsl(i - 1) + 1
            endif
            k = k + 1
            !Add the original point into tempverts
            tempverts(k,:) = L(i,:)
            !get the block index
            call get_index ( (tempverts(k,:) + tempverts(k-1,:)) / 2.0, tempids(k-1,:), outside )
        enddo
        nsa = k
        
        !allocate the new line and store segment vectors and block ids
        if(allocated(dsl)) deallocate(dsl)
        if(allocated(bsl)) deallocate(bsl)
        
        allocate(dsl(nsa-1,3), bsl(nsa-1,3), stat = test1)
        !if(test1) call out_of_mem ( 'segment_lineset  _2' )
        dsl = tempverts(2:nsa,:) - tempverts(1:nsa-1,:)
        bsl = tempids
        deallocate(tempverts, tempids)
        Return
    end subroutine split_line
    
    !Split up a line segment defined by endpoints a and b based on a
    !regular grid definition.  Return the resulting line and block
    !ids for each segment in xps and ids.
    subroutine split_segment (a, b, ns, xps, ids, outside)
        real(kind=8), dimension ( 3 ), intent (in)  :: a, b     !Endpoints of line
        real(kind=8), pointer :: xps(:,:) !segment endpoints
        integer, pointer :: ids(:,:) !block indexes
        integer, intent (out) :: ns       !number of points along line ab
        logical, intent (out) :: outside !have we stepped outside the grid
        
        !INTERNALS
        integer,                   dimension ( 3 ) :: ida, idb !blocks a and b are in
        real(kind=8),              dimension ( 3 ) :: dab !b - a
        real(kind=8), allocatable, dimension ( : ) :: t !parametric line intervals
        real(kind=8) :: x !intersection point, stores x, y or z values
        integer      :: i, j, k
        
        !First determine which blocks a and b are in.  Note that if a or b are beyond
        !the grid, the problem exists outside of this routine.  It will be detected and
        !the program terminated!
        outside = .False.
        call get_index (a, ida, outside)
        call get_index (b, idb, outside)
        do i = 1, 3
            if(ida(i) < 1 .or. ida(i) > n(i) .or. &
               idb(i) < 1 .or. idb(i) > n(i))then
                !call gen_error ( 'line beyond grid', 'split_line' )
                outside = .True. ; Return
            endif
        enddo
        dab = b - a !setup parametric version of line ab
        
        !calculate number of planes between a and b in x, y and z directions by
        !using the block indices.  This will be a maximum on potential intersections.
        !Include the endpoints ( + 2 )
        ns = sum(abs(ida - idb)) + 2
        allocate(t(ns), stat = test1)
        !if(test1) call out_of_mem ( 'split_line  _1' )
        
        !calculate values in t along each coordinate.  Add location 0 to the
        !start and location 1 to the end.  If line endpoints match a plane
        !exactly, ignore those points.
        j = 1
        t(1) = 0.0
        do k = 1, 3
            if(ida(k) /= idb(k))then !kth-direction
                do i = min(ida(k),idb(k)),max(ida(k),idb(k)) - 1
                    j = j + 1
                    x = ( o(k) - sss(k)/2.0 ) + real(i) * sss(k)
                    t(j) = ( x - a(k) ) / ( b(k) - a(k) )
                    if(abs(t(j)) <= eps .or. abs(t(j) - 1.0) <= eps)then
                        j = j - 1
                        t(ns-1) = -1.0
                        ns = ns - 1
                    endif
                enddo
            endif
        enddo
        t(ns) = 1.0
        
        !Now sort values in t and check for repeated intersections.
        call QUICK_SORT(t, 2, ns-1)
        
        !Check for repeated intersection points.
        i = 0
        do
            i = i + 1
            if(i == ns - 1) exit
            if(abs(t(i) - t(i+1)) <= eps)then
                t(i:ns-1) = t(i+1:ns)
                t(ns) = -1.0
                ns = ns - 1
                i = i - 1
            endif
        enddo
        
        !Allocate memory for line segment set.  Deallocate memory that
        !possibly already exists.
        if(associated(xps)) deallocate(xps)
        if(associated(ids)) deallocate(ids)
        allocate(xps(ns,3), ids(ns-1,3), stat = test1)
        !if(test1) call out_of_mem ( 'split_line  _2' )
        
        !Calculate intersections and block indexes
        do i = 1, ns
            xps(i,:) = dab * t(i) + a !parametric equation of line
            if(i > 1)then
                call get_index ( (xps(i,:) + xps(i-1,:)) / 2.0, ids(i-1,:), outside )
            endif
        enddo
        deallocate( t )
        Return
    end subroutine split_segment
    
    !Calculate squared anisotropic distance between two points.  It is up
    !to the user to ensure that the rotation matrix R is applicatble to
    !the line ab for locally varying anisotropy fields (LVAF's)
    real(kind=8) function dsqrd ( a, b, R ) result ( d )
        real(kind=8), dimension ( 3 ), intent (in) :: a, b !Endpoints
        real(kind=8), dimension (:,:), intent (in) :: R    !Anisotropic rotation matrix
        real(kind=8), dimension ( 3 )      :: ba   !Vector b - a
        ba = matmul ( R, b - a )
        d  = dot_product ( ba, ba )
    end function dsqrd
    !same function but for a vector b - a
    real(kind=8) function dsqrd_vect ( ba, R ) result ( d )
        real(kind=8), dimension ( 3 ), intent (in) :: ba !Vector b - a
        real(kind=8), dimension (:,:), intent (in) :: R  !Anisotropic rotation matrix
        real(kind=8), dimension ( 3 )      :: v  !copy of ba
        v = matmul ( R, ba )
        d = dot_product ( v, v )
    end function dsqrd_vect
    !same function for isotropic distance
    real(kind=8) function dsqrd_iso ( a, b ) result ( d )
        real(kind=8), dimension ( 3 ), intent (in) :: a, b
        real(kind=8), dimension ( 3 ) :: ba
        ba = b - a
        d = dot_product( ba, ba )
    end function dsqrd_iso
    !same function for isotropic vector
    real(kind=8) function dsqrd_iso_vect ( ba ) result ( d )
        real(kind=8), dimension ( 3 ), intent (in) :: ba
        d = dot_product( ba, ba )
    end function dsqrd_iso_vect
    !same function for a line L
    real(kind=8) function dsqrt_line ( n, L, outside ) result ( d )
        integer, intent (in) :: n !number of vertices including endpoints
        real(kind=8), dimension (n,3), intent (in) :: L   !line
        logical, intent (out) :: outside
        !INTERNALS
        real(kind=8), dimension (3,3) :: R  !Rotation matrix
        real(kind=8), dimension ( 3 ) :: angles
        real(kind=8), dimension ( 2 ) :: ratios
        integer :: i, j
        
        d  = 0.0
        call split_line ( L, n, outside )
        if(outside)then
            d = 1.0e21 ; Return
        endif
        do i = 1, n - 1 !Loop over the large segments
            do j = nsl(i), nsl(i+1) - 1 !Loop over the small segments
                !set angles, anisotropic ratios and rotation matrix
                !for this small segment.
                
                angles = grid( bsl(j,1), bsl(j,2), bsl(j,3), 1:3 )
                ratios = grid( bsl(j,1), bsl(j,2), bsl(j,3), 4:5 )
                !call fix_angles  ( angles )
                
                
                call get_rmatrix ( angles, ratios, R )
                !Accumulate the distance
                d = d + dsqrt( dsqrd_vect ( dsl(j,:), R ) )
            enddo
        enddo
        Return
    end function dsqrt_line
    
    !Convert three angles of anisotropy to those that are usable.  Note that
    !the vector of angles will be updated and returned
    subroutine fix_angles ( angles )
        real(kind=8), dimension ( 3 ), intent (inout) :: angles
        !1st - angle between the major axis of anisotropy and the
        !      E-W axis. Note: Counter clockwise is positive.
        !2nd - angle between major axis and the horizontal plane.
        !      (The dip of the ellipsoid measured positive down)
        !3rd - Angle of rotation of minor axis about the major axis
        !      of the ellipsoid.
        if(angles(1) >= 0.0 .and. angles(1) < 270.0) then
            angles(1) = ( 90.0 - angles(1)) * DEG2RAD ; else
            angles(1) = (450.0 - angles(1)) * DEG2RAD ; endif
        angles(2) = -1.0 * angles(2) * DEG2RAD
        angles(3) =  1.0 * angles(3) * DEG2RAD
        Return
    end subroutine fix_angles
    
    !Setup an anisotropic rotation matrix from angles and anisotropy ratios.
    !NOTE: the angles should have been adjusted (fix_angles) prior to calling
    !      this routine.
    subroutine get_rmatrix (angles, ratios, R)
        real(kind=8), dimension ( 3 ), intent (in)    :: angles
        real(kind=8), dimension ( 2 ), intent (in)    :: ratios
        real(kind=8), dimension (3,3), intent (inout) :: R
        real(kind=8) :: sina, sinb, sint, cosa, cosb, cost
        real(kind=8), dimension ( 3 ) :: a
        real(kind=8) :: r1, r2
        
        ! Get the required sines and cosines:
        a = real(angles, 8)
        sina = dsin(a(1)) ; sinb = dsin(a(2)) ; sint = dsin(a(3))
        cosa = dcos(a(1)) ; cosb = dcos(a(2)) ; cost = dcos(a(3))

        ! Construct the rotation matrix:
        r1 = 1.0 / ratios(1)
        r2 = 1.0 / ratios(2)
        R(1,:) = (/ cosb * cosa , cosb * sina , -sinb /)
        R(2,:) = (/ r1 * (-cost * sina + sint * sinb * cosa) , &
                    r1 * ( cost * cosa + sint * sinb * sina) , &
                    r1 * ( sint * cosb) /)
        R(3,:) = (/ r2 * ( sint * sina + cost * sinb * cosa) , &
                    r2 * (-sint * cosa + cost * sinb * sina) , &
                    r2 * ( cost * cosb) /)
        Return
    end subroutine get_rmatrix
    
!    !Determine a step parameter t using backtracking line search
!    !given descent direction D and gradient G.
!    real(kind=8) function backtrack ( D, G, L0, n, verts, Fn ) result ( t )
!        real(kind=8), parameter :: phi = 1.618033988749895D+00
!
!        real, parameter :: alpha = 0.01
!        real, parameter :: beta  = 0.9
!        integer, intent (in) :: n  !number of vertices
!        real(kind=8), intent (in) :: L0 !initial length of line
!        real(kind=8), dimension ( : ),   intent (in) :: D, G  !descent direction and gradient
!        real,         dimension ( n, 3), intent (in) :: verts !line vertices
!        real(kind=8), intent (inout) :: Fn
!        !INTERNALS
!        integer :: i
!        logical :: outside
!        real(kind=8) :: C0
!        real, dimension (n,3) :: Vn
!
!        t = 1.0
!        Vn = verts
!        C0 = dot_product(G, D)
!        do i = 1, n-2
!            Vn(i+1,:) = verts(i+1,:) - t * D((i-1)*3+1:i*3)
!        enddo
!        Fn = dsqrt_line ( n, Vn, outside )
!        do while ( Fn > L0 + alpha * t * C0 )
!            t = beta * t
!            do i = 1, n-2
!                Vn(i+1,:) = verts(i+1,:) - t * D((i-1)*3+1:i*3)
!            enddo
!            Fn = dsqrt_line ( n, Vn, outside )
!        enddo
!    end function backtrack
    
    !Perform a golden section line search
    real(kind=8) function goldsection ( D, L0, n, Line, f ) result ( u )
        real(kind=8), parameter :: phi  = 1.618033988749895D+00
        real(kind=8), parameter :: beta = 0.9
        integer, intent (in) :: n  !number of vertices
        real(kind=8), intent (in) :: L0 !initial length of line
        real(kind=8), dimension ( : ),   intent (in) :: D  !descent direction
        real(kind=8), dimension ( n, 3), intent (in) :: Line !line vertices
        real(kind=8), intent (inout) :: f
        !INTERNALS
        real(kind=8), dimension (n,3) :: tline
        real(kind=8) :: fa, fb, fc, df
        real(kind=8) :: ua, ub, uc
        integer :: i
        logical :: outside, left
        
        fc = L0 ; uc = 0.0
        tline = Line
        !Find two points in +/- direction from Line along step D
        ua = max(1.0, 1.0 / maxval(abs(D))) ; outside = .True.
        do while (outside)
            ua = ua * beta
            do i = 1, n-2
                tline(i+1,:) = Line(i+1,:) + ua * D((i-1)*3+1:i*3)
            enddo
            fa = dsqrt_line ( n, tline, outside )
        enddo
        if(fa < L0)then
            f = fa ; u = ua ; return
        endif
        ub = -max(1.0, 1.0 / maxval(abs(D))) ; outside = .True.
        do while (outside)
            ub = ub * beta
            do i = 1, n-2
                tline(i+1,:) = Line(i+1,:) + ub * D((i-1)*3+1:i*3)
            enddo
            fb = dsqrt_line ( n, tline, outside )
        enddo
        if(fb < L0)then
            f = fb ; u = ub ; return
        endif
        
        !Use golden section search to find u
        df = 1.0
        do while ( df > sqrt(epsilon(1.0)) )
            !Propose a new location
            if(ua - uc > uc - ub)then
                u = uc + (ua - uc) / phi
                left = .False.
            else
                u = uc + (ub - uc) / phi
                left = .True.
            endif
            !Evaluate the function at u
            do i = 1, n-2
                tline(i+1,:) = Line(i+1,:) + u * D((i-1)*3+1:i*3)
            enddo
            f = dsqrt_line ( n, tline, outside )
            !Choose new positions
            if(f <= fc)then
                if(left)then
                    fa = fc ; ua = uc ; else
                    fb = fc ; ub = uc ; endif
                df = abs(f - fc)
                fc = f  ; uc = u
            else
                if(left)then
                    fb = f ; ub = u ; else
                    fa = f ; ua = u ; endif
                df = abs(f - fc)
            endif
        enddo
        u = uc ; f = fc
    end function goldsection
    
!Debug support function

    !Write a path defined by L out to the debug file with initial
    !distance D0 and final Dk
    subroutine dbg_print_path ( L, D0, Dk )
        real(kind=8), dimension (:,:), intent (in) :: L
        real(kind=8), intent (in) :: D0, Dk
        integer :: i
        write(dbg,'(a)') 'Start new path'
        write(dbg,'(a,2f10.3)') 'Initial and final distance: ',D0,Dk
        write(dbg,'(a)') 'Vertices: x, y, z'
        do i = 1, ubound(L,1)
            write(dbg,'(3f10.3)') L(i,:)
        enddo
        write(dbg,'(a)') 'End'
        write(dbg,'(a)')
        Return
    end subroutine dbg_print_path

end module class_lvaf2
