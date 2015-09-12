module grid_info
    implicit none
    
    integer, public, dimension ( 3 ) :: n,temp_n    !grid size (x,y,z)
    real*8,    public, dimension ( 3 ) :: o, sss,temp_o,temp_sss !origine and block size (x,y,z)
    real*8,    public, allocatable, dimension (:,:,:,:) :: grid,temp_grid !grid information
    
  contains
      integer function block_count (a, b)
        real(kind=8), dimension ( 3 ), intent(in) :: a, b
        integer, dimension ( 3 ) :: ida, idb
        logical :: outside
        call get_index ( a, ida, outside )
        if(.not.outside)then
            call get_index ( b, idb, outside )
            if(.not.outside)then
                block_count = sum(abs(ida - idb))
            else
                block_count = -1
            endif
        else
            block_count = -1
        endif
    end function block_count
    
    !Return the index of a gridblock for point p
    subroutine get_index (p, ind, outside)
        real(kind=8), dimension ( 3 ), intent (in)     :: p !point of interest
        integer,      dimension ( 3 ), intent (in out) :: ind  !resulting index
        logical, intent(out) :: outside
        integer :: i
        outside = .False.
        do i = 1, 3
            if(n(i) > 1)then
                ind(i) = int((p(i) - o(i)) / max(1.0,sss(i)) + 1.5)
                if(ind(i) < 1 .or. ind(i) > n(i))then
                    outside = .True.
                endif
            else
                ind(i) = 1
            endif
        enddo
        Return
    end subroutine get_index
    
    
end module grid_info
