module sort
    implicit none
    public :: QUICK_SORT, SWAP

  contains

    recursive  subroutine  QUICK_SORT (A, Left, Right)
    !! This subroutine implements Quicksort by C. A. R. Hoare.
    !! INCOMING: A     = an array whose elements are to be sorted;
    !!           Left  = a pointer to the left-hand end of the partition to be sorted;
    !!           Right = a pointer to the right-hand end of the partition to be sorted;
    !! OUTGOING: A     = an array whose elements A(Left) through A(Right) are sorted.

    real(kind=8), dimension ( : ), intent (in out) :: A
    integer, intent (in) :: Left, Right
    real(kind=8)    :: Ref         !! To hold the Reference Element.
    integer :: L, R        !! L = a pointer used when searching
                            !! from the left.
                            !! R = a pointer, used when searching
                            !! from the right.

    if (Right <= Left) then                        !! The partition is empty or contains one element.
        return                                      !! There's no work to do.
    else if (Right - Left == 1) then               !! There are 2 elements in the partition.
        if (A(Left) > A(Right)) then
            call SWAP (A(Left), A(Right))
        end  if
        return
    end  if

    L = Left
    R = Right + 1

    call SWAP (A(Left), A((Left + Right)/2) )      !! Bring the Reference Element to the
                                                                                    !! left end of the partition.
    Ref = A(Left)                                  !! Select the Reference Element.

    !! Partition the elements into three groups.
    do
        if (L >= R) then
            exit
        end  if
        do L = L + 1, R - 1                         !! Scan from the left for an element
            if (A(L) > Ref) then                     !! that is larger than the Reference.
                exit
            end  if
        end  do

        do R = R - 1, L, -1                         !! Scan from the right for an element that
            if (A(R) <= Ref) then                    !! is less than or equal to the Reference.
                exit
            end  if
        end  do

        if (L < R) then                             !! Swap two elements that are in the
                                                                                    !! wrong partitions.
            call SWAP (A(L), A(R))
        end  if
    end  do
                                                                                    !! Partitioning is complete.
    if (Left < R) then                             !! Swap the Reference Element into its
                                                                                    !! final position R in the array.
        A(Left) = A(R)
        A(R)    = Ref
    end  if
    !! At this point, A(R) is in its correct position in the list.  Elements A(Left) to A(R-1) are
    !! less than or equal to A(R), and elements A(R+1) to A(Right) are greater then A(R).

    call QUICK_SORT (A, Left, R-1)                 !! Partition the left segment.
                                                                                    !! The element at R is in its correct position.
    call QUICK_SORT (A, R+1, Right)                !! Partition the right segment.

    end  subroutine  QUICK_SORT

    !! This subroutine swaps the element Left_Element with Right_Element.
    subroutine SWAP (Left_Element, Right_Element)
    real(kind=8), intent (in out) :: Left_Element, Right_Element
    real(kind=8) :: Temp
    Temp     = Left_Element
    Left_Element   = Right_Element
    Right_Element  = Temp
    end subroutine SWAP

end module sort