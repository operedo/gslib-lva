subroutine get_landmark_pts(ndmax,nd,nx,ny,nz)
      use        geostat
      use        global !to get some global variables
      use        graph     !for using graph theory representation
      use        graph_vars
      !use        dfport
      use        ifport

    integer test,cnt,nodesfl,sys_call_out,j,ndmax,nx,ny,nz
    real*4 LWORK(2)
    nodesfl=9698
    


    write(*,*) 
    write(*,*) 
    write(*,*) '**********************************************************'
    write(*,*) '*Cal. grid to landmark pnts for ISOMAP multidimen scaling*'
    write(*,*) '**********************************************************'
    write(*,*) 
    


    !write out the nodes to calculate distances for
    open(nodesfl,file='nodes2cal.out',status='unknown')
    write(nodesfl,*) xyzland
    write(nodesfl,'(I10)') (landpts)
    close(nodesfl)
    
    start_time=secnds(0.0)
    if(MDS_opt==2) then
       if(idijkstra==1)then
          write(*,*)'Starting multi-thread Boost Dijkstra...'
          sys_call_out = systemqq('rm dist_cpp.out* && ./Boost_dijkstra_openmp && cat dist_cpp.out_* > dist_cpp.out') !calculate the distances to the landmark points
       else
          write(*,*)'Starting single-thread Boost Dijkstra...'
          sys_call_out = systemqq('./Boost_dijkstra') !calculate the distances to the landmark points
       endif
    end if
    sum_time = secnds(start_time)
Write(*,*) 'calculating Dijkstra with Boost library',sum_time
 

    !read in the distances from the file:

    start_time=secnds(0.0)
    open(nodesfl,file='dist_cpp.out',status='unknown')
    
    allocate(coord_ISOMAP(NODES,xyzland),stat = test)
    if(test.ne.0)then
          write(*,*)'ERROR: Allocation of coord_ISOMAP failed due to', &
                ' insufficient memory.'
          stop
    end if

    do j=1,xyzland
        read(nodesfl,*) coord_ISOMAP(1:NODES,j)
    end do
    
    close(nodesfl)
    sum_time = secnds(start_time)
Write(*,*) 'reading dist_cpp.out into coord_ISOMAP',sum_time
    
    
end subroutine get_landmark_pts

subroutine get_landmark_pts_onlymem(ndmax,nd,nx,ny,nz,NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array)
      use        geostat
      use        global !to get some global variables
      use        graph     !for using graph theory representation
      use        graph_vars
      !use        dfport
      use        ifport
      use        boostdijkstra

    integer test,cnt,nodesfl,sys_call_out,j,ndmax,nx,ny,nz
    real*4 LWORK(2)

    integer, intent(in) :: NODES_LENGTH
    integer, intent(in) :: GRID_OUT_LENGTH
    !integer, intent(in) :: NODES2CAL_LENGTH
    integer, intent(in) :: cur_edge_node_array1(GRID_OUT_LENGTH)
    integer, intent(in) :: cur_edge_node_array2(GRID_OUT_LENGTH)
    real*8, intent(in) :: edge_dist_array(GRID_OUT_LENGTH)
    !integer, intent(in) :: nodes2cal_array(GRID_OUT_LENGTH)
     

    nodesfl=9698
    


    write(*,*) 
    write(*,*) 
    write(*,*) '**********************************************************'
    write(*,*) '*Cal. grid to landmark pnts for ISOMAP multidimen scaling*'
    write(*,*) '**********************************************************'
    write(*,*) 
   
    sum_time = secnds(start_time)
    allocate(coord_ISOMAP(NODES,xyzland),stat = test)
    if(test.ne.0)then
          write(*,*)'ERROR: Allocation of coord_ISOMAP failed due to', &
                ' insufficient memory.'
          stop
    end if

    !coord_ISOMAP(1,1)=1.0
    !coord_ISOMAP(2,1)=2.0
    !coord_ISOMAP(1,19)=3.0
    !coord_ISOMAP(2,19)=4.0
    !coord_ISOMAP(1,23)=5.0
    !coord_ISOMAP(2,23)=6.0
    !write(*,*) size(cur_edge_node_array1(:)),size(cur_edge_node_array2(:)),size(edge_dist_array(:)) 
    !write(*,*) coord_ISOMAP(1,1),coord_ISOMAP(2,1),coord_ISOMAP(1,19),coord_ISOMAP(2,19),coord_ISOMAP(1,23),coord_ISOMAP(2,23)  
    call dijkstra(NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array,xyzland,landpts,coord_ISOMAP)
    !write(*,*) coord_ISOMAP(1,1),coord_ISOMAP(2,1),coord_ISOMAP(1,19),coord_ISOMAP(2,19),coord_ISOMAP(1,23),coord_ISOMAP(2,23)  
!stop

!    !write out the nodes to calculate distances for
!    open(nodesfl,file='nodes2cal.out',status='unknown')
!    write(nodesfl,*) xyzland
!    write(nodesfl,'(I10)') (landpts)
!    close(nodesfl)
!    
!    start_time=secnds(0.0)
!    if(MDS_opt==2) then
!       if(idijkstra==1)then
!          write(*,*)'Starting multi-thread Boost Dijkstra...'
!          sys_call_out = systemqq('rm dist_cpp.out* && ./Boost_dijkstra_openmp && cat dist_cpp.out_* > dist_cpp.out') !calculate the distances to the landmark points
!       else
!          write(*,*)'Starting single-thread Boost Dijkstra...'
!          sys_call_out = systemqq('./Boost_dijkstra') !calculate the distances to the landmark points
!       endif
!    end if
!    sum_time = secnds(start_time)
!Write(*,*) 'calculating Dijkstra with Boost library',sum_time
! 
!
!    !read in the distances from the file:
!
!    start_time=secnds(0.0)
!    open(nodesfl,file='dist_cpp.out',status='unknown')
!    
!    allocate(coord_ISOMAP(NODES,xyzland),stat = test)
!    if(test.ne.0)then
!          write(*,*)'ERROR: Allocation of coord_ISOMAP failed due to', &
!                ' insufficient memory.'
!          stop
!    end if
!
!    do j=1,xyzland
!        read(nodesfl,*) coord_ISOMAP(1:NODES,j)
!    end do
!    
!    close(nodesfl)
    sum_time = secnds(start_time)
Write(*,*) 'reading dist_cpp.out into coord_ISOMAP',sum_time
    
    
end subroutine get_landmark_pts_onlymem

subroutine get_landmark_pts_onlymem_xyzland(ndmax,nd,nx,ny,nz,NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array,xyzland2,landpts2)
      use        geostat
      use        global !to get some global variables
      use        graph     !for using graph theory representation
      use        graph_vars
      !use        dfport
      use        ifport
      use        boostdijkstra

    integer test,cnt,nodesfl,sys_call_out,j,ndmax,nx,ny,nz
    real*4 LWORK(2)

    integer, intent(in) :: NODES_LENGTH
    integer, intent(in) :: GRID_OUT_LENGTH
    !integer, intent(in) :: NODES2CAL_LENGTH
    integer, intent(in) :: cur_edge_node_array1(GRID_OUT_LENGTH)
    integer, intent(in) :: cur_edge_node_array2(GRID_OUT_LENGTH)
    real*8, intent(in) :: edge_dist_array(GRID_OUT_LENGTH)
    !integer, intent(in) :: nodes2cal_array(GRID_OUT_LENGTH)
    integer, intent(in) :: xyzland2
    integer, intent(in) :: landpts2(xyzland2)

    nodesfl=9698
    


    write(*,*) 
    write(*,*) 
    write(*,*) '**********************************************************'
    write(*,*) '*Cal. grid to landmark pnts for ISOMAP multidimen scaling*'
    write(*,*) '**********************************************************'
    write(*,*) NODES,xyzland2 
   
    sum_time = secnds(start_time)
    allocate(coord_ISOMAP(NODES,xyzland2),stat = test)
    if(test.ne.0)then
          write(*,*)'ERROR: Allocation of coord_ISOMAP failed due to', &
                ' insufficient memory.'
          stop
    end if

    call dijkstra(NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array,xyzland2,landpts2,coord_ISOMAP)

    sum_time = secnds(start_time)
Write(*,*) 'reading dist_cpp.out into coord_ISOMAP',sum_time
    
    
end subroutine get_landmark_pts_onlymem_xyzland



subroutine MDS_ISOMAP(ndmax,nd,nx,ny,nz)    
      use        geostat
      use        global    !to get some global variables 
      use        graph     !for using graph theory representation
      use        graph_vars
    integer i,j,stress_nodes,nx,ny,nz,ndmax
    real*8 D1,vec_temp(xyzland),mds_diff
    real*8, allocatable, dimension (:,:) :: A1,B1
    integer test,cnt

! variables for dgemm mkl
    real*8, allocatable, dimension (:,:) :: coord_ISOMAP_out
    integer M,N,K
    real*8 ALPHA,BETA

    open(3676,file='grid_dists.out',status="unknown")
     

    !need to make the subB matrix:
    
    !subB = -.5*(D.^2 - sum(D'.^2)'*ones(1,nl)/nl - ones(N,1)*sum(D.^2)/N+sum(sum(D.^2))/(N*nl));
    !subB = -.5*(D.^2 -            A              -          B           +             D     );
    
    allocate(A1(NODES,1))
    allocate(B1(1,xyzland))
    
    if(cal_stress==1) then
        allocate(temp_coord(NODES,xyzland))
        temp_coord=coord_ISOMAP
    end if    
    
    start_time=secnds(0.0)
    coord_ISOMAP=coord_ISOMAP**2
    sum_time = secnds(start_time)
Write(*,*) 'calculating coord_ISOMAP**2',sum_time
    
    start_time=secnds(0.0)
    do i=1,NODES
      A1(i,1)=sum(coord_ISOMAP(i,:))/xyzland
    end do
    sum_time = secnds(start_time)
Write(*,*) 'calculating A1(:,1)',sum_time

    start_time=secnds(0.0)
    do i=1,xyzland
      B1(1,i)=sum(coord_ISOMAP(:,i))/NODES
    end do
    sum_time = secnds(start_time)
Write(*,*) 'calculating B1(:,1)',sum_time

    start_time=secnds(0.0)
    D1=sum(coord_ISOMAP)/(NODES*xyzland)
    sum_time = secnds(start_time)
Write(*,*) 'calculating D1',sum_time
    
    start_time=secnds(0.0)
    do j=1,NODES
    do i=1,xyzland
        coord_ISOMAP(j,i)= coord_ISOMAP(j,i) - A1(j,1) - B1(1,i) + D1
    end do
    end do
    sum_time = secnds(start_time)
Write(*,*) 'calculating updated coord_ISOMAP',sum_time
    
    coord_ISOMAP=-0.5*coord_ISOMAP

    !subB=-0.5*(transpose(d2lanmark) - A1 - B1 + D1 )

    if(allocated(B1)) deallocate(B1)
    allocate(subB2(xyzland,xyzland))
      
    start_time=secnds(0.0)
    !subB2=matmul(transpose(coord_ISOMAP),coord_ISOMAP)

M=xyzland
K=NODES
N=xyzland
ALPHA=1.0
BETA=0.0
CALL DGEMM('T','N',M,N,K,ALPHA,coord_ISOMAP,K,coord_ISOMAP,K,BETA,subB2,M)




    sum_time = secnds(start_time)
Write(*,*) 'calculating subB2',sum_time

    start_time=secnds(0.0)
    call eig ( subB2, evalues, vectors, ubound(subB2,1) )
    sum_time = secnds(start_time)
Write(*,*) 'calculating the eigenvalues and vectors',sum_time
    write(T_eigs,'(3I12,f19.5)')  NODES, edges, ubound(subB2,1), secnds(start_time)

!Write(*,*) 'calculating the coord of all grid points'
    
    start_time=secnds(0.0)
    
    !need to reverse array order
    !use prexisting arrays as temp arrays
    subB2=-999 !array for vectors
    A1=-999 !array for values
    subB2(1:xyzland,1:xyzland)=vectors
    A1(1:xyzland,1)=evalues
    
    do i=1,xyzland
        evalues(i)=A1(xyzland-i+1,1)
        vectors(:,i) = subB2(1:xyzland,xyzland-i+1)
    end do

    !find the smallest eigenvalue that is pos
    do i=1,xyzland
        if(evalues(i)<=0) then
            if(dim==-1) then !make the max possible
                dim=i-1
            else !set to user defined value if it will be pos def
                if(dim>i) dim=i
            end if
            exit
        end if
    end do
    if(dim==xyzland) dim=xyzland-1 !discard the smallest if there are none neg, as per ISOMAP (see their paper)
    if(dim==0 .or. dim==-1) dim=xyzland-1 !discard the smallest if there are none neg, as per ISOMAP (see their paper)
    write(*,*) 'number of dimensions to acutally use:', int(dim)

    allocate(val(xyzland,xyzland))
     
    !start_time=secnds(0.0)    


!!this replaces the ellements in a1 with a1*a2, essentially the matrix mult "matmul(coord_ISOMAP,vectors(:,1:dim))" but with less memory

!!$omp parallel default(shared) private(i,j,vec_temp)
!!$omp do schedule(static)
!do i=1,NODES
!do j=1,dim  
!    vec_temp(j) = sum(  coord_ISOMAP(i,:)*vectors(:,j)  )
!end do
!!now dont need the row in a1
!coord_ISOMAP(i,1:dim) = vec_temp
!end do
!!$omp end do nowait
!!$omp end parallel
!!dgemm/sgemm -> blas

M=NODES
K=xyzland
N=dim
ALPHA=1.0
BETA=0.0
allocate(coord_ISOMAP_out(M,N))
CALL DGEMM('N','N',M,N,K,ALPHA,coord_ISOMAP,M,vectors,K,BETA,coord_ISOMAP_out,M)
coord_ISOMAP(:,1:dim) = coord_ISOMAP_out
deallocate(coord_ISOMAP_out)

!write(*,*) coord_ISOMAP(:,1)


if(allocated(vectors)) deallocate(vectors)
      
    !need the eigenvalues as an array
    val=0
    do i=1,xyzland
    val(i,i) = 1/sqrt(evalues(i))
    end do
      
      !now do vec1*val**(-1)
    if(allocated(evalues)) deallocate(evalues)
       do i=1,dim
        coord_ISOMAP(:,i)=coord_ISOMAP(:,i)*val(i,i)
       end do
       
    A1=0
    do i=1,dim
      A1(i,1)=1/sqrt(val(i,i))
    end do
    
    if(allocated(val)) deallocate(val)
    
    do i=1,NODES
      coord_ISOMAP(i,1:dim)=coord_ISOMAP(i,1:dim)*A1(1:dim,1)
    end do
    
    !write out the coords?
!    open(47861, file='coords.out',status='unknown')
!    write(47861,*) 'coordinates of the points'
!    write(47861,*) dim
!    do i=1,dim
!    write(47861,'(a,I5)') 'dim ',i
!    end do
!    
!    do i=1,NODES
!     write(47861,'(3000f18.8)') coord_ISOMAP(i,1:dim)
!    end do
    
    if(allocated(subB2)) deallocate(subB2)
    if(allocated(A1)) deallocate(A1)

    sum_time = secnds(start_time)
Write(*,*) 'calculating the coord of all grid points',sum_time


    if(cal_stress==0) return    !return if you do not want to compare original distances with isomap distances


MDS_diff=1

    
!need to calculate the stress


stress = 0
sum_dist = 0

stress_nodes  = xyzland !just landmark

    do j=1,stress_nodes 
    cur_node = landpts(j)
    do i=1,NODES
    !what is the after MDS dist - pre MDS
    !dist between node i and j after MDS
        post_MDS_dist = sqrt(  sum ( (   coord_ISOMAP(cur_node,1:dim) -  coord_ISOMAP(i,1:dim)    )**2   )   )
        stress = stress + (   temp_coord(i,j) - post_MDS_dist   )**2
        sum_dist = sum_dist + temp_coord(i,j)**2
    end do
    end do
    
    stress=stress/sum_dist

    write(*,*) 'the stress is ', stress
    open(9367, file = 'stress.out', status='unknown')
    write(9367,*) xyzland,stress, nx*ny*nz




end subroutine MDS_ISOMAP


subroutine MDS_ISOMAP_xyzland(ndmax,nd,nx,ny,nz,xyzland2)    
      use        geostat
      use        global    !to get some global variables 
      use        graph     !for using graph theory representation
      use        graph_vars
    integer i,j,stress_nodes,nx,ny,nz,ndmax,xyzland2
    real*8 D1,vec_temp(xyzland),mds_diff
    real*8, allocatable, dimension (:,:) :: A1,B1
    integer test,cnt

! variables for dgemm mkl
    real*8, allocatable, dimension (:,:) :: coord_ISOMAP_out
    integer M,N,K
    real*8 ALPHA,BETA

    open(3676,file='grid_dists.out',status="unknown")
     

    !need to make the subB matrix:
    
    !subB = -.5*(D.^2 - sum(D'.^2)'*ones(1,nl)/nl - ones(N,1)*sum(D.^2)/N+sum(sum(D.^2))/(N*nl));
    !subB = -.5*(D.^2 -            A              -          B           +             D     );

    write(*,*) xyzland2

    if(allocated(evalues)) deallocate(evalues)
    if(allocated(vectors)) deallocate(vectors)
    allocate(evalues(xyzland2))
    allocate(vectors(xyzland2,xyzland2))


    
    allocate(A1(NODES,1))
    allocate(B1(1,xyzland2))
    
    if(cal_stress==1) then
        allocate(temp_coord(NODES,xyzland2))
        temp_coord=coord_ISOMAP
    end if    
    
    start_time=secnds(0.0)
    coord_ISOMAP=coord_ISOMAP**2
    sum_time = secnds(start_time)
Write(*,*) 'calculating coord_ISOMAP**2',sum_time
    
    start_time=secnds(0.0)
    do i=1,NODES
      A1(i,1)=sum(coord_ISOMAP(i,:))/xyzland2
    end do
    sum_time = secnds(start_time)
Write(*,*) 'calculating A1(:,1)',sum_time

    start_time=secnds(0.0)
    do i=1,xyzland2
      B1(1,i)=sum(coord_ISOMAP(:,i))/NODES
    end do
    sum_time = secnds(start_time)
Write(*,*) 'calculating B1(:,1)',sum_time

    start_time=secnds(0.0)
    D1=sum(coord_ISOMAP)/(NODES*xyzland2)
    sum_time = secnds(start_time)
Write(*,*) 'calculating D1',sum_time
    
    start_time=secnds(0.0)
    do j=1,NODES
    do i=1,xyzland2
        coord_ISOMAP(j,i)= coord_ISOMAP(j,i) - A1(j,1) - B1(1,i) + D1
    end do
    end do
    sum_time = secnds(start_time)
Write(*,*) 'calculating updated coord_ISOMAP',sum_time
    
    coord_ISOMAP=-0.5*coord_ISOMAP

    !subB=-0.5*(transpose(d2lanmark) - A1 - B1 + D1 )

    if(allocated(B1)) deallocate(B1)
    allocate(subB2(xyzland2,xyzland2))
      
    start_time=secnds(0.0)
    !subB2=matmul(transpose(coord_ISOMAP),coord_ISOMAP)

    M=xyzland2
    K=NODES
    N=xyzland2
    ALPHA=1.0
    BETA=0.0
    CALL DGEMM('T','N',M,N,K,ALPHA,coord_ISOMAP,K,coord_ISOMAP,K,BETA,subB2,M)




    sum_time = secnds(start_time)
Write(*,*) 'calculating subB2',sum_time

    start_time=secnds(0.0)
    call eig ( subB2, evalues, vectors, ubound(subB2,1) )
    sum_time = secnds(start_time)
Write(*,*) 'calculating the eigenvalues and vectors',sum_time
    write(T_eigs,'(3I12,f19.5)')  NODES, edges, ubound(subB2,1), secnds(start_time)

!Write(*,*) 'calculating the coord of all grid points'
    
    start_time=secnds(0.0)
    
    !need to reverse array order
    !use prexisting arrays as temp arrays
    subB2=-999 !array for vectors
    A1=-999 !array for values
    subB2(1:xyzland2,1:xyzland2)=vectors
    A1(1:xyzland2,1)=evalues
    
    do i=1,xyzland2
        evalues(i)=A1(xyzland2-i+1,1)
        vectors(:,i) = subB2(1:xyzland2,xyzland2-i+1)
    end do

    !find the smallest eigenvalue that is pos
    do i=1,xyzland2
        if(evalues(i)<=0) then
            if(dim==-1) then !make the max possible
                dim=i-1
            else !set to user defined value if it will be pos def
                if(dim>i) dim=i
            end if
            exit
        end if
    end do
    !if(dim==xyzland2) dim=xyzland2-1 !discard the smallest if there are none neg, as per ISOMAP (see their paper)
    dim=xyzland2-1 !discard the smallest if there are none neg, as per ISOMAP (see their paper)
    if(dim==0 .or. dim==-1) dim=xyzland2-1 !discard the smallest if there are none neg, as per ISOMAP (see their paper)
    write(*,*) 'number of dimensions to acutally use:', int(dim)

    allocate(val(xyzland2,xyzland2))
     
    !start_time=secnds(0.0)    


!!this replaces the ellements in a1 with a1*a2, essentially the matrix mult "matmul(coord_ISOMAP,vectors(:,1:dim))" but with less memory

!!$omp parallel default(shared) private(i,j,vec_temp)
!!$omp do schedule(static)
!do i=1,NODES
!do j=1,dim  
!    vec_temp(j) = sum(  coord_ISOMAP(i,:)*vectors(:,j)  )
!end do
!!now dont need the row in a1
!coord_ISOMAP(i,1:dim) = vec_temp
!end do
!!$omp end do nowait
!!$omp end parallel
!!dgemm/sgemm -> blas

M=NODES
K=xyzland2
N=dim
ALPHA=1.0
BETA=0.0
allocate(coord_ISOMAP_out(M,N))
CALL DGEMM('N','N',M,N,K,ALPHA,coord_ISOMAP,M,vectors,K,BETA,coord_ISOMAP_out,M)
coord_ISOMAP(:,1:dim) = coord_ISOMAP_out
deallocate(coord_ISOMAP_out)

!write(*,*) coord_ISOMAP(:,1)


if(allocated(vectors)) deallocate(vectors)
      
    !need the eigenvalues as an array
    val=0
    do i=1,xyzland2
    val(i,i) = 1/sqrt(evalues(i))
    end do
      
      !now do vec1*val**(-1)
    if(allocated(evalues)) deallocate(evalues)
       do i=1,dim
        coord_ISOMAP(:,i)=coord_ISOMAP(:,i)*val(i,i)
       end do
       
    A1=0
    do i=1,dim
      A1(i,1)=1/sqrt(val(i,i))
    end do
    
    if(allocated(val)) deallocate(val)
    
    do i=1,NODES
      coord_ISOMAP(i,1:dim)=coord_ISOMAP(i,1:dim)*A1(1:dim,1)
    end do
    
    !write out the coords?
!    open(47861, file='coords.out',status='unknown')
!    write(47861,*) 'coordinates of the points'
!    write(47861,*) dim
!    do i=1,dim
!    write(47861,'(a,I5)') 'dim ',i
!    end do
!    
!    do i=1,NODES
!     write(47861,'(3000f18.8)') coord_ISOMAP(i,1:dim)
!    end do
    
    if(allocated(subB2)) deallocate(subB2)
    if(allocated(A1)) deallocate(A1)

    sum_time = secnds(start_time)
Write(*,*) 'calculating the coord of all grid points',sum_time


    if(cal_stress==0) return    !return if you do not want to compare original distances with isomap distances


MDS_diff=1

    
!need to calculate the stress


stress = 0
sum_dist = 0

stress_nodes  = xyzland2 !just landmark

    do j=1,stress_nodes 
    cur_node = landpts(j)
    do i=1,NODES
    !what is the after MDS dist - pre MDS
    !dist between node i and j after MDS
        post_MDS_dist = sqrt(  sum ( (   coord_ISOMAP(cur_node,1:dim) -  coord_ISOMAP(i,1:dim)    )**2   )   )
        stress = stress + (   temp_coord(i,j) - post_MDS_dist   )**2
        sum_dist = sum_dist + temp_coord(i,j)**2
    end do
    end do
    
    stress=stress/sum_dist

    write(*,*) 'the stress is ', stress
    open(9367, file = 'stress.out', status='unknown')
    write(9367,*) xyzland2,stress, nx*ny*nz




end subroutine MDS_ISOMAP_xyzland


