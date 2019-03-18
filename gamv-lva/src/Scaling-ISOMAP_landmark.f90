subroutine get_landmark_pts(ndmax,nd)
      use geostat
      use global !to get some global variables for the robust solver
      use graph     !for using graph theory representation
      use graph_vars
      use ifport
      
    integer sys_call_out,j,i,k,ii,jj,kk,ndmax,no_cpp,indk,indj,indi,ind

    integer test111,xspace,yspace,zspace,cnt,nodesfl
    !real*4 LWORK(2)
    real(kind=sp) LWORK(2)
    nodesfl=9698
    
    xyzland=xland*yland*zland           
     
     
      allocate(landpts(xyzland),stat = test111)
            if(test111.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
      allocate(evalues(xyzland),stat = test111)
            if(test111.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
      allocate(vectors(xyzland,xyzland),stat = test111)
            if(test111.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
           
       !   !xyzland=int(real(xyzland)*(6.0/12.0))
       !   print *,'NODES*xyzland=',(NODES*xyzland)

       !   !if(NODE*xyzland.lt.2**26)then 
       !   allocate(coord_ISOMAP(NODES,xyzland),stat = test111) 

       ! print *,'size(coord_ISOMAP)=',size(coord_ISOMAP)
       ! !do i=1,NODES
       ! !do j=1,xyzland
       ! !   coord_ISOMAP(i,j)=i+j;
       ! !end do
       ! !end do
       ! !call sleep(10)


       ! if(test111.ne.0)then
       !       write(*,*)'ERROR: Allocation of coord_ISOMAP failed due to', &
       !             ' insufficient memory.'
       !       stop
       ! end if
       ! deallocate(coord_ISOMAP)

        !else
        !   write(*,*)'Using out-of-core computations (overflow in memory usage)'
        !end if


    !stop

    write(*,*) 
    write(*,*) 
    write(*,*) '**********************************************************'
    write(*,*) '*Cal. grid to landmark pnts for ISOMAP multidimen scaling*'
    write(*,*) '**********************************************************'
    write(*,*) 
    



    !where are our landmark points?  xyzland
    !evenly space them out

    print *, nx,ny,nz

    xspace=int(real(nx/xland))
    yspace=int(real(ny/yland))
    zspace=int(real(nz/zland))
    no_cpp=0
    cnt=0

    do k=1,zland
    indk=zspace*k    - 0.5*zspace + 1
    do j=1,yland
    indj=yspace*j    - 0.5*yspace + 1
    do i=1,xland
      cnt=cnt+1    
      indi=xspace*i - 0.5*xspace + 1
      ind=get_node(indi,indj,indk,nx,ny,nz)
      landpts(cnt)=ind
      
      
      
      
    end do
    end do
    end do



    !write out the nodes to calculate distances for
    open(nodesfl,file='nodes2cal.out',status='unknown')
    write(nodesfl,*) xyzland
    write(nodesfl,'(I10)') (landpts)
    close(nodesfl)
    
    if(MDS_opt==2) then
       if(idijkstra==1)then
          write(*,*)'Starting multi-thread Boost Dijkstra...'
          sys_call_out = systemqq('rm dist_cpp.out* && ./Boost_dijkstra_openmp && cat dist_cpp.out_* > dist_cpp.out') !calculate the distances to the landmark points
       else
          write(*,*)'Starting single-thread Boost Dijkstra...'
          sys_call_out = systemqq('./Boost_dijkstra') !calculate the distances to the landmark points
       endif
    end if
    
    !read in the distances from the file:
    

          !if(NODE*xyzland.lt.2**26)then 
    
    open(nodesfl,file='dist_cpp.out',status='unknown')
    
        allocate(coord_ISOMAP(NODES,xyzland),stat = test111) 
        if(test111.ne.0)then
              write(*,*)'ERROR: Allocation of coord_ISOMAP failed due to', &
                    ' insufficient memory.'
              stop
        end if



    do j=1,xyzland
    do i=1,NODES
        read(nodesfl,*) coord_ISOMAP(i,j)
    end do
    end do
    close(nodesfl)
   
    !end if
 
end subroutine get_landmark_pts








subroutine MDS_ISOMAP(ndmax,nd,nx,ny,nz)    
      use        geostat
      use        global !to get some global variables for the robust solver
      use        graph     !for using graph theory representation
      use        graph_vars

    integer i,j,k,ii,jj,kk,stress_nodes,ndmax

    !real*8 D1,vec_temp(xyzland),start_time,sum_time,t_eigs
    real(kind=dp) D1,vec_temp(xyzland),start_time,sum_time,t_eigs
    !real*8, allocatable, dimension (:,:) :: A1,B1
    real(kind=dp), allocatable, dimension (:,:) :: A1,B1

    integer test111,cnt
    
    !to do: include the matrix_mult program


    !print *,'IN'
    
    !open(3676,file='grid_dists.out',status="unknown")
    !open(3676,file='dist_cpp.out',status="unknown")

    !print *,'OUT'
     

    !need to make the subB matrix:
    !subB = -.5*(D.^2 - sum(D'.^2)'*ones(1,nl)/nl - ones(N,1)*sum(D.^2)/N+sum(sum(D.^2))/(N*nl));
    !subB = -.5*(D.^2 -            A              -          B           +             D     );
    

    allocate(A1(NODES,1))
    allocate(B1(1,xyzland))
    coord_ISOMAP=coord_ISOMAP**2
    do i=1,NODES
      A1(i,1)=sum(coord_ISOMAP(i,:))/xyzland
    end do
    !do i=2,xyzland
    !  A1(:,i)=A1(:,1)
    !end do
    do i=1,xyzland
      B1(1,i)=sum(coord_ISOMAP(:,i))/NODES
    end do
    !do i=2,NODES
    !  B1(i,:)=B1(1,:)
    !end do
    
    D1=sum(coord_ISOMAP)/(NODES*xyzland)
!!$omp parallel default(shared) private(i,j)
!!$omp do schedule(static) 
    do j=1,NODES
       do i=1,xyzland
          coord_ISOMAP(j,i)= coord_ISOMAP(j,i) - A1(j,1) - B1(1,i) + D1
       end do
    end do
!!$omp end do 
!!$omp end parallel

    coord_ISOMAP=-0.5*coord_ISOMAP

    !subB=-0.5*(transpose(d2lanmark) - A1 - B1 + D1 )
    !d2lanmark=d2lanmark**(0.5)
    
    if(allocated(B1)) deallocate(B1)
    allocate(subB2(xyzland,xyzland))
!     
      
    start_time=secnds(0.0)
    subB2=matmul(transpose(coord_ISOMAP),coord_ISOMAP)
    sum_time = secnds(start_time)

    Write(*,*) 'calculating the eigenvalues and vectors'

    start_time=secnds(0.0)
    call eig ( subB2, evalues, vectors, ubound(subB2,1) )
    write(T_eigs,'(3I12,f19.5)')  NODES, edges, ubound(subB2,1), secnds(start_time)

    Write(*,*) 'calculating the coord of all grid points'
    
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
                dim=i
            else !set to user defined value if it will be pos def
                if(dim>i) dim=i
            end if
            exit
        end if
    end do

    if(dim==xyzland) dim=xyzland-1 !discard the smallest if there are none neg, as per ISOMAP (see their paper)
    if(dim==0 .or. dim==-1) dim=xyzland-1 !discard the smallest if there are none neg, as per ISOMAP (see their paper)
    write(*,*) 'number of dimensions to acutally use:', int(dim)
    
    !make the vec array

    allocate(val(xyzland,xyzland))
    !      
    start_time=secnds(0.0)    
    !need a temp array for matmul or else get a stack overflow, strange b/c this way should take more memory

    allocate(subB(NODES,dim))
    !write(*,*)'matmul(coord_ISOMAP,vectors(:,1:dim))'
    subB=matmul(coord_ISOMAP,vectors(:,1:dim))


    !write(T_matrix,'(3I12,f19.5)')  NODES, edges, dim, sum_time + secnds(start_time)
    coord_ISOMAP(:,1:dim)=subB

    !!this replaces the ellements in a1 with a1*a2, essentially the matrix mult "matmul(coord_ISOMAP,vectors(:,1:dim))" but with less memory
!$omp parallel default(shared) private(i,j,vec_temp)
!$omp do schedule(static)
    do i=1,NODES
      !write(*,*) i,'/',NODES,': sum(coord_ISOMAP(',i,')*vectors(:,...))'
      do j=1,dim   
        !write(*,*) i,'/',NODES,', ',j,'/',dim,': sum(coord_ISOMAP(',i,')*vectors(:,',j,'))'
        vec_temp(j) = sum(  coord_ISOMAP(i,:)*vectors(:,j)  )
      end do
      !now dont need the row in a1
      coord_ISOMAP(i,1:dim) = vec_temp
    end do
!$omp end do nowait
!$omp end parallel

    deallocate(subB)

    !if(allocated(subB)) deallocate(subB)
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
     !vec now has the coords
    if(allocated(A1)) deallocate(A1)

    write(*,*) 'exiting MDS_ISOMAP'

    !write(*,*) 'writing modified coord_ISOMAP'
   
    !do j=1,xyzland
    !do i=1,NODES
    !    write(*,*) coord_ISOMAP(i,j)
    !end do
    !end do
 
    return    

end subroutine MDS_ISOMAP


