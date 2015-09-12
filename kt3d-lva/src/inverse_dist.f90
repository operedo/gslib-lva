subroutine inv_d(outfl,ndmax,ndmin)
      use        geostat
      use        global   
      use        graph     !for using graph theory representation
      use        graph_vars
      use        ifport
      implicit none

    character outfl*512
    real*8, allocatable, dimension (:) :: in_est
    integer ind55,ndmax,ndmin,i,j,k

    open(987,file=outfl,status='unknown')
    write(987,*) 'inverse dist estimate with power', power
    write(987,*) '1'
    write(987,*) 'est'

    call get_inv_dists(ndmax)  !to get the distances
    
    coord_ISOMAP=coord_ISOMAP**power !get the weights from the distances
    rad_inv_dis=rad_inv_dis**power !make the trimming rad a power, then trim based on the weigth
    do j=1,ndmax
    do i=1,NODES
        if(  coord_ISOMAP(i,j) <rad_inv_dis) coord_ISOMAP(i,j) = 0  !will give a weight = 0
    end do
    end do    
    
    !loop over all the grid nodes and get the estiamtes
    allocate(in_est(NODES))  
    in_est=0
    do i=1,NODES
        if(sum(coord_ISOMAP(i,:)) ==0) then
            in_est(i)=-999 !check if no close data
        end if
        coord_ISOMAP(i,:) = coord_ISOMAP(i,:) / sum(coord_ISOMAP(i,:)) !standardize the weights
        coord_ISOMAP(i,:) = coord_ISOMAP(i,:)*vr(coord_ISOMAP_ind(i,:))      !wt*data
        if(in_est(i) /=-999) in_est(i)=sum(coord_ISOMAP(i,:))   !est
        !did we find enough data:
        ind55=ndmin
        do j=1,ndmax
            if(coord_ISOMAP_ind(i,j) < 0 ) then
                ind55=j-1
                exit
            end if
        end do
        if(ind55 < ndmin ) in_est(i) = -999 
    end do
    
    in_est(data_ind(1:nd)) = vr(1:nd)!reset the data values to the data values
    write(987,'(f25.8)')   in_est   !need to write out the inv dis est

stop   !done
end subroutine inv_d


subroutine get_inv_dists(ndmax)
      use        geostat
      use        global 
      use        graph     !for using graph theory representation
      use        graph_vars
      use        ifport
    implicit none
    
    real*8 max_dist
    integer test,cnt,nodesfl,ind55,ndmax,sys_call_out,i,j,k,ii,jj,kk
    real*4 LWORK(2)
    nodesfl=9698
  
    !write out the nodes to calculate distances for, for inv_dist this is all the data locations, get dist from data to all nodes
    open(nodesfl,file='nodes2cal.out',status='unknown')
    write(nodesfl,*) nd
    write(nodesfl,'(I10)') (data_ind)
    close(nodesfl)
    
    if(MDS_opt==2) sys_call_out = systemqq('Boost_dijkstra.exe') !calculate the distances to the landmark points
    
    !read in the distances from the file:
    !use the coord_ISOMAP array, resize it to hold the distances to data
    if(allocated(coord_ISOMAP)) deallocate(coord_ISOMAP)
    
    !only need the distances within a search radius, for each location, keep the ndmax, and the indicator for which point that is:
    allocate(coord_ISOMAP(NODES,ndmax),stat = test) 
    if(test.ne.0)then
          write(*,*)'ERROR: Allocation for inverse distance failed due to insufficient memory.'
          stop
    end if

    allocate(coord_ISOMAP_ind(NODES,ndmax),stat = test)
    if(test.ne.0)then
          write(*,*)'ERROR: Allocation for inverse distance failed due to insufficient memory.'
          stop
    end if

    
    allocate(temp_coord_ISOMAP(NODES),stat = test)
    if(test.ne.0)then
          write(*,*)'ERROR: Allocation for inverse distance failed due to insufficient memory.'
          stop
    end if

    coord_ISOMAP=1e21
    coord_ISOMAP_ind=1
    open(nodesfl,file='dist_cpp.out',status='unknown')
    
    do j=1,nd
        read(nodesfl,*) temp_coord_ISOMAP
        
        !for this data location, loop over the grid, all locatoins where the distance>search, store it as close:
        do ii=1,NODES
            if(temp_coord_ISOMAP(ii) < rad_inv_dis ) then
                max_dist = 0.9e21
                do jj=1,ndmax
                    !is the proposed distance shorter than what is currently sotred...
                    !what is the largest distance here:
                    if(coord_ISOMAP(ii,jj) > max_dist) then
                        max_dist = coord_ISOMAP(ii,jj) 
                        ind55 = jj
                    end if
                end do

                !if the proposed dist is shoerter, then store
                if(coord_ISOMAP(ii,ind55) > temp_coord_ISOMAP(ii) ) then
                  coord_ISOMAP(ii,ind55) = temp_coord_ISOMAP(ii)
                  coord_ISOMAP_ind(ii,ind55) = j
                end if
                
            end if
        
        end do
        
    end do
    close(nodesfl)
    
end subroutine get_inv_dists



