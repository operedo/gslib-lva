

!-------------------------------------------------------------------
!-------------------------------------------------------------------


module graph_vars
    implicit none
    integer, allocatable, dimension (:) :: KF  !FLIST,
!    real*8, allocatable, dimension (:) :: DFLIST
    integer NODES
    integer*8 EDGES
    real*8 refine
    integer n_close_data,ldbg2,debug_dist,graph_offset,data_node,LANDNODES,DATANODES
end module graph_vars


module graph
  use global
  implicit none
  
  contains
  
subroutine set_graph(MDS_opt)
    !initializes the graph, determines nodes and links between nodes
    !if there is enough memory a template array is used to speed up
    !the distance calcualtion
    
    !this program also writes out the distances to a grid.out file so that
    !the c++ program can be used to cal all of the distances (Dijkstra with BOOST)
    use graph_vars
    use global
    use grid_info
    use class_lvaf2
    
    include  'sgs.inc'
    
    integer mds_opt,lpsout,kkk,kkkk,ind,cnt,no_cpp
    real*8 andist,adist
    integer i,j,k,ii,jj,kk,ix,iy,iz,nx2,ny2,nz2,indx,indy,indz
    
    real*8 sang1,sang2,sang3,sanis1,sanis2,na,ndb
    integer isrot
    
    
    integer old_max_seg,inflag,cur_node,CUR_EDGE,Pmin(3),Pmax(3),use_template
    real*8 p1(3),p2(3),junk
    real*8 xsiz2,ysiz2,zsiz2,ymn2,zmn2,xmn2,dist,ang(3),rat(2),old_ang(3),old_rat(2),testd
    integer xtest,ytest,ztest    
    real(kind=8), dimension ( 3 ) :: ba !Vector b - a
    real(kind=8), dimension (3,3) :: R  !Anisotropic rotation matrix
    real(kind=8), dimension ( 3 ) :: v  !copy of ba
    integer, allocatable, dimension (:,:,:) :: redund_temp !is the offset redundant
    real*8, allocatable :: temp_graph_LLE(:,:,:,:) !is the offset redundant
    integer, allocatable :: cnt_LLE(:,:) !is the offset redundant
    integer sizeG(10),curNODE_G
    real*8 LLE_node(4)

   !make some temprary files
        ldbg2=987
        open(ldbg2,file='grid.out',status='unknown')
    
    nx2=nx*refine;ny2=ny*refine;nz2=nz*refine
    xsiz2=xsiz/refine; ysiz2=ysiz/refine; zsiz2=zsiz/refine
    xmn2=xmn/refine; ymn2=ymn/refine; zmn2=zmn/refine
    allocate(redund_temp(-graph_offset:graph_offset,-graph_offset:graph_offset,-graph_offset:graph_offset))
    n_close_data=0
    
add_paths_land=0   
add_paths_data=0 

    LANDNODES=0
    DATANODES=0
    
    if(add_paths_land ==1 ) LANDNODES=xyzland*nx2*ny2*nz2
    if(add_paths_data ==1 ) DATANODES=ndata*nx2*ny2*nz2
    
    if(MDS_opt ==3) then
        deallocate(redund_temp)
        NODES=nx2*ny2*nz2
        return 
    end if
    
    NODES=nx2*ny2*nz2
        
    old_max_seg=max_seg !have to use the straightline dist, so set the max_seg=0 for john's min dist solution
    max_seg=0
    write(*,*) ; write(*,*)  ; write(*,*) 
    write(*,*) 'calcualting the distances between nodes for the graph theory implementation'
    write(*,*) ; write(*,*)  ; write(*,*) 
    
    !lookup way for the distances
    !need a normalizing constant for each offset, depending on number of refinements.
    
    i=graph_offset
    ii=graph_offset
    if(nz2==1) ii=1

    use_template=1 !make =0 if template allocation fails
    allocate(norm_D(-i:i,-i:i,-ii:ii),stat = test)
    if(test /= 0 ) then
        use_template=0
        write(*,*) '!!!!WARRNING!!!'
        write(*,*) 'could not allocate an array big enought for distance offsets'
        write(*,*) 'this will slow down the calcualtion of the STRAIGHTLINE distances in the graph'
        if(allocated(template_D)) deallocate(template_D)
        if(allocated(norm_D)) deallocate(norm_D)
    end if
    
    allocate(template_D(-i:i,-i:i,-ii:ii,nx2,ny2,nz2),stat = test)
    if(test /= 0 ) then
        use_template=0
        write(*,*) '!!!!WARRNING!!!'
        write(*,*) 'could not allocate an array big enought for distance offsets'
        write(*,*) 'this will slow down the calcualtion of the STRAIGHTLINE distances in the graph'
        if(allocated(template_D)) deallocate(template_D)
        if(allocated(norm_D)) deallocate(norm_D)
    end if
    
    if(use_template==1) then
        norm_D=0
        kkk=(graph_offset)   
        kkkk= (graph_offset)   
        if(nz2==1) kkkk=0

            !now loop over all possible offsets for each grid node
        old_ang=0
        old_rat=0
        do kk=1,nz2
        do jj=1,ny2
        do ii=1,nx2
        !get the grid parameters at this location
            ang = grid(ii,jj,kk,1:3 )
            rat = grid(ii,jj,kk, 4:5 )
            call get_rmatrix ( ang, rat, R )
                do k=-1*kkkk,1*kkkk
                do j=-1*kkk,1*kkk
                do i=-1*kkk,1*kkk
                    !what is the direction?
                    ba(1)=i*xsiz ; ba(2) = j*xsiz ; ba(3) = k*xsiz
                    v = matmul (  R,ba )
                    template_D(i,j,k,ii,jj,kk) = ( dot_product ( v, v ) )**(.5) !* norm_D(i,j,k)
                end do
                end do
                end do
                old_ang=ang
                old_rat=rat
        end do
        end do
        end do
    end if !if using the template to cal distances
    
!now need a template if the offset is is more than one
!some of the links in the offset>1 are redundant
!will make the overall graph smaller if we ignore these links
redund_temp=0
do k=-1*kkk,1*kkk
do j=-1*kkk,1*kkk
do i=-1*kkk,1*kkk
    !is this offset redundant with one of the first offsets?
    do ind=2,graph_offset
    !check all original offsets
    do kk=-1,1
    do jj=-1,1
    do ii=-1,1
        if (ii*ind ==i .and. jj*ind ==j .and. kk*ind==k) redund_temp(i,j,k)=1
    end do
    end do
    end do
    
    end do
end do
end do
end do
redund_temp(0,0,0)=1
    !count number of edges we need
    do k=1,nz2
    do j=1,ny2
    do i=1,nx2
        !which nodes are we attached to (if 2D)(check 26 possible nodes)
        do kk=-kkkk,kkkk
        do jj=-kkk,kkk
        do ii=-kkk,kkk
            !are we inside the model still
            ix=i+ii
            iy=j+jj
            iz=k+kk
            if(ix>0 .and. ix<=nx2 .and.   iy>0 .and. iy<=ny2 .and.   iz>0 .and. iz<=nz2 .and. redund_temp(ii,jj,kk)==0) EDGES=EDGES+1 !we are in the model
        end do
        end do
        end do
    end do
    end do
    end do
    
    EDGES=EDGES



    write(ldbg2,'(I10,I10)') NODES,EDGES+LANDNODES+DATANODES !/2 !THIS FILE NEEDS NUMBER OF EDGES AND NODES
    
    allocate(KF(nx2*ny2*nz2+1+ndata),stat = test)
        if(test.ne.0)then
              write(*,'(a,I10,a)') 'ERROR: Allocation failed your graph is too large and requires ',EDGES,' edges'
              stop
        end if
    
    KF=0
    cur_node=0
    CUR_EDGE=0
    
    do k=1,nz2
    do j=1,ny2
    do i=1,nx2
        cur_node=cur_node+1 !our base node
        if(cur_node /= get_node(i,j,k,nx2,ny2,nz2) ) then !just a check
            write(*,*) get_node(i,j,k,nx2,ny2,nz2)
            write(*,*)
        end if
        p1(1) = xmn2 + real(i-1)*xsiz2
        p1(2) = ymn2 + real(j-1)*ysiz2
        p1(3) = zmn2 + real(k-1)*zsiz2
    
        !which nodes are we attached to (if 3D)(check 26 possible nodes)
        do kk=-kkkk,kkkk
        do jj=-kkk,kkk
        do ii=-kkk,kkk
            !are we inside the model still
            ix=i+ii; iy=j+jj; iz=k+kk
            if(ix>0 .and. ix<=nx2 .and.   iy>0 .and. iy<=ny2 .and.   iz>0 .and. iz<=nz2 .and. redund_temp(ii,jj,kk)==0) then  !we are in the model
                    CUR_EDGE=CUR_EDGE+1
                    p2(1) = xmn2 + real(ix-1)*xsiz2
                    p2(2) = ymn2 + real(iy-1)*ysiz2
                    p2(3) = zmn2 + real(iz-1)*zsiz2
                    edge_NODE=get_node(ix,iy,iz,nx2,ny2,nz2)
                    KF(cur_node+1)=CUR_EDGE                   
                    if(use_template==1) then
                        Pmin(1)=i ; Pmin(2)=j ; Pmin(3)=k ; Pmax=Pmin
                        if(Pmin(1)>ix) Pmin(1)=ix
                        if(Pmin(2)>iy) Pmin(2)=iy
                        if(Pmin(3)>iz) Pmin(3)=iz
                        if(Pmax(1)<ix) Pmax(1)=ix
                        if(Pmax(2)<iy) Pmax(2)=iy
                        if(Pmax(3)<iz) Pmax(3)=iz
                        cnt=sum(Pmax-Pmin)
                        edge_dist= andist(p1,p2,Pmin,Pmax,int(cnt+3),xsiz2,ysiz2,zsiz2,ymn2,zmn2,xmn2,ii,jj,kk)
                    else
                        edge_dist=adist(p1,p2,junk)
                        
                    end if
 !                   to write out the graph: 
                     write(ldbg2,'(I8,xx,I8,xx,g15.9)' ) cur_node,edge_NODE,edge_dist
                
            end if !check if we are in the model still
        end do
        end do
        end do
    end do
    end do
    end do

    no_cpp=0
    cnt=0
    
      

    write(*,*) ; write(*,*)  ; write(*,*) 
    write(*,'(a,I10,a,I10,a)') 'DONE calcualting the distances between ',NODES, ' nodes, with ', EDGES,' edges'
    write(*,*) ; write(*,*)  ; write(*,*) 
    
    max_seg=old_max_seg !reset the max segment to the user defined value
    close(ldbg2)

if(use_template==1) then
    !dont need these anymore
    deallocate(template_D)
    deallocate(KF)
    if(allocated(redund_temp)) deallocate(redund_temp)
    if(allocated(norm_D)) deallocate(norm_D)
    if(allocated(template_D)) deallocate(template_D)
end if

deallocate(grid)
    
end subroutine set_graph




integer function get_node(i,j,k,nx2,ny2,nz2)
!gets the node index given the model location and the model parameters
integer i,j,k,nx2,ny2,nz2,lpsout
    get_node=i+(j-1)*nx2 + (k-1)*nx2*ny2
end function get_node



end module graph


 
real*8 function adist(a,b,temp_D)
!calculate the anisotropic distance
    use class_lvaf2
    use global
    implicit none
    real(kind=8), dimension ( 3 ) :: a, b
    real(kind=8) :: temp_D
    
    n_seg=0 !number of times line has been segmented
    adist = min_dsqrd ( a, b )
end function adist
!-------------------------------------------------------------------
!-------------------------------------------------------------------

 
real*8 function andist(a,b,Pmin,Pmax,inter,xsiz2,ysiz2,zsiz2,ymn2,zmn2,xmn2,ii,jj,kk)
!calculate the anisotropic distance (jeff)
    use class_lvaf2
    use global
    use grid_info
    USE IFPORT
    
    implicit none
    real(kind=8), dimension ( 3 ) :: a, b
    real(kind=8) :: temp_D
    integer ii,jj,kk !the offsets for this particular line
    integer Pmin(3),Pmax(3),ind(3),inter  !min and max loc of the interseactions, inter=max number of intersections
    integer cnt,i,j,k
    real*8 line(3),Tline(inter+2),locc,xsiz2,ysiz2,zsiz2,ymn2,zmn2,xmn2,d,length,p(3)
    logical outside
  
    
    !get eqn of the line from point a to be (starting at b?
    line=(b-a)
    !how many intersections do we need to consider?
    !check all X interseactions (from Pmin -> Pmax)
    cnt=1
    Tline(1)=0
        !get the value of t along the line define in (line)
        !get the actual x,y locations for this interseaction 
    if(line(1) /=0) then   
    do i=Pmin(1),Pmax(1)
        cnt=cnt+1
        locc=xmn2 + real(i-1)*xsiz2  + xsiz2/2 !actually want the block intersection, no the block center
        Tline(cnt)=(b(1)-locc)/line(1)
        if(Tline(cnt) <= 0 .or. Tline(cnt) >= 1 ) cnt=cnt-1
    end do
    end if
    
    if(line(2) /=0) then  
    do j=Pmin(2),Pmax(2)
        cnt=cnt+1
        locc=ymn2 + real(j-1)*ysiz2  + ysiz2/2
        Tline(cnt)=(b(2)-locc)/line(2)
        if(Tline(cnt) <= 0 .or. Tline(cnt) >= 1 ) cnt=cnt-1
    end do
    end if
    
    if(line(3) /=0) then  
    do k=Pmin(3),Pmax(3)
        cnt=cnt+1
        locc=zmn2 + real(k-1)*zsiz2  + zsiz2/2
        Tline(cnt)=(b(3)-locc)/line(3)
        if(Tline(cnt) <= 0 .or. Tline(cnt) >= 1 ) cnt=cnt-1
    end do
    end if
    cnt=cnt+1
    Tline(cnt)=1
    Call SORTQQ (LOC(Tline), cnt, SRT$REAL8)
    
    ! for all positive intersections along the line T split it and calcualte the distance
    ! in that block

    andist=0
    do i=2,cnt !there are cnt+1 segments of the line
        !how long is the segment in this block:
        length= Tline(i) - Tline(i-1)
        if(length>0) then 
            !which block are we in? Need xyz loc for center of segment:
            p=a+line*Tline(i-1) +line*length/2
            call get_index (p, ind, outside)
            !get the aniso dist of this segment from the template lookup
            andist=andist+template_D(ii,jj,kk,ind(1),ind(2),ind(3) ) * length !*norm_D(ii,jj,kk) !length is the % of the line, template has the actual line distance
        end if
        
        
    
    end do
    
    
    
    
    
end function andist
!-------------------------------------------------------------------
!-------------------------------------------------------------------
