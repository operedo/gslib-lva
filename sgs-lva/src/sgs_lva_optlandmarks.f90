


!-----------------------------------------------------------------------
!
!this version of the program uses the graph centered on the data and gets all distances from the data to the grid,
!rather than visiting each node and then getting the close data
! Title of Paper Submitted: Programs for Kriging and Sequential Gaussian Simulation with Locally Varying Anisotropy Using Non-Euclidean Distances
!
! Authors of paper:
! J. B. BOISVERT
! C. V. DEUTSCH


!-----------------------------------------------------------------------
!
!                     KRIGING WITH LOCALLY VARYING ANISOTROPY
!                     *****************************************
! This program uses the optimized distance between points and kriges a grid.
! 
!
!
! AUTHOR: Jeff Boisvert                             DATE: 2009
! MODIFIED FROM:
! AUTHOR (GSLIB): Clayton V. Deutsch                DATE: 1989-1999
!-----------------------------------------------------------------------


! Module to declare dynamic arrays in multiple subroutines:
!
!
   module geostat
      integer,allocatable :: nisb(:),ixsbtosr(:),iysbtosr(:),izsbtosr(:),path(:)
      real*8,allocatable    :: x(:),y(:),z(:),vr(:),ve(:),dh(:),tmp(:), &
               close(:),close2(:),xa(:),ya(:),za(:),vra(:,:),vea(:),xdb(:),ydb(:), &
               zdb(:),cut(:),cdf(:), temp_x(:), temp_y(:), temp_z(:), &
               temp_xdb(:),temp_ydb(:),temp_zdb(:),dist_data(:,:), &
               lambda(:),grid2grid(:,:),grid2grid_MDS(:,:),grid2grid_LLE(:,:)
      real*8 asum,mmm,radsqd
      real*8,allocatable  :: r(:),rr(:),s(:),a(:),bb(:,:),jsum(:),isum(:),aaa(:,:),atemp(:,:)
      integer ng_nodes,grid_NODES,debug_node,cnt_lambda(100),minE,ID1,ID2
      integer idijkstra
      
   end module
!
      program main
!-----------------------------------------------------------------------
!
!             Kriging (SK,OK,KT) of a 3-D Rectangular Grid
!             ********************************************
!
! The program is executed with no command line arguments.  The user
! will be prompted for the name of a parameter file.
!
!
!
! AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
!-----------------------------------------------------------------------
      use geostat
      use global
      use graph_vars
      implicit none
      include  'sgs.inc'
      character outfl*512
      real*8 t_total,sang1,sang2,sang3,sanis1,sanis2
      integer maxdis,maxsby,maxsbx,maxsbz,isrot,na,ndb
      
      data_node=-1
      n_indef=0; 
      na_total=0
      open(T_eigs,file='time_eig.out',status="unknown")
      open(T_matrix,file='time_matrix.out',status="unknown")
      open(T_total,file='time_total.out',status="unknown")
      open(T_search,file='time_search.out',status="unknown")
      
      
      
      open(time_out,file='time.out',status="unknown")
      open(time_out2,file='time2.out',status="unknown")
      
      total=secnds(real(0.0,4))
!
! Read the parameters, the data, and open the output files:

      call readparm(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,outfl)
!
! Call sgs to simulate the grid:
      call sgs(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,outfl)
!
! Finished:
!
      close(ldbg)
      close(lout)


      !write out the time to run kriging
      elapsed=secnds(total)
      write(T_total,'(4I12,f19.5)') NODES, edges, dim, nd,elapsed !subtract time to search

      write(time_out,'(25f16.8)') real(ndmax),real(aa(1)),real(elapsed),real(n_indef),real(na_total)/real(na_cnt),real(n_close_data)/2

      
      write(*,9998) VERSION
 9998 format(/' SGS_LVA Version: ',f5.3, ' Finished'/)
      stop

      end program main 
 
      subroutine readparm(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,outfl)
!-----------------------------------------------------------------------
!
!                  Initialization and Read Parameters
!                  **********************************
!
! The input parameters and data are read in from their files. Some quick
! error checking is performed and the statistics of all the variables
! being considered are written to standard output.
!
!
!
!-----------------------------------------------------------------------
!      use       msflib
      use       geostat
      use       grid_info ! for grid 
      use       global    ! for some global variables
      use       graph     !for using graph theory representation
      use       graph_vars
      implicit none
      include  'sgs.inc'
      integer mv,lin,gin,no_hash,i,j,k,idhl,ixl,iyl,izl,ivrl,iextv,maxcut,ia1g,ia2g,ia3g,ir1,ir2,jj,locx,locy,locz,nvari,maxdis,maxsam,maxsbx,maxsby,maxsbz
      integer maxsb,mxsx,mxsxy,maxdat,ind_test,isrot,na,ndb
      real*8 xmax1,ymax1,zmax1,xmax2,ymax2,zmax2,xx,yy,zz,sanis1,sanis2,sang1,sang2,sang3,aa1,aa2,av,ss,vrt
      parameter(MV=200)
      real*8      var(MV)
      character datafl*512,jackfl*512,extfl*512,outfl*512,dbgfl*512,str*512,title*80,grid_fl*512
      logical   testfl,inflag
      real(kind=8), dimension ( 3 ) :: angles
      integer seg_opt
   
! FORTRAN Units:
!
      lin   = 1
      ldbg  = 3
      lout  = 4
      gin   = 5 !file for the grid
      lext  = 7
      ljack = 8
      
      LLE_opt=1
      no_hash=0
!
! Note VERSION number:
!
      write(*,9999) VERSION
 9999 format(/' SGS_LVA Version: ',f5.3/)
!
! Get the name of the parameter file - try the default name if no input:
!
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'sgs_lva.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'sgs_lva.par            ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
!
! Find Start of Parameters:
!
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
!
! Read Input Parameters:
!
      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)
      idhl=0
      read(lin,*,err=98) ixl,iyl,izl,ivrl
      iextv=0
      write(*,'(a20,6(I4))') ' columns = ',idhl,ixl,iyl,izl,ivrl

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

koption=0
!
! This is an undocumented feature to have sgs construct an IK-type
! distribution:
!
      iktype = 0
      if(koption.lt.0) then
            iktype  = 1
            koption = -koption
      end if
      if(iktype.eq.1) then

            read(lin,*,err=98) ncut
            write(*,*) ' number of cutoffs = ',ncut
!
! Find the needed parameter:
!
            MAXCUT = ncut
!
! Allocate the needed memory:
!21
            allocate(cut(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to', &
                             ' insufficient memory.'
                        stop
                  end if
!22
            allocate(cdf(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to', &
                             ' insufficient memory.'
                        stop
                  end if
!
            read(lin,*,err=98) (cut(i),i=1,ncut)
            write(*,*) ' cutoffs = ',(cut(i),i=1,ncut)

      end if

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)
      
      
      read(lin,*,err=98) nreal
      write(*,*) ' n realizations:',nreal

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' nx, xmn, xsiz = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' ny, ymn, ysiz = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' nz, zmn, zsiz = ',nz,zmn,zsiz
      
      read(lin,*,err=98) seedd
      write(*,*) ' random number seed:',seedd

      
      nxdis=1
      nydis=1
      nzdis=1
      
      read(lin,'(a)',err=98) grid_fl
      call chknam(grid_fl,512)
      write(*,*) ' grid file to use:',grid_fl(1:40)
      
      read(lin,*,err=98) ia1g,ia2g,ia3g,ir1,ir2 !read in col for angles (1-3) and ratios (1-2)
      write(*,'(a20,5(I4))') ' columns for LVA grid = ',ia1g,ia2g,ia3g,ir1,ir2

      read(lin,*,err=98) n(1),o(1),sss(1)
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' ngx, xgmn, xgsiz = ',n(1),o(1),sss(1)

      read(lin,*,err=98) n(2),o(2),sss(2)
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' ngy, ygmn, ygsiz = ',n(2),o(2),sss(2)

      read(lin,*,err=98) n(3),o(3),sss(3)
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' ngz, zgmn, zgsiz = ',n(3),o(3),sss(3)

      refine=1
      
      read(lin,*,err=98) graph_offset
      if(graph_offset/=int(graph_offset) .or. graph_offset<1) then
          write(*,*) 'ERROR WITH THE graph_offset VALUE, SETTING graph_offset TO 1'
          graph_offset=1
      end if
      write(*,*) ' offset parameter',graph_offset
      
      read(lin,*,err=98) MDS_opt
      
      
      cal_stress=0
      if(MDS_opt <0) then
      
      write(*,*) '****will calculate stress******'
      cal_stress = 1
      MDS_opt=-MDS_opt
      end if
      
      
      write(*,*) ' Use MDS:',MDS_opt
      
      read(lin,*,err=98) xland,yland,zland
      write(*,*) ' spacing of landmark points:',xland,yland,zland
      
      read(lin,*,err=98) dim
      write(*,*) ' max number of dimens to use:',dim
      
      !quick check on grid coverage:
      xmax1=xmn+nx*xsiz;ymax1=ymn+ny*ysiz;zmax1=zmn+nz*zsiz
      xmax2=o(1)+n(1)*sss(1);ymax2=o(2)+n(2)*sss(2);zmax2=o(3)+n(3)*sss(3)
    !  if(xmax1>xmax2 .or. ymax1>ymax2 .or. zmax1>zmax2 .or. o(1)>xmn .or. o(2)>ymn .or. o(3)>zmn) &
    !  write(*,*)'**WARNING** LVA GRID MAY NOT COVER KRIGING AREA, MAY CAUSE ERROR'
   

      !read in the grid and store the necessary angles etc.
     inquire(file=grid_fl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the grid file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            stop
      endif
     open(gin,file=grid_fl,status='OLD')

     read(gin,*)
      read(gin,*,err=99)       nvari
      do i=1,nvari
            read(gin,*)
      end do
      
       allocate(grid(nx,ny,nz,5))


      if(nx /=n(1) .or. ny /= n(2) .or. nz /= n(3) .or. xmn /= o(1) .or. ymn /= o(2) .or. zmn /= o(3) .or. xsiz /= sss(1) .or. ysiz /= sss(2) .or. zsiz /= sss(3)  ) then
      write(*,*) 
      write(*,*) 
      write(*,*) '*******WARNING*******'
      write(*,*) 'simulation grid does not match LVA grid'
      write(*,*) 'will regrid the LVA grid to match the simulation grid'
      write(*,*) 
      write(*,*) 
      temp_o = o
      temp_sss = sss
      temp_n = n
      allocate(temp_grid(n(1),n(2),n(3),5))
      n(1)=nx ; n(2)=ny ; n(3)=nz
      o(1)=xmn ; o(2)=ymn ; o(3)=zmn
      sss(1) = xsiz ;sss(2) = ysiz ;sss(3) = zsiz
      
      
     
         do k=1,temp_n(3)
          do j=1,temp_n(2)
           do i=1,temp_n(1)
            read(gin,*,err=98)  (var(jj),jj=1,nvari)
            temp_grid(i,j,k,1)=var(ia1g)
            temp_grid(i,j,k,2)=var(ia2g)
            temp_grid(i,j,k,3)=var(ia3g)
            temp_grid(i,j,k,4)=var(ir1)
            temp_grid(i,j,k,5)=var(ir2)
           end do
          end do
         end do
         
         do k=1,n(3)
          do j=1,n(2)
           do i=1,n(1)
           
           !where are we in the model (xx,yy,zz)
           
            xx = xmn + real(i-1)*xsiz
            yy = ymn + real(j-1)*ysiz
            zz = zmn + real(k-1)*zsiz
            !where is this in the LVA grid
            call getindx(temp_n(1),temp_o(1),temp_sss(1),xx,locx,inflag)
            call getindx(temp_n(2),temp_o(2),temp_sss(2),yy,locy,inflag)
            call getindx(temp_n(3),temp_o(3),temp_sss(3),zz,locz,inflag)
            
            !fill in the grid
            grid(i,j,k,:) = temp_grid(locx,locy,locz,:)
            
           end do
          end do
         end do
        deallocate(temp_grid)
     
     else !read in per normal
     
     
     do k=1,n(3)
      do j=1,n(2)
       do i=1,n(1)
        read(gin,*,err=98)  (var(jj),jj=1,nvari)
        grid(i,j,k,1)=var(ia1g)
        grid(i,j,k,2)=var(ia2g)
        grid(i,j,k,3)=var(ia3g)
        grid(i,j,k,4)=var(ir1)
        grid(i,j,k,5)=var(ir2)
       end do
      end do
     end do
     
     end if !loop over regrid
     
     do k=1,n(3)
      do j=1,n(2)
       do i=1,n(1)
        !if necessary, fix the three angles read in
        if(k==1 .and. j==1 .and. i==1  .and. grid(i,j,k,4)==0) then
        write(*,*) 'WARNING***** you cannot have an anisotropy ratio =0, it will be reset to 1:1****'
        end if
        
        if(k==1 .and. j==1 .and. i==1  .and. grid(i,j,k,5)==0) then
        write(*,*) 'WARNING***** you cannot have an anisotropy ratio =0, it will be reset to 1:1****'
        end if
        
        if(grid(i,j,k,4)==0) grid(i,j,k,4)=1
        if(grid(i,j,k,5)==0) grid(i,j,k,5)=1
        
        angles(1) = grid(i,j,k,1) ; angles(2) = grid(i,j,k,2) ; angles(3) = grid(i,j,k,3) 
        call fix_angles  ( angles )
        grid(i,j,k,1) = angles(1) ; grid(i,j,k,2) = angles(2) ; grid(i,j,k,3) = angles(3) ; 
       end do
      end do
      
     end do
     close (gin)
     assign_nodes=1

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' ndmin,ndmax = ',ndmin,ndmax
      noct=0

      read(lin,*,err=98) radius
      write(*,'(a,f18.4)') ' search radii = ',radius
      
      read(lin,*,err=98) d_tree
      write(*,'(a,I4)') ' search dim to use = ',d_tree
      
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius / radius
      sanis2 = radius / radius
      sang1=0
      sang2=0
      sang3=0

      read(lin,*,err=98) power,skmean
      backspace(lin)
      read(lin,*,err=98) ktype,skmean


      write(*,'(a,I4,xx,f8.4)') ' ktype, skmean =',ktype,skmean
      
      inv_dist=.false.
            
      read(lin,*,err=98) nst(1),c0(1)
      write(*,'(a,I4,xx,f8.4)') ' nst, c0 = ',nst(1),c0(1)

      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/, &
                  ' The c or a values can be set to zero')
            stop
      endif

      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),aa(i)
            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            write(*,'(a,I4,xx,f8.4,xx,f8.4,xx,f8.4,xx,f8.4)') ' it,cc,ang[1,2,3]; ',it(i),cc(i), &
                        ang1(i),ang2(i),ang3(i)
            write(*,'(a,xx,f8.4,xx,f8.4,xx,f8.4)') ' a1 a2 a3: ',aa(i),aa1,aa2
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
      end do

      read(lin,*,err=98) idijkstra
      write(*,'(a,I4)') ' flag to use multi-thread Boost Dijkstra = ',idijkstra

      close(lin)
!
! Find the needed parameters:
!
      MAXDIS = nxdis*nydis*nzdis
      MAXSAM = ndmax + 1
      MAXEQ = MAXSAM + MAXDT + 2
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2.00)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
!
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2.00)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
!
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2.00)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
!
      MAXSB = MAXSBX*MAXSBY*MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX  = 2 * MAXSBX
!
! Allocate the needed memory:
!1
      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                       ' insufficient memory.'
                  stop
            end if
!2
      allocate(ixsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!3
      allocate(iysbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!4
      allocate(izsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!13
      allocate(xa(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if

      allocate(ya(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!15
      allocate(za(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!16
      allocate(vra(MAXSAM,nreal),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!17
      allocate(vea(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!18
      allocate(xdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!19
      allocate(ydb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!20
      allocate(zdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!23
      allocate(r(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
      allocate(r_tmp(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if       
            
!24
      allocate(rr(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!25
      allocate(s(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!26
      allocate(a(MAXEQ * MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!
! Perform some quick error checking:
!
      if(ndmax.gt.MAXSAM) stop 'ndmax is too big - modify .inc file'
      if(ktype.eq.3.and.iextv.le.0) stop 'must have external variable'
      if(ixl.le.0.and.nx.gt.1) write(*,*) ' WARNING: ixl=0 and nx>1 ! '
      if(iyl.le.0.and.ny.gt.1) write(*,*) ' WARNING: iyl=0 and ny>1 ! '
      if(izl.le.0.and.nz.gt.1) write(*,*) ' WARNING: izl=0 and nz>1 ! '
!
! Check to make sure the data file exists, then either read in the
! data or write an error message and stop:
!
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR data file ',datafl,' does not exist!'
            stop
      endif
!
! The data file exists so open the file and read in the header
! information. Initialize the storage that will be used to summarize
! the data found in the file:
!
      title(1:42) = 'Sequential Gaussian Simulation with LVA: '
      open(lin,file=datafl,status='OLD')
      read(lin,*)
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,*)
      end do
      MAXDAT = 0 ; ndata=0
 22   read(lin,*,end=33,err=99) (var(j),j=1,nvari)
      if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) go to 22
      MAXDAT = MAXDAT + 1
      ndata=ndata+1
      go to 22
 33   continue
 
!
! Allocate the needed memory:
!5
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!6
      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!7
      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!8
      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!9
      allocate(ve(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!10
      allocate(dh(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!11
      allocate(tmp(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!12
      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!
!
 allocate(data_ind(MAXDAT))

      rewind(lin)
      read(lin,'(a)') title(42:80)
      read(lin,*,err=99)       nvari
      nd = 0
      av = 0.0
      ss = 0.0
      do i=1,nvari
            read(lin,'(a40)',err=99) str
      end do
!
! Some tests on column numbers:
!
      if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.ivrl.gt.nvari) &
           then
            write(*,*) 'There are only ',nvari,' columns in input data'
            write(*,*) '  your specification is out of range'
            stop
      end if
!
! Read all the data until the end of the file:

    allocate(x_loc(MAXDAT))
    allocate(y_loc(MAXDAT))
    allocate(z_loc(MAXDAT))


!
 2    read(lin,*,end=3,err=99) (var(j),j=1,nvari)
      vrt = var(ivrl)
      if(vrt.lt.tmin.or.vrt.ge.tmax) go to 2
      !if we have a data at this node already go to 2:
      
      
        call getindx(nx,xmn,xsiz,var(ixl),xind,inflag)
         if(not(inflag)) go to 2
        call getindx(ny,ymn,ysiz,var(iyl),yind,inflag)
         if(not(inflag)) go to 2
        zind=1
        if(izl>0) call getindx(nz,zmn,zsiz,var(izl),zind,inflag)
         if(not(inflag)) go to 2
        ind_test=get_node(xind,yind,zind,nx,ny,nz)
         do j=1,nd !check data so far
            if(data_ind(j)== ind_test) then !skip it
                go to 2
            end if  
         end do

      nd = nd + 1
      if(nd.gt.MAXDAT) then
            write(*,*) ' ERROR: Exceeded available memory for data'
            stop
      end if
!
! Establish the location of this datum:
!
      if(idhl.le.0) then
            dh(nd) = -99
      else
            dh(nd) = var(idhl)
      endif
      if(ixl.le.0) then
            x(nd) = xmn
      else
            x(nd) = var(ixl)
      endif
      if(iyl.le.0) then
            y(nd) = ymn
      else
            y(nd) = var(iyl)
      endif
      if(izl.le.0) then
            z(nd) = zmn
      else
            z(nd) = var(izl)
      endif

!
! Establish the external drift variable (if needed):
!
      vr(nd) = vrt
      av     = av + vrt
      ss     = ss + vrt*vrt
    
    !assign this data to nodes  
    call getindx(nx,xmn,xsiz,x(nd),xind,inflag)
    if(inflag) x(nd)=xmn+(xind-1)*xsiz
    
    call getindx(ny,ymn,ysiz,y(nd),yind,inflag)
    if(inflag) y(nd)=ymn+(yind-1)*ysiz
    
    call getindx(nz,zmn,zsiz,z(nd),zind,inflag)
    if(inflag) z(nd)=zmn+(zind-1)*zsiz
    
    data_ind(nd)=get_node(xind,yind,zind,nx,ny,nz)
    
    
    x_loc(nd) = xmn + real(xind-1)*xsiz
    y_loc(nd) = ymn + real(yind-1)*ysiz
    z_loc(nd) = zmn + real(zind-1)*zsiz

      go to 2
 3    close(lin)
 
!
! Compute the averages and variances as an error check for the user:
!
      av = av / max(real(nd),1.0)
      ss =(ss / max(real(nd),1.0)) - av * av
      write(*,*) 'Data for SGS: Variable number ',ivrl
      write(*,*) '  Number   = ',nd
      write(*,*) '  Average  = ',av
      write(*,*) '  Variance = ',ss
      if(nd.lt.1) then
            write(*,*) ' ERROR: there are no data'
            stop
      end if
!
! Open the debugging and output files:
!
    open(ldbg,file=dbgfl,status='UNKNOWN')
    open(lout,file=outfl,status='UNKNOWN')
      write(lout,'(a80)') title

           write(lout,201) 1,nx,ny,nz,nreal
           write(lout,'(a)') 'Simulated Values'

 201  format(5(1x,i4))

      if(iktype.eq.1) then
            if(koption.eq.0) then
                  write(lout,201) ncut,nx,ny,nz
            else
                  write(lout,201) ncut+1
            end if
            do i=1,ncut
                  write(lout,104) i,cut(i)
 104              format('Threshold: ',i2,' = ',f12.5)
            end do
            if(koption.eq.1) write(lout,105)
 105        format('true value')
      end if
!
! Open the external drift file if needed and position it at the
! first grid node in the file:
!
! Set up for cross validation:
!
      if(koption.eq.1) then
            jackfl = datafl
            idhlj  = idhl
            ixlj   = ixl
            iylj   = iyl
            izlj   = izl
            ivrlj  = ivrl
            iextvj = iextv
            open(ljack,file=jackfl,status="unknown")
            
            
              read(ljack,*)
              read(ljack,*,err=99)       nvarij
              do i=1,nvarij
                    read(ljack,*)
              end do
      end if
!
! Open the file with the jackknife data?
!

! Finished here:
!
      return
!
! Error in an Input File Somewhere:
!
 96   stop 'ERROR in jackknife file!'
 97   stop 'ERROR in external drift file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end subroutine readparm



      subroutine sgs(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,outfl)
!-----------------------------------------------------------------------
!
!                Simulate a 3-D Grid of Rectangular Blocks
!                **************************************
!
! This subroutine estimates point values of one variable by
! simple kriging and draw a stochastic value for the variance (line 1481)
!
!
!
! PROGRAM NOTES:
!
!   1. The data and parameters are passed in common blocks defined
!      in sgs.inc.  Local storage is allocated in the subroutine
!      for kriging matrices, i.e.,
!         - xa,ya,za,vra   arrays for data within search neighborhood
!         - a,r,rr,s       kriging arrays
!         - xdb,ydb,zdb    relative position of discretization points
!         - cbb            block covariance
!   2. The kriged value and the kriging variance is written to Fortran
!      unit number "lout".
!
!
!
!
! Original:  A.G. Journel and C. Lemmer                             1981
! Revisions: A.G. Journel and C. Kostov                             1984
!-----------------------------------------------------------------------
      use        geostat
      !use        dfport
      use        global !to get some global variables
      use        graph     !for using graph theory representation
      use        graph_vars
      use        kdtree2_module
      use        random2
      use        srch       !for the exhaustive serach
      implicit none
      include    'sgs.inc'
      real pmx
      integer mv,lin,gin,no_hash,i,j,idhl,ixl,iyl,izl,ivrl,iextv,maxcut,ia1g,ia2g,ia3g,ir1,ir2,jj,locx,locy,locz,nvari,maxdis,maxsam,maxsbx,maxsby,maxsbz
      integer maxsb,mxsx,mxsxy,maxdat,ind_test,isrot,na,ndb,no_cpp,cnt,indk,indi,ind,nsec,mdt,ix,iy,iz,indj,nk,nxy,nxyz,nloop,irepo,ivarg,istart,is,ist,nclose,ireal
      integer index2,indexmneq,ind1,ind2,kadim,ksdim,nrhs,nv,ie,ising,ierr,index,neq
      real*8 xmax1,ymax1,zmax1,xmax2,ymax2,zmax2,xx,yy,zz,sanis1,sanis2,sang1,sang2,sang3,aa1,aa2,av,ss,vrt,xdis,ydis,xloc,yloc,zloc,zdis,xk,vk,xkmae
      real*8 xkmse,ddh,cmax,secj,extest,estv,true

      character outfl*512
      

    integer, allocatable :: useablemap(:)

      real*8     cbb
      real*8       var(20),tdist(ndata)
      logical    first,fircon,accept,inflag
      data       fircon/.true./
      integer    pt_ID(2)
      character*2000 str
      integer(4) sys_call_out
      
      !declarations for search tree:
    real(kdkind), dimension(:,:), allocatable :: my_array
    real(kdkind), allocatable :: query_vec(:)
    type(kdtree2), pointer :: tree  ! this is how you declare a tree in your main program
    integer :: k
    type(kdtree2_result),allocatable :: results(:), resultsb(:)
    integer   :: nnbrute, rind
    real*8      :: t0, t1, sps, avgnum, maxdeviation
    real(kdkind) :: rv 
    integer, parameter  :: nnn = 5
    integer, parameter  :: nr2 = 5
    real*8 r2array(5)
    data r2array / 1.0e-4,1.0e-3,5.0e-3,1.0e-2,2.0e-2 / 
    integer   :: nnarray(5)
    data nnarray / 1, 5, 10, 25, 500/ 
    real*8 rand_num
    real*4, allocatable :: sim(:,:)
    real*8, allocatable :: temp_order(:,:)
    integer, allocatable :: order(:)
  

    integer :: NODES_LENGTH 
    integer :: GRID_OUT_LENGTH 
    integer :: NODES2CAL_LENGTH 
    integer, allocatable :: cur_edge_node_array1(:)
    integer, allocatable :: cur_edge_node_array2(:)
    real*8, allocatable :: edge_dist_array(:)
    integer, allocatable :: nodes2cal_array(:)

    real :: rndtest
    real*8 :: sumdiff, sumdiffsq, maxdiff, mindiff
    integer :: STATS_SAMPLE
 
    allocate(results(ndmax))


    refine=refine
    
    !open a second debuging file with some interesting results
    ldbg2 = 999
    open(ldbg2,file='indef_mat.out',status='unknown')
    write(ldbg2,*) '1=indefinate matrix, 2=poorly conditioned matrix'
    write(ldbg2,*) 1
    write(ldbg2,*) 'indicator only works with robust solver'
    
    PMX    = 999.0
    c1=MAXSBX
    c2=MAXSBY
    c3=MAXSBZ
    
    elapsed=secnds(total)
    write(*,'(f10.4,a)')      elapsed, 's - time to read pars'
    write(time_out2,'(f10.4,a)')      elapsed, 's - time to read pars'

    !where are our landmark points?  xyzland
    !evenly space them out
    xyzland=xland*yland*zland   
     
    allocate(landpts(xyzland),stat = test)
    if(test.ne.0)then
        write(*,*)'ERROR: Allocation failed due to', &
        ' insufficient memory.'
        stop
    end if
    allocate(evalues(xyzland),stat = test)
    if(test.ne.0)then
        write(*,*)'ERROR: Allocation failed due to', &
        ' insufficient memory.'
        stop
    end if
    allocate(vectors(xyzland,xyzland),stat = test)
    if(test.ne.0)then
        write(*,*)'ERROR: Allocation failed due to', &
        ' insufficient memory.'
        stop
    end if

    xspace=int(real(nx/xland))
    yspace=int(real(ny/yland))
    zspace=int(real(nz/zland))
    no_cpp=0
    cnt=0
    do k=1,zland
    indk=zspace*k - 0.5*zspace + 1
    do j=1,yland
    indj=yspace*j - 0.5*yspace + 1
    do i=1,xland
      cnt=cnt+1    
      indi=xspace*i - 0.5*xspace + 1
      ind=get_node(indi,indj,indk,nx,ny,nz)
      landpts(cnt)=ind
    end do
    end do
    end do

!*************************************************************************************
! Set up the graph for the input grid, calculates edge lengths
    !call set_graph(MDS_opt) 
    call set_graph_onlymem(MDS_opt,NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array) 
!*************************************************************************************
  

    elapsed=secnds(total)
    write(*,'(f10.4,a)')      elapsed, 's - time to cal edge lengths'
    write(time_out2,'(f10.4,a)')      elapsed, 's - time to cal edge lengths'

    nsec = 2
    mdt = 1
    if(ktype.eq.3) mdt = mdt + 1
    if(ktype.eq.0) mdt = 0
    if(ktype.eq.2) mdt = 0

! In all cases the offsets are relative to the lower left corner.
! This is done for rescaling the drift terms in the kriging matrix.

    if(nxdis.lt.1) nxdis = 1
    if(nydis.lt.1) nydis = 1
    if(nzdis.lt.1) nzdis = 1
    ndb=1
    xdis = xsiz  / max(real(nxdis),1.0)
    ydis = ysiz  / max(real(nydis),1.0)
    zdis = zsiz  / max(real(nzdis),1.0)
    i    = 0
    xloc = -0.5*(xsiz+xdis)
    do ix =1,nxdis
          xloc = xloc + xdis
          yloc = -0.5*(ysiz+ydis)
          do iy=1,nydis
                yloc = yloc + ydis
                zloc = -0.5*(zsiz+zdis)
                do iz=1,nzdis
                      zloc = zloc + zdis
                      i = i+1
                      xdb(i) = xloc + 0.5*xsiz
                      ydb(i) = yloc + 0.5*ysiz
                      zdb(i) = zloc + 0.5*zsiz
                end do
          end do
    end do
!
! Initialize accumulators:
!
    nk    = 0
    xk    = 0.0
    vk    = 0.0
    xkmae = 0.0
    xkmse = 0.0
! Report on progress from time to time:
!
    if(koption.eq.0) then
          nxy   = nx*ny
          nxyz  = nx*ny*nz
          nloop = nxyz
          irepo = max(1,min((nxyz/10),10000))
    else
          nloop = nx*ny*nz
          irepo = max(1,min((nd/10),10000))
    end if
    ddh = 0.0
      
    elapsed=secnds(total)
    write(*,'(f10.4,a)')      elapsed, 's - setup for ISOMAP'
!    write(time_out2,'(f10.4,a)')      elapsed, 's - time to cal edge lengths'


!
!this routine will fix the distances according to landmark (ISOMAP) multi dim scaling
!
    if(MDS_opt==2 .or. MDS_opt==3) then

        call init_genrand(seedd)  !intial the random generator
        do i=1,10000
            rand_num = grnd()
        end do
        
        allocate(temp_order(nx*ny*nz,2))
        !use sim array as a temp array for the order
        do i=1,nx*ny*nz
            temp_order(i,1) = grnd()
            temp_order(i,2) = dble(i)
        end do
    
        !call sortem(1,nx*ny*nz,temp_order(:,1),1,temp_order(:,2))   !now have the simulation order
        call sortem2dp(1,nx*ny*nz,temp_order(:,1),1,temp_order(:,2))   !now have the simulation order
    
        allocate(order(nx*ny*nz))
        order = int(temp_order(:,2))
        


        call get_landmark_pts_onlymem(ndmax,nd,nx,ny,nz,NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array) !use c++ program to get distances to landmark poitns
        call MDS_ISOMAP(ndmax,nd,nx,ny,nz) !do the multidimensinoal scaling

        
        STATS_SAMPLE = 50000
        do i=1,STATS_SAMPLE
           temp_order(i,1)=sqrt ( sum( ( coord_ISOMAP(order(nxyz-i),1:dim)-coord_ISOMAP(order(nxyz-i-1),1:dim) ) **2    )    )
           !write(*,*) order(nxyz-i),order(nxyz-i-1),temp_order(i,1)  ,coord_ISOMAP(order(nxyz-i),1:dim)  
        end do

        do NODES2CAL_LENGTH=1000,100,-50
        !NODES2CAL_LENGTH=950

        if(allocated(coord_ISOMAP)) deallocate(coord_ISOMAP)
        if(allocated(nodes2cal_array)) deallocate(nodes2cal_array)

        allocate(nodes2cal_array(NODES2CAL_LENGTH))
        nodes2cal_array = order(1:NODES2CAL_LENGTH) 
        
        call get_landmark_pts_onlymem_xyzland(ndmax,nd,nx,ny,nz,NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array,NODES2CAL_LENGTH,nodes2cal_array) !use c++ program to get distances to landmark poitns
        !call MDS_ISOMAP(ndmax,nd,nx,ny,nz) !do the multidimensinoal scaling
        call MDS_ISOMAP_xyzland(ndmax,nd,nx,ny,nz,NODES2CAL_LENGTH) !do the multidimensinoal scaling

        sumdiff=0.0
        sumdiffsq=0.0
        maxdiff=-1000000000.0
        mindiff= 1000000000.0
        do i=1,STATS_SAMPLE
           temp_order(i,2)=sqrt ( sum( ( coord_ISOMAP(order(nxyz-i),1:dim)-coord_ISOMAP(order(nxyz-i-1),1:dim) ) **2    )    )
           !write(*,*) order(nxyz-i),order(nxyz-i-1),temp_order(i,2)   
           sumdiff = sumdiff + temp_order(i,2)/temp_order(i,1)
           !sumdiffsq = sumdiffsq + (temp_order(i,2)/temp_order(i,1))*(temp_order(i,2)/temp_order(i,1)) 
           maxdiff = MAX(maxdiff,temp_order(i,2)/temp_order(i,1))
           mindiff = MIN(mindiff,temp_order(i,2)/temp_order(i,1))
        end do
        write(*,*) 'STATS: ',STATS_SAMPLE,NODES2CAL_LENGTH,(sumdiff/STATS_SAMPLE),mindiff,maxdiff

        end do

        deallocate(temp_order)

    end if

    elapsed=secnds(total)
    write(*,'(f10.4,a)')      elapsed, 's - time to finish multidimensional scaling'
!!    write(time_out2,'(f10.4,a)')      elapsed, 's - time to finish multidimensional scaling'
!
!!    start_time=secnds(0.0)
!!what is the max C?
!    ivarg = 1
!    cmax   = c0(ivarg)
!    istart = 1 + (ivarg-1)*MAXNST          
!    do is=1,nst(ivarg)
!          ist = istart + is - 1
!          if(it(ist).eq.4) then
!                cmax = cmax + PMX
!          else
!                cmax = cmax + cc(ist)
!          endif
!    end do
!
!
!!set up the kdtree for searching:
!    if(d_tree>dim .or. d_tree<1) then
!      d_tree=dim
!      write(*,'(a,I5)') 'actual number of dim to use in the search: ',d_tree
!    end if
!    
!    
!    nclose = min(nd,ndmax)
!  
!!set up the logical array for searching the tree:
!    allocate(is_usable(nx*ny*nz), useablemap(nx*ny*nz))
!
!!now set all the data locations to usable:
!    is_usable=.false.
!    is_usable(data_ind(1:nd)) = .true.
!
!
!    tree => kdtree2_create(nx,ny,nz,NODES,sort=.true.,rearrange=.true.)  ! this is how you create a tree. 
! 
!    
!!Build a map array so that the is_usable array can be updated
!!in the case where the kdtree rearranges the data
!    if(tree%rearrange)then
!       do i=1, nx*ny*nz
!           useablemap(tree%ind(i)) = i
!       enddo
!    endif
! 
!    sum_time=0
!
!! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID
!    write(*,*)
!    write(*,*) 'Preparing the simulation '
!
!
!!need to make an order for the simulation:
!    
!    call init_genrand(seedd)  !intial the random generator
!    do i=1,10000
!        rand_num = grnd()
!    end do
!    
!    allocate(temp_order(nx*ny*nz,2))
!    !use sim array as a temp array for the order
!    do i=1,nx*ny*nz
!        temp_order(i,1) = grnd()
!        temp_order(i,2) = dble(i)
!    end do
!
!    !call sortem(1,nx*ny*nz,temp_order(:,1),1,temp_order(:,2))   !now have the simulation order
!    call sortem2dp(1,nx*ny*nz,temp_order(:,1),1,temp_order(:,2))   !now have the simulation order
!
!    allocate(order(nx*ny*nz))
!    order = int(temp_order(:,2))
!    
!    deallocate(temp_order)
!
!
!    allocate(sim(nx*ny*nz,nreal),stat = test)
!    if(test.ne.0)then
!    
!          write(*,*)'ERROR: Allocation of sim array failed due to', &
!                ' insufficient memory.'
!          do 
!            nreal = nreal-1
!            if(allocated(sim)) deallocate(sim)
!            allocate(sim(nx*ny*nz,nreal),stat = test)
!            if(test == 0 .or. nreal==0) then
!                exit
!            end if           
!          end do
!          nreal=nreal-1 !leave some memory for onther opperations
!          deallocate(sim)
!          allocate(sim(nx*ny*nz,nreal),stat = test)
!          if(nreal <= 0 ) then
!              write(*,'(a)' ) 'ERROR: your grid is too large, cannot allocate for even a single realization'
!              write(*,'(a)' ) 'STOPPING'
!              stop
!          end if
!          write(*,'(a,I5,a)') 'ERROR: can only simulate ', nreal, ' realizations'
!    end if
!    
!    allocate(est(nreal) )
!    sim=-999  !reset the sim array
!    
!    do ireal=1,nreal
!       sim(data_ind(1:nd),ireal) = vr(1:nd)
!    end do
!
!    use_kd_tree=.true.  ! if =.false. will start with an exhaustive search, and used kd when it is more efficient
!    !use_kd_tree=.false.  ! if =.false. will start with an exhaustive search, and used kd when it is more efficient
!    !can set this =.true. and will always use the kd tree.
!    
!    allocate(ex_results_d(nclose))
!    allocate(ex_results(nclose))
!    
!    allocate(ex_sim_array(nloop+nd))
!    n_searched=0    
!    ex_sim_array=-999
!    ex_sim_array(1:nd) = data_ind(1:nd)
!    
!!    sum_time = secnds(start_time)
!!Write(*,*) 'setting arrays and variables for simulation',sum_time
!    elapsed=secnds(total)
!    write(*,'(f10.4,a)')      elapsed, 's - time to setting arrays and variables for simulation'
!    
!    write(*,*)
!    write(*,*) 'Working on the simulation '
!    
!!    start_time=secnds(0.0)
!    do index2=1,nloop
!!    do index2=1,10001
!
!      if((int(index2/irepo)*irepo).eq.index2) then
!         write(*,103) index2
! 103  format('   currently on estimate ',i9)
!      end if      
!!
!! Where are we making an estimate?
!!
!      if(koption.eq.0) then
!
!            !write(*,*)'koption.eq.0'
!
!            index = order(index2)
!            iz   = int((index-1)/nxy) + 1
!            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
!            ix   = index - (iz-1)*nxy - (iy-1)*nx
!            xloc = xmn + real(ix-1)*xsiz
!            yloc = ymn + real(iy-1)*ysiz
!            zloc = zmn + real(iz-1)*zsiz
!            !write(*,*)index2,index,ix,iy,iz,xloc,yloc,zloc
!      else
!        !make cross validation
!            !write(*,*)'koption.ne.0'
!            read(ljack,*,err=96,end=2) (var(i),i=1,nvarij)
!            ddh  = 0.0
!            xloc = xmn
!            yloc = ymn
!            zloc = zmn
!            true = UNEST
!            secj = UNEST
!            if(idhlj.gt.0)  ddh    = var(idhlj)
!            if(ixlj.gt.0)   xloc   = var(ixlj)
!            if(iylj.gt.0)   yloc   = var(iylj)
!            if(izlj.gt.0)   zloc   = var(izlj)
!            if(ivrlj.gt.0)  true   = var(ivrlj)
!            if(iextvj.gt.0) extest = var(iextvj)
!            if(true.lt.tmin.or.true.ge.tmax) true = UNEST
!            !index in now ...
!            call getindx(nx,xmn,xsiz,xloc,ix,inflag)
!            call getindx(ny,ymn,ysiz,yloc,iy,inflag)
!            call getindx(nz,zmn,zsiz,zloc,iz,inflag)
!
!            index = ix + (iy-1)*nx + (iz-1) * nx*ny 
!      end if
!      
!      !get the close data in an array 'results'
!      nclose = min(nd,ndmax)
!     
!      !the tree is only efficient if it is fairly full, do 1000 with an exhaustive search and then test tree,
!      !once the tree is efficient only use the tree
!      if(n_searched == 1000 .and. not(use_kd_tree)) then
!              
!!          n_searched=0
!!          time_temp=secnds(real(0.0,4))
!!          do i=1,10
!!          call kdtree2_n_nearest_around_point(real(coord_ISOMAP(index,1:d_tree)),d_tree,tp=tree,idxin=-1,correltime=-1,nn=nclose, results=results)
!!          end do
!!          kd_time=secnds(time_temp)
!!          
!!          
!!          !do exhaustive search to test for time
!!          time_temp=secnds(real(0.0,4))
!!          do i=1,10
!!          call exhaustive_srch(real(coord_ISOMAP(index,1:d_tree)),d_tree,nclose, ex_results, ex_results_d,nd+index2)
!!          end do
!!          exhaustive_srch_time=secnds(time_temp)
!!          
!!          if(kd_time<exhaustive_srch_time) use_kd_tree = .true.  !from now on use the kd tree
!         
!      else
!          n_searched=n_searched+1
!      end if
!     
!     
!      !do the search, for the first few locations it is faster to do an exhaustive search
!      if(use_kd_tree) then
!         call kdtree2_n_nearest_around_point(real(coord_ISOMAP(index,1:d_tree)),d_tree,tp=tree,idxin=-1,correltime=-1,nn=nclose, results=results)
!         !write(*,*) 'calculating neighbours (kdtree) for node ',index 
!      else
!         call exhaustive_srch(real(coord_ISOMAP(index,1:d_tree)),d_tree,nclose, ex_results,ex_results_d,nd+index2)
!         
!         !fill the results array:
!
!         
!         results(:).dis = ex_results_d(1:nclose)
!         results(:).idx = ex_results(1:nclose)
!         ex_sim_array(nd+index2) = index
!         
!         !write(*,*) 'calculating neighbours (exhaustive) for node ',index 
!      end if
!      tree%REARRANGED_IS_USABLE(useablemap(index)) = .true.  !need to tell the tree that this point is now simulated:
!
!      !do we need to trim out some of the points found?
!      do i=1,min(nd,ndmax)
!         if(results(i).dis>radsqd) then
!            nclose=nclose-1
!         end if
!      end do
!
!!      !how many eqns for kriging?
!!      neq=nclose
!!      if(ktype==1) neq=neq+1
!!      !
!!      ! Initialize the main kriging matrix:
!!      a = 0.0
!!      !
!!      ! Fill in the LHS kriging matrix:
!!      !
!!
!!      !need to calculate the covariances:
!!      do i=1,nclose
!!         ind1=results(i).idx
!!         do j=i,nclose
!!            ind2=results(j).idx
!!            dist=sqrt ( sum( ( coord_ISOMAP(ind1,1:dim)-coord_ISOMAP(ind2,1:dim) ) **2    )    )
!!            call cova3_1D(dist,1,nst,MAXNST,c0,it,cc,aa,cmax,a(neq*(i-1)+j)) 
!!            a(neq*(j-1)+i) = a(neq*(i-1)+j)
!!         end do
!!      end do
!!      !
!!      ! Fill in the OK unbiasedness portion of the matrix (if not doing SK):
!!      !
!!      if(ktype==1) then
!!         do i=1,nclose
!!            a(neq*(i-1)+nclose+1) = dble(cmax)
!!            a(neq*nclose+i)       = dble(cmax)
!!         end do
!!      endif
!!      !
!!      ! Fill in the RHS kriging matrix (r):
!!      !
!!      r=cmax !so that if we are doing OK the last term is set to 1 (sum of weights=1)
!!      do i=1,nclose
!!          !cal the distance first then the covariance
!!          ind1 = results(i).idx
!!          vra(i,:)=sim(ind1,:) 
!!          if(d_tree/=dim) then
!!              dist=results(i).dis + sum( ( coord_ISOMAP(ind1,d_tree+1:dim)-coord_ISOMAP(index,d_tree+1:dim) ) **2    ) !dist in tree is only for the first few dims
!!          else
!!              dist=results(i).dis 
!!          end if
!!          dist=sqrt(dist)
!!          call cova3_1D(dist,1,nst,MAXNST,c0,it,cc,aa,cmax,r(i)) 
!!      end do
!!      !
!!      ! Copy the right hand side to compute the kriging variance later:
!!      !
!!      do k=1,neq
!!            rr(k) = r(k)
!!      end do
!!      kadim = neq * neq
!!      ksdim = neq
!!      nrhs  = 1
!!      nv    = 1
!!!
!!! Write out the kriging Matrix if Seriously Debugging:
!!!
!!      if(idbg.eq.3) then
!!            write(ldbg,*) 'Estimating node index : ',ix,iy,iz
!!            is = 1 - neq
!!            do i=1,neq
!!                  is = 1 + (i-1)*neq
!!                  ie = is + neq - 1
!!                  write(ldbg,100) i,r(i),(a(j),j=is,ie)
!! 100              format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
!!            end do
!!      endif
!!      
!!      estv=-999
!!      !
!!      ! Solve the kriging system:
!!      !
!!      if(nclose>0 ) call ktsol(neq,nrhs,nv,a,r,s,ising,MAXEQ)
!!
!!      ! Compute the solution:
!!      if(nclose == 1 ) then
!!         ising=0
!!         s(1) = r(1) ;
!!      end if
!!      if(ising.ne.0) then
!!            if(idbg.ge.3) write(ldbg,*) ' Singular Matrix ',ix,iy,iz
!!            est  = UNEST
!!            estv = UNEST
!!      else if(nclose<ndmin) then
!!            if(idbg.ge.3) write(ldbg,*) ' not enough data to estimate '
!!            est  = UNEST
!!            estv = UNEST
!!      else
!!            est  = 0.0
!!            estv = real(cmax)
!!            do j=1,neq   
!!            
!!                      estv = estv - real(s(j))*rr(j)
!!                      if(j.le.nclose) then
!!                            if(ktype.eq.0) then
!!                                  est(:) = est(:) + real(s(j))*(vra(j,:)-skmean)
!!                            else
!!                                  est(:) = est(:) + real(s(j))*vra(j,:)
!!                            endif
!!                      endif
!!            end do
!!            if(ktype.eq.0.or.ktype.eq.2) est = est + skmean
!!            nk   = nk + 1
!!
!!            !
!!            ! Write the kriging weights and data if debugging level is above 2:
!!            !
!!            if(idbg.ge.2) then
!!                  write(ldbg,*) '       '
!!                  write(ldbg,*) 'BLOCK: ',ix,iy,iz,' at ',xloc,yloc,zloc
!!                  write(ldbg,*) '       '
!!                  if(ktype.ne.0)  &
!!                  write(ldbg,*) '  Lagrange : ',s(na+1)*unbias
!!                  write(ldbg,*) '  BLOCK EST: x,y,z,vr,wt '
!!                  do i=1,na
!!                        xa(i) = xa(i) + xloc - 0.5*xsiz
!!                        ya(i) = ya(i) + yloc - 0.5*ysiz
!!                        za(i) = za(i) + zloc - 0.5*zsiz
!!                        write(ldbg,'(50f12.3)') xa(i),ya(i),za(i), &
!!                                              vra(i,:),s(i)
!!                  end do
!!                  write(ldbg,*) '  estimate, variance  ',est,estv
!!            endif
!!      endif
!
!!
!! END OF MAIN KRIGING LOOP:
!!
! 1    continue
! 
! !need to draw from the gaussian distribution with this est and var:  (est,estv)
!
!!      if(nclose<ndmin .or. nclose == 0) then
!!         est = 0.0  !the global mean
!!         estv = 1.0
!!      end if
!!      
!!      !if (abs(estv) <= 1e-10 ) estv = 0.0
!!      if (abs(estv) <= 1e-7 ) estv = 0.0
!!
!!      estv=sqrt(estv)     
!!               
!!      do ireal=1,nreal
!!         pp = grnd()
!!         call gauinv(pp,xpp,ierr)
!!         sim(index,ireal) = xpp * estv + est(ireal)
!!      end do
!
!      end do
! 2    continue
!      if(koption.gt.0) close(ljack)
!!    sum_time = secnds(start_time)
!!Write(*,*) 'simulating points',sum_time
!    elapsed=secnds(total)
!    write(*,'(f10.4,a)')      elapsed, 's - time of simulation'
!!
!! Write statistics of kriged values:
!!
! 
!      if(nk.gt.0.and.idbg.gt.0) then
!            xk    = xk/real(nk)
!            vk    = vk/real(nk) - xk*xk
!            xkmae = xkmae/real(nk)
!            xkmse = xkmse/real(nk)
!            write(ldbg,105) nk,xk,vk
!            write(*,   105) nk,xk,vk
! 105        format(/,'Estimated   ',i8,' blocks ',/, &
!                    '  average   ',g14.8,/,'  variance  ',g14.8,/)
!            if(koption.ne.0) then
!                  write(*,106) xkmae,xkmse
! 106              format(/,'  mean error',g14.8,/,'  mean sqd e',g14.8)
!            end if
!      endif
!
!!
!! All finished the kriging:
!!
!      do i=1,nreal
!         write(lout,'(f10.4)') sim(:,i)
!      end do
!
!      deallocate(order)
!      deallocate(coord_ISOMAP)
!    elapsed=secnds(total)
!    write(*,'(f10.4,a)')      elapsed, 's - time of write out simulations'

      return
 96   stop 'ERROR in jackknife file!'
      end subroutine sgs



    subroutine fix_angles ( angles )
        parameter(DEG2RAD=3.141592654/180.0,EPSLON=1.e-20)
    
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
    

subroutine makepar
!----------------------------------------------------------------------
!                       WRITE A PARAMETER FILE                         
!----------------------------------------------------------------------
    integer, parameter :: lun=99
    open(lun,file='sgs_lva.par',status='UNKNOWN')
    write(lun,10)
10  format(  &
    '',/, &
    '                  Parameters for sgs_lva',/, &
    '                  *************************',/, &
    '',/, &
    'START OF PARAMETERS:',/, &
    'nsdata.out         -file with data',/, &
    '1 2 0 3            -columns for X,Y,Z,var,sec var',/, &
    '-998   1.0e21      -trimming limits',/, &
    '0                  -debugging level: 0,1,2,3',/, &
    'sgs.dbg            -file for debugging output',/, &
    'sgs.out            -file for output',/, &
    '10                 -nreal',/, &
    '100 0.5 1          -nx,xmn,xsiz (ESTIMATION GRID see NOTE1)',/, &
    '100 0.5 1          -ny,ymn,ysiz (ESTIMATION GRID see NOTE1)',/, &
    '1   0.5 1          -nz,zmn,zsiz (ESTIMATION GRID see NOTE1)',/, &
    '32146              -random number seed',/, &
    'grid_LVA.out       -file containing the locally varying anisotropy grid (LVAG)',/, &
    '1 2 3 4 5          -LVA columns for ang1, ang2, ang3, aniso ratio min/max, aniso ratio vert/max',/, &
    '100 0.5 1          -nx,xmn,xsiz (LVA GRID see NOTE1)',/, &
    '100 0.5 1          -ny,ymn,ysiz (LVA GRID see NOTE1)',/, &
    '1   0.5 1          -nz,zmn,zsiz (LVA GRID see NOTE1)',/, &
    '1                  -noffsets for graph (number of offsets described below, see NOTE2)',/, &
    '2                  -use MDS? 2=L-ISOMAP 3=read dist from  grid_cpp.out  and use L-ISOMAP (see NOTE3)',/, &
    '10 10 1            -number of landmark points in x,y,z (evenly spaced grid nodes)',/, &
    '-1                 -max number of dimensions to use (set -1 to use max, see NOTE4)',/, &
    '2 30               -min, max nodes for simulation',/, &
    '1000               -maximum search radii (a 1D isotropic distance in q dimensions)',/, &
    '-1                 -maximum number of dimensions to use in search (-1 uses all dimensions see NOTE5)',/, &
    '0     0            -0=SK,1=OK, simple kriging mean',/, &
    '1 0.00             -# of nested structures, nugget effect (1D variogram)',/, &
    '2 1 35             -it,cc,range (see NOTE6)',/, &
    '0                  -use multi-thread Boost Dijkstra (0=no, 1=yes)',/, &
    '',/, &
    '',/, &
    '',/, &
    'NOTE1: The estimation grid and the LVA grid can be identical.',/, &
    '       If the grids are different the LVA GRID is resampled to the ESTIMATION GRID.',/, &
    '       The estimation grid affects how the multidimensional scaling occurs',/, &
    '       therefore, if the estimation grid is different it needs to be accounted for',/, &
    '       in the variogram calculation.',/, &
    '',/, &
    'NOTE2: This program calls on the Boost_dijkstra.exe program to calculate the shortest paths.',/, &
    '       When calculating the shortest paths there is the option to connect nodes separated by more than one grid block.',/, &
    '       This allows for more flexible paths.  Keeping offset=1 connects each node to the nearest 8 (in 2D) and only',/, &
    '       allows paths with 0,45,90,135 degree increments.  Using more offsets is recommended as paths will be smoother and',/, &
    '       shorter but this requires more memory and may be infeasible for large 3D models',/, &
    '',/, &
    'NOTE3: The Dijkstra program can be CPU intensive.  If the program is run once a file  grid_cpp.out  is',/, &
    '       created with all the distances between nodes and landmark points required.  Subsequent runs of the',/, &
    '       program can simply read in this file rather than call the Dijkstra algorithm again.',/, &
    '       If the number of offsets, the grid or the configuration of landmark points are changed the Dijkstra program',/, &
    '       MUST be run again to calculate the new UPDATED distances.',/, &
    '',/, &
    'NOTE4: Retaining fewer dimensions than the max will reduce memory requirements a bit and will speed up the program.',/, &
    '       The most important dimensions are retained.',/, &
    '',/, &
    'NOTE5: Using fewer than the maximum dimensions can speed up the program.  The most important dimensions are retained.',/, &
    '       This may reduce the memory requirements of the program (i.e. the kdtree) if using large grids and many landmark points.',/, &
    '       Set to -1 unless CPU time is an issue.  The number of dimensions searched to determine the nearest data',/, &
    '       can be different than the number of dimensions used in the calculation of the distance between locations.',/, &
    '       Using fewer data in the search can speed up the program, but may result in some error when finding close data.',/, &
    '',/, &
    'NOTE6: Variogram must be positive definate in q dimensions.  The exponential variogram is recomended:',/, &
    '       1-spherical variogram',/, &
    '       2-exponenital variogram',/, &
    '       3-gaussian variogram',/, &
    '       4-power model',/, &
    '       5-hole effect model',/, &
    '',/)

    close(lun)
end subroutine makepar

