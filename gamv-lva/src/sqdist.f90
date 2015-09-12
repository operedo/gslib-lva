      !real(kind=dp) function sqdist(index1,index2,opt)
      subroutine sqdist(index1,index2,opt,h)
      use kind
!      use LVA_dist
      use global

!-----------------------------------------------------------------------
!
!    Squared Anisotropic Distance Calculation Given Matrix Indicator
!    ***************************************************************

!this program modified by Jeff Boisvert April 2007.
!
!
!
!program now takes x,y,z coordinates and an anisotropy field.  Optimization
!of the path ia performed and the distance is returned.  If optimization is wanted (option=0)
!then the true distance is returned not the squared distance.  If the optimization is ignored (option>0)
!then the squared distance is returned.  This is because the superbloc search is set up for squared
!distance, but the optimization (from john manchuk) returns the true distance.
!
! INPUT VARIABLES:
!
!   x1,y1,z1         Coordinates of first point
!   x2,y2,z2         Coordinates of second point
!   variables 'maxrot, ind and rotmat' are not used
!   pt_ID           contains the point id for the two points.  if it is a -1 then it is not a data point but a grid location
!                   if it is a data to data distance it will be stored in the hash table
!   opt             if opt=0 then optimize the distance and return the TRUE distance (not the squared distacnce)
!                   if opt=1 then return the straight line squared distacnce, this is required for the super block search
!
! OUTPUT VARIABLES:
!
!   sqdist           The squared distance accounting for the locally varying anisotropy
!
!
!  EXTERNAL REFERENCES
!  uses the grid defintion as read in kt3d_aniso.f90
!
!-----------------------------------------------------------------------
      !use global
      !real*8 cont,dx,dy,dz,p1(3),p2(3),sqdist2,temp_D
      real(kind=dp) cont,dx,dy,dz,p1(3),p2(3),sqdist2,temp_D,sumtmp,h
      integer opt, position(2),i
      integer :: pt_ID(2)
      integer key(2)
      

      !sqdist= dsqrt(   dble(  sum( ( coord_ISOMAP(index1,1:dim)- &
      !coord_ISOMAP(index2,1:dim) ) **2    )    )  )!dist in tree is only for the first few dims 

      sumtmp=0.0
      do i=1,dim
         sumtmp=sumtmp + ( real(coord_ISOMAP(index1,i))-&
                real(coord_ISOMAP(index2,i)))**2 
      end do
      !write(*,*)'sumtmp=',sumtmp 

      !sqdist= sqrt(   real(  sum( ( coord_ISOMAP(index1,1:dim)- &
      !coord_ISOMAP(index2,1:dim) ) **2    )    )  )!dist in tree is only for the first few dims 

      !sqdist = sqrt(real(sumtmp))
      h = sqrt(real(sumtmp))
      
      !write(*,*)'sqrt(real(sumtmp))=',sqdist 

      end
