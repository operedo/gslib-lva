module srch
  use global
  implicit none

contains 

!kdtree2_n_nearest_around_point(real   (coord_ISOMAP(index,1:d_tree)),d_tree,nn=nclose, results=results)
  subroutine exhaustive_srch           (vecin,dim,nn,results,results_d,n2srch)
    ! Find the 'nn' vectors nearest to  'vecin',
    ! returing results in results(:)
    !n2srch has the number of previously simulated nodes to search in the global array ex_sim_array
    
    !real*8                           :: vecin(dim),d,largest_d,results_d(nn),results_temp(nn)
    integer, intent(in)            :: dim,n2srch
    real, intent(in)               :: vecin(dim)
    integer, intent(inout)            :: nn,results(nn)
    real*8, intent(inout)            :: results_d(nn)
    integer                        :: idxin, correltime, i,j,k,ind,tyty(1),largest
    real*8         :: results_temp(nn),d,largest_d 
    logical        ::          skip

   !do an exhaustive search
   
   results_d=1.0e21
   
   !get the largest dist:
   largest   = 1
   largest_d = 1.0e21
   results_temp=-999
   
   do i=1,n2srch

      !get the distance to this node
      ind = ex_sim_array(i)
      if(ind == -999) exit !(we are done now)   
      

      d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
      
      if(d<largest_d) then  !if d==0.0 then we are at a simualtion location and skip
      
          !do we have this location already?
          skip=.false.
          do k=1,nn
              !if(results_temp(k) ==real(ind)) then 
              if(results_temp(k) ==dble(ind)) then 
                  skip = .true.
                  exit
              end if
          end do
      
          if(not(skip) ) then 
              results_d(largest) = d
              !results_temp(largest) = real(ind)
              results_temp(largest) = dble(ind)
              
              !get the new largest value
              tyty = MAXLOC(results_d)
              largest=tyty(1)
              largest_d=results_d(largest)
          end if
      end if
   end do

   !call sortem(1,nn,results_d,1,results_temp) 
   call sortem2dp(1,nn,results_d,1,results_temp)


 
   if(minval(results_temp) == -999) then !need to get new nclose
      do i=1,nn
          if(results_temp(i) == -999) then
              nn=i-1
              exit
          end if
      end do
   end if


   results=int(results_temp)


end subroutine exhaustive_srch


end module srch
