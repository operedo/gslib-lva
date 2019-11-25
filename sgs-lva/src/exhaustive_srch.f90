module srch
  use global
!  implicit none

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
    
    real               :: vecin_tmp(dim),ds


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
      
      !vecin_tmp =coord_ISOMAP(ind,1:dim)-vecin 
      !ds = SNRM2(dim,vecin_tmp,1)
      !d = dble(ds)
      !d = d*d
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

subroutine exhaustive_srch_trans(vecin,dim,nn,results,results_d,n2srch)
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
    
    real               :: vecin_tmp(dim),ds


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
      
      !vecin_tmp =coord_ISOMAP(ind,1:dim)-vecin 
      !ds = SNRM2(dim,vecin_tmp,1)
      !d = dble(ds)
      !d = d*d
      !d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
      d= dble( sum( ( coord_ISOMAP_trans(1:dim,ind)-vecin ) **2    ))
      
      if(d<largest_d) then  !if d==0.0 then we are at a simualtion location and skip
      
          !!do we have this location already?
          !skip=.false.
          !do k=1,nn
          !    !if(results_temp(k) ==real(ind)) then 
          !    if(results_temp(k) ==dble(ind)) then
          !        write(*,*)'ERROR' 
          !        skip = .true.
          !        exit
          !    end if
          !end do
      
          !if(not(skip) ) then 
              results_d(largest) = d
              !results_temp(largest) = real(ind)
              results_temp(largest) = dble(ind)
              
              !get the new largest value
              tyty = MAXLOC(results_d)
              largest=tyty(1)
              largest_d=results_d(largest)
          !end if
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


end subroutine exhaustive_srch_trans

subroutine exhaustive_srch_trans_local(vecin,dim,nn,results,results_d,n2srch,nndd,tid,ex_sim_array_local,nsize)
    ! Find the 'nn' vectors nearest to  'vecin',
    ! returing results in results(:)
    !n2srch has the number of previously simulated nodes to search in the global array ex_sim_array
    
    !real*8                           :: vecin(dim),d,largest_d,results_d(nn),results_temp(nn)
    integer, intent(in)            :: dim,n2srch,nndd,tid,nsize
    real, intent(in)               :: vecin(dim)
    integer, intent(inout)            :: nn,results(nn),ex_sim_array_local(nsize)
    real*8, intent(inout)            :: results_d(nn)
    integer                        :: idxin, correltime, i,j,k,ind,tyty(1),largest
    real*8         :: results_temp(nn),d,largest_d 
    logical        ::          skip
    
    real               :: vecin_tmp(dim),ds

   !if(tid.eq.2)write(*,*)'[AFTER] ex_sim_array_local(',n2srch-1,')=',ex_sim_array_local(n2srch-1)
   !do an exhaustive search
   
   results_d=1.0e21
   
   !get the largest dist:
   largest   = 1
   largest_d = 1.0e21
   results_temp=-999
   
   do i=1,n2srch

      !get the distance to this node
      ind = ex_sim_array_local(i)
      !if(i.gt.nndd)write(*,*)'-->',i,ind
      if(ind == -999) exit !(we are done now)   
      
      !vecin_tmp =coord_ISOMAP(ind,1:dim)-vecin 
      !ds = SNRM2(dim,vecin_tmp,1)
      !d = dble(ds)
      !d = d*d
      !d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
      d= dble( sum( ( coord_ISOMAP_trans(1:dim,ind)-vecin ) **2    ))
      
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


end subroutine exhaustive_srch_trans_local


subroutine exhaustive_srch_trans_opt01(vecin,dim,nn,results,results_d,n2srch)
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
    
    real               :: vecin_tmp(dim),ds


   !do an exhaustive search
   
   results_d=1.0e21
   
   !get the largest dist:
   largest   = 1
   largest_d = 1.0e21
   results_temp=-999
   
   do i=1,n2srch-1,2

      !get the distance to this node
      ind = ex_sim_array(i)
      if(ind .ne. -999) then !(we are done now)   
      
      !vecin_tmp =coord_ISOMAP(ind,1:dim)-vecin 
      !ds = SNRM2(dim,vecin_tmp,1)
      !d = dble(ds)
      !d = d*d
      !d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
      d= dble( sum( ( coord_ISOMAP_trans(1:dim,ind)-vecin ) **2    ))
      
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

      end if

      !get the distance to this node
      ind = ex_sim_array(i+1)
      if(ind == -999) exit !(we are done now)   
      
      !vecin_tmp =coord_ISOMAP(ind,1:dim)-vecin 
      !ds = SNRM2(dim,vecin_tmp,1)
      !d = dble(ds)
      !d = d*d
      !d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
      d= dble( sum( ( coord_ISOMAP_trans(1:dim,ind)-vecin ) **2    ))
      
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

   do j=i,n2srch

      !get the distance to this node
      ind = ex_sim_array(j)
      if(ind == -999) exit !(we are done now)   
      
      !vecin_tmp =coord_ISOMAP(ind,1:dim)-vecin 
      !ds = SNRM2(dim,vecin_tmp,1)
      !d = dble(ds)
      !d = d*d
      !d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
      d= dble( sum( ( coord_ISOMAP_trans(1:dim,ind)-vecin ) **2    ))
      
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


end subroutine exhaustive_srch_trans_opt01




subroutine exhaustive_srch_opt01           (vecin,dim,nn,results,results_d,n2srch)
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
    
    real               :: vecin_tmp(dim),ds


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
      
      !vecin_tmp =coord_ISOMAP_trans(1:dim,ind)-vecin 
      !ds = SNRM2(dim,vecin_tmp,1)
      !d = dble(ds)
      !d = d*d
      !d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
      d= dble( sum( ( coord_ISOMAP_trans(1:dim,ind)-vecin ) **2    ))
      
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


end subroutine exhaustive_srch_opt01


!subroutine exhaustive_srch_omp(vecin,dim,nn,results,results_d,n2srch)
!    ! Find the 'nn' vectors nearest to  'vecin',
!    ! returing results in results(:)
!    !n2srch has the number of previously simulated nodes to search in the global array ex_sim_array
!    
!    !real*8                           :: vecin(dim),d,largest_d,results_d(nn),results_temp(nn)
!    integer, intent(in)            :: dim,n2srch
!    real, intent(in)               :: vecin(dim)
!    integer, intent(inout)            :: nn,results(nn)
!    real*8, intent(inout)            :: results_d(nn)
!    integer                        :: idxin, correltime, i,j,k,ind,tyty(1),largest
!    real*8         :: results_temp(nn),d,largest_d 
!    logical        ::          skip
!    
!    real               :: vecin_tmp(dim),ds
!
!    integer, allocatable :: results_omp(:)
!    real*8, allocatable :: results_temp_omp(:),results_d_omp(:)
!    integer :: thread_id,num_threads
!
!    INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
!
!   !do an exhaustive search
!   
!   !results_d=1.0e21
!   results_d_omp=1.0e21
!   
!   !get the largest dist:
!   largest   = 1
!   largest_d = 1.0e21
!   !results_temp=-999
!   results_temp_omp=-999
!   
!!$omp parallel private(thread_id,num_threads)
!   thread_id = omp_get_thread_num()
!   num_threads = omp_get_num_threads()
!   write(*,*) 'INFO: ',thread_id,num_threads 
!!$omp end parallel
!   allocate(results_d_omp(nn*num_threads))
!   allocate(results_omp(nn*num_threads))
!   allocate(results_temp_omp(nn*num_threads))
!
!!$omp parallel default(shared) private(thread_id,i,ind,d,skip,k,results_d,results_tmp,tyty,largest,largest_d)
!   thread_id = omp_get_thread_num()
!   num_threads = omp_get_num_threads()
!!$omp do schedule(static)
!   do i=1,n2srch
!
!      !get the distance to this node
!      ind = ex_sim_array(i)
!      !if(ind == -999) exit !(we are done now)   
!      if(ind .ne. -999) then
!      
!      !vecin_tmp =coord_ISOMAP(ind,1:dim)-vecin 
!      !ds = SNRM2(dim,vecin_tmp,1)
!      !d = dble(ds)
!      !d = d*d
!      d= dble( sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    ))
!      
!      if(d<largest_d) then  !if d==0.0 then we are at a simualtion location and skip
!          !do we have this location already?
!          skip=.false.
!          do while(skip.eq..false.)
!          !do k=1,nn
!              write(*,*) 'thread:',thread_id,' ',i,' ',k 
!              !if(results_temp(k) ==real(ind)) then 
!              !if(results_temp(k) ==dble(ind)) then 
!              if(results_temp_omp(k+thread_id*nn) ==dble(ind)) then 
!                  skip = .true.
!                  !exit
!              end if
!          end do
!      
!          if(not(skip) ) then 
!              !results_d(largest) = d
!              results_d_omp(largest + thread_id*nn) = d
!              !results_temp(largest) = dble(ind)
!              results_temp_omp(largest + thread_id*nn) = dble(ind)
!!              
!!              !get the new largest value
!!              !tyty = MAXLOC(results_d)
!!              tyty = MAXLOC(results_d_omp((thread_id*nn+1):((thread_id+1)*nn)))
!!              largest=tyty(1)
!!              !largest_d=results_d(largest)
!!              largest_d=results_d_omp(largest + thread_id*nn)
!          end if
!      end if
!      end if
!   end do
!!$omp end do nowait
!!$omp end parallel
!
!   !call sortem(1,nn,results_d,1,results_temp) 
!   !call sortem2dp(1,nn,results_d,1,results_temp)
!   !call sortem2dp(1,(nn*num_threads),results_d_omp,1,results_temp_omp)
!
!
!   write(*,*) results_temp_omp
! 
!!   if(minval(results_temp) == -999) then !need to get new nclose
!!      do i=1,nn
!!          if(results_temp(i) == -999) then
!!              nn=i-1
!!              exit
!!          end if
!!      end do
!!   end if
!!
!!
!!   results=int(results_temp)
!   deallocate(results_d_omp)
!   deallocate(results_omp)
!
!
!end subroutine exhaustive_srch_omp
!

end module srch
