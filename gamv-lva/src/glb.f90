module global
    use kind

    implicit none
    
    real(kind=dp), public :: F_squared
    
    real(kind=dp), allocatable, dimension (:,:,:,:,:,:) :: template_D !template to speed up dist cals
    real(kind=dp), allocatable, dimension (:,:,:) :: norm_D !template of lengths of offsets
    real(kind=dp), allocatable, dimension (:) ::  r_dist
    integer, allocatable, dimension (:) :: IPIV
    integer, allocatable :: landpts(:)
    
    real(kind=dp), allocatable, dimension (:,:) :: data_dist ! store the distance between all points
    real(kind=dp), save, allocatable, dimension (:,:) :: a_JM,d2lanmark,vec,vec1,val,subB2,subB22,coord_ISOMAP,vec_test
    real(kind=dp), allocatable, dimension (:,:) :: subB
    real(kind=dp), allocatable, dimension (:) :: r_JM
    real(kind=dp), allocatable, dimension (:) :: x_JM,r_tmp
    integer, allocatable, dimension (:) :: bad_data
    real(kind=dp) :: bbb = 0.5 !for John Manchuk's robust solver    
    real(kind=dp), pointer, dimension (:) :: hash_dist ! hash table to store the distances
    integer, pointer, dimension (:,:) :: keys ! hash table to store the keys
    integer ndata,no_hash,which_solver
    integer :: nbad_matrix !count the number of bad matricies
    integer :: ncol !count the number of colisions for the hash tabel
    integer :: n_seg ! number of times optimization segments the line
    integer :: max_seg ! if greater than -1, number of times to segment the line
    integer :: n_hash_table !count number of entries in the hash table
    integer T_size ! the size of the hash table
    integer max_hash !most number of entries in the hash table
    real(kind=dp)      :: elapsed, total, min_dist_tol !for timers
    integer, parameter :: time_out= 11
    logical err_beta
    character*40 error_message
    real(kind=dp) min_eigen
    integer n_indef,na_total,na_cnt,xyzland,dim
    integer ng_nodes,grid_NODES,debug_node,cnt_lambda(100),MDS_opt,minE,ID1,ID2,xland,yland,zland

    integer nx,ny,nz,xind,yind,zind
    real(kind=dp) xmn,ymn,zmn,xsiz,ysiz,zsiz
    real(kind=dp),allocatable  :: evalues(:), vectors(:,:),bb(:,:),jsum(:),isum(:),aaa(:,:),atemp(:,:)
    
end module global
