module global
!just contains some global variables, mostly for trouble shooting
    implicit none
    
    real*8, public :: F_squared
    real*8, allocatable, dimension (:,:,:,:,:,:) :: template_D !template to speed up dist cals
    real*8, allocatable, dimension (:,:,:) :: norm_D !template of lengths of offsets
    real*8, allocatable, dimension (:,:) :: data_cov,data_cross ! store the distance between all points
    real*8, allocatable, dimension (:,:) :: a_JM,d2lanmark,vec,vec1,val,subB2,subB22,vec_test
    real*8, allocatable, dimension (:,:) :: coord_ISOMAP,temp_coord,coord_LLE,anticline_ex_dij_orij
    integer, allocatable, dimension (:,:) :: coord_ISOMAP_ind
    real*8, allocatable, dimension (:,:) :: subB
    real*8, allocatable :: evalues(:), vectors(:,:),x_loc(:),y_loc(:),z_loc(:),temp_coord_ISOMAP(:),dist_4_cross(:)
    integer xland,yland,zland,xspace,yspace,zspace,add_paths_data,add_paths_land
    
    
    real*8, allocatable, dimension (:) :: r_JM, r_dist,dist_temp
    real*8, allocatable, dimension (:) :: x_JM,r_tmp,cross_ind
    integer, allocatable, dimension (:) :: bad_data,IPIV,data_ind
    integer, allocatable :: landpts(:)
    real*8 :: bbb = 0.5 !for John Manchuk's robust solver    
    real*8, pointer, dimension (:) :: hash_dist ! hash table to store the distances
    integer, pointer, dimension (:,:) :: keys ! hash table to store the keys
    integer ndata,cal_stress,cur_node
    real*8 post_MDS_dist,sum_dist
    integer :: nbad_matrix !count the number of bad matricies
    integer :: ncol !count the number of colisions for the hash tabel
    integer :: n_seg ! number of times optimization segments the line
    integer :: max_seg ! if greater than -1, number of times to segment the line
    integer :: n_hash_table !count number of entries in the hash table
    integer T_size ! the size of the hash table
    integer max_hash !most number of entries in the hash table
    real*8  ::  min_dist_tol,stress !for timers
    real*4 :: elapsed, total, start_time,sum_time
    integer, parameter :: time_out= 11,time_out2= 13, T_eigs=145, T_matrix=146, T_dijk=147,T_search=148
    logical err_beta,inv_dist
    character*40 error_message
    real*8 min_eigen
    integer n_indef,na_total,na_cnt,good,graph_opt,ldistout,all_dist
    integer  indef_mat,c1,c2,c3,assign_nodes,xind,yind,zind
    integer xyzland
    real*8 :: DIAG1,dist,pp1(20),pp2(20)
    integer dim,d_tree
    real*8 edge_dist
    integer edge_NODE
    real*8 power,rad_inv_dis
    integer nd,DHID_col,cross_cnt

end module global