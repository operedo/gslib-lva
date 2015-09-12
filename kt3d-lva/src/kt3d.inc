!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 2003, Statios Software and Services Incorporated.  All %
! rights reserved.                                                     %
!                                                                      %
! This program has been modified from the one distributed in 1996 (see %
! below).  This version is also distributed in the hope that it will   %
! be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
! code may be redistributed without restriction; however, this code is %
! for one developer only. Each developer or user of this source code   %
! must purchase a separate copy from Statios.                          %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
!                  Kriging of a 3-D Rectangular Grid
!                  *********************************
!
! This file controls both the static dimensioning and allocates most of
! the storage in common blocks.
!
!   MAXSBX    maximum super block nodes in X direction
!   MAXSBY    maximum super block nodes in Y direction
!   MAXSBZ    maximum super block nodes in Z direction
!   MAXDAT    maximum number of data points
!   MAXDT     maximum number of drift terms
!   MAXSAM    maximum number of data points to use in one kriging system
!   MAXDIS    maximum number of discretization points per block
!   MAXNST    maximum number of nested structures
!
!-----------------------------------------------------------------------
!
! User Adjustable Parameters:
!
    integer MAXDT,MAXNST,MAXROT
      parameter(MAXDT  =   9, MAXNST =   4)
      real*8 UNEST,EPSLON,VERSION
!
! Fixed Parameters:
!
      parameter(UNEST = -999.0, EPSLON = 0.000001, MAXROT = MAXNST + 1, &
               VERSION = 1.000)
!
! Static Array Dimensioning:
!
      integer idrif(MAXDT),nst(1),it(MAXNST),iktype,ncut,test
      integer nx,ny,nz,ndmin,noct,nxdis,nydis,nzdis,itrend,ktype
      integer idbg,koption,iout,iext,lext,iedtve,iextve,ljack,idhlj
      integer ixlj,iylj,izlj,ivrlj,iextvj,nvarij,MAXEQ
      integer na,ndb,bv,unbias,isrot
      real*8  c0(1),cc(MAXNST),aa(MAXNST),ang1(MAXNST), &
              ang2(MAXNST),ang3(MAXNST),anis1(MAXNST),anis2(MAXNST)
      real*8  rotmat(MAXROT,3,3),tmin,tmax, xmn,ymn,zmn,ysiz,xsiz,zsiz
      real*8  radius,skmean
      real*8  sang1,sang2,sang3,sanis1,sanis2
      
!
! The data and other input variables:
!
      common /datcom/ tmin,tmax,nx,ny,nz,xmn,ymn,zmn,test, &
                     xsiz,ysiz,zsiz,ndmax,ndmin,radius,noct,nxdis, &
                     nydis,nzdis,idrif,itrend,ktype,skmean,koption, &
                     idbg,ldbg,iout,lout,iext,lext,iextve,ljack, &
                     idhlj,ixlj,iylj,izlj,ivrlj,iextvj,nvarij,MAXEQ
!
! Variogram Parameters:
!
      common /vargdt/ nst,it,c0,cc,aa,ang1,ang2,ang3,anis1,anis2,rotmat
!
! Search variables and data for kriging:
!
      common /srccom/ sang1,sang2,sang3,sanis1,sanis2,isrot,na,ndb, &
                     bv,unbias,iktype,ncut
