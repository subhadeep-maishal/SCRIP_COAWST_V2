!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module assigns the masking for each model and calls 
!     the SCRIP Package routines    
!
!-----------------------------------------------------------------------
!     This is a wrapper around SCRIP package developed at USGS, 
!     Woods Hole Science And Marine Center. This routine calls the 
!     SCRIP's package routine that is the driver for computing the 
!     addresses and weights for interpolating between two grids on 
!     a sphere.
!
!---- Written by John C. Warner-------------------------------------- 
!-----         Tarandeep S. Kalra -------------------------------------
!--------------Date: 08/20/2015----------------------------------------
!***********************************************************************

      module create_masks

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types 
      use scripwrap_mod ! defines input file variable names for SCRIP-COAWST wrapper
      use scrip         ! Module that calls the scrip package subroutin

      implicit none

!     Variables that are used in all the subroutines below
      character(char_len),    save :: grid1_file, grid2_file  
      character(char_len),    save :: interp_file1, interp_file2   
      character(char_len),    save :: map1_name, map2_name 
      character(char_len),    save :: mo_string, mw_string, ma_string
      integer(kind=int_kind), save :: i, j, c, N
      integer(kind=int_kind), save :: mo, mw, ma, nx, ny, grdsize
      integer(kind=int_kind), save :: ncstat, nc_file_id

      contains

      subroutine ocn2wav_mask()

      implicit none 

! Allocate the source mask 
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%src_mask(nx,ny))   
      end do 
! Allocate the destination mask 
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        allocate(ngrd_sw(mw)%dst_mask(nx,ny))
      end do 
!  Save the ocean grid mask to src_mask 
      do mo=1,Ngrids_roms
        do j=1,ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,j)=ngrd_rm(mo)%mask_rho_o(i,j)
          end do 
        end do 
      end do 
!  Save the wave grid mask to dst_mask 
      do mw=1,Ngrids_swan
        do j=1,ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%dst_mask(i,j)=ngrd_sw(mw)%mask_rho_w(i,j)
          end do
        end do 
      end do 
!  Determine destination masks for the swan grids using ocean as the source.
      if ((Ngrids_swan>0).and.(Ngrids_roms>0)) then 
        do mo=1,Ngrids_roms
          if (Ngrids_roms>mo) then ! mask out child portion of this ocn grid
            do j=ngrd_rm(mo)%jstr_o,ngrd_rm(mo)%jend_o
              do i=ngrd_rm(mo)%istr_o,ngrd_rm(mo)%iend_o
                ngrd_rm(mo)%src_mask(i,j)=0
              end do 
            end do 
          end if 
        end do 
      end if 
!  Send the data to SCRIP routines
      do mw=1,Ngrids_swan
        do mo=1,Ngrids_roms 
          write(mo_string,'(i1)')mo
          write(mw_string,'(i1)')mw
          interp_file1='ocn'//trim(mo_string)//'_to_'//'wav'//
     &                       trim(mw_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='ROMS to SWAN distwgt Mapping'
          map2_name='unknown' 
                 
          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", roms_grids(mo)
          write(stdout,*)"The dst grid is: ", swan_coord(mw)

          grid1_file=roms_grids(mo)
          grid2_file=swan_coord(mw)

          call scrip_package(grid1_file, grid2_file, 
     &                             ngrd_rm(mo)%xi_size,
     &                             ngrd_rm(mo)%eta_size,
     &                             ngrd_rm(mo)%lon_rho_o, 
     &                             ngrd_rm(mo)%lat_rho_o, 
     &                             ngrd_rm(mo)%x_full_grid, 
     &                             ngrd_rm(mo)%y_full_grid, 
     &                             ngrd_sw(mw)%Numx_swan,
     &                             ngrd_sw(mw)%Numy_swan,
     &                             ngrd_sw(mw)%lon_rho_w,
     &                             ngrd_sw(mw)%lat_rho_w,
     &                             ngrd_sw(mw)%x_full_grid, 
     &                             ngrd_sw(mw)%y_full_grid, 
     &                              interp_file1, interp_file2, 
     &                                    map1_name, map2_name, 
     &                                    ngrd_rm(mo)%src_mask,
     &                                    ngrd_sw(mw)%dst_mask)
        end do 
      end do
! DOUBLE CHECK THESE LINES ARE OKAY 
       do mo=1,Ngrids_roms
          deallocate(ngrd_rm(mo)%src_mask)
         end do
         do mw=1,Ngrids_swan
          deallocate(ngrd_sw(mw)%dst_mask)
       end do
      end subroutine ocn2wav_mask

!======================================================================

      subroutine wav2ocn_mask()
      implicit none 

!     local variables 
      character(char_len) :: grid1_file, grid2_file  
      character(char_len) :: interp_file1, interp_file2   
      character(char_len) :: map1_name, map2_name 
      character(char_len) :: mo_string, mw_string
      integer(int_kind)   :: i, j, c, N
      integer(int_kind)   :: mo, mw, nx, ny, grdsize
      integer(int_kind)   :: ncstat, nc_file_id

! Allocate the source mask 
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        allocate(ngrd_sw(mw)%src_mask(nx,ny))
      end do 
! Allocate the destination mask 
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%dst_mask(nx,ny))   
      end do 
!  Save the wave grid mask to dst_mask 
      do mw=1,Ngrids_swan
        do j=1,ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%src_mask(i,j)=ngrd_sw(mw)%mask_rho_w(i,j)
          end do
        end do 
      end do 
!  Save the ocean grid mask to src_mask 
      do mo=1,Ngrids_roms
        do j=1,ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%dst_mask(i,j)=ngrd_rm(mo)%mask_rho_o(i,j)
          end do 
        end do 
      end do 
!  Determine destination masks for the ROMS grids using SWAN as the source.
      if ((Ngrids_swan>0).and.(Ngrids_roms>0)) then 
        do mw=1,Ngrids_swan
          if (Ngrids_swan>mw) then ! mask out child portion of this ocn grid
            do j=ngrd_sw(mw)%jstr_w,ngrd_sw(mw)%jend_w
              do i=ngrd_sw(mw)%istr_w,ngrd_sw(mw)%iend_w
                ngrd_sw(mw)%src_mask(i,j)=0
              end do 
            end do 
          end if 
        end do 
      end if 
!  Send the data to SCRIP routines
      do mo=1,Ngrids_roms
        do mw=1,Ngrids_swan 
          write(mw_string,'(i1)')mw
          write(mo_string,'(i1)')mo
          interp_file1='wav'//trim(mw_string)//'_to_'//'ocn'//
     &                       trim(mo_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='SWAN to ROMS distwgt Mapping'
          map2_name='unknown' 
                 
          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", swan_coord(mw)
          write(stdout,*)"The dst grid is: ", roms_grids(mo)

          grid1_file=swan_coord(mw)
          grid2_file=roms_grids(mo)

          call scrip_package(grid1_file, grid2_file, 
     &                             ngrd_sw(mw)%Numx_swan,
     &                             ngrd_sw(mw)%Numy_swan,
     &                             ngrd_sw(mw)%lon_rho_w,
     &                             ngrd_sw(mw)%lat_rho_w,
     &                             ngrd_sw(mw)%x_full_grid, 
     &                             ngrd_sw(mw)%y_full_grid, 
     &                             ngrd_rm(mo)%xi_size,
     &                             ngrd_rm(mo)%eta_size,
     &                             ngrd_rm(mo)%lon_rho_o, 
     &                             ngrd_rm(mo)%lat_rho_o, 
     &                             ngrd_rm(mo)%x_full_grid, 
     &                             ngrd_rm(mo)%y_full_grid, 
     &                              interp_file1, interp_file2, 
     &                                    map1_name, map2_name, 
     &                                    ngrd_sw(mw)%src_mask,
     &                                    ngrd_rm(mo)%dst_mask)
        end do 
      end do
! DOUBLE CHECK THESE LINES ARE OKAY 
       do mw=1,Ngrids_swan
          deallocate(ngrd_sw(mw)%src_mask)
       end do
       do mo=1,Ngrids_roms
          deallocate(ngrd_rm(mo)%dst_mask)
       end do

      end subroutine wav2ocn_mask
          
!======================================================================

      subroutine ocn2atm_mask()
      implicit none 

!     local variables 
      character(char_len) :: grid1_file, grid2_file  
      character(char_len) :: interp_file1, interp_file2   
      character(char_len) :: map1_name, map2_name 
      character(char_len) :: mo_string, mw_string
      integer(int_kind)   :: i, j, c, N
      integer(int_kind)   :: mo, ma, nx, ny, grdsize
      integer(int_kind)   :: ncstat, nc_file_id

! Allocate the source mask 
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%src_mask(nx,ny))   
      end do 
! Allocate the destination mask 
      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        allocate(ngrd_wr(ma)%dst_mask(nx,ny))
      end do 
!  Save the ocean grid mask to src_mask 
      do mo=1,Ngrids_roms
        do j=1,ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,j)=ngrd_rm(mo)%mask_rho_o(i,j)
          end do 
        end do 
      end do 
!  Save the wrf grid mask to dst_mask 
      do ma=1,Ngrids_wrf
        do j=1,ngrd_wr(ma)%sn_size
          do i=1,ngrd_wr(ma)%we_size
            ngrd_wr(ma)%dst_mask(i,j)=ngrd_wr(ma)%mask_rho_a(i,j)
          end do
        end do 
      end do 
!  Determine destination masks for the WRF grids using ROMS as the source.
      if ((Ngrids_roms>0).and.(Ngrids_wrf>0)) then 
        do mo=1,Ngrids_roms
          if (Ngrids_roms>mo) then ! mask out child portion of this ocn grid
            do j=ngrd_rm(mo)%jstr_o,ngrd_rm(mo)%jend_o
              do i=ngrd_rm(mo)%istr_o,ngrd_rm(mo)%iend_o
                ngrd_rm(mo)%src_mask(i,j)=0
              end do 
            end do 
          end if 
        end do 
      end if 
!  Send the data to SCRIP routines
      do ma=1,Ngrids_wrf
        do mo=1,Ngrids_roms
          write(ma_string,'(i1)')ma
          write(mo_string,'(i1)')mo
          interp_file1='ocn'//trim(mo_string)//'_to_'//'atm'//
     &                       trim(ma_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='ROMS to WRF distwgt Mapping'
          map2_name='unknown' 
                 
          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", roms_grids(mo)
          write(stdout,*)"The dst grid is: ", wrf_grids(ma)

          grid1_file=roms_grids(mo)
          grid2_file=wrf_grids(ma)

          call scrip_package(grid1_file, grid2_file, 
     &                             ngrd_rm(mo)%xi_size,
     &                             ngrd_rm(mo)%eta_size,
     &                             ngrd_rm(mo)%lon_rho_o,
     &                             ngrd_rm(mo)%lat_rho_o,
     &                             ngrd_rm(mo)%x_full_grid, 
     &                             ngrd_rm(mo)%y_full_grid, 
     &                             ngrd_wr(ma)%we_size,
     &                             ngrd_wr(ma)%sn_size,
     &                             ngrd_wr(ma)%lon_rho_a, 
     &                             ngrd_wr(ma)%lat_rho_a, 
     &                             ngrd_wr(ma)%x_full_grid, 
     &                             ngrd_wr(ma)%y_full_grid, 
     &                              interp_file1, interp_file2, 
     &                                    map1_name, map2_name, 
     &                                    ngrd_rm(mo)%src_mask,
     &                                    ngrd_wr(ma)%dst_mask)
        end do 
      end do
! DOUBLE CHECK THESE LINES ARE OKAY 
       do mo=1,Ngrids_roms
          deallocate(ngrd_rm(mo)%src_mask)
       end do
       do ma=1,Ngrids_wrf
          deallocate(ngrd_wr(ma)%dst_mask)
       end do

      end subroutine ocn2atm_mask

!======================================================================

      subroutine atm2ocn_mask()
      implicit none 

!     local variables 
      character(char_len) :: grid1_file, grid2_file  
      character(char_len) :: interp_file1, interp_file2   
      character(char_len) :: map1_name, map2_name 
      character(char_len) :: mo_string, mw_string
      integer(int_kind)   :: i, j, c, N
      integer(int_kind)   :: mo, ma, nx, ny, grdsize
      integer(int_kind)   :: ncstat, nc_file_id

! Allocate the source mask 
      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        allocate(ngrd_wr(ma)%src_mask(nx,ny))
      end do 
! Allocate the destination mask 
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%dst_mask(nx,ny))   
      end do 
!  Save the wrf grid mask to src_mask 
      do ma=1,Ngrids_wrf
        do j=1,ngrd_wr(ma)%sn_size
          do i=1,ngrd_wr(ma)%we_size
            ngrd_wr(ma)%src_mask(i,j)=ngrd_wr(ma)%mask_rho_a(i,j)
          end do 
        end do 
      end do 
!  Save the grid mask to dst mask
      do mo=1,Ngrids_roms
        do j=1,ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%dst_mask(i,j)=ngrd_rm(mo)%mask_rho_o(i,j)
          end do
        end do 
      end do 
!  Determine destination masks for the WRF grids using ROMS as the source.
      if ((Ngrids_wrf>0).and.(Ngrids_roms>0)) then 
        do ma=1,Ngrids_wrf
          if (Ngrids_wrf>ma) then ! mask out child portion of this ocn grid
            do j=ngrd_wr(ma)%jstr_a,ngrd_wr(ma)%jend_a
              do i=ngrd_wr(ma)%istr_a,ngrd_wr(ma)%iend_a
                ngrd_wr(ma)%src_mask(i,j)=0
              end do 
            end do 
          end if 
        end do 
      end if 
!  Send the data to SCRIP routines
      do mo=1,Ngrids_roms
        do ma=1,Ngrids_wrf
          write(mo_string,'(i1)')mo
          write(ma_string,'(i1)')ma
          interp_file1='atm'//trim(ma_string)//'_to_'//'ocn'//
     &                       trim(mo_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='WRF to ROMS distwgt Mapping'
          map2_name='unknown' 
                 
          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", wrf_grids(ma)
          write(stdout,*)"The dst grid is: ", roms_grids(mo)

          grid1_file=wrf_grids(ma)
          grid2_file=roms_grids(mo)

          call scrip_package(grid1_file, grid2_file, 
     &                             ngrd_wr(ma)%we_size,
     &                             ngrd_wr(ma)%sn_size,
     &                             ngrd_wr(ma)%lon_rho_a, 
     &                             ngrd_wr(ma)%lat_rho_a, 
     &                             ngrd_wr(ma)%x_full_grid, 
     &                             ngrd_wr(ma)%y_full_grid, 
     &                             ngrd_rm(mo)%xi_size,
     &                             ngrd_rm(mo)%eta_size,
     &                             ngrd_rm(mo)%lon_rho_o,
     &                             ngrd_rm(mo)%lat_rho_o,
     &                             ngrd_rm(mo)%x_full_grid, 
     &                             ngrd_rm(mo)%y_full_grid, 
     &                              interp_file1, interp_file2, 
     &                                    map1_name, map2_name, 
     &                                    ngrd_wr(ma)%src_mask,
     &                                    ngrd_rm(mo)%dst_mask)
        end do 
      end do
! DOUBLE CHECK THESE LINES ARE OKAY 
       do ma=1,Ngrids_wrf
          deallocate(ngrd_wr(ma)%src_mask)
       end do
       do mo=1,Ngrids_roms
          deallocate(ngrd_rm(mo)%dst_mask)
       end do

      end subroutine atm2ocn_mask
          
      subroutine atm2wav_mask()
      implicit none 

!     local variables 
      character(char_len) :: grid1_file, grid2_file  
      character(char_len) :: interp_file1, interp_file2   
      character(char_len) :: map1_name, map2_name 
      character(char_len) :: mo_string, mw_string
      integer(int_kind)   :: i, j, c, N
      integer(int_kind)   :: mo, ma, nx, ny, grdsize
      integer(int_kind)   :: ncstat, nc_file_id

! Allocate the source mask 
      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        allocate(ngrd_wr(ma)%src_mask(nx,ny))
      end do 
! Allocate the destination mask 
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        allocate(ngrd_sw(mw)%dst_mask(nx,ny))   
      end do 
!  Save the wrf grid mask to src_mask 
      do ma=1,Ngrids_wrf
        do j=1,ngrd_wr(ma)%sn_size
          do i=1,ngrd_wr(ma)%we_size
            ngrd_wr(ma)%src_mask(i,j)=ngrd_wr(ma)%mask_rho_a(i,j)
          end do 
        end do 
      end do 
!  Save the swan grid mask to dst mask
      do mw=1,Ngrids_swan
        do j=1,ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%dst_mask(i,j)=ngrd_sw(mw)%mask_rho_w(i,j)
          end do
        end do 
      end do 
!  Determine destination masks for the SWAN grids using WRF as the source.
      if ((Ngrids_wrf>0).and.(Ngrids_swan>0)) then 
        do ma=1,Ngrids_wrf
          if (Ngrids_wrf>ma) then ! mask out child portion of this ocn grid
            do j=ngrd_wr(ma)%jstr_a,ngrd_wr(ma)%jend_a
              do i=ngrd_wr(ma)%istr_a,ngrd_wr(ma)%iend_a
                ngrd_wr(ma)%src_mask(i,j)=0
              end do 
            end do 
          end if 
        end do 
      end if 
!  Send the data to SCRIP routines
      do mw=1,Ngrids_swan
        do ma=1,Ngrids_wrf
          write(mw_string,'(i1)')mw
          write(ma_string,'(i1)')ma
          interp_file1='wav'//trim(ma_string)//'_to_'//'atm'//
     &                       trim(mw_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='WRF to SWAN distwgt Mapping'
          map2_name='unknown' 
                 
          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", wrf_grids(ma)
          write(stdout,*)"The dst grid is: ", swan_coord(mw)

          grid1_file=wrf_grids(ma)
          grid2_file=swan_coord(mw)

          call scrip_package(grid1_file, grid2_file, 
     &                             ngrd_wr(ma)%we_size,
     &                             ngrd_wr(ma)%sn_size,
     &                             ngrd_wr(ma)%lon_rho_a, 
     &                             ngrd_wr(ma)%lat_rho_a, 
     &                             ngrd_wr(ma)%x_full_grid, 
     &                             ngrd_wr(ma)%y_full_grid, 
     &                             ngrd_sw(mw)%Numx_swan,
     &                             ngrd_sw(mw)%Numy_swan,
     &                             ngrd_sw(mw)%lon_rho_w,
     &                             ngrd_sw(mw)%lat_rho_w,
     &                             ngrd_sw(mw)%x_full_grid, 
     &                             ngrd_sw(mw)%y_full_grid, 
     &                              interp_file1, interp_file2, 
     &                                    map1_name, map2_name, 
     &                                    ngrd_wr(ma)%src_mask,
     &                                    ngrd_sw(mw)%dst_mask)
        end do 
      end do
! DOUBLE CHECK THESE LINES ARE OKAY 
       do ma=1,Ngrids_wrf
          deallocate(ngrd_rm(ma)%src_mask)
       end do
       do mw=1,Ngrids_swan
          deallocate(ngrd_sw(mw)%dst_mask)
       end do

      end subroutine atm2wav_mask
          
      subroutine wav2atm_mask()
      implicit none 

!     local variables 
      character(char_len) :: grid1_file, grid2_file  
      character(char_len) :: interp_file1, interp_file2   
      character(char_len) :: map1_name, map2_name 
      character(char_len) :: mo_string, mw_string
      integer(int_kind)   :: i, j, c, N
      integer(int_kind)   :: mo, ma, nx, ny, grdsize
      integer(int_kind)   :: ncstat, nc_file_id

! Allocate the source mask 
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        allocate(ngrd_sw(mw)%src_mask(nx,ny))   
      end do 
! Allocate the destination mask 
      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        allocate(ngrd_wr(ma)%dst_mask(nx,ny))
      end do 
!  Save the wave grid mask to src_mask 
      do mw=1,Ngrids_swan
        do j=1,ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%src_mask(i,j)=ngrd_sw(mw)%mask_rho_w(i,j)
          end do
        end do 
!  Save the wrf grid mask to dst mask
      do ma=1,Ngrids_wrf
        do j=1,ngrd_wr(ma)%sn_size
          do i=1,ngrd_wr(ma)%we_size
            ngrd_wr(ma)%dst_mask(i,j)=ngrd_wr(ma)%mask_rho_a(i,j)
          end do 
        end do 
      end do 
      end do 
!  Determine destination masks for the WRF grids using SWAN as the source.
      if ((Ngrids_wrf>0).and.(Ngrids_swan>0)) then 
        do mw=1,Ngrids_swan
          if (Ngrids_swan>mw) then ! mask out child portion of this ocn grid
            do j=ngrd_sw(mw)%jstr_w,ngrd_sw(mw)%jend_w
              do i=ngrd_sw(mw)%istr_w,ngrd_sw(mw)%iend_w
                ngrd_sw(mw)%src_mask(i,j)=0
              end do 
            end do 
          end if 
        end do 
      end if 
!  Send the data to SCRIP routines
      do ma=1,Ngrids_wrf
        do mw=1,Ngrids_swan
          write(ma_string,'(i1)')ma
          write(mw_string,'(i1)')mw
          interp_file1='atm'//trim(ma_string)//'_to_'//'wav'//
     &                       trim(mw_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='SWAN to WRF distwgt Mapping'
          map2_name='unknown' 
                 
          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", swan_coord(mw)
          write(stdout,*)"The dst grid is: ", wrf_grids(ma)

          grid1_file=swan_coord(mw)
          grid2_file=wrf_grids(ma)

          call scrip_package(grid1_file, grid2_file, 
     &                             ngrd_sw(mw)%Numx_swan,
     &                             ngrd_sw(mw)%Numy_swan,
     &                             ngrd_sw(mw)%lon_rho_w,
     &                             ngrd_sw(mw)%lat_rho_w,
     &                             ngrd_sw(mw)%x_full_grid, 
     &                             ngrd_sw(mw)%y_full_grid, 
     &                             ngrd_wr(ma)%we_size,
     &                             ngrd_wr(ma)%sn_size,
     &                             ngrd_wr(ma)%lon_rho_a, 
     &                             ngrd_wr(ma)%lat_rho_a, 
     &                             ngrd_wr(ma)%x_full_grid, 
     &                             ngrd_wr(ma)%y_full_grid, 
     &                              interp_file1, interp_file2, 
     &                                    map1_name, map2_name, 
     &                                    ngrd_sw(mw)%src_mask,
     &                                    ngrd_wr(ma)%dst_mask)
        end do 
      end do
! DOUBLE CHECK THESE LINES ARE OKAY 
       do mw=1,Ngrids_swan
          deallocate(ngrd_sw(mw)%src_mask)
       end do
       do ma=1,Ngrids_wrf
          deallocate(ngrd_wr(ma)%dst_mask)
       end do

      end subroutine wav2atm_mask

      end module create_masks 
