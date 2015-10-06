!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains routines reading ROMS grids  
!
!-----------------------------------------------------------------------
!     This is a wrapper around SCRIP package developed at USGS, 
!     Woods Hole Science And Marine Center. This routine calls the 
!     SCRIP's package routine that is the driver for computing the 
!     addresses and weights for interpolating between two grids on 
!     a sphere.
!
!---- Rewritten by John C. Warner-------------------------------------- 
!-----         Tarandeep S. Kalra -------------------------------------
!--------------Date: 08/04/2015----------------------------------------
!***********************************************************************

      module read_roms

!-----------------------------------------------------------------------
      use kinds_mod      ! defines common data types 
      use scripwrap_mod  ! defines input file variable names for SCRIP-COAWST wrapper
      use create_fullgrid! call the subroutine to create full grid of psi points 
      use netcdf_mod     !netCDF include file and a netcdf error handling routine.
 
      contains 

      subroutine load_roms_grid()

      implicit none

!      include 'netcdf.inc'

      integer(int_kind) :: i, j, iunit
      integer(int_kind) :: nx, ny, mo 
      integer(int_kind)  :: ncstat, nc_file_id, nc_grdsize_id, 
     &                      nc_grdlat_id, nc_grdlon_id, 
     &                      nc_grdcrnrlat_id, nc_grdcrnrlon_id,
     &                      nc_grdmsk_id

      real(dbl_kind) :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind) :: dist1, dist_max
      real(dbl_kind), allocatable :: lon_psi_o(:,:), lat_psi_o(:,:)

      allocate(ngrd_rm(Ngrids_roms))

      do mo=1,Ngrids_roms

!     Open the file. 
        ncstat=nf_open(roms_grids(mo),nf_nowrite,nc_file_id)
        call netcdf_error_handler(ncstat)

!     Read dimension id 
        ncstat=nf_inq_dimid(nc_file_id,'xi_rho',nc_grdsize_id) 
        call netcdf_error_handler(ncstat)
!     Get the grid size in each direction 
        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,                  
     &                            ngrd_rm(mo)%xi_size)
        call netcdf_error_handler(ncstat)

        ncstat=nf_inq_dimid(nc_file_id,'eta_rho',nc_grdsize_id) 
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,
     &                           ngrd_rm(mo)%eta_size)
        call netcdf_error_handler(ncstat)

!     Read variable id 
!        ncstat = nf_inq_varid(nc_file_id, 'spherical',checksphere_id)
!        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_varid(nc_file_id,'lon_rho',nc_grdlon_id)
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_varid(nc_file_id,'lat_rho',nc_grdlat_id)
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_varid(nc_file_id,'mask_rho',nc_grdmsk_id)
        call netcdf_error_handler(ncstat)

        ncstat=nf_inq_varid(nc_file_id,'lon_psi',nc_grdcrnrlon_id)
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_varid(nc_file_id,'lat_psi',nc_grdcrnrlat_id)
        call netcdf_error_handler(ncstat)

!     Allocate arrays
        allocate(ngrd_rm(mo)%
     &           lon_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%
     &           lat_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%
     &           mask_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%x_full_grid(ngrd_rm(mo)%xi_size+1,
     &                                   ngrd_rm(mo)%eta_size+1))
        allocate(ngrd_rm(mo)%y_full_grid(ngrd_rm(mo)%xi_size+1
     &                                  ,ngrd_rm(mo)%eta_size+1))

!     Get variables 
!       ncstat = nf_get_var_double(nc_file_id, checksphere_id,spherical)
!       call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_double(nc_file_id, nc_grdlon_id, 
     &                                     ngrd_rm(mo)%lon_rho_o)
        call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_double(nc_file_id, nc_grdlat_id, 
     &                                     ngrd_rm(mo)%lat_rho_o)
        call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_int(nc_file_id, nc_grdmsk_id, 
     &                                     ngrd_rm(mo)%mask_rho_o)
        call netcdf_error_handler(ncstat)

!       Making lon_psi and lat_psi as local arrays 
!       Local arrays need to deallocated, no_psi_pts=no_rho_pts-1
        allocate(lon_psi_o(ngrd_rm(mo)%xi_size-1,
     &                     ngrd_rm(mo)%eta_size-1))
        allocate(lat_psi_o(ngrd_rm(mo)%xi_size-1,
     &                     ngrd_rm(mo)%eta_size-1))

        ncstat=nf_get_var_double(nc_file_id, nc_grdcrnrlon_id,
     &                                       lon_psi_o)
        call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_double(nc_file_id, nc_grdcrnrlat_id,
     &                                       lat_psi_o)
        call netcdf_error_handler(ncstat)
 

        call create_extra_rho_grid(ngrd_rm(mo)%xi_size-1,
     &                             ngrd_rm(mo)%eta_size-1,
     &                             lon_psi_o, lat_psi_o,
     &                             ngrd_rm(mo)%x_full_grid, 
     &                             ngrd_rm(mo)%y_full_grid)

        deallocate(lon_psi_o, lat_psi_o) 
!
      end do

!   Find child grid with respect to parent grid 
! 
      do mo=1,Ngrids_roms-1
        allocate(ngrd_rm(mo)%istr_o,ngrd_rm(mo)%jstr_o,
     &           ngrd_rm(mo)%iend_o,ngrd_rm(mo)%jend_o)

        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size 
        dist_max=10e6
!       First the (2,2) point 
        xx2=ngrd_rm(mo+1)%lon_rho_o(2,2)
        yy2=ngrd_rm(mo+1)%lat_rho_o(2,2)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_rm(mo)%lon_rho_o(i,j)-xx2)**2+
     &                 (ngrd_rm(mo)%lat_rho_o(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_rm(mo)%istr_o=i
              ngrd_rm(mo)%jstr_o=j
            endif
          enddo
        enddo
        dist_max=10e6
!       The second last point 
        xxend_1=ngrd_rm(mo+1)%lon_rho_o(ngrd_rm(mo+1)%xi_size-1,        
     &                                ngrd_rm(mo+1)%eta_size-1 )
        yyend_1=ngrd_rm(mo+1)%lat_rho_o(ngrd_rm(mo+1)%xi_size-1,
     &                                ngrd_rm(mo+1)%eta_size-1 )
        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_rm(mo)%lon_rho_o(i,j)-xxend_1)**2+
     &                 (ngrd_rm(mo)%lat_rho_o(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_rm(mo)%iend_o=i
              ngrd_rm(mo)%jend_o=j
            endif
          enddo
        enddo
         print*,"ROMS starting i & j index of parent w.r.t child grid"
         print*,ngrd_rm(mo)%istr_o, ngrd_rm(mo)%jstr_o
         print*,"ending i & j index of parent w.r.t child grid--"
         print*,ngrd_rm(mo)%iend_o, ngrd_rm(mo)%jend_o
      enddo



      end subroutine load_roms_grid
!======================================================================

      end module read_roms 



