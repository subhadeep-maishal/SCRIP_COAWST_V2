!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module helps create the masks between each model   
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
!--------------Date: 08/20/2015----------------------------------------
!***********************************************************************

      module create_masks

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types 
      use scripwrap_mod ! defines input file variable names 

      implicit none

      contains

!      PRIVATE :: create_scrip_roms_infile 
!      PUBLIC  :: ocn2wav_mask

      subroutine ocn2wav_mask()
      implicit none 

!     local variables 
      integer(int_kind) :: i, j, c, N
      integer(int_kind) :: mo, mw, nx, ny, grdsize

      real(dbl_kind) :: area_src_box, area_total
      real(dbl_kind), allocatable :: zz(:)
      real(dbl_kind), allocatable :: grd(:)

      real(dbl_kind) :: x, y
      real(dbl_kind) :: x_src(4), y_src(4)
      real(dbl_kind) :: px, py

! Allocate the destination mask 
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        allocate(ngrd_sw(mw)%dst_mask(nx,ny))
      end do 

! Allocate the destination mask 
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%src_mask(nx,ny))   
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
            ngrd_sw(mw)%dst_mask(i,j)=0
          end do
        end do 
      end do 
!
!  Determine destination masks for the swan grids using ocean as the source.
      if ((Ngrids_swan>0).and.(Ngrids_roms>0)) then 
        do mo=1,Ngrids_roms
          if (Ngrids_roms>mo) then ! mask out child portion of this ocn grid
            do j=1,ngrd_rm(mo)%eta_size
              do i=1,ngrd_rm(mo)%xi_size
                ngrd_rm(mo)%src_mask(i,j)=0
              end do 
            end do 
          end if 
!
!          do mw=1,Ngrids_swan
!            do j=1,ngrd_sw(mw)%Numy_swan(mw)
!              do i=1,ngrd_sw(mw)%Numx_swan(mw)
!                px=ngrd_sw(mw)%lon_rho_w(i,j)
!                py=ngrd_sw(mw)%lat_rho_w(i,j)

! find out if dst point lies in/on/out of the src grd
! if it lies in or on the src grd then find the nearest
! src grd point from that dst point and make its mask = mask of src point 
           
!  Compute vertices of the source box 
!                x_src(1)=ngrd_rm(mo)%lon_rho_o(1,1)
!                y_src(1)=ngrd_rm(mo)%lat_rho_o(1,1)
             
!                x_src(2)=ngrd_rm(mo)%lon_rho_o(ngrd_rm(mo)%xi_size,1)
!                y_src(2)=ngrd_rm(mo)%lat_rho_o(ngrd_rm(mo)%xi_size,1)
            
!                x_src(3)=ngrd_rm(mo)%lon_rho_o(ngrd_rm(mo)%xi_size,     &
!     &                   ngrd_rm(mo)%eta_size)
!                y_src(3)=ngrd_rm(mo)%lat_rho_o(ngrd_rm(mo)%xi_size,     &
!     &                   ngrd_rm(mo)%eta_size)

!                x_src(4)=ngrd_rm(mo)%lon_rho_o(1,ngrd_rm(mo)%eta_size)
!                y_src(4)=ngrd_rm(mo)%lat_rho_o(1,ngrd_rm(mo)%eta_size)
            
!  Get the area of the polygon formed by the ocean grid 
!                area_src_box=polyarea(x_src, y_src, 4)
                 
!  Check if the wave grid has a point inside/on/outside ocean grid
!                area_total=area_total(x_src, y_src, px, py)

!                if (area_src_box>=area_total) then 
!                  c=1
! This point will get its masking from the nearest source point 
!                  ngrd_rm(mw)%dst_mask(i,j)=1
!                else
!                  c=0
!                  ngrd_rm(mw)%dst_mask(i,j)=0
!                end if
!              end do 
!            end do 
!          end do 
!
        end do 
      end if 


      call create_scrip_roms_infile()
      end subroutine ocn2wav_mask

      subroutine create_scrip_roms_infile() 
!     This should generate a src.nc file 
      use kinds_mod     ! defines common data types 
      use scripwrap_mod ! defines input file variable names 
      use netcdf_mod    !netCDF include file and a netcdf error handling routine.

      implicit none 

      integer(int_kind) :: i, j, iunit
      integer(int_kind) :: nx, ny, mw
      integer(int_kind)  :: ncstat, nc_file_id, nc_grdsize_id,
     &         nc_grdlat_id, nc_grdlon_id, nc_grdmsk_id

      real(dbl_kind), allocatable :: lat_psi(:,:), lon_psi(:,:) 
      do mw=1,Ngrids_roms

        allocate(lat_psi
       ! Open the file. 
        ncstat=nf_open(roms_grids(mw),nf_nowrite,nc_file_id)
        call netcdf_error_handler(ncstat)

!     Read dimension id 
!        ncstat=nf_inq_dimid(nc_file_id,'xi_rho',nc_grdsize_id)
!        call netcdf_error_handler(ncstat)
!     Get the grid size in each direction 
!        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,
!     &                            ngrd_rm(mw)%xi_size)
!        call netcdf_error_handler(ncstat)

!        ncstat=nf_inq_dimid(nc_file_id,'eta_rho',nc_grdsize_id)
!        call netcdf_error_handler(ncstat)
!        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,
!     &                           ngrd_rm(mw)%eta_size)
!        call netcdf_error_handler(ncstat)

!       HOW TO KNOW THE SIZE 
!        ncstat = nf_inq_varid(nc_file_id,'spherical',checksphere_id)
!        call netcdf_error_handler(ncstat)

!     Read variable id 
!        ncstat = nf_get_var_char(nc_file_id,checksphere_id,'spherical')
!        call netcdf_error_handler(ncstat)

!         if ('spherical'.eq.0.or.'spherical'.eq.'F'.or.                &
!     &                          'spherical'.eq.'f') then           
!            ncstat=nf_inq_varid(nc_file_id,'x_psi',nc_grdlon_id)
!            call netcdf_error_handler(ncstat)
!            ncstat=nf_inq_varid(nc_file_id,'y_psi',nc_grdlat_id)
!            call netcdf_error_handler(ncstat)
!          else   
        ncstat=nf_inq_varid(nc_file_id,'lon_psi',nc_grdlon_id)
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_varid(nc_file_id,'lat_psi',nc_grdlat_id)
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_varid(nc_file_id,'lon_rho',nc_grdlon_id)
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_varid(nc_file_id,'lat_rho',nc_grdlon_id)
        call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_double(nc_file_id,nc_grdlon_id,               &
     &                                     ngrd_rm(mw)%lon_psi)
        call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_double(nc_file_id,nc_grdlat_id,               &
     &                                     ngrd_rm(mw)%lat_psi)
        call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_double(nc_file_id, nc_grdlon_id,              &
     &                                     ngrd_rm(mw)%lon_rho_o)
        call netcdf_error_handler(ncstat)
        ncstat=nf_get_var_double(nc_file_id, nc_grdlat_id,              &
     &                                     ngrd_rm(mw)%lat_rho_o)
        call netcdf_error_handler(ncstat)

!       send lat_psi, lon_psi
        nx_psi = ngrd_rm(mw)%xi_size+1
        ny_psi = ngrd_rm(mw)%eta_size+1
        call create_extra_rho_grid(nx_psi,ny_psi,lon_psi,lat_psi,       &
                                   x_full_grid,y_full_grid)
          
!        

      end do 
      end subroutine create_scrip_roms_infile() 

      subroutine create_extra_rho_grid(nx_psi,ny_psi,x_rho,y_rho,       &
                                        x_full_grid,y_full_grid)

       implicit none 

       integer(int_kind), intent(in) :: nx_psi, ny_psi
       real(dbl_kind), intent(in)    :: x_rho(nx_psi,ny_psi)
       real(dbl_kind), intent(in)    :: y_rho(nx_psi,ny_psi)
       real(dbl_kind), intent(out)   :: x_full_grid(nx_psi+1,ny_psi+1)
       real(dbl_kind), intent(out)   :: y_full_grid(nx_psi+1,ny_psi+1)
!      local variables
       integer :: i, j
       integer :: iend, jend

       jend=ny_psi
       iend=nx_psi 

!      new grid should have (iend,jend)
       do j=2,jend
         do i=2,iend
           !copy interior values 
           x_full_grid(i,j)=x_rho(i-1,j-1) 
         end do
       end do  
       do i=2,iend-1
           !extrapolate South boundary
         x_full_grid(i,1)=x_rho(i,1)-(x_rho(i,2)-x_rho(i,1))
       end do 
       do i=2,iend-1
           !extrapolate North boundary
         x_full_grid(i,jend)=x_rho(i,jend)-(x_rho(i,jend)-              &
     &                                                x_rho(i,jend-1))
       end do 
       do j=2,jend-1
           !extrapolate West boundary
         x_full_grid(1,j)=x_rho(1,j)-(x_rho(2,j)-x_rho(i,1))
       end do 
       do j=2,jend-1
           !extrapolate East boundary
         x_full_grid(iend,j)=x_rho(iend,j)-(x_rho(iend,j)-        &
                                                x_rho(iend-1,j))
       end do

!         NW      NE
!            _____
!           |     |
!           |     |
!           |_ _ _|
!
!         SW      SE
!
       ! Now extrapolate to four corners 
       ! SouthEast 
       x_full_grid(1,1)=x_full_grid(1,2)-                               &
     &                   (x_full_grid(1,3)-x_full_grid(1,2))
       ! NorthWest 
       x_full_grid(1,jend)=x_full_grid(1,jend-1)+(x_full_grid(1,jend-1)-&
     &                        x_full_grid(1,jend-2))
       ! SouthWest 
       x_full_grid(iend,1)=x_full_grid(iend,2)-(x_full_grid(iend,3)-    &
     &                        x_full_grid(iend,2));
       ! NorthEast 
       x_full_grid(iend,jend)=x_full_grid(iend,jend-1)+                 &
     &              (x_full_grid(iend,jend)-x_full_grid(iend,jend-1))

!-------------------------------------------------------------------------
       do j=2,jend
         do i=2,iend
           !copy interior values 
           y_full_grid(i,j)=y_rho(i,j) 
         end do
       end do  

       do i=2,iend-1
           !extrapolate South boundary
           y_full_grid(i,1)=y_rho(i,1)-(y_rho(i,2)-y_rho(i,1))
       end do 
       do i=2,iend-1
           !extrapolate North boundary
         y_full_grid(i,jend)=y_rho(i,jend)-(y_rho(i,jend)-              &
     &                                                y_rho(i,jend-1))
       end do 
       do j=2,jend-1
           !extrapolate West boundary
         y_full_grid(1,j)=y_rho(1,j)-(y_rho(2,j)-y_rho(i,1))
       end do 
       do j=2,jend-1
           !extrapolate East boundary
         y_full_grid(iend,j)=y_rho(iend,j)-(y_rho(iend,j)-              &
                                                y_rho(iend-1,j))
       end do

!         NW      NE
!            _____
!           |     |
!           |     |
!           |_ _ _|
!
!         SW      SE
!
       ! Now extrapolate to four corners 
       ! SouthEast 
       y_full_grid(1,1)=y_full_grid(1,2)-                               &
     &                   (y_full_grid(1,3)-y_full_grid(1,2))
       ! NorthWest 
       y_full_grid(1,jend)=y_full_grid(1,jend-1)+(y_full_grid(1,jend-1)-&
     &                     y_full_grid(1,jend-2))
       ! SouthWest 
       y_full_grid(iend,1)=y_full_grid(iend,2)-(y_full_grid(iend,3)-    &
     &                     y_full_grid(iend,2));
       ! NorthEast 
       y_full_grid(iend,jend)=y_full_grid(iend,jend-1)+                 &
     &              (y_full_grid(iend,jend-1)-y_full_grid(iend,jend-2))
 
        end do 

        end subroutine create_scrip_roms_infile
      end module create_masks 
!        ncstat=nf_inq_varid(nc_file_id,'mask_rho',nc_grdmsk_id)
!        call netcdf_error_handler(ncstat)

!     load the roms grid
!%netcdf_load(roms_grid_file)
!        x_psi=ncread(roms_grid_file,'lon_psi');
!        y_psi=ncread(roms_grid_file,'lat_psi');
!        %mask_rho=ncread(roms_grid_file,'mask_rho');
!        mask_rho=mask;
!        spherical=ncread(roms_grid_file,'spherical');

![LP,MP]=size(h);
!%gridsize=LP*MP;

!if ((spherical=='F')||(spherical=='f'))||(spherical==0)
!%  lon_rho=x_rho;
!%  lat_rho=y_rho;
!%  lon_psi=x_psi;
!%  lat_psi=y_psi;
!   x_psi=ncread(roms_grid_file,'x_psi');
!   y_psi=ncread(roms_grid_file,'y_psi');
!   x_rho=ncread(roms_grid_file,'x_rho');
!   y_rho=ncread(roms_grid_file,'y_rho');
!   projection='mercator';
!   m_proj(projection);
!   [lon_rho, lat_rho] = m_xy2ll(x_rho/6371000, y_rho/6371000);   % Degrees.
!   [lon_psi, lat_psi] = m_xy2ll(x_psi/6371000, y_psi/6371000);   % Degrees.
!else
!  lon_psi=ncread(roms_grid_file,'lon_psi');
!  lat_psi=ncread(roms_grid_file,'lat_psi');
!  lon_rho=ncread(roms_grid_file,'lon_rho');
!  lat_rho=ncread(roms_grid_file,'lat_rho');
!end

!%create a full set of psi points
![x_full_grid,y_full_grid]=create_extra_rho_grid(lon_psi,lat_psi);

!%create the srcip netcdf grid file
!nc = netcdf.create(out_file, 'clobber');

!%% Global attributes:
!disp(' ## Defining Global Attributes...')

!netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','Scrip file');
!netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',['Created by ', mfilename ', on ', datestr(now)]);

!L = LP-1;       % The psi dimension.
!M = MP-1;       % The psi dimension.

!disp(' ## Defining Dimensions...')

!!grid_size = netcdf.defDim(nc,'grid_size',LP*MP);
!grid_corners = netcdf.defDim(nc,'grid_corners',4);
!grid_rank = netcdf.defDim(nc,'grid_rank',2);

!disp(' ## Defining Variables')

!v1 = netcdf.defVar(nc,'grid_dims','short',grid_rank);
!netcdf.putAtt(nc,v1,'long_name','grid dimensions')
!netcdf.putAtt(nc,v1,'units','---')
!
!v2 = netcdf.defVar(nc,'grid_imask','short',grid_size);
!netcdf.putAtt(nc,v2,'long_name','grid masking')
!netcdf.putAtt(nc,v2,'units','---')
!
!v3 = netcdf.defVar(nc,'grid_center_lat','double',grid_size);
!netcdf.putAtt(nc,v3,'units','radians')
!
!v4 = netcdf.defVar(nc,'grid_center_lon','double',grid_size);
!netcdf.putAtt(nc,v4,'units','radians')
!
!v5 = netcdf.defVar(nc,'grid_corner_lat','double',[grid_corners grid_size]);
!netcdf.putAtt(nc,v5,'units','radians')
!
!v6 = netcdf.defVar(nc,'grid_corner_lon','double',[grid_corners grid_size]);
!netcdf.putAtt(nc,v6,'units','radians')

!netcdf.endDef(nc)
!
!%now fill that netcdf file
!netcdf.putVar(nc,v1, [MP, LP]);

!disp('step 0/4, filling grid mask')
!%get grid mask
!count=0;
!for jj=1:MP
!  for ii=1:LP
!    count=count+1;
!    ztemp(count)=mask_rho(ii,jj);
!  end
!end
!netcdf.putVar(nc,v2,ztemp);
!clear ztemp

!scale=pi/180;

!%get grid centers
!disp('step 1/4, filling grid lat centers')
!ztemp=reshape(lat_rho,MP*LP,1);
!netcdf.putVar(nc,v3,ztemp*scale);
!clear ztemp

!disp('step 2/4, filling grid lon centers')
!ztemp=reshape(lon_rho,MP*LP,1);
!!netcdf.putVar(nc,v4,ztemp*scale);
clear ztemp

!%get grid corners, counterclockwise
!disp('step 3/4, filling grid lat corners')
!c1=reshape(y_full_grid(1:LP,1:MP),MP*LP,1);
!c2=reshape(y_full_grid(2:LP+1,1:MP),MP*LP,1);
!c3=reshape(y_full_grid(2:LP+1,2:MP+1),MP*LP,1);
!c4=reshape(y_full_grid(1:LP,2:MP+1),MP*LP,1);
!ztemp(:,:)=[c1 c2 c3 c4].';
!netcdf.putVar(nc,v5,ztemp*scale);
!clear ztemp

!disp('step 4/4, filling grid lon corners')
!c1=reshape(x_full_grid(1:LP,1:MP),MP*LP,1);
!c2=reshape(x_full_grid(2:LP+1,1:MP),MP*LP,1);
!c3=reshape(x_full_grid(2:LP+1,2:MP+1),MP*LP,1);
!c4=reshape(x_full_grid(1:LP,2:MP+1),MP*LP,1);
!ztemp(:,:)=[c1 c2 c3 c4].';
!netcdf.putVar(nc,v6,ztemp*scale);
!clear ztemp
!
!%close the file.
!netcdf.close(nc)



!======================================================================

!*****************************START OF FUNCTIONS************************* 
!      real(dbl_kind) function polyarea(x, y, N)
!      integer(int_kind), intent(in) :: i, j, N 
!      real(dbl_kind), intent(in)    :: x(N), y(N), area
!
!      area = 0.0_r8
!
!      j=N-1      
!      do i=1,N       
!        if (i==(N)) then 
!          j=1  
!        else 
!          j=i+1
!        end if 
!        area = area+( x(j)+x(i))*(y(j)-y(i) )
!      end do 
!
!      area=area*0.5_r8 
!      return 
!      end function polyarea 
!
!      real(dbl_kind) function area_total(x_src, y_src, px, py, N)
!      integer(int_kind), intent(in):: N 
!      real(dbl_kind), intent(in)   :: px, py
!      real(dbl_kind), intent(in)   :: x_src(N), y_src(N) 
!      integer(int_kind)            :: i, j 
!      real(dbl_kind)               :: x(4), y(4), area(4)
!      real(dbl_kind)               :: area_tot
 
!      area_total=0.0_r8
!      do j=1,N
!        do i=1,N
!          if (i==j) then 
!            x(i)=px 
!            y(i)=py
!          else 
!            x(i)=x_src(i)
!            y(i)=y_src(i)
!          end if 
!        end do 
!        area(j)=polyarea(x, y, 4)
!      end do 

!      do j=1,N
!        area_total=area(j)+area_total
!      end do 

!      end function area_total

!======================================================================


!***********************************************************************

