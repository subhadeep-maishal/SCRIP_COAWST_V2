!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module declares routines the global parameters for 
!     SCRIP_COAWST package 
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

      module read_swan

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types 
      use scripwrap_mod ! defines input file variable names 

      implicit none

      contains

      subroutine load_swan_grid()

      type get_grid
        integer(int_kind), allocatable :: istr_w(:), jstr_w(:),
     &                                    iend_w(:), jend_w(:)
        real(dbl_kind), allocatable :: xx(:,:)
        real(dbl_kind), allocatable :: lon_rho_w(:,:), lat_rho_w(:,:)
        real(dbl_kind) :: mask_rho_w
      end type get_grid

      type(get_grid), allocatable :: ngrd(:)

      allocate(ngrd(Ngrids_swan))

      do mw=1,Ngrids_swan
! Notice that the indices are reversed in order for these files to be read 
         nx=Numy_swan(mw)
         ny=Numx_swan(mw)

         allocate(ngrd(mw)%xx(nx,ny))
         open(unit=iunit,file=bath_file(mw),status='old')
         do i=1,nx 
           read(iunit,*) (ngrd(mw)%xx(i,j),j=1,ny)
           do j=1,ny  
             ngrd(mw)%mask_rho_w=ngrd(mw)%xx(i,j)
             if(ngrd(mw)%xx(i,j)==9999)then
               ngrd(mw)%mask_rho_w=0
             endif 
            end do                 
         end do
         close(iunit)
      end do 
!----------------------------------------------------------------------
!     Work on the SWAN_COORD input file 
!----------------------------------------------------------------------
      do mw=1,Ngrids_swan
        nx=Numx_swan(mw)
        ny=Numy_swan(mw)

        allocate(ngrd(mw)%lon_rho_w(nx,ny))
        allocate(ngrd(mw)%lat_rho_w(nx,ny))

        open(unit=iunit,file=swan_coord(mw),status='old')
        n=0
        do
          read(iunit,*,end=1)
          n=n+1
        end do
1       rewind(iunit)

        grdsize=n/2

        allocate(grd(n),zz(grdsize))

        do i=1,n
          read(iunit,*) bf
          grd(i)=bf
        end do
        do i=1,grdsize
          zz(i)=grd(i)
        end do

!  every point will have lon_rho,lat_rho
        do j=1,ny
          do i=1,nx
            ngrd(mw)%lon_rho_w(i,j)=zz(nx*(j-1)+i)
          end do
        end do

        do i=grdsize+1,n
          zz(i-grdsize)=grd(i)
        end do

        do j=1,ny
          do i=1,nx
            ngrd(mw)%lat_rho_w(i,j)=zz(nx*(j-1)+i)

          end do
        end do

        deallocate(grd, zz)

      end do

      do mw=1,Ngrids_swan-1
        allocate(ngrd(mw)%istr_w(Ngrids_swan-1),
     &           ngrd(mw)%jstr_w(Ngrids_swan-1),
     &           ngrd(mw)%iend_w(Ngrids_swan-1),
     &           ngrd(mw)%jend_w(Ngrids_swan-1))

        nx=Numx_swan(mw)
        ny=Numy_swan(mw)
        dist_max=10e6
        xx2=ngrd(mw+1)%lon_rho_w(2,2)
        yy2=ngrd(mw+1)%lat_rho_w(2,2)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd(mw)%lon_rho_w(i,j)-xx2)**2+
     &                 (ngrd(mw)%lat_rho_w(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd(mw)%istr_w(mw)=i
              ngrd(mw)%jstr_w(mw)=j
            endif
          enddo
        enddo

        dist_max=10e6
        xxend_1=
     &        ngrd(mw+1)%lon_rho_w(Numx_swan(mw+1)-1,Numy_swan(mw+1)-1)
        yyend_1= 
     &        ngrd(mw+1)%lat_rho_w(Numx_swan(mw+1)-1,Numy_swan(mw+1)-1)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd(mw)%lon_rho_w(i,j)-xxend_1)**2+
     &                 (ngrd(mw)%lat_rho_w(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd(mw)%iend_w(mw)=i
              ngrd(mw)%jend_w(mw)=j
            endif
          enddo
        enddo
         print*,"SWAN starting i & j index of parent w.r.t child grid--"
         print*,ngrd(mw)%istr_w(mw), ngrd(mw)%jstr_w(mw)
         print*,"ending i & j index of parent w.r.t child grid--"
         print*,ngrd(mw)%iend_w(mw), ngrd(mw)%jend_w(mw)
      enddo


      end subroutine load_swan_grid
!======================================================================


      end module read_swan

!***********************************************************************

