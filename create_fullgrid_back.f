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
      module create_fullgrid
!-----------------------------------------------------------------------

      use kinds_mod    ! defines data types

      contains
 
      subroutine create_extra_rho_grid(nx_psi, ny_psi, x_psi, y_psi,
     &                                   x_full_psi, y_full_psi)

      implicit none

      integer(int_kind), intent(in) :: nx_psi, ny_psi
      real(dbl_kind), intent(in)    :: x_psi(nx_psi,ny_psi)
      real(dbl_kind), intent(in)    :: y_psi(nx_psi,ny_psi)
      real(dbl_kind), intent(out)   :: x_full_psi(nx_psi+2,ny_psi+2)
      real(dbl_kind), intent(out)   :: y_full_psi(nx_psi+2,ny_psi+2)

      real(dbl_kind) :: x_full_psi1(nx_psi+2,ny_psi+2)
      real(dbl_kind) :: y_full_psi1(nx_psi+2,ny_psi+2)
      ! local variables
      integer :: i, j
      integer :: iend, jend

      iend=nx_psi+2
      jend=ny_psi+2
                   
      ! new grid should have (iend,jend)
      do j=2,jend-1
        do i=2,iend-1
         !copy interior values
          x_full_psi(i,j)=x_psi(i-1,j-1)
          y_full_psi(i,j)=y_psi(i-1,j-1)
        end do
      end do
      do i=2,iend-1
          !extrapolate South boundary
        x_full_psi(i,1)=x_psi(i-1,1)-(x_psi(i-1,2)-x_psi(i-1,1))
        y_full_psi(i,1)=y_psi(i-1,1)-(y_psi(i-1,2)-y_psi(i-1,1))
      end do
      do i=2,iend-1
          !extrapolate North boundary
        x_full_psi(i,jend)=x_psi(i-1,jend-2)+(x_psi(i-1,jend-2)-       
     &                                     x_psi(i-1,jend-3))
      end do
      do j=2,jend-1
          !extrapolate West boundary
        x_full_psi(1,j)=x_psi(1,j-1)-(x_psi(2,j-1)-x_psi(1,j-1))
      end do
      do j=2,jend-1
          !extrapolate East boundary
        x_full_psi(iend,j)=x_psi(iend-2,j-1)+(x_psi(iend-2,j-1)-       
     &                                      x_psi(iend-3,j-1))
      end do
!______________________________________________________________________
!         NW      NE
!            _____
!           |     |
!           |     |
!           |_ _ _|
!         SW      SE
!_____________________________________________________________________
      ! Now extrapolate to four corners 
      ! SouthEast 
      x_full_psi(1,1)=x_full_psi(1,2)-(x_full_psi(1,3)
     &                                 -x_full_psi(1,2))
      ! NorthWest 
      x_full_psi(1,jend)=x_full_psi(1,jend-1)+(x_full_psi(1,jend-1)-
     &                                          x_full_psi(1,jend-2))
      ! SouthWest 
      x_full_psi(iend,1)=x_full_psi(iend,2)-(x_full_psi(iend,3)-    
     &                                        x_full_psi(iend,2))
       ! NorthEast 
      x_full_psi(iend,jend)=x_full_psi(iend,jend-1)+                 
     &             (x_full_psi(iend,jend-1)-x_full_psi(iend,jend-2))

!-------------------------------------------------------------------------
      do j=2,jend-1
        do i=2,iend-1
          !copy interior values 
          y_full_psi(i,j)=y_psi(i-1,j-1)
        end do
      end do

      do i=2,iend-1
           !extrapolate South boundary
          y_full_psi(i,1)=y_psi(i-1,1)-(y_psi(i-1,2)-y_psi(i-1,1))
      end do
      do i=2,iend-1
          !extrapolate North boundary
        y_full_psi(i,jend)=y_psi(i-1,jend-2)+(y_psi(i-1,jend-2)-      
     &                                       y_psi(i-1,jend-3))
      end do
      do j=2,jend-1
          !extrapolate West boundary
        y_full_psi(1,j)=y_psi(1,j-1)-(y_psi(2,j-1)-y_psi(1,j-1))
      end do
      do j=2,jend-1
          !extrapolate East boundary
        y_full_psi(iend,j)=y_psi(iend-2,j-1)+(y_psi(iend-2,j-1)-       
     &                                       y_psi(iend-3,j-1))
      end do
!______________________________________________________________________
!         NW      NE
!            _____
!           |     |
!           |     |
!           |_ _ _|
!         SW      SE
!_____________________________________________________________________
      ! Now extrapolate to four corners 
      ! SouthEast 
      y_full_psi(1,1)=y_full_psi(1,2)-                               
     &               (y_full_psi(1,3)-y_full_psi(1,2))
       ! NorthWest 
      y_full_psi(1,jend)=y_full_psi(1,jend-1)+(y_full_psi(1,jend-1)-
     &                                          y_full_psi(1,jend-2))
       ! SouthWest 
      y_full_psi(iend,1)=y_full_psi(iend,2)-(y_full_psi(iend,3)-    
     &                                        y_full_psi(iend,2))
      ! NorthEast 
      y_full_psi(iend,jend)=y_full_psi(iend,jend-1)+                 
     &              (y_full_psi(iend,jend-1)-y_full_psi(iend,jend-2))

      do j=1,jend
        do i=1,iend
          write(60,*)x_full_psi(i,j),y_full_psi(i,j)
         end do 
       end do 
       do j=1,jend-2
         do i=1,iend-2
           write(61,*)x_psi(i,j),y_psi(i,j)
          end do 
        end do 
      end subroutine create_extra_rho_grid
!======================================================================

!***********************************************************************
      end module create_fullgrid 
