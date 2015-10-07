!======================================================================
!     This is a wrapper around SCRIP package developed at USGS, 
!     Woods Hole Science And Marine Center. This routine calls the 
!     SCRIP's package routine that is the driver for computing the 
!     addresses and weights for interpolating between two grids on 
!     a sphere.
!
!---- Written by John C. Warner-------------------------------------- 
!-----         Tarandeep S. Kalra -------------------------------------
!--------------Date: 08/04/2015----------------------------------------
!======================================================================

      program coawst_scrip 

      use kinds_mod
      use constants
      use iounits 
      use scripwrap_mod 
      use read_swan 
      use read_roms 
      use read_wrf
      use create_masks

      implicit none 

      open(unit=iunit,file='input_file',status='old',form='formatted')
      call read_inputs()
      close(iunit)

      call load_swan_grid()
      call load_roms_grid()
      call load_wrf_grid()
!      call ocn2wav_mask() 
!      call wav2ocn_mask() 
!      call wrf2ocn_mask() 
!       call ocn2wrf_mask() 
!      call wav2wrf_mask() 
       call wrf2wav_mask() 
!      call scrip_package()
       
      end program coawst_scrip

!======================================================================
      subroutine read_inputs()
      use kinds_mod
      use iounits 
      use scripwrap_mod
 
      implicit none 

!     local variables
      integer(int_kind) :: i

!     Allocate input variables
      namelist /inputs/ Ngrids_roms, Ngrids_swan, Ngrids_wrf, 
     &                  roms_grids, swan_coord, 
     &                  Numx_swan_temp, Numy_swan_temp, cartesian, 
     &                  cartesian, 
     &                  bath_file, wrf_grids

      write(stdout,*)"================================================"
      write(stdout,*) ' Read input_file for SCRIP_COAWST Wrapper '
      write(stdout,*)"================================================"
      read (iunit, inputs)
     
      write(stdout,*) "Ngrid_roms=",Ngrids_roms
      write(stdout,*) "Ngrid_swan=",Ngrids_swan
      write(stdout,*) "Ngrid_wrf =",Ngrids_wrf
      do i = 1,Ngrids_roms 
        write(stdout,*) "User input roms grid",i,"=", roms_grids(i)
      end do 
      do i = 1,Ngrids_swan
        write(stdout,*)"User input swan grid",i,"=", swan_coord(i)
        write(stdout,*)"Input bathymtry file",i,"=", bath_file(i)
      end do 
      do i = 1,Ngrids_wrf
        write(stdout,*)"User input WRF grid ", i,"=", wrf_grids(i)
      end do 
      write(stdout,*)"Cartesian input In meters=", cartesian(1)
      write(stdout,*)"Cartesian input In meters=", cartesian(2)
      write(stdout,*)"================================================"
!     end do 


      end subroutine read_inputs 
 

