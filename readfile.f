      program test 

        integer :: iunit,nx,ny
        integer, allocatable::xx(:,:),yy(:,:)
        iunit = 10 
        nx = 4 
        ny = 3
        allocate(xx(nx,ny),yy(nx,ny)) 

        open(unit=iunit,file='test_file',status='old')
        do j=1,ny
         do i=1,nx
          read(iunit,*)xx(i,j)!xx(i,j),yy(i,j) 
          write(6,*),xx(i,j)!,yy(i,j)
         end do
        end do 
        close(iunit)
!        n=0
!        do
!          read(iunit,*,end=1)
!          n=n+1
!        end do 
!        print*, "Number of line", n
!1       rewind(iunit)

!        allocate(zz(n))  
!        do i=1,n
!          read(iunit,*) bf
!          zz(i)=bf 
!        end do 
!        close(iunit)       


!       write(6,*),bath_file(mw)
!        read(5,*) bf 
!        close(iunit)       
!        zz=ones(size(bath{mw}))
!        zz(bath{mw}==9999)=0
!        mask_rho_w{mw}=zz
!        nx = Numx_swan(mw)
!        ny = Numy_swan(mw)

!        open(unit=iunit,file=swan_coord(mw),status='old')
!        write(6,*),swan_coord(mw)

      end program 

