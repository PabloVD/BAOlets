! MODULE
module constants

implicit none

integer, parameter :: FIELDMAX = 3893
integer, parameter :: CENTERMAX = 1599
real(8), parameter :: radstep = 1.0
integer, parameter :: nsteps = 300
real(8), parameter :: pi = acos(-1.0)
real(8), parameter :: Rmin = 20, Rmax = 225
real(8), parameter :: smin = 3, smax = 60
integer, parameter :: Rbins = 125, sbins = 20
real(8), parameter :: bandwidth = 10.0 ! Epanechnikov kernel

end module constants

! MAIN PROGRAM
program baolet

use constants
implicit none

integer i, j, k, stat, kmin, kmax, FIELDMAX_new, CENTERMAX_new
real(8), dimension(FIELDMAX) :: x_main, y_main, z_main
real(8), dimension(CENTERMAX) :: x_center, y_center, z_center
real(8), dimension(CENTERMAX,FIELDMAX) :: dist
real(8), dimension(0:nsteps) :: rho, rvals
real(8) Rstep, sstep, rho_norm, kernel
real(8), dimension(Rbins) :: R
real(8), dimension(sbins) :: s
real(8), dimension(Rbins,sbins) :: B_Rs


open (unit=20, file="main_default.dat")
open (unit=30, file="lrg_default.dat")

do i=1, FIELDMAX

    read (20,*,iostat=stat) x_main(i), y_main(i), z_main(i)
        !if (stat /= 0) exit

enddo

do i=1, CENTERMAX

    read (30,*,iostat=stat) x_center(i), y_center(i), z_center(i)
        !if (stat /= 0) exit

enddo

close(20)
close(30)

do i=0, nsteps
    rho(i) = 0.0
    rvals(i) = i*radstep
enddo


! Loop over centers
do i=1, CENTERMAX
    !write(*,*) "Working in centre no.", i
    ! Loop over field catalogue
    do j=1, FIELDMAX
        ! Calculate distance
        dist(i,j)=sqrt( (x_center(i)-x_main(j))**2 + (y_center(i)-y_main(j))**2 + (z_center(i)-z_main(j))**2 )
        !if (dist(i,j) > (radstep*nsteps + bandwidth)) continue

        kmin = int((dist(i,j) - bandwidth)/radstep)
        if (kmin<0) then
            kmin=0
        endif

        kmax = int((dist(i,j) + bandwidth)/radstep)
        if (kmax>nsteps) then
            kmax=nsteps
        endif

        do k=kmin, kmax
            ! Epanechnikov smoothing kernel
            !if (abs(rvals(k)-dist(i,j))/bandwidth<=1) then
            kernel=(0.75/bandwidth)*(1.-abs(rvals(k)-dist(i,j))/bandwidth)**2
            rho(k) = rho(k) + kernel
            !endif
        enddo

    enddo
enddo

! Normalize and write the profile
open (unit=40, file="prof.dat")

! Normalization differs from Arnalte code in a factor of 2
rho_norm = 2.*pi*CENTERMAX*FIELDMAX/1e9

do k=1, nsteps
    rho(k)=rho(k)/(rho_norm*rvals(k)**2)
    write(40,*) rvals(k), rho(k)
enddo

close (40)

end program baolet
