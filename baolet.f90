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
integer, parameter :: Rbins = 205, sbins = 57

contains

    function spline(x)

        real(8) x, spline
        spline = (1./12.)*( abs(x-2.)**3. - (4.*abs(x-1.)**3.) + (6.*abs(x)**3.) - (4.*abs(x+1.)**3.) + abs(x+2.)**3. )
        return

    end function spline

    function wavelet(x,R,s)

        real(8) wavelet, x, R, s
        if ( x<(R - 2.*s) .or. x>(R + 2.*s) ) then
            wavelet = 0
        else
            wavelet = (1./(4*pi))*( 2*spline( 2*(x-R)/s ) - spline( (x-R)/s ) )
        endif
        return

    end function wavelet

    function baolet_norm(R,s)

        real(8) baolet_norm, R, s, sum, f
        integer i

        sum = 0

        do i=1, nsteps
        f=(wavelet(i*radstep,R,s)/(i*radstep))**2
        sum = sum + f*radstep
        enddo

        baolet_norm = 1/sqrt(4*pi*sum)

        return

    end function baolet_norm

    function baolet_norm2(R,s)

        real(8) baolet_norm2, R, s

        baolet_norm2 = (((32*s**2)*(15*R**3 + 95*s*R**2 + 90*R*s**2 - 64*s**3) &
        +15*((R - 2*s)**5)*log(R - 2*s) - 15*((R + 2*s)**5)*log(R + 2*s))/(90*s**6))**(-1/2)

        return

    end function baolet_norm2

end module constants

! MAIN PROGRAM
program baolet

use constants
implicit none

integer i, j, k
real(8), dimension(FIELDMAX) :: x_main, y_main, z_main
real(8), dimension(CENTERMAX) :: x_center, y_center, z_center
real(8), dimension(CENTERMAX,FIELDMAX) :: dist
real(8), dimension(0:nsteps) :: rho, rvals
real(8) Rstep, sstep, rho_norm
real(8), dimension(Rbins) :: R
real(8), dimension(sbins) :: s
real(8), dimension(Rbins,sbins) :: B_Rs


open (unit=20, file="main_default.dat")
open (unit=30, file="lrg_default.dat")

do i=1, FIELDMAX

    read (20,*) x_main(i), y_main(i), z_main(i)

enddo

do i=1, CENTERMAX

    read (30,*) x_center(i), y_center(i), z_center(i)

enddo

close(20)
close(30)


! Read the profile
open (unit=40, file="prof.dat")

do k=0, nsteps-1
    read(40,*) rvals(k), rho(k)
enddo

close (40)

Rstep = (Rmax - Rmin)/Rbins
sstep = (smax - smin)/sbins

do k=1, Rbins
    R(k)=Rmin + k*Rstep
enddo

do k=1, sbins
    s(k)=smin + k*sstep
enddo

open (unit=50,file="Brs.dat")

! Calculate B(R,s)
do i=1, Rbins
    write(*,*) "Working in R =", R(i)
    do j=1, sbins
        B_Rs(i,j) = 0
        do k=1, nsteps
            B_Rs(i,j) = B_Rs(i,j) + rho(k)*wavelet(rvals(k),R(i),s(j))*radstep
        enddo
        !B_Rs(i,j) = B_Rs(i,j) - 0.5*rho(1)*wavelet(rvals(0),R(i),s(j))*radstep
        !B_Rs(i,j) = B_Rs(i,j) - 0.5*rho(nsteps-1)*wavelet(rvals(nsteps-1),R(i),s(j))*radstep
        !write(51,*) R(i), s(j), baolet_norm2(R(i),s(j))
        if (s(j)>0.5*R(i)) then
            write(50,*) R(i), s(j), -100
        else
            write(50,*) R(i), s(j), 10*B_Rs(i,j)*baolet_norm2(R(i),s(j))
        endif

    enddo
    write(50,*)
enddo

close(50)

end program baolet
