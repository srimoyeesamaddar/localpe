!!!
!!! Authors:  Douglas P. Drob and John Emmert, NRL 7632
!!!
!!! =============================================
!!! Altitude to Geopotential Height
!!! =============================================

real(8) function hgt2gph(theta,z)

  implicit none
  real(8),intent(in)  :: theta
  real(8),intent(in)  :: z

  real(8)  :: xa,xb,xc
  real(8)  :: x0,x1,x2,x3,x4

  real(8),parameter :: epsqr = 6.69437999014e-3
  real(8),parameter :: ra = 6378137.0
  real(8),parameter :: rb = 6335439.32729283
  real(8),parameter :: rc = 6393052.93059087

  real(8),parameter :: c1 = 2.71365980406076e-10
  real(8),parameter :: c2 = 40683298295332.0
  real(8),parameter :: c3 = 5.37532939282856e24
  real(8),parameter :: c4 = 66063097365.0104
  real(8),parameter :: c5 = 22021032455.0035
  real(8),parameter    :: deg2rad = 0.017453292519943295d0                                                                        
  xa = sin(theta*deg2rad)**2                                                                                                     
  xb = cos(theta*deg2rad)**2 

!  xa = sind(theta)**2
!  xb = cosd(theta)**2

  x0 = 1.0d0/sqrt(-epsqr*xa + 1.0d0)
  x1 = xb * (ra * x0 + z) ** 2
  x2 = xa * (rb * x0 + z) ** 2
  x3 = x1 + x2
  x4 = 1.0d0/x3

  hgt2gph = -c1 * x1 + rc - (-c2*x4*(c4*x2*x4 - c5) + c2)/sqrt(x3)
  return

end function hgt2gph

!!! =============================================
!!! Geopotential Height to Altitude
!!! =============================================

real(8) function gph2hgt(theta,gph)

  implicit none
  real(8),intent(in)  :: theta
  real(8),intent(in)  :: gph

  real(8),external  :: hgt2gph
  real(8)           :: x,dx,y,dydz
  integer           :: n
  integer,parameter :: maxn = 10
  real(8),parameter :: epsilon = 0.0005

  x = gph
  n = 0
  dx = epsilon + epsilon
  do while ((abs(dx) .gt. epsilon) .and. (n .lt. 10))
     y = hgt2gph(theta,x)
     dydz = (hgt2gph(theta,x+dx) - y)/dx
     dx = (gph - y)/dydz
     x = x + dx
     n = n + 1
  end do

  gph2hgt = x

end function gph2hgt
