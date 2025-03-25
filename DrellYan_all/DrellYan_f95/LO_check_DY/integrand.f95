function flo2_LO(yy, vwgt)
  implicit none
  use parameters
  double precision :: flo2_LO, vwgt, rs, xa, xb, rsp, xinvmass, scale, xmuf, xmur, xmu2, AL, sig, xnorm, wgt
  double precision, parameter :: pi = 3.14159265358979d0
  double precision, parameter :: hbarc2 = 0.3894d12 ! in fb
  double precision, dimension(10) :: yy
  double precision, dimension(-6:6) :: f1, f2
  double precision, dimension(15) :: xl
  double precision, dimension(0:3) :: p1, p2, p3, p4, q, xp1, xp2
  double precision, dimension(0:3,1:4) :: p
  double precision, dimension(1:2) :: c, Born
  double precision, dimension(2,-2:0) :: coef
  double precision, dimension(-2:0) :: SumI
  character(len=50) :: name
  integer :: leg, ipass
  double precision :: eps, xlow, xhigh, xcut
  common /pdfname/ name
  common /leg_choice/ leg
  common /energy/ s
  common /distribution/ xq
  common /renor_scale/ scale
  common /usedalpha/ AL, ge
  external :: Born_uU2eE

  rs = sqrt(s)
  xa = yy(1)
  xb = yy(2)

  rsp = sqrt(xa * xb * s)

  ipass = 0
  eps = 1.0d0
  xlow = xq - eps
  xhigh = xq + eps

  xcut = xq - 10.0d0

  call kinvar2(yy, xinvmass, p1, p2, p3, p4)
  scale = xinvmass

  if (scale >= xlow .and. scale <= xhigh) then
    xmuf = scale
    xmur = scale
    xmu2 = xmuf**2

    call pdf(xa, xmuf, f1)
    call pdf(xb, xmuf, f2)
    call setlum(f1, f2, xl)
    AL = alphasPDF(xmur)

    sig = xl(1) * Born_uU2eE(0, p1, p2, p3, p4)

    xnorm = hbarc2 / (16d0 * pi * xa * xb * s)
    wgt = xnorm * sig * vwgt
    flo2_LO = wgt / vwgt / 2d0 / eps
  else
    flo2_LO = 0d0
  end if

  return
end function flo2_LO

