module ref_sol_2d_hel
  use Iso_C_Binding
  implicit none
contains
  !---------------------------------------
  function bessel_3_fnu(zr, zi, fnu) result(res)
    implicit none

    real(c_double),intent(in) :: zr
    real(c_double),intent(in) :: zi
    real(c_double),intent(in) :: fnu

    complex(c_double_complex) :: res(3)
    real(8) :: cyr(3)
    real(8) :: cyi(3)
    real(8) :: cyr_tmp
    real(8) :: cyi_tmp
    integer, parameter :: n = 3
    integer, parameter :: kode = 1
    integer :: nz, ierr

    if(dble(floor(fnu)) .eq. fnu) then
       res = dcmplx(0d0, 0d0)

       if(fnu .ge. 1d0) then
          call ZBESJ(ZR, ZI, FNU - 1d0, KODE, N, CYR, CYI, NZ, IERR)
       else if(fnu .eq. 0d0) then
          call ZBESJ(ZR, ZI, 0d0, KODE, 2, CYR(2), CYI(2), NZ, IERR)
          cyr(1) = 0d0
          cyi(1) = 0d0
          cyr(3) = 2d0*cyr(3)
          cyi(3) = 2d0*cyi(3)
       else if(fnu .eq. -1d0) then
          call ZBESJ(ZR, ZI, 0d0, KODE, N, CYR, CYI, NZ, IERR)
          cyr_tmp = cyr(1)
          cyi_tmp = cyi(1)
          cyr(1) = cyr(3)
          cyi(1) = cyi(3)
          cyr(2) = -cyr(2)
          cyi(2) = -cyi(2)
          cyr(3) = cyr_tmp
          cyi(3) = cyi_tmp
       else if(fnu .le. -2d0) then
          call ZBESJ(ZR, ZI, -(FNU + 1d0), KODE, N, CYR, CYI, NZ, IERR)
          CYR(1) = CYR(1)*(-1d0)**(-fnu-1d0)
          CYI(1) = CYI(1)*(-1d0)**(-fnu-1d0)
          CYR(2) = CYR(2)*(-1d0)**(-fnu)
          CYI(2) = CYI(2)*(-1d0)**(-fnu)
          CYR(3) = CYR(3)*(-1d0)**(-fnu+1d0)
          CYI(3) = CYI(3)*(-1d0)**(-fnu+1d0)
          cyr_tmp = CYR(1)
          cyi_tmp = CYI(1)
          CYR(1) = CYR(3)
          CYI(1) = CYI(3)
          CYR(3) = cyr_tmp
          CYI(3) = cyi_tmp
       else
          error stop
       end if
       if(ierr .ne. 0) then
          write(*,*) 'in zbesj, ierr, nz', ierr, nz
          error stop
       end if
       res(:) = dcmplx(cyr(:), cyi(:))
    else
       write(*,*) 'in special function, fnu MUST be integer'
       error stop
    end if

  end function bessel_3_fnu
  !---------------------------------------
  function hankel_3_fnu(zr, zi, fnu) result(res)
    implicit none

    real(c_double),intent(in) :: zr
    real(c_double),intent(in) :: zi
    real(c_double),intent(in) :: fnu

    complex(c_double_complex) :: res(3)
    real(8) :: cyr(3)
    real(8) :: cyi(3)
    real(8) :: cyr_tmp
    real(8) :: cyi_tmp
    integer, parameter :: n = 3
    integer, parameter :: kode = 1
    integer, parameter :: m = 1
    integer :: nz, ierr

    if(dble(floor(fnu)) .eq. fnu) then
       res = dcmplx(0d0, 0d0)

       if(fnu .ge. 1d0) then
          call ZBESH(ZR, ZI, FNU - 1d0, KODE, M, N, CYR, CYI, NZ, IERR)
       else if(fnu .eq. 0d0) then
          call ZBESH(ZR, ZI, 0d0, KODE, M, 2, CYR(2), CYI(2), NZ, IERR)
          cyr(1) = 0d0
          cyi(1) = 0d0
          cyr(3) = 2d0*cyr(3)
          cyi(3) = 2d0*cyi(3)
       else if(fnu .eq. -1d0) then
          call ZBESH(ZR, ZI, 0d0, KODE, M, N, CYR, CYI, NZ, IERR)
          cyr_tmp = cyr(1)
          cyi_tmp = cyi(1)
          cyr(1) = cyr(3)
          cyi(1) = cyi(3)
          cyr(2) = -cyr(2)
          cyi(2) = -cyi(2)
          cyr(3) = cyr_tmp
          cyi(3) = cyi_tmp
       else if(fnu .le. -2d0) then
          call ZBESH(ZR, ZI, -(FNU + 1d0), KODE, M, N, CYR, CYI, NZ, IERR)
          CYR(1) = CYR(1)*(-1d0)**(-fnu-1d0)
          CYI(1) = CYI(1)*(-1d0)**(-fnu-1d0)
          CYR(2) = CYR(2)*(-1d0)**(-fnu)
          CYI(2) = CYI(2)*(-1d0)**(-fnu)
          CYR(3) = CYR(3)*(-1d0)**(-fnu+1d0)
          CYI(3) = CYI(3)*(-1d0)**(-fnu+1d0)
          cyr_tmp = CYR(1)
          cyi_tmp = CYI(1)
          CYR(1) = CYR(3)
          CYI(1) = CYI(3)
          CYR(3) = cyr_tmp
          CYI(3) = cyi_tmp
       else
          error stop
       end if
       if(ierr .ne. 0) then
          write(*,*) 'in zbesh, ierr, nz', ierr, nz
          error stop
       end if
       res(:) = dcmplx(cyr(:), cyi(:))
    else
       write(*,*) 'in special function, fnu MUST be integer'
       error stop
    end if

  end function hankel_3_fnu
  !---------------------------------------
  function bessel_inc_wave(k_1, xm1, xm2, rad, fnu) result(res) bind(c)
    implicit none

    complex(c_double_complex),intent(in) :: k_1
    real(c_double),intent(in) :: xm1
    real(c_double),intent(in) :: xm2
    real(c_double),intent(in) :: rad
    real(c_double),intent(in) :: fnu

    complex(c_double_complex) :: res
    real(8) :: zr, zi
    complex*16 :: bessel(3)
    complex*16,parameter :: iunit = dcmplx(0d0, 1d0)

    res = dcmplx(0d0, 0d0)

    zr = dble(k_1*rad)
    zi = dimag(k_1*rad)

    bessel(:) = bessel_3_fnu(zr, zi, fnu)

    res = bessel(2)*exp(iunit*fnu*atan2(xm2, xm1))

  end function bessel_inc_wave
  !---------------------------------------
  function make_sol_mie_series_diri(k_1, xm1, xm2, rad, fnu) result(res) bind(c)
    implicit none

    complex(c_double_complex),intent(in) :: k_1
    real(c_double),intent(in) :: xm1
    real(c_double),intent(in) :: xm2
    real(c_double),intent(in) :: rad
    real(c_double),intent(in) :: fnu

    complex(c_double_complex) :: res
    real(8) :: zr, zi
    complex*16 :: hankel(3)
    complex*16 :: bessel(3)
    complex*16,parameter :: iunit = dcmplx(0d0, 1d0)

    res = dcmplx(0d0, 0d0)

    zr = dble(k_1*rad)
    zi = dimag(k_1*rad)

    bessel(:) = bessel_3_fnu(zr, zi, fnu)
    hankel(:) = hankel_3_fnu(zr, zi, fnu)

    res = 0.5d0*k_1&
         & *(-bessel(2)/hankel(2)*(hankel(1) - hankel(3)) + bessel(1) - bessel(3))&
         & *exp(iunit*fnu*atan2(xm2, xm1))

  end function make_sol_mie_series_diri
  !---------------------------------------
  function make_sol_mie_series_neum(k_1, xm1, xm2, rad, fnu) result(res) bind(c)
    implicit none

    complex(c_double_complex),intent(in) :: k_1
    real(c_double),intent(in) :: xm1
    real(c_double),intent(in) :: xm2
    real(c_double),intent(in) :: rad
    real(c_double),intent(in) :: fnu

    complex(c_double_complex) :: res
    real(8) :: zr, zi
    complex*16 :: hankel(3)
    complex*16 :: bessel(3)
    complex*16,parameter :: iunit = dcmplx(0d0, 1d0)

    res = dcmplx(0d0, 0d0)

    zr = dble(k_1*rad)
    zi = dimag(k_1*rad)

    bessel(:) = bessel_3_fnu(zr, zi, fnu)
    hankel(:) = hankel_3_fnu(zr, zi, fnu)

    res = (-(bessel(1) - bessel(3))/(hankel(1) - hankel(3))*hankel(2) + bessel(2))&
         & *exp(iunit*fnu*atan2(xm2, xm1))

  end function make_sol_mie_series_neum
  !---------------------------------------
  subroutine make_sol_mie_series_transmission(resu, resq, k_1, k_2, dielec_1, dielec_2, xm1, xm2, rad, fnu) bind(c)
    implicit none

    complex(c_double_complex),intent(out) :: resu
    complex(c_double_complex),intent(out) :: resq
    complex(c_double_complex),intent(in) :: k_1
    complex(c_double_complex),intent(in) :: k_2
    complex(c_double_complex),intent(in) :: dielec_1
    complex(c_double_complex),intent(in) :: dielec_2
    real(c_double),intent(in) :: xm1
    real(c_double),intent(in) :: xm2
    real(c_double),intent(in) :: rad
    real(c_double),intent(in) :: fnu

    real(8) :: zr, zi
    complex*16 :: hankel(3)
    complex*16 :: bessel(3)
    complex*16 :: bessel_inner(3)
    complex*16 :: bb
    complex*16,parameter :: iunit = dcmplx(0d0, 1d0)

    resu = dcmplx(0d0, 0d0)
    resq = dcmplx(0d0, 0d0)

    zr = dble(k_1*rad)
    zi = dimag(k_1*rad)
    bessel(:) = bessel_3_fnu(zr, zi, fnu)
    hankel(:) = hankel_3_fnu(zr, zi, fnu)

    zr = dble(k_2*rad)
    zi = dimag(k_2*rad)
    bessel_inner(:) = bessel_3_fnu(zr, zi, fnu)

    bb = k_1/dielec_1*(bessel(2)/hankel(2)*(hankel(1) - hankel(3)) - (bessel(1) - bessel(3)))/(bessel_inner(2)/hankel(2)*k_1/dielec_1*(hankel(1) - hankel(3)) - k_2/dielec_2*(bessel_inner(1) - bessel_inner(3)))

    resu = bb*bessel_inner(2)*exp(iunit*fnu*atan2(xm2, xm1))
    resq = bb*0.5d0*k_2/dielec_2*(bessel_inner(1) - bessel_inner(3))*exp(iunit*fnu*atan2(xm2, xm1))

  end subroutine make_sol_mie_series_transmission
  !---------------------------------------
end module ref_sol_2d_hel
