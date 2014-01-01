! file: spectrum.f90
! author: Bozhin Karaivanov
! date: 28.07.2013
! last modified: 31.12.2013
!
! This module contains subroutines for calculation
! of power spectrum, described in arXiv:1307.4938
! Eqs (3) and (4)

module spectrum

   public :: fourier, spec
   double precision, parameter, private :: pi = 4.d0*datan(1.d0)
   double precision, parameter, private :: twopi = 8.d0*datan(1.d0)
   public :: spectrum_init, fourier1, fourier2, spectrum_finalize
   double complex, allocatable, dimension(:,:), private :: zexponents
   logical, private :: is_initialized

contains

   subroutine fourier (T, N, st, sf)
   ! This subroutine performs the discrete Fourier transformation.
   ! See Eq. (3)
   !
   ! Dummy variables:
   ! T = time length: t = 0, ..., T-1, f = 0, ..., T-1
   ! N = Lattice size
   ! st = input vector of size T: time series of states s_x(t)
   !      of site x for moment t
   ! sf = output vector of size T; the Fourier-tranformed image of st
      implicit none
      integer, intent (in) :: T, N
      integer, dimension(T, N), intent (in) :: st
      double complex, dimension(T, N), intent(out) :: sf
      double complex :: zT
      double precision :: dT, tpdt
      double complex :: phase, zxp
      integer :: f, i, x
      
      dT = dfloat(T)
      zT = (1.d0, 0.d0) / dcmplx(dT, 0.d0)
      tpdt = twopi / dT
      ! Initialization
      do x = 1, N
         do f = 1, T
            sf(f,x) = (0.d0, 0.d0)
         enddo
      enddo
      ! Calculation
      do i = 1, T
         do f = 1, T
            phase = dcmplx(0.d0, -tpdt*dfloat((f-1)*(i-1)))
            zxp = cdexp (phase) * zT ! zT here means normaization
            do x = 1, N
               sf(f, x) = sf(f, x) + dcmplx(st(i, x),0) * zxp
            enddo
         enddo
      enddo
   end subroutine fourier
   
   subroutine spec(T, N, sf, s)
   ! This subroutine calculates the spectrum distribution.
   ! See Eq. (4)
      implicit none
      integer, intent (in) :: T
      integer, intent (in) :: N
      double complex, dimension(T, N), intent(in) :: sf
      double precision, dimension(T), intent(out) :: s
      !double precision :: amp
      integer :: f, x
      
      do f = 1, T
         s(f) = 0.d0
         do x = 1, N
            s(f) = s(f) + dble(dconjg(sf(f,x))*sf(f,x))
         enddo
      enddo
   end subroutine spec
   
   subroutine spectrum_init(T)
      implicit none
      integer, intent (in) :: T
      double complex :: zT
      double precision :: dT, tpdt
      double complex :: phase
      integer :: f, i
      is_initialized = .false.
      if (T.lt.1) return
      allocate (zexponents(T,T))
      dT = dfloat(T)
      zT = (1.d0, 0.d0) / dcmplx(dT, 0.d0)
      tpdt = twopi / dT
      ! Calculation
      do i = 0, T-1
         do f = 0, T-1
            phase = dcmplx(0.d0, -tpdt*dfloat(f * i))
            zexponents(i+1,f+1) = cdexp (phase) * zT ! zT here means normaization
         enddo
      enddo
      is_initialized = .true.
   end subroutine spectrum_init
   
   subroutine spectrum_finalize()
      deallocate (zexponents)
   end subroutine spectrum_finalize
   
   subroutine fourier1(T, N, st, sf)
      implicit none
      integer, intent (in) :: T, N
      integer, dimension(T, N), intent (in) :: st
      double complex, dimension(T, N), intent(out) :: sf
      double complex :: zst
      integer :: f, i, x
      ! Initialization
      do x = 1, N
         do f = 1, T
            sf(f,x) = (0.d0, 0.d0)
         enddo
      enddo
      if (.not. is_initialized) return
      ! Calculation
      do i = 1, T
         do x = 1, N
            zst = dcmplx(st(i, x),0)
            do f = 1, T
               sf(f, x) = sf(f, x) + zst * zexponents(i, f)
            enddo
         enddo
      enddo
   end subroutine fourier1
   
   subroutine fourier2(T, N, st, sf)
      ! This subroutine must be used only if st(i,j) takes values 0 or 1
      implicit none
      integer, intent (in) :: T, N
      integer, dimension(T, N), intent (in) :: st
      double complex, dimension(T, N), intent(out) :: sf
      integer :: f, i, x
      ! Initialization
      do x = 1, N
         do f = 1, T
            sf(f,x) = (0.d0, 0.d0)
         enddo
      enddo
      if (.not. is_initialized) return
      ! Calculation
      do i = 1, T
         do x = 1, N
            if (st(i,x).gt.0) then
               do f = 1, T
                  sf(f, x) = sf(f, x) + zexponents(i, f)
               enddo
            endif
         enddo
      enddo
   end subroutine fourier2

end module spectrum

