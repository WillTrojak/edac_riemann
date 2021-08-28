module rusanov
   use precision
   implicit none
   
   private
   
   public :: riemann_flux


contains


   function riemann_flux(func, ql, qr, n) result(f)
      use constants, only : ed_zeta
      use flux, only : edac_flux_1d
      use lambda
      use transform, only : transform_to, transform_from
      implicit none

      procedure(wave_speed), pointer :: func
      
      real(kind=fptype), intent(in) :: ql(:), qr(:), n(:)

      real(kind=fptype) :: f(size(ql))

      real(kind=fptype) :: lambda, ft(size(ql))
      real(kind=fptype) :: qtl(size(ql)), qtr(size(ql))
      real(kind=fptype) :: ftl(size(ql)), ftr(size(qr))

      qtl = transform_to(n, ql, 1)
      qtr = transform_to(n, qr, 1)

      lambda = func(qtl, qtr)

      ftl = edac_flux_1d(qtl)
      ftr = edac_flux_1d(qtr)
      
      ft = 0.5*(ftr + ftl) - 0.5*lambda*(qtr - qtl)

      f = transform_from(n, ft, 1)
      
   end function riemann_flux
     
     
end module rusanov
