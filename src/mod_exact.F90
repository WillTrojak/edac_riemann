module exact
   use precision
   implicit none

   private

   public :: riemann_flux


contains


   function riemann_flux(ql, qr, n) result(f)
      use constants, only : ed_zeta
      use flux, only : edac_flux_1d
      use transform, only : transform_from, transform_to
      implicit none

      
      real(kind=fptype), intent(in) :: ql(:), qr(:), n(:)

      real(kind=fptype) :: f(size(ql))

      real(kind=fptype) :: ft(size(ql))
      real(kind=fptype) :: qtl(size(ql)), qtr(size(ql))
      
      
      qtl = transform_to(n, ql, 1)
      qtr = transform_to(n, qr, 1)
      
   end function riemann_flux


   subroutine newton_parts(ql, qr, us, ps, dpl, dpr, gdpl, gdpr)
      use constants, only : ed_zeta
      implicit none

      real(kind=fptype), intent(in) :: ql(:), qr(:), us, ps

      real(kind=fptype), intent(out) :: dpl, dpr, gdpl, gdpr

      real(kind=fptype) :: dl, dr, ds
      
      dl = sqrt(0.25*ql(2)*ql(2) + ed_zeta + ql(1))
      dr = sqrt(0.25*qr(2)*qr(2) + ed_zeta + qr(1))
      ds = sqrt(0.25*us*us + ed_zeta + ps)

      if(1.5*ql(2) - dl .le. 1.5*us - ds) then
         ! Rarefaction
         dpl = 0.25(ql(2)*ql(2) - us*us)
         gdpl = -0.5*us
      else
         ! Shock
         dpl = (us - ql(2))*(-0.5*ql(2) - dl)
         gdpl = -(0.5*ul + dl)
      endif

      if(1.5*us + ds .le. 1.5*qr(2) + dr) then
         ! Rarefaction
         dpr = 0.75*(qr(2)*qr(2) - us*us)
         gdpr = -1.5*us
      else
         ! Shock
         dpr = (qu(2) - us)*(-0.5*qr(2) + dr
         gdpr = 0.5*qr(2) - dr
      endif         
      
   end subroutine newton_parts
   
   
end module exact
