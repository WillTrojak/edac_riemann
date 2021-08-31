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
      
      real(kind=fptype) :: ps, us, psl, psr
      real(kind=fptype) :: dpl, dpr, gdpl, gdpr

      integer(kind=itype) :: k
      integer(kind=itype), parameter :: k_max = 10

      qtl = transform_to(n, ql, 1)
      qtr = transform_to(n, qr, 1)
      
      ! Initialise solution
      ps = 0.5*(ql(1) + qr(1))
      us = 0.5*(ql(2) + qr(2))
      psl = ps; psr = ps

      ! Interate
      do k=1,k_max
         call newton_parts(qtl, qtr, us, psl, psr, dpl, dpr, gdpl, gdpr)
         us = us - (ql(1) - qr(1) + dpr + dpl)/(gdpr + gdpl)

         psl = dpl + ql(1)
         psr = qr(1) - dpr
         print *, us, dpl, dpr, gdpl, gdpr
      enddo
      ps = qr(1) - dpr
      ps = ql(1) + dpl

      print *,'Presure', qr(1) - dpr,ql(1) + dpl
      ft(1) = us*(ps + ed_zeta)
      ft(2) = us*us + ps

      f = transform_from(n, ft, 1)

   end function riemann_flux


   subroutine newton_parts(ql, qr, us, psl, psr, dpl, dpr, gdpl, gdpr)
      use constants, only : ed_zeta
      implicit none

      real(kind=fptype), intent(in) :: ql(:), qr(:), us, psl, psr

      real(kind=fptype), intent(out) :: dpl, dpr, gdpl, gdpr

      real(kind=fptype) :: dl, dr, dsr, dsl
      
      dl = sqrt(0.25*ql(2)*ql(2) + ed_zeta + ql(1))
      dr = sqrt(0.25*qr(2)*qr(2) + ed_zeta + qr(1))
      dsl = sqrt(0.25*ql(2)*ql(2) + ed_zeta + psl)
      dsr = sqrt(0.25*qr(2)*qr(2) + ed_zeta + psr)

      if(1.5*ql(2) - dl .le. 1.5*us - dsl) then
         ! Rarefaction
         print *,'lrare'
         dpl = 0.25*(ql(2)*ql(2) - us*us)
         gdpl = -0.5*us
      else
         ! Shock
         print *,'lshock'
         dpl = (us - ql(2))*(-0.5*ql(2) - dl)
         gdpl = -0.5*ql(2) - dl
      endif

      if(1.5*us + dsr .le. 1.5*qr(2) + dr) then
         ! Rarefaction
         print *,'rrare'
         dpr = -0.25*(qr(2)*qr(2) - us*us)
         gdpr = 0.5*us
      else
         ! Shock
         print *,'rshock'
         dpr = (qr(2) - us)*(-0.5*qr(2) + dr)
         gdpr = 0.5*qr(2) - dr
      endif
      
   end subroutine newton_parts
   
   
end module exact

! Left
! dP/du = (P + z)/(u/2 - d) = (d + u/2)(d - u/2)/(u/2 - d) = -d - u/2
! dP/du = - u/2 - sqrt(u*u/4 + P + z)
! P = - u*u/4 - z
! => dP/du = -u/2
! dP/ud = -u/2
!
! 
! P + z = d^2 - u^2/4 = (d - u/2)(d + u/2)
! -(P + z)/(d + u/2) = u/2 - d
!
! Right
! dP/du = (P + z)/(d + u/2) = (d - u/2)
! P = -z - u^2/4
