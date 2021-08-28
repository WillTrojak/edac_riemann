module lambda
   use precision
   implicit none

   private

   public :: davis, wave_speed


   abstract interface
      function wave_speed(ql, qr) result(lambda)
         use precision
         real(kind=fptype), intent(in) :: ql(:), qr(:)
         real(kind=fptype) :: lambda
      end function wave_speed
   end interface


   ! All functions assume ql and qr are transformed to [1,0,...]^T
contains


   function davis(ql, qr) result(lambda)
      use constants, only : ed_zeta
      implicit none

      real(kind=fptype), intent(in) :: ql(:), qr(:)

      real(kind=fptype) :: lambda

      real(kind=fptype) sl, sr

      sl = 1.5*abs(ql(2)) + sqrt(0.25*ql(2)*ql(2) + ql(1) + ed_zeta)
      sr = 1.5*abs(qr(2)) + sqrt(0.25*qr(2)*qr(2) + qr(1) + ed_zeta)

      lambda = max(sl, sr)
      
   end function davis

   
end module lambda
