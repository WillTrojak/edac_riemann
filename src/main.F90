program main
   use constants
   use lambda
   use precision
   use rusanov, only : rusanov_flux => riemann_flux
   use exact, only : exact_flux => riemann_flux
   implicit none
   
   procedure(wave_speed), pointer :: lam_func

   real(kind=fptype) :: ql(2), qr(2), n(1)
   real(kind=fptype) :: f(2), fe(2)

   integer(kind=itype) :: i
   
   lam_func => davis

   ql = [2., 1.]
   qr = [1., 1.]
   n = [1.]

   call init_constants(10._fptype)

   f = rusanov_flux(lam_func, ql, qr, n)
   fe = exact_flux(ql, qr, n)

   do i=1,size(ql)
      print *, ql(i), qr(i), f(i), fe(i)
   enddo
   
end program main
