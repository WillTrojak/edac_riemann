module flux
   use precision
   implicit none

   private

   public :: edac_flux_1d


contains

  
  function edac_flux_1d(q) result(f)
     use constants, only : ed_zeta
     implicit none

     real(kind=fptype), intent(in) :: q(:)

     real(kind=fptype) :: f(size(q))

     integer(kind=itype) :: i
     
     f(1) = q(2)*(q(1) + ed_zeta)
     do i=2,size(q)
        f(i) = q(2)*q(i)
     enddo
     f(2) = f(2) + q(1)
     
  end function edac_flux_1d

  
end module flux
