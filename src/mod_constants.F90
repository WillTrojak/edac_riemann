module constants
   use precision
   implicit none

   private 

   real(kind=fptype), public, save :: ed_zeta

   public :: init_constants

   
contains


   subroutine init_constants(zeta)
      implicit none

      real(kind=fptype), intent(in) :: zeta
      
      ed_zeta = zeta
      
   end subroutine init_constants

    
end module constants
