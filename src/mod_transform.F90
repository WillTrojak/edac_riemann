module transform
   use precision
   implicit none

   private

   real(kind=fptype), parameter :: ttol = 0.99e0_fptype


   public :: transform_from, transform_to


contains


   function transform_to(n, u, offset) result(t)
      implicit none
      real(kind=fptype), intent(in) :: n(:), u(:)

      integer(kind=itype), value :: offset

      real(kind=fptype) :: t(size(u))

      select case(size(n))
      case(1)
         t = transform_to_1d(n, u, offset)
      case(2)
         t = transform_to_2d(n, u, offset)
      case(3)
         t = transform_to_3d(n, u, offset)
      end select

   end function transform_to


   function transform_to_1d(n, u, offset) result(t)
      implicit none

      real(kind=fptype), intent(in) :: n(1), u(:)

      real(kind=fptype) :: t(size(u))

      integer(kind=itype), value :: offset

      t = u

      t(offset + 1) = n(1)*u(offset + 1)

   end function transform_to_1d


   function transform_to_2d(n, u, offset) result(t)
      implicit none

      real(kind=fptype), intent(in) :: n(2), u(:)

      real(kind=fptype) :: t(size(u))

      integer(kind=itype), value :: offset

      t = u

      t(offset + 1) =  n(1)*u(offset + 1) + n(2)*u(offset + 2)
      t(offset + 2) = -n(2)*u(offset + 1) + n(1)*u(offset + 2)

   end function transform_to_2d


   function transform_to_3d(n, u, offset) result(t)
      implicit none

      real(kind=fptype), intent(in) :: n(3), u(:)

      real(kind=fptype) :: t(size(u))

      integer(kind=itype), value :: offset

      integer(kind=itype) :: os
      real(kind=fptype) :: h

      os = offset
      t = u

      if(abs(n(1)) .lt. ttol) then
         h = 1./(1 + n(1))

         t(os + 1) =  n(1)*u(os + 1) + n(2)*u(os + 2) + n(3)*u(os + 3);
         t(os + 2) = -n(2)*u(os + 1) + (n(1) + h*n(3)*n(3))*u(os + 2) - h*n(2)*n(3)*u(os + 3)
         t(os + 3) = -n(3)*u(os + 1) - h*n(2)*n(3)*u(os + 2) + (n(1) + h*n(2)*n(2))*u(os + 3)

      elseif(abs(n(2)) .lt. abs(n(3))) then
         h = 1./(1 - n(2))

         t(os + 1) = n(1)*u(os + 1) + n(2)*u(os + 2) + n(3)*u(os + 3)
         t(os + 2) =  (1 - h*n(1)*n(1))*u(os + 1) + n(1)*u(os + 2) - h*n(1)*n(3)*u(os + 3)
         t(os + 3) = -h*n(1)*n(3)*u(os + 1) + n(3)*u(os + 2) + (1 - h*n(3)*n(3))*u(os + 3)

      else
         h = 1./(1 - n(3))

         t(os + 1) = n(1)*u(os + 1) + n(2)*u(os + 2) + n(3)*u(os + 3)
         t(os + 2) = -h*n(1)*n(2)*u(os + 1) + (1 - h*n(2)*n(2))*u(os + 2) + n(2)*u(os + 3)
         t(os + 3) = (1 - h*n(1)*n(1))*u(os + 1) - h*n(1)*n(2)*u(os + 2) + n(1)*u(os + 3)
      endif

      return
   end function transform_to_3d


   function transform_from(n, t, offset) result(u)
      implicit none
      real(kind=fptype), intent(in) :: n(:), t(:)

      integer(kind=itype), value :: offset

      real(kind=fptype) :: u(size(t))

      select case(size(n))
      case(1)
         u = transform_from_1d(n, t, offset)
      case(2)
         u = transform_from_2d(n, t, offset)
      case(3)
         u = transform_from_3d(n, t, offset)
      end select

   end function transform_from


   function transform_from_1d(n, t, offset) result(u)
      implicit none

      real(kind=fptype), intent(in) :: n(1), t(:)

      integer(kind=itype), value :: offset

      real(kind=fptype) :: u(size(t))

      u = t

      u(offset + 1) = -n(1)*t(offset + 1)

   end function transform_from_1d


   function transform_from_2d(n, t, offset) result(u)
      implicit none

      real(kind=fptype), intent(in) :: n(2), t(:)

      integer(kind=itype), value :: offset

      real(kind=fptype) :: u(size(t))

      u = t

      u(offset + 1) = n(1)*t(offset + 1) - n(2)*t(offset + 2)
      u(offset + 2) = n(2)*t(offset + 1) + n(1)*t(offset + 2)

   end function transform_from_2d


   function transform_from_3d(n, t, offset) result(u)
      implicit none

      real(kind=fptype), intent(in) :: n(3), t(:)

      integer(kind=itype), value :: offset

      real(kind=fptype) :: u(size(t))

      integer(kind=itype) :: os
      real(kind=fptype) :: h

      os = offset

      u = t

      if(abs(n(1)) .lt. ttol) then
         h = 1./(1 + n(1))

         u(os + 1) =  n(1)*t(os + 1) - n(2)*t(os + 2) - n(3)*t(os + 3)
         u(os + 2) =  n(2)*t(os + 1) + (n(1) + h*n(3)*n(3))*t(os + 2) - h*n(2)*n(3)*t(os + 3)
         u(os + 3) =  n(3)*t(os + 1) - h*n(2)*n(3)*t(os + 2) + (n(1) + h*n(2)*n(2))*t(os + 3)
      elseif(abs(n(2)) .lt. abs(n(3))) then
         h = 1./(1 - n(2))

         u(os + 1) = n(1)*t(os + 1) +  (1 - h*n(1)*n(1))*t(os + 2) - h*n(1)*n(3)*t(os + 3)
         u(os + 2) = n(2)*t(os + 1) + n(1)*t(os + 2) + n(3)*t(os + 3)
         u(os + 3) = n(3)*t(os + 1) - h*n(1)*n(3)*t(os + 2) + (1 - h*n(3)*n(3))*t(os + 3)
      else
         h = 1./(1 - n(3))

         u(os + 1) = n(1)*t(os + 1) - h*n(1)*n(2)*t(os + 2) + (1 - h*n(1)*n(1))*t(os + 3)
         u(os + 2) = n(2)*t(os + 1) + (1 - h*n(2)*n(2))*t(os + 2) - h*n(1)*n(2)*t(os + 3)
         u(os + 3) = n(3)*t(os + 1) + n(2)*t(os + 2) + n(1)*t(os + 3)
      endif

   end function transform_from_3d


end module transform
