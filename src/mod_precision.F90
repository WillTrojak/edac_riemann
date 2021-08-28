module precision
   use, intrinsic :: iso_fortran_env
   implicit none
   
   private

   integer, parameter, public :: i8 = int8
   integer, parameter, public :: i16 = int16
   integer, parameter, public :: i32 = int32
   integer, parameter, public :: i64 = int64
   integer, parameter, public :: itype = i32

   integer, parameter, public :: fp32 = real32
   integer, parameter, public :: fp64 = real64
   !integer, parameter, public :: fp80 = real80
   integer, parameter, public :: fptype = fp64

   integer, parameter, public :: wl = 4

end module precision
