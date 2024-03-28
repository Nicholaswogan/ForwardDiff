add_library(forwarddiff
  forwarddiff_const.f90
  forwarddiff_dual.f90
  forwarddiff_derivative.f90
  forwarddiff.f90
)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(forwarddiff PRIVATE -Wunused -Wimplicit-interface -fimplicit-none)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(forwarddiff PRIVATE -fcheck=all,no-array-temps)
  endif()
endif()