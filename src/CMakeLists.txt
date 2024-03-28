add_library(forwarddif 
  forwarddif_const.f90
  forwarddif_dual.f90
  forwarddif_derivative.f90
  forwarddif.f90
)

# add_executable(test_forwarddif test_forwarddif.f90)
# target_link_libraries(test_forwarddif forwarddif)

# add_library(dnad dnad.F90)
# add_executable(test_forwarddif test_forwarddif.f90)
# target_link_libraries(test_forwarddif dnad)