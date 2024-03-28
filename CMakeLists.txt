cmake_minimum_required(VERSION "3.14")

project(FORWARDDIF LANGUAGES Fortran)

add_library(forwarddif forwarddif.f90)

add_executable(test_forwarddif test_forwarddif.f90)
target_link_libraries(test_forwarddif forwarddif)
