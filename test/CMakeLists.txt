add_executable(test_forwarddiff test_forwarddiff.f90)
target_link_libraries(test_forwarddiff forwarddiff)
target_include_directories(test_forwarddiff PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(test_sparse test_sparse.f90)
target_link_libraries(test_sparse forwarddiff)
target_include_directories(test_sparse PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

find_program(FYPP fypp REQUIRED)
include(${PROJECT_SOURCE_DIR}/cmake/fypp.cmake)

set(fppFiles fypp_example.fypp)
fypp_f90("" "${fppFiles}" outFiles)
add_executable(fypp_example ${outFiles})
target_link_libraries(fypp_example forwarddiff)
target_include_directories(fypp_example PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})