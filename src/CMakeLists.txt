configure_file(mudirac.in.hpp ${CMAKE_CURRENT_SOURCE_DIR}/mudirac.hpp)
add_link_options(-lgsl -lgslcblas -lm)
add_executable(mudirac mudirac.cpp)
target_link_libraries(mudirac mudiraclib)

