cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(OBD)
set(OBD_VERSION_MAJOR 0)
set(OBD_VERSION_MINOR 3)
configure_file("${PROJECT_SOURCE_DIR}/OBDConfig.h.in"
               "${PROJECT_BINARY_DIR}/OBDConfig.h")

include_directories("${PROJECT_BINARY_DIR}")


find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

find_package(Doxygen)

add_subdirectory(bikegui)
add_subdirectory(whipple)

include(CPack)

enable_testing()
