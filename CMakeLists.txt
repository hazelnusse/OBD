# Option to build shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

# Option for building all of OBD's tests
option(obd_build_tests "Build all of OBD's tests." ON)

# Option for building the bicycle GUI
option(obd_build_gui "Build Qt4 gui tool." ON)

# Option for building documentation using Doxygen
option(obd_build_docs "Build Doxygen documentation." OFF)

# We need to call enable_testing() before any subdirectories are added
if (obd_build_tests)
    enable_testing()
endif()

# Add custom FindXXX.cmake files cmake/Modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(OBD CXX C)
cmake_minimum_required(VERSION 2.8)

#set(OBD_VERSION_MAJOR 0)
#set(OBD_VERSION_MINOR 3)

# Include custom cmake functions
#include(cmake/internal_utils.cmake)
#config_compiler_and_linker()  # Defined in internal_utils.cmake.

#configure_file("${PROJECT_SOURCE_DIR}/OBDConfig.h.in"
#               "${PROJECT_BINARY_DIR}/OBDConfig.h")
#
#include_directories("${PROJECT_BINARY_DIR}")

# Find GSL
find_package(GSL REQUIRED)

# Directories for each bicycle model
add_subdirectory(whipple)

if (obd_build_docs)
    find_package(Doxygen)
endif()

if (obd_build_gui)
    find_package(Qt4 REQUIRED)
    include(${QT_USE_FILE})
    add_subdirectory(bikegui)
endif()

if (obd_build_tests)
    add_subdirectory(gtest)
endif()

#include(CPack)
