# Option to build shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

# Option for building all of OBD's tests
option(obd_build_tests "Build all of OBD's tests." ON)

# Option for building the bicycle GUI
option(obd_build_gui "Build Qt4 gui tool." ON)

# Option for building documentation using Doxygen
option(obd_build_docs "Build Doxygen documentation." OFF)

# Option for building code for deriving equations of motion
option(obd_build_derivation "Build EOM derivation (requires GiNaC)" ON)

# We need to call enable_testing() before any subdirectories are added
if (obd_build_tests)
    enable_testing()
endif()

# Add custom FindXXX.cmake files cmake/Modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(OBD CXX C)
cmake_minimum_required(VERSION 2.8)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Set g++ compiler flags if GNU C++ compiler is being used
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

set(OBD_VERSION_MAJOR 0)
set(OBD_VERSION_MINOR 3)

# If we find git installed, set OBD_VERSION_COMMIT to most recent SHA1
find_package(Git)
if(GIT_FOUND)
    # TODO:  check that ${PROJECT_SOURCE_DIR} has a .git folder to ensure it is
    # a git repo rather than just a source tar-ball.
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=oneline
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                    OUTPUT_VARIABLE OBD_VERSION_COMMIT)

    # Keep everything up to the first space, add quotes
    string(REGEX REPLACE "([^ ]+).*"
           "\"\\1\"" OBD_VERSION_COMMIT
           ${OBD_VERSION_COMMIT})
else()
    set(OBD_VERSION_COMMIT "\"\"")
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/OBDConfig.h.in"
               "${PROJECT_BINARY_DIR}/include/OBDConfig.h")

include_directories("${PROJECT_BINARY_DIR}/include")

# Find GiNaC
if (obd_build_derivation)
    find_package(Ginac REQUIRED)
    find_package(Cln REQUIRED)
endif()

# Find GSL
find_package(GSL REQUIRED)

# Directories for each bicycle model
add_subdirectory(whipple)

if (obd_build_docs)
    find_package(Doxygen)
    if (DOXYGEN_FOUND STREQUAL "NO")
        message(FATAL_ERROR "Doxygen not found. Please get a copy http://www.doxygen.org")
    endif()
    # prepare doxygen configuration file
    configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in
                   ${PROJECT_BINARY_DIR}/Doxyfile)
    add_custom_target(doc
                      ${DOXYGEN_EXECUTABLE}
                      ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

if (obd_build_gui)
    find_package(Qt4 REQUIRED)
    add_subdirectory(bikegui)
endif()

if (obd_build_tests)
    add_subdirectory(gtest)
endif()
