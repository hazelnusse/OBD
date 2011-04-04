# We need to be able to find the gtest header files
include_directories(${PROJECT_SOURCE_DIR}/gtest)
include_directories(${PROJECT_SOURCE_DIR}/gtest/include)

# For now, lets just build static libraries and disable thread support
add_definitions(-DGTEST_CREATE_SHARED_LIBRARY=0 -DGTEST_HAS_PTHREAD=0)

# Add the two Google Test libraries, built statically
add_library(gtest STATIC src/gtest-all.cc)
add_library(gtest_main STATIC src/gtest_main.cc)

# Link gtest_main to gtest
target_link_libraries(gtest_main gtest)
