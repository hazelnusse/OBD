# We need to be able to find the gtest header files
include_directories(${PROJECT_SOURCE_DIR}/gtest)
include_directories(${PROJECT_SOURCE_DIR}/gtest/include)

# Eigenvalue test-gtest
#add_executable(whipple-eigs-unittest whipple-eigs-unittest.cpp)
#target_link_libraries(whipple-eigs-unittest
#    gtest_main
#    whipple
#    ${GSL_LIBRARIES}
#    ${GSLCBLAS_LIBRARIES})
#
## Add test to testing framework
#add_test(NAME whipple-eigs-unittest
#         COMMAND whipple-eigs-unittest)
