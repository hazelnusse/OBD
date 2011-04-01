add_subdirectory(common)

# Make sure we can find the header files
include_directories("${PROJECT_SOURCE_DIR}/whipple/common")

# Add tests subdirectory if tests are enabled
if (obd_build_tests)
    add_subdirectory(tests)
endif()

# Executables which demonstrate the class
add_executable(whippleeig whippleeig.cpp)
target_link_libraries(whippleeig
    whipple
    ${GSL_LIBRARIES}
    ${GSLCBLAS_LIBRARIES})

add_executable(whipplesim whipplesim.cpp)
target_link_libraries(whipplesim
    whipple
    ${GSL_LIBRARIES}
    ${GSLCBLAS_LIBRARIES})

add_executable(whipplesteady whipplesteady.cpp)
target_link_libraries(whipplesteady
    whipple
    ${GSL_LIBRARIES}
    ${GSLCBLAS_LIBRARIES})

add_executable(whipplesteadyreview whipplesteady_review.cpp)
target_link_libraries(whipplesteadyreview
    whipple
    ${GSL_LIBRARIES}
    ${GSLCBLAS_LIBRARIES})
