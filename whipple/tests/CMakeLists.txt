
add_executable(testeigs testeigs.cpp)
target_link_libraries(testeigs
    whipple
    ${GSL_LIBRARIES}
    ${GSLCBLAS_LIBRARIES})

add_test(EigenvalueTest testeigs)
