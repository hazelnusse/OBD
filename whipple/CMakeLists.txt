# Make sure we can find the header files
include_directories("${PROJECT_SOURCE_DIR}/whipple/include")

add_subdirectory(src)

# Add tests subdirectory if tests are enabled
#if (obd_build_tests)
#    add_subdirectory(tests)
#endif()
