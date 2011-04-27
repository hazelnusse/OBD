add_library(gyrostat GyroStat.cpp)
add_library(whipple whipple.cpp)
target_link_libraries(whipple gyrostat)
