include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

include_directories(${PROJECT_SOURCE_DIR}/whipple/common)
set(obdgui_sources main.cpp mainwindow.cpp parameters.cpp)
set(obdgui_headers mainwindow.h parameters.h)

qt4_wrap_cpp(obdgui_headers_moc ${obdgui_headers})

add_executable(obdgui
               ${obdgui_sources}
               ${obdgui_headers_moc})
target_link_libraries(obdgui
               whipple
               ${GSL_LIBRARIES}
               ${GSLCBLAS_LIBRARIES}
               ${QT_LIBRARIES})
