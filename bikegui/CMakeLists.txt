include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})
set(obdgui_sources main.cpp mainwindow.cpp)
set(obdgui_headers mainwindow.h)

qt4_wrap_cpp(obdgui_headers_moc ${obdgui_headers})

add_executable(obdgui
               ${obdgui_sources}
               ${obdgui_headers_moc})
target_link_libraries(obdgui ${QT_LIBRARIES})
