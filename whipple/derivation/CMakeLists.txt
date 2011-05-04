add_executable(eoms eoms.cpp)
target_link_libraries(eoms ${GINAC_LIBS} ${CLN_LIBS})
