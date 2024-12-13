FILE(GLOB CPP *.cpp)

add_library(aggcol ${CPP})
target_link_libraries(aggcol PUBLIC ${EXTERNAL_LIBS})
install(TARGETS aggcol LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
