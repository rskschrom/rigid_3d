add_executable(two_collide two_collide.cpp)
target_link_libraries(two_collide aggcol)
install(TARGETS two_collide RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/../runs)
