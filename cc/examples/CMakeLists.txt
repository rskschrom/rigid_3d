#add_executable(two_collide two_collide.cpp)
#target_link_libraries(two_collide aggcol)

add_executable(freefall_particle freefall_particle.cpp)
target_link_libraries(freefall_particle aggcol)

add_executable(create_particle create_particle.cpp)
target_link_libraries(create_particle aggcol)

#install(TARGETS two_collide RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/../runs)
install(TARGETS freefall_particle RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/../runs)
install(TARGETS create_particle RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/../runs)
