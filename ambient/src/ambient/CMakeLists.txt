# add_definitions(-DMPICH_IGNORE_CXX_SEEK)

 add_library (ambient ambient.cpp)

 install(TARGETS ambient COMPONENT libraries
                         DESTINATION lib
                         EXPORT AmbientTargets)
 export(TARGETS ambient APPEND FILE ${CMAKE_BINARY_DIR}/AmbientTargets.cmake)

