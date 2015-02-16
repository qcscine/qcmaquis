add_library (ambient ambient.cpp)

install(TARGETS ambient COMPONENT libraries
                        DESTINATION lib
                        EXPORT AmbientTargets)
export(TARGETS ambient APPEND FILE ${CMAKE_BINARY_DIR}/AmbientTargets.cmake)
