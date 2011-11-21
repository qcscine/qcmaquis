
add_library(mtm STATIC mtm_scheduler.cpp)

# Install / export
list(APPEND MAQUIS_TYPES_LIBRARIES mtm)

export(TARGETS mtm APPEND FILE ${CMAKE_BINARY_DIR}/MaquisTypesTargets.cmake)
install(TARGETS mtm COMPONENT libraries
        DESTINATION lib
        EXPORT MaquisTypesTargets)

