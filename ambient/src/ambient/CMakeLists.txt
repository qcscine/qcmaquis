# add_definitions(-DMPICH_IGNORE_CXX_SEEK)
# Note: regenerate the preprocessor includes using the line below
# g++ ambient/interface/pp/kernel_inliner.pp.hpp.template.h -P -E -I/apps/eiger/boost_1_46_1/include/ -I. > ambient/interface/pp/kernel_inliner.pp.hpp

 add_library (ambient ambient.cpp)

 install(TARGETS ambient COMPONENT libraries
                         DESTINATION lib
                         EXPORT MaquisAmbientTargets)
 export(TARGETS ambient APPEND FILE ${CMAKE_BINARY_DIR}/MaquisAmbientTargets.cmake)

