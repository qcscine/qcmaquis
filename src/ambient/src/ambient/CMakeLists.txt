#add_definitions(-DMPICH_IGNORE_CXX_SEEK)
add_library (ambient 
             auxiliary.cpp
             packets/packet.cpp
             packets/packet_t.cpp
             groups/packet_manager.cpp
             groups/multirank.cpp
             groups/group.cpp
             core/p_object.cpp
             core/memblock.cpp
             core/select.cpp
             core/scope_context.cpp
             core/workgroup_context.cpp
             core/operation/operation.cpp
             core/layout.cpp
             model.cpp
             utils/zout.cpp
             ambient.cpp)


 install(TARGETS ambient COMPONENT libraries
                         DESTINATION lib
                         EXPORT MaquisAmbientTargets)
 export(TARGETS ambient APPEND FILE ${CMAKE_BINARY_DIR}/MaquisAmbientTargets.cmake)

