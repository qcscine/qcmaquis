# add_definitions(-DMPICH_IGNORE_CXX_SEEK)
# Note: regenerate the preprocessor includes using the line below
# g++ ambient/models/operation/pp/operation.pp.hpp.template.h -P -E -I/apps/eiger/boost_1_46_1/include/ -I. > ambient/models/operation/pp/operation.pp.hpp

add_library (ambient
             models/v_model.cpp
             models/object.cpp
             models/revision.cpp
             models/reduction.cpp
             models/layout.cpp
             models/modifier.cpp
             models/operation/operation.cpp
             models/operation/select.cpp
             controllers/v_controller.cpp
             controllers/context.cpp
             channels/mpi_channel.cpp
             channels/packets/packet_t.cpp
             channels/packets/packet.cpp
             channels/groups/group.cpp
             channels/groups/multirank.cpp
             utils/hashmap.cpp
             utils/tasklist.cpp)

 install(TARGETS ambient COMPONENT libraries
                         DESTINATION lib
                         EXPORT MaquisAmbientTargets)
 export(TARGETS ambient APPEND FILE ${CMAKE_BINARY_DIR}/MaquisAmbientTargets.cmake)

