# Create a library called "ambient" which includes the source files. 
# The extension is already found. Any number of sources could be listed here.
add_definitions(-g -O0)
include_directories (${DMRG_SOURCE_DIR}) 
add_library (ambient 
             auxiliary.cpp
             packets/packet.cpp
             packets/packet_t.cpp
             groups/packet_manager.cpp
             groups/multirank.cpp
             groups/group.cpp
             core/p_profile.cpp
             core/workgroup.cpp
             core/select.cpp
             core/smp.cpp
             core/operation.cpp
             core/layout.cpp
             ambient.cpp)
