# Create a library called "ambient" which includes the source files. 
# The extension is already found. Any number of sources could be listed here.
enable_mpi()
add_library (ambient 
             auxiliary.cpp
             packets/packet.cpp
             packets/packet_t.cpp
             groups/packet_manager.cpp
             groups/multirank.cpp
             groups/group.cpp
             core/p_profile.cpp
             core/memblock.cpp
             core/select.cpp
             core/scope_context.cpp
             core/workgroup_context.cpp
             core/operation/operation.cpp
             core/layout.cpp
             ambient.cpp)
