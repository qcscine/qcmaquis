# Make sure the linker can find the vli library once it is built. 
link_directories (${PROJECT_BINARY_DIR}/src/vli) 
#set(order_list 1 2 3 4 5 6 7 8 9 10 11 12 13 14)
set(order_list 6)

if(VLI_COMPILE_GPU)
    add_definitions(-DVLI_USE_GPU -DVLI_SPLIT_PARAM=${VLI_SPLIT_PARAM} )
    set(LIBRAIRIES vli accelerator_gpu ${CUDA_CUDART_LIBRARY} ${CUDA_LIBRARIES} ${OpenMP_CXX_FLAGS} ${GMP_LIBRARY_CXX})
else(VLI_COMPILE_GPU)
    set(LIBRAIRIES vli ${OpenMP_CXX_FLAGS} ${GMP_LIBRARY_CXX})
endif(VLI_COMPILE_GPU)

foreach(order ${order_list})
        add_executable(main_order_${order} main.cpp)
        SET_TARGET_PROPERTIES(main_order_${order} PROPERTIES COMPILE_FLAGS "-D__ORDER__=${order} -w")
        target_link_libraries(main_order_${order} ${LIBRAIRIES})
endforeach(order)

#the data programm
#add_executable(data data.cpp)
#SET_TARGET_PROPERTIES(data PROPERTIES COMPILE_FLAGS "-D__ORDER__=14 -w -std=gnu++0x ")
#target_link_libraries(data)
