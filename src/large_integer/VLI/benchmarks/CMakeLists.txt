
if(VLI_COMPILE_GPU)
    include_directories(${CUDA_INCLUDE_DIRS})
    add_definitions(-DVLI_USE_GPU -DVLI_SPLIT_PARAM=${VLI_SPLIT_PARAM} )
endif(VLI_COMPILE_GPU)

add_definitions (-DBOOST_TEST_DYN_LINK -Wno-unused-variable)

set(size_vli 2 3 4)
set(order_poly   2 3 4 5)
#set(order_poly 10)
set(variable_poly 1 2 3 4 )
set(benchmarks benchmarks_max_order_each benchmarks_max_order_combined)

foreach(bench ${benchmarks})
    foreach(order ${order_poly})
        foreach(size ${size_vli})
            foreach(variable ${variable_poly})
              add_executable(${bench}_OrderPoly_${order}_SizeVli_${size}_NumVariable_${variable}_benchmark ${bench}.cpp )
              SET_TARGET_PROPERTIES(${bench}_OrderPoly_${order}_SizeVli_${size}_NumVariable_${variable}_benchmark PROPERTIES COMPILE_FLAGS "-DVLI_SIZE=${size} -DORDER_POLY=${order} -DVARIABLE_POLY=${variable}")
   
              if(VLI_COMPILE_GPU)
                       target_link_libraries(${bench}_OrderPoly_${order}_SizeVli_${size}_NumVariable_${variable}_benchmark ${Boost_LIBRARIES} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${PROJECT_BINARY_DIR}/src/vli/libvli_gpu.a ${CUDA_CUDART_LIBRARY} ${CUDA_LIBRARIES} ${OpenMP_CXX_FLAGS})
              endif(VLI_COMPILE_GPU)
   
              target_link_libraries(${bench}_OrderPoly_${order}_SizeVli_${size}_NumVariable_${variable}_benchmark ${Boost_LIBRARIES} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${OpenMP_CXX_FLAGS})
              add_test(${bench}_OrderPoly_${order}_SizeVli_${size}_NumVariable_${variable}_benchmark ${bench}_OrderPoly_${order}_SizeVli_${size}_NumVariable_${variable}_benchmark)
            endforeach(variable)
        endforeach(size)
    endforeach(order)
endforeach(bench)

