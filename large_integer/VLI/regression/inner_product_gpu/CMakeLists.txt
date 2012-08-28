# Make sure the linker can find the vli library once it is built. 
link_directories (${PROJECT_BINARY_DIR}/src/vli) 

if(VLI_COMPILE_GPU)
    include_directories(${CUDA_INCLUDE_DIRS})
    add_definitions(-DVLI_USE_GPU -DVLI_SPLIT_PARAM=${VLI_SPLIT_PARAM} )

add_definitions (-DBOOST_TEST_DYN_LINK -Wno-unused-variable)

set(inner_product inner_product_max_order_each inner_product_max_order_combined)

set(test_sizes 2 3 4)

foreach(test ${inner_product})
    foreach(size ${test_sizes})
        add_executable(${test}_${size}_test ${test}.cpp)
        SET_TARGET_PROPERTIES(${test}_${size}_test PROPERTIES COMPILE_FLAGS "-DVLI_SIZE=${size}")
        target_link_libraries(${test}_${size}_test ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} vli accelerator_gpu ${CUDA_CUDART_LIBRARY} ${CUDA_LIBRARIES} ${OpenMP_CXX_FLAGS})
        add_test(${test}_${size}_test ${test}_${size}_test)
    endforeach(size)
endforeach(test)

endif(VLI_COMPILE_GPU)
