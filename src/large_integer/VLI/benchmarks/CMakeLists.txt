
add_definitions (-DBOOST_TEST_DYN_LINK -Wno-unused-variable)

set(size_vli 2 3 4)
set(order_poly 11 12 13 14 15)

foreach(order ${order_poly})
    foreach(size ${size_vli})
        add_executable(OrderPoly_${order}_SizeVli_${size}_benchmark benchmark.cpp )
        SET_TARGET_PROPERTIES(OrderPoly_${order}_SizeVli_${size}_benchmark PROPERTIES COMPILE_FLAGS "-DVLI_SIZE=${size} -DORDER_POLY=${order}")
        target_link_libraries(OrderPoly_${order}_SizeVli_${size}_benchmark ${Boost_LIBRARIES} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${OpenMP_CXX_FLAGS})
        add_test(OrderPoly_${order}_SizeVli_${size}_benchmark OrderPoly_${order}_SizeVli_${size}_benchmark)
    endforeach(size)
endforeach(order)
