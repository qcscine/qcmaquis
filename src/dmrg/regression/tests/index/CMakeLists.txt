
add_executable(index index.cpp)
target_link_libraries(index ${MAQUIS_Boost_LIBRARIES})

add_test(index index)
