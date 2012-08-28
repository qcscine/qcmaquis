
# Make sure the linker can find the vli library once it is built. 
link_directories (${PROJECT_BINARY_DIR}/src/vli) 

set(VLI_TEST_LIBRARIES ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} vli)
add_subdirectory (vli)
add_subdirectory (polynomial)
add_subdirectory (inner_product_gpu)
