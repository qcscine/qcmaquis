
set(VLI_TEST_LIBRARIES ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} ${PROJECT_BINARY_DIR}/src/vli/libvli.a)
add_subdirectory (vli)
add_subdirectory (polynomial)
add_subdirectory (inner_product_gpu)
