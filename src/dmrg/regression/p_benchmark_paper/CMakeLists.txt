add_executable (p_BenchOne.paper BenchOne.cpp) 
add_executable (p_BenchTwo.paper BenchTwo.cpp) 

target_link_libraries (p_BenchOne.paper ambient acml_mp gfortran -fopenmp)
target_link_libraries (p_BenchTwo.paper ambient acml_mp gfortran -fopenmp)
