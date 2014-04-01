
if(USE_AMBIENT)
  compile_with_ambient()
  list(INSERT LIBRARIES 0 ambient ${MPI_CXX_LIBRARIES})
endif(USE_AMBIENT) 

#set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
add_definitions(-DHAVE_ALPS_HDF5)

set(tests block_matrix serialize) #thing to break the block_matrix tests by group, constructor, operators, free functions ...
set(types_of_matrices DENSE PDENSE) #the MT matrice ??

if(USE_AMBIENT)
    set(types_of_matrices DENSE PDENSE) #the MT matrice ??
else()
    set(types_of_matrices DENSE) #the MT matrice ??
endif(USE_AMBIENT)

foreach(test ${tests})
    foreach(type ${types_of_matrices})
        add_executable(${test}_${type}_test ${test}.cpp)
        target_link_libraries(${test}_${type}_test ${DMRG_LIBRARIES})
        SET_TARGET_PROPERTIES(${test}_${type}_test PROPERTIES COMPILE_FLAGS "-D${type}") 
        if(USE_AMBIENT)
            target_link_libraries(${test}_${type}_test ${LIBRARIES})
        endif(USE_AMBIENT)
        add_test(${test}_${type}_test ${test}_${type}_test)
    endforeach(type)
endforeach(test)

