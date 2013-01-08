
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)  

get_symm_files(APPEND DMRG_MODELS_SOURCES "models/factory/factory_{SYMM}.cpp")

if(ENABLE_ALPS_MODELS)
  add_definitions(-DENABLE_ALPS_MODELS)
  get_symm_files(APPEND DMRG_MODELS_SOURCES "models/alps/symm_impl/model_{SYMM}.cpp")
endif(ENABLE_ALPS_MODELS)

if(ENABLE_LL_MODELS)
  list(APPEND DMRG_MODELS_SOURCES "models/ll/ll_models.cpp")
endif(ENABLE_LL_MODELS)


# *** Libraries
add_library(dmrg_utils STATIC utils/utils.cpp utils/stream_storage.cpp sim/symm_factory.cpp)
add_library(dmrg_models STATIC ${DMRG_MODELS_SOURCES})


# *** Install

install(TARGETS dmrg_utils COMPONENT libraries
                           DESTINATION lib)

install(TARGETS dmrg_models COMPONENT libraries
                            DESTINATION lib)
