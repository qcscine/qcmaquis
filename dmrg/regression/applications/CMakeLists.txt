if(BUILD_DMRG)
  add_subdirectory(dmrg)
endif(BUILD_DMRG)

if(BUILD_TIME_EVOLUTION)
 add_subdirectory(dmrg_tevol)
endif(BUILD_TIME_EVOLUTION)

#if(BUILD_MULTIGRID)
#  add_subdirectory(multigrid)
#endif(BUILD_MULTIGRID)
