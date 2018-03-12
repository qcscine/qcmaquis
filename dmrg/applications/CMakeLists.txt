
option(BUILD_DMRG "Build DMRG Apps" ON)
option(BUILD_ALPS_APPS "Build ALPS applications (requires NU1 symmetry)" ON)
option(BUILD_TOOLS "Build tools" ON)
option(BUILD_MPI_APPS "Build mpi apps" OFF)
option(BUILD_CUSTOMAPP "Build custom applications" OFF)
option(BUILD_OLD_TARGETS "Enable targets with old names" OFF)


if(BUILD_DMRG)
  add_subdirectory(dmrg)
  add_subdirectory(measure)
  add_subdirectory(init)
  add_subdirectory(cideas)
else(BUILD_DMRG)
  add_subdirectory(dmrg EXCLUDE_FROM_ALL)
  add_subdirectory(measure EXCLUDE_FROM_ALL)
  add_subdirectory(init EXCLUDE_FROM_ALL)
endif(BUILD_DMRG)

if(BUILD_TOOLS)
  add_subdirectory(tools)
endif(BUILD_TOOLS)

if(BUILD_MPI_APPS)
  add_subdirectory(pdmrg)
  add_subdirectory(pcorr)
endif(BUILD_MPI_APPS)

if(BUILD_CUSTOMAPP)
  add_subdirectory(strcorr)
endif(BUILD_CUSTOMAPP)

if(DMRG_HAS_NU1 AND BUILD_ALPS_APPS)
  add_subdirectory(alps)
endif(DMRG_HAS_NU1 AND BUILD_ALPS_APPS)
