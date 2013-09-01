
option(BUILD_DMRG "Build DMRG Apps" ON)
option(BUILD_TIME_EVOLUTION "Build DMRG-Time evolution Apps" ON)
option(BUILD_TOOLS "Build tools" ON)
option(BUILD_DIAG "Build Diagonalization App" OFF)
option(BUILD_MULTIGRID "Build DMRG-Multigrid Apps" OFF)
option(BUILD_CUSTOMAPP "Build custom applications" OFF)
option(BUILD_EXAMPLES "Build example applications" OFF)
option(BUILD_OLD_TARGETS "Enable targets with old names" OFF)


if(BUILD_DMRG)
  add_subdirectory(dmrg)
  add_subdirectory(init)
else(BUILD_DMRG)
  add_subdirectory(dmrg EXCLUDE_FROM_ALL)
  add_subdirectory(init EXCLUDE_FROM_ALL)
endif(BUILD_DMRG)

if(BUILD_TIME_EVOLUTION)
  add_subdirectory(time_evolution)
else(BUILD_TIME_EVOLUTION)
  add_subdirectory(time_evolution EXCLUDE_FROM_ALL)
endif(BUILD_TIME_EVOLUTION)

if(BUILD_TOOLS)
  add_subdirectory(tools)
endif(BUILD_TOOLS)

if(BUILD_DIAG)
  add_subdirectory(diag)
else(BUILD_DIAG)
  add_subdirectory(diag EXCLUDE_FROM_ALL)
endif(BUILD_DIAG)

if(BUILD_MULTIGRID)
  add_subdirectory(multigrid)
else(BUILD_MULTIGRID)
  add_subdirectory(multigrid EXCLUDE_FROM_ALL)
endif(BUILD_MULTIGRID)

if(BUILD_CUSTOMAPP)
  add_subdirectory(strcorr)
endif(BUILD_CUSTOMAPP)


if(BUILD_EXAMPLES)
  add_subdirectory(simple_main/dmrg_u1_current)
  add_subdirectory(simple_main/dmrg_current)
endif(BUILD_EXAMPLES)
