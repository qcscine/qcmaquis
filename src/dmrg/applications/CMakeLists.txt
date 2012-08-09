
option(BUILD_DMRG "Build DMRG Apps" ON)
option(BUILD_TIME_EVOLUTION "Build DMRG-Time evolution Apps" ON)
option(BUILD_DIAG "Build Diagonalization App" ON)
option(BUILD_MULTIGRID "Build DMRG-Multigrid Apps" OFF)
option(BUILD_CUSTOMAPP "Build custom applications" OFF)
option(BUILD_OLD_TARGETS "Enable targets with old names" OFF)


if(BUILD_DMRG)
  add_subdirectory(dmrg)
else(BUILD_DMRG)
  add_subdirectory(dmrg EXCLUDE_FROM_ALL)
endif(BUILD_DMRG)

if(BUILD_TIME_EVOLUTION)
  add_subdirectory(time_evolution)
else(BUILD_TIME_EVOLUTION)
  add_subdirectory(time_evolution EXCLUDE_FROM_ALL)
endif(BUILD_TIME_EVOLUTION)

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
  add_subdirectory(misc)
endif(BUILD_CUSTOMAPP)

