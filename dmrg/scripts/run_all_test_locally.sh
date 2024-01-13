#!/bin/env bash

cd ..

printf "\n=== Running Electronic ===\n"
DIRECTORY=build/elec
cmake -B $DIRECTORY -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON
cmake --build $DIRECTORY

printf "=== Running Relativistic ===\n"
DIRECTORY=build/rel
cmake -B $DIRECTORY -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG" -DCMAKE_BUILD_TYPE=Release
cmake --build $DIRECTORY

printf "=== Running preBO ===\n"
DIRECTORY=build/preBO
cmake -B $DIRECTORY -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1" -DCMAKE_BUILD_TYPE=Release -DBUILD_PREBO=ON -DDMRG_NUMSYMM=5
cmake --build $DIRECTORY

printf "=== Running Vibrational ===\n"
DIRECTORY=build/vibrational
cmake -B $DIRECTORY -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="NU1;NONE" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRATIONAL=ON -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON -DBUILD_MPS2CI=ON -DDMRG_NUMSYMM=5 -DDMRG_ORDERNONE=12
cmake --build $DIRECTORY

printf "=== Running Vibronic ===\n"
DIRECTORY=build/vibronic
cmake --build -B $DIRECTORY -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1PG;SU2U1PG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRONIC=ON
cmake --build $DIRECTORY

printf "=== Running TD Conventional ===\n"
DIRECTORY=build/td-conv
cmake -B $DIRECTORY -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_TRANSCORRELATED_DMRG=ON
cmake --build $DIRECTORY

printf "=== Running TD All ===\n"
DIRECTORY=build/td-all
cmake -B $DIRECTORY -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_PREBO=ON -DBUILD_VIBRONIC=ON -DBUILD_MPS2CI=ON -DBUILD_MPS_OVERLAP=ON -DBUILD_TRANSCORRELATED_DMRG=ON
cmake --build $DIRECTORY

ctest --test-dir build/elec
ctest --test-dir build/rel
ctest --test-dir build/preBO
ctest --test-dir build/vibronic
ctest --test-dir build/vibrational
ctest --test-dir build/td-conv
ctest --test-dir build/td-all
