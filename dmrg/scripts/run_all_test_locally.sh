#!/bin/env bash

cd ..

printf "\n=== Running Electronic ===\n"
cmake -B build/elec -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON
cmake --build build/elec

printf "=== Running Relativistic ===\n"
cmake -B build/rel -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG" -DCMAKE_BUILD_TYPE=Release
cmake --build build/rel

printf "=== Running preBO ===\n"
cmake -B build/preBO -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1" -DCMAKE_BUILD_TYPE=Release -DBUILD_PREBO=ON -DDMRG_NUMSYMM=5
cmake --build build/preBO

printf "=== Running Vibronic ===\n"
cmake -B build/vibronic -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="NU1;NONE" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRATIONAL=ON -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON -DBUILD_MPS2CI=ON -DDMRG_NUMSYMM=5 -DDMRG_ORDERNONE=12
cmake --build build/vibronic

printf "=== Running Vibrational ===\n"
cmake --build -B build/vibrational -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1PG;SU2U1PG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRONIC=ON
cmake --build build/vibrational

printf "=== Running TD Conventional ===\n"
cmake -B build/td-conv -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_TRANSCORRELATED_DMRG=ON
cmake --build build/td-conv

printf "=== Running TD All ===\n"
cmake -B build/td-all -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_PREBO=ON -DBUILD_VIBRONIC=ON -DBUILD_MPS2CI=ON -DBUILD_MPS_OVERLAP=ON -DBUILD_TRANSCORRELATED_DMRG=ON
cmake --build build/td-all

ctest --test-dir build/elec
ctest --test-dir build/rel
ctest --test-dir build/preBO
ctest --test-dir build/vibronic
ctest --test-dir build/vibrational
ctest --test-dir build/td-conv
ctest --test-dir build/td-all
