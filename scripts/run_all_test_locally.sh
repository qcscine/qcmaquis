#!/bin/env bash

export CC=gcc CXX=g++

src_dir=".."

printf "\n=== Running Electronic ===\n"
build_dir=build/elec
cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON
cmake --build $build_dir

printf "=== Running Relativistic ===\n"
build_dir=build/rel
cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG" -DCMAKE_BUILD_TYPE=Release
cmake --build $build_dir

printf "=== Running preBO ===\n"
build_dir=build/preBO
cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1" -DCMAKE_BUILD_TYPE=Release -DBUILD_PREBO=ON -DDMRG_NUMSYMM=5
cmake --build $build_dir

printf "=== Running Vibrational ===\n"
build_dir=build/vibrational
cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="NU1;NONE" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRATIONAL=ON -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON -DBUILD_MPS2CI=ON -DDMRG_NUMSYMM=5 -DDMRG_ORDERNONE=12
cmake --build $build_dir

printf "=== Running Vibronic ===\n"
build_dir=build/vibronic
cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1PG;SU2U1PG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRONIC=ON
cmake --build $build_dir

printf "=== Running TD Conventional ===\n"
build_dir=build/td-conv
cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_TRANSCORRELATED_DMRG=ON
cmake --build $build_dir

printf "=== Running TD All ===\n"
build_dir=build/td-all
cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_PREBO=ON -DBUILD_VIBRONIC=ON -DBUILD_MPS2CI=ON -DBUILD_MPS_OVERLAP=ON -DBUILD_TRANSCORRELATED_DMRG=ON
cmake --build $build_dir

ctest --test-dir build/elec
ctest --test-dir build/rel
ctest --test-dir build/preBO
ctest --test-dir build/vibronic
ctest --test-dir build/vibrational
ctest --test-dir build/td-conv
ctest --test-dir build/td-all
