src_dir=".."

all: configure-all build-all test-all

configure-all: configure-elec configure-rel configure-prebo configure-vibrational configure-vibrational configure-tdconv configure-tdall
build-all: build-elec build-rel build-prebo build-vibrational build-vibrational build-tdconv build-tdall
test-all: test-elec test-rel test-prebo test-vibrational test-vibrational test-tdconv test-tdall

configure-elec:
	printf "\n=== Running Electronic ===\n"
	build_dir=build/elec
	cmake  -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON

build-elec:
  build_dir=build/elec
  cmake --build $build_dir

configure-rel:
  printf "=== Running Relativistic ===\n"
  build_dir=build/rel
  cmake  -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG" -DCMAKE_BUILD_TYPE=Release

build-rel:
  build_dir=build/rel
  cmake --build $build_dir

configure-prebo:
  printf "=== Running preBO ===\n"
  build_dir=build/preBO
  cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1" -DCMAKE_BUILD_TYPE=Release -DBUILD_PREBO=ON -DDMRG_NUMSYMM=5

build-prebo:
  build_dir=build/preBO
  cmake --build $build_dir

configure-vibrational:
  printf "=== Running Vibrational ===\n"
  build_dir=build/vibrational
  cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="NU1;NONE" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRATIONAL=ON -DBUILD_DMRG_FEAST=ON -DBUILD_SRCAS=ON -DBUILD_MPS2CI=ON -DDMRG_NUMSYMM=5 -DDMRG_ORDERNONE=12

build-vibrational:
  build_dir=build/vibrational
  cmake --build $build_dir

configure-vibronic:
  printf "=== Running Vibronic ===\n"
  build_dir=build/vibronic
  cmake -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1PG;SU2U1PG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_VIBRONIC=ON

build-vibronic:
  build_dir=build/vibronic
  cmake --build $build_dir

configure-tdconv:
  printf "=== Running TD Conventional ===\n"
  build_dir=build/td-conv
  cmake  -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_TRANSCORRELATED_DMRG=ON

build-tdconv:
  build_dir=build/td-conv
  cmake --build $build_dir

configure-tdall:
  printf "=== Running TD All ===\n"
  build_dir=build/td-all
  cmake  -S $src_dir -B $build_dir -DQCMAQUIS_TESTS=ON -DBUILD_SYMMETRIES="TwoU1;TwoU1PG;SU2U1;SU2U1PG;U1DG;NU1;U1" -DCMAKE_BUILD_TYPE=Release -DBUILD_DMRG_EVOLVE=ON -DBUILD_PREBO=ON -DBUILD_VIBRONIC=ON -DBUILD_MPS2CI=ON -DBUILD_MPS_OVERLAP=ON -DBUILD_TRANSCORRELATED_DMRG=ON

build-tdall:
  build_dir=build/td-all
  cmake --build $build_dir

test-elec:
  ctest --test-dir build/elec

test-rel:
  ctest --test-dir build/rel

test-prebo:
  ctest --test-dir build/preBO

test-vibrational:
  ctest --test-dir build/vibrational
  
test-vibronic:
  ctest --test-dir build/vibronic

test-tdconv:
  ctest --test-dir build/td-conv

test-tdall:
  ctest --test-dir build/td-all
