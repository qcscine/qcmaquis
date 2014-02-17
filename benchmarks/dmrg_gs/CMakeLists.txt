
macro(add_dmrg_bench BASEBENCH BOND_DIMENSION)
  set(bench_name dmrg_gs.bench.${BASEBENCH}.${BOND_DIMENSION})
  add_test(${bench_name} dmrg ${BASEBENCH}/parms.${BOND_DIMENSION} ${BASEBENCH}/model)
  set_property(TEST ${bench_name} PROPERTY LABELS performance)
endmacro(add_dmrg_bench)



# micro.L4
add_dmrg_bench(micro.L4  50)

# spinless.L6
add_dmrg_bench(spinless.L6  6000)

# spinless.L8
add_dmrg_bench(spinless.L8  6000)
add_dmrg_bench(spinless.L8  8000)
add_dmrg_bench(spinless.L8  10000)

# spinful.L8
add_dmrg_bench(spinful.L8  6000)
add_dmrg_bench(spinful.L8  8000)
add_dmrg_bench(spinful.L8  10000)
