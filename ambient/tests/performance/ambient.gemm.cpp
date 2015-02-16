#include "utils/testing.hpp"

TEST_CASE( "Matrix multiplication performance measured", "[ambient::gemm]" )
{
    measurement::timer gtime("total"); gtime.begin();
    measurement params;

    size_t x = params.num_cols();
    size_t y = params.num_rows();

    matrix<double> A(x, y);
    matrix<double> B(x, y);
    matrix<double> C(x, y);
    matrix<double> C_orig(x, y);

    generate(A);
    generate(B);
    ambient::sync();

    printf("ambient::gemm strassen...\n");
    ambient::numeric::gemm_strassen(std::move(A), std::move(B), std::move(C)); 
    measurement::timer time("ambient::gemm_strassen"); time.begin();
    ambient::sync();
    time.end();

    printf("ambient::gemm...\n");
    ambient::numeric::gemm(A, B, C_orig); 
    measurement::timer time_orig("ambient::gemm"); time_orig.begin();
    ambient::sync();
    time_orig.end();

    params.report(gflops::gemm, time_orig.get_time());
    params.report(gflops::gemm, time.get_time());

    gtime.end();
    std::cout << "Global time: " << gtime.get_time() << "\n";

    REQUIRE((C == C_orig));
}
