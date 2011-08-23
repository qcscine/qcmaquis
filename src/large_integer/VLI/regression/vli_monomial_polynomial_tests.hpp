
#define MAXORDER 4

BOOST_AUTO_TEST_CASE_TEMPLATE(monome_polynome, Vli, vli_types)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    //Init CPU
    Vli a;
    fill_random(a);    
    monomial<Vli> ma(a);
    polynomial_cpu<Vli, MAXORDER> pa;    
    fill_poly_random(pa);

    //Init GPU
    vli_gpu<typename Vli::value_type,Vli::size> a_gpu(a);
    monomial<vli_gpu<typename Vli::value_type,Vli::size> > magpu(a);
    polynomial_gpu<vli_gpu<typename Vli::value_type,Vli::size>, MAXORDER> pagpu(pa);

    pa = pa * ma; 
    pagpu = pagpu * magpu; 
 
    BOOST_CHECK_EQUAL(pa, pagpu);
}
