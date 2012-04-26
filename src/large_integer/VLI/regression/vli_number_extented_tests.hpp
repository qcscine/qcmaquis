
BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign_extented_CB, Vli, vli_extented_type)
{
    Vli a,b,c;

    std::size_t size (Vli::size);   

    for(std::size_t i(0); i < size-1; ++i) 
        a[i] = 0xffffffffffffffff;  
    b[0] = 1; 
    a+=b; 
    
    c[size-1] = 1;  
   
    BOOST_CHECK_EQUAL(a,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_extented_BB, Vli, vli_extented_type)
{
    Vli a,b,c;

    std::size_t size (Vli::size);   
    a[size-1] = 0x1;  
    b[0] = 0x1; 
    a-=b; 
    
    for(std::size_t i(0); i < size-1; ++i) 
        c[i] = 0xffffffffffffffff;  
   
    BOOST_CHECK_EQUAL(a,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( positive_multiply_64bits_assign_extented_BB, Vli, vli_extented_type)
{
    Vli a,c;
    unsigned long int b=2;

    std::size_t size (Vli::size);   

    for(std::size_t i(0); i < size-1; ++i) 
        a[i] = 0xffffffffffffffff;

    a*=b; 
    
    c[0] = 0xfffffffffffffffe;
 
    for(std::size_t i(1); i < size-1; ++i) 
        c[i] = 0xffffffffffffffff;  
     
    c[size-1] = 1;
   
    BOOST_CHECK_EQUAL(a,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( negative_multiply_64bits_assign_extented_BB, Vli, vli_extented_type)
{
    Vli a,c;
    long int b=-2;
    long int d=-1;

    std::size_t size (Vli::size);   

    for(std::size_t i(0); i < size-1; ++i) 
        a[i] = 0xffffffffffffffff;

    a*=b; 
    
    c[0] = 0x2;
 
    for(std::size_t i(1); i < size; ++i) 
        c[i] = 0;  
     
    c[size-1] = 0xfffffffffffffffe;  
    BOOST_CHECK_EQUAL(a,c);

    a*=d;

    c[0] = 0xfffffffffffffffe;
 
    for(std::size_t i(1); i < size-1; ++i) 
        c[i] = 0xffffffffffffffff;  
     
    c[size-1] = 1;

    BOOST_CHECK_EQUAL(a,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiply_nbits_assign_extented_BB, Vli, vli_extented_type)
{
    Vli a,c,b;
      
    b[0]=2;    

    std::size_t size (Vli::size);   

    for(std::size_t i(0); i < size-1; ++i) 
        a[i] = 0xffffffffffffffff;

    a*=b; 

    c[0] = 0xfffffffffffffffe;
 
    for(std::size_t i(1); i < size-1; ++i) 
        c[i] = 0xffffffffffffffff;  
     
    c[size-1] = 1;


    BOOST_CHECK_EQUAL(a,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ext_multiplies_positive_negative_int_gmp, Vli,vli_extented_type  )
{
    Vli a;
    fill_random(a,Vli::size/2);
    int b = rnd_valid_int<Vli>();
    
    b =-b;
    
    mpz_class agmp(a.get_str());
    
    a*=b;
    agmp *= b;
    
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_negative_numbers_gmp_192_192_192_ex, Vli, vli_extented_type )
{
    Vli a;
    Vli b;
    fill_random(a,Vli::size/2);
    fill_random(b,Vli::size/2); 
    b.negate();
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    a*=b;
    mpz_class cgmp = agmp * bgmp;
    
    
    BOOST_CHECK_EQUAL(a.get_str(),cgmp.get_str());
}
