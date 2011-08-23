

BOOST_AUTO_TEST_CASE_TEMPLATE ( constructor, Vli, vli_types )
{
    monomial<Vli> ma;
    monomial<Vli> mb(0,0);
    BOOST_CHECK_EQUAL(ma,mb);

    monomial<Vli> mc(1,2);
    monomial<Vli> md(2,1);

    BOOST_CHECK_EQUAL(mc == md, false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( copy_constructor, Vli, vli_types )
{
    monomial<Vli> ma;
    monomial<Vli> mb(8,5);
    monomial<Vli> mc(mb);
    BOOST_CHECK_EQUAL(mb,mc);
    BOOST_CHECK_EQUAL(ma == mc,false);

    mc.j_exp_ = 1;
    BOOST_CHECK_EQUAL(mb == mc,false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( multiply_by_vli, Vli, vli_types )
{
    Vli a;
    fill_random(a,Vli::size-1);
    Vli a_orig(a);
    monomial<Vli> ma(2,3);
    monomial<Vli> ma_orig(ma);

    monomial<Vli> mb = ma * a;
    monomial<Vli> mc = a * ma;
    BOOST_CHECK_EQUAL(ma, ma_orig);

    ma *= a;
    BOOST_CHECK_EQUAL(ma,mb);
    BOOST_CHECK_EQUAL(ma,mc);

    BOOST_CHECK_EQUAL(ma == ma_orig, false);
    BOOST_CHECK_EQUAL(ma.coeff_,a);

    BOOST_CHECK_EQUAL(a, a_orig);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( multiply_by_int, Vli, vli_types )
{
    int a = 10;
    int a_orig(a);
    monomial<Vli> ma(5,7);
    monomial<Vli> ma_orig(ma);

    monomial<Vli> mb = ma * a;
    monomial<Vli> mc = a* ma;
    BOOST_CHECK_EQUAL(ma, ma_orig);

    ma *= a;
    BOOST_CHECK_EQUAL(ma,mb);
    BOOST_CHECK_EQUAL(ma,mc);

    BOOST_CHECK_EQUAL(ma == ma_orig, false);
    BOOST_CHECK_EQUAL(ma.coeff_,Vli(a));
    
    BOOST_CHECK_EQUAL(a, a_orig);
}
