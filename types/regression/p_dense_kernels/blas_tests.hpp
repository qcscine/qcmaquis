

BOOST_AUTO_TEST_CASE_TEMPLATE( minus, T, test_types)
{
    ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

    pMatrix pA(T::value,T::value);
    pMatrix pB(T::value,T::value);
    pMatrix pC(T::value,T::value);


    sMatrix sA(T::value,T::value);
    sMatrix sB(T::value,T::value);
    sMatrix sC(T::value,T::value);

    pA.set_init(ambient::random_i<typename T::dbl>);
    pB.set_init(ambient::random_i<typename T::dbl>);
    pC.set_init(ambient::null_i<typename T::dbl>);

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast
    sC = maquis::traits::matrix_cast<sMatrix>(pC); // playout is inside the cast
 
    sC = sA + sB; 
    pC = pA + pB; 

    BOOST_CHECK(pC==sC); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

BOOST_AUTO_TEST_CASE_TEMPLATE( dgemm, T, test_types)
{
    ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

    pMatrix pA(T::value,T::value);
    pMatrix pB(T::value,T::value);
    pMatrix pC(T::value,T::value);


    sMatrix sA(T::value,T::value);
    sMatrix sB(T::value,T::value);
    sMatrix sC(T::value,T::value);

    pA.set_init(ambient::random_i<typename T::dbl>);
    pB.set_init(ambient::random_i<typename T::dbl>);
    pC.set_init(ambient::null_i<typename T::dbl>);

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast
    sC = maquis::traits::matrix_cast<sMatrix>(pC); // playout is inside the cast
 
    gemm(sA,sB,sC); 
    gemm(pA,pB,pC); 

    BOOST_CHECK(pC==sC); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

