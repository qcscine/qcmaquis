/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include "utils/io.hpp" // has to be first include because of impi
#include <iostream>

#include "maquis_dmrg.h"


// Test 1: H2 with d=3 Angstrom, singlet, 6-31G* basis set, FCI, integrals generated by Blueberry
BOOST_AUTO_TEST_CASE( PreBO_Test1 )
{

    DmrgParameters p;

    const maquis::integral_map<double, chem::Hamiltonian::PreBO> integrals {
        {{-1,-1,-1,-1,-1,-1,-1,-1},0.1763923530387111},
        {{0,0,0,0,-1,-1,-1,-1},-0.6802711028685324},
        {{0,0,0,2,-1,-1,-1,-1},0.07526840630272968},
        {{0,1,0,1,-1,-1,-1,-1},-0.6543632168967157},
        {{0,1,0,3,-1,-1,-1,-1},0.08728167494843309},
        {{0,2,0,0,-1,-1,-1,-1},0.07526840630272949},
        {{0,2,0,2,-1,-1,-1,-1},0.2283902278635599},
        {{0,3,0,1,-1,-1,-1,-1},0.08728167494843292},
        {{0,3,0,3,-1,-1,-1,-1},0.3423882983918999},
        {{0,0,0,0,0,0,0,0},0.1842790692925557},
        {{0,0,0,0,0,0,0,2},-0.03763420337488414},
        {{0,0,0,0,0,1,0,1},0.09872952891658415},
        {{0,0,0,0,0,1,0,3},-0.03477123411947585},
        {{0,0,0,0,0,2,0,0},-0.03763420337488413},
        {{0,0,0,0,0,2,0,2},0.03655689640589711},
        {{0,0,0,0,0,3,0,1},-0.03477123411947588},
        {{0,0,0,0,0,3,0,3},0.03327236099371499},
        {{0,0,0,1,0,0,0,1},0.09872952891658418},
        {{0,0,0,1,0,0,0,3},-0.03477123411947584},
        {{0,0,0,1,0,1,0,0},0.1866311043452651},
        {{0,0,0,1,0,1,0,2},-0.04128027869732104},
        {{0,0,0,1,0,2,0,1},-0.04481688352860439},
        {{0,0,0,1,0,2,0,3},0.03664246246971625},
        {{0,0,0,1,0,3,0,0},-0.03920603585723025},
        {{0,0,0,1,0,3,0,2},0.03658798479823142},
        {{0,0,0,2,0,0,0,0},-0.03763420337488416},
        {{0,0,0,2,0,0,0,2},0.03655689640589712},
        {{0,0,0,2,0,1,0,1},-0.04481688352860438},
        {{0,0,0,2,0,1,0,3},0.03664246246971625},
        {{0,0,0,2,0,2,0,0},0.1947198053616306},
        {{0,0,0,2,0,2,0,2},-0.05831200436614524},
        {{0,0,0,2,0,3,0,1},0.1068127293548811},
        {{0,0,0,2,0,3,0,3},-0.05322835824445955},
        {{0,0,0,3,0,0,0,1},-0.03477123411947588},
        {{0,0,0,3,0,0,0,3},0.03327236099371499},
        {{0,0,0,3,0,1,0,0},-0.03920603585723025},
        {{0,0,0,3,0,1,0,2},0.03658798479823138},
        {{0,0,0,3,0,2,0,1},0.1068127293548811},
        {{0,0,0,3,0,2,0,3},-0.05322835824445957},
        {{0,0,0,3,0,3,0,0},0.1906570567143694},
        {{0,0,0,3,0,3,0,2},-0.05745793841483472},
        {{0,1,0,0,0,0,0,1},0.1866311043452651},
        {{0,1,0,0,0,0,0,3},-0.03920603585723022},
        {{0,1,0,0,0,1,0,0},0.09872952891658415},
        {{0,1,0,0,0,1,0,2},-0.0448168835286044},
        {{0,1,0,0,0,2,0,1},-0.04128027869732109},
        {{0,1,0,0,0,2,0,3},0.03658798479823136},
        {{0,1,0,0,0,3,0,0},-0.03477123411947588},
        {{0,1,0,0,0,3,0,2},0.03664246246971627},
        {{0,1,0,1,0,0,0,0},0.09872952891658415},
        {{0,1,0,1,0,0,0,2},-0.04481688352860439},
        {{0,1,0,1,0,1,0,1},0.189844745754033},
        {{0,1,0,1,0,1,0,3},-0.04219754853151676},
        {{0,1,0,1,0,2,0,0},-0.04481688352860438},
        {{0,1,0,1,0,2,0,2},0.04106123456023925},
        {{0,1,0,1,0,3,0,1},-0.04219754853151678},
        {{0,1,0,1,0,3,0,3},0.03750279013385537},
        {{0,1,0,2,0,0,0,1},-0.04128027869732104},
        {{0,1,0,2,0,0,0,3},0.03658798479823134},
        {{0,1,0,2,0,1,0,0},-0.04481688352860438},
        {{0,1,0,2,0,1,0,2},0.04106123456023927},
        {{0,1,0,2,0,2,0,1},0.1984667878135254},
        {{0,1,0,2,0,2,0,3},-0.06104527396609304},
        {{0,1,0,2,0,3,0,0},0.1068127293548811},
        {{0,1,0,2,0,3,0,2},-0.06414087787820752},
        {{0,1,0,3,0,0,0,0},-0.03477123411947587},
        {{0,1,0,3,0,0,0,2},0.0366424624697162},
        {{0,1,0,3,0,1,0,1},-0.04219754853151677},
        {{0,1,0,3,0,1,0,3},0.03750279013385537},
        {{0,1,0,3,0,2,0,0},0.1068127293548812},
        {{0,1,0,3,0,2,0,2},-0.06414087787820755},
        {{0,1,0,3,0,3,0,1},0.1950618399114621},
        {{0,1,0,3,0,3,0,3},-0.05920528748005911},
        {{0,2,0,0,0,0,0,0},-0.03763420337488416},
        {{0,2,0,0,0,0,0,2},0.1947198053616306},
        {{0,2,0,0,0,1,0,1},-0.0448168835286044},
        {{0,2,0,0,0,1,0,3},0.1068127293548811},
        {{0,2,0,0,0,2,0,0},0.03655689640589711},
        {{0,2,0,0,0,2,0,2},-0.05831200436614524},
        {{0,2,0,0,0,3,0,1},0.03664246246971621},
        {{0,2,0,0,0,3,0,3},-0.05322835824445955},
        {{0,2,0,1,0,0,0,1},-0.04481688352860439},
        {{0,2,0,1,0,0,0,3},0.1068127293548811},
        {{0,2,0,1,0,1,0,0},-0.04128027869732105},
        {{0,2,0,1,0,1,0,2},0.1984667878135254},
        {{0,2,0,1,0,2,0,1},0.04106123456023928},
        {{0,2,0,1,0,2,0,3},-0.06414087787820746},
        {{0,2,0,1,0,3,0,0},0.03658798479823136},
        {{0,2,0,1,0,3,0,2},-0.06104527396609301},
        {{0,2,0,2,0,0,0,0},0.03655689640589711},
        {{0,2,0,2,0,0,0,2},-0.05831200436614522},
        {{0,2,0,2,0,1,0,1},0.04106123456023928},
        {{0,2,0,2,0,1,0,3},-0.06414087787820748},
        {{0,2,0,2,0,2,0,0},-0.0583120043661453},
        {{0,2,0,2,0,2,0,2},0.2278333714168219},
        {{0,2,0,2,0,3,0,1},-0.06414087787820751},
        {{0,2,0,2,0,3,0,3},0.13577048665041},
        {{0,2,0,3,0,0,0,1},0.0366424624697162},
        {{0,2,0,3,0,0,0,3},-0.05322835824445953},
        {{0,2,0,3,0,1,0,0},0.03658798479823135},
        {{0,2,0,3,0,1,0,2},-0.061045273966093},
        {{0,2,0,3,0,2,0,1},-0.06414087787820749},
        {{0,2,0,3,0,2,0,3},0.13577048665041},
        {{0,2,0,3,0,3,0,0},-0.05745793841483472},
        {{0,2,0,3,0,3,0,2},0.2231085965042993},
        {{0,3,0,0,0,0,0,1},-0.03920603585723025},
        {{0,3,0,0,0,0,0,3},0.1906570567143694},
        {{0,3,0,0,0,1,0,0},-0.03477123411947589},
        {{0,3,0,0,0,1,0,2},0.1068127293548812},
        {{0,3,0,0,0,2,0,1},0.03658798479823135},
        {{0,3,0,0,0,2,0,3},-0.05745793841483469},
        {{0,3,0,0,0,3,0,0},0.03327236099371495},
        {{0,3,0,0,0,3,0,2},-0.05322835824445952},
        {{0,3,0,1,0,0,0,0},-0.03477123411947589},
        {{0,3,0,1,0,0,0,2},0.1068127293548812},
        {{0,3,0,1,0,1,0,1},-0.04219754853151678},
        {{0,3,0,1,0,1,0,3},0.1950618399114621},
        {{0,3,0,1,0,2,0,0},0.03664246246971622},
        {{0,3,0,1,0,2,0,2},-0.06414087787820745},
        {{0,3,0,1,0,3,0,1},0.03750279013385537},
        {{0,3,0,1,0,3,0,3},-0.059205287480059},
        {{0,3,0,2,0,0,0,1},0.03658798479823136},
        {{0,3,0,2,0,0,0,3},-0.05745793841483467},
        {{0,3,0,2,0,1,0,0},0.03664246246971622},
        {{0,3,0,2,0,1,0,2},-0.06414087787820748},
        {{0,3,0,2,0,2,0,1},-0.06104527396609303},
        {{0,3,0,2,0,2,0,3},0.2231085965042992},
        {{0,3,0,2,0,3,0,0},-0.05322835824445951},
        {{0,3,0,2,0,3,0,2},0.1357704866504099},
        {{0,3,0,3,0,0,0,0},0.03327236099371494},
        {{0,3,0,3,0,0,0,2},-0.0532283582444595},
        {{0,3,0,3,0,1,0,1},0.03750279013385537},
        {{0,3,0,3,0,1,0,3},-0.05920528748005899},
        {{0,3,0,3,0,2,0,0},-0.05322835824445953},
        {{0,3,0,3,0,2,0,2},0.1357704866504099},
        {{0,3,0,3,0,3,0,1},-0.05920528748005899},
        {{0,3,0,3,0,3,0,3},0.2196551979193185}};


    std::vector<std::vector<std::vector<double>>> RDM
    {{{1.13103300179178,-2.52211470905664e-08,-0.0947313411126442,-9.1603632625322e-09},
    {-2.52211470905664e-08,0.855019531445115,3.53874389777801e-08,-0.0716856100438375},
    {-0.0947313411126442,3.53874389777801e-08,0.00793626405671728,-2.03230542773261e-09},
    {-9.1603632625322e-09,-0.0716856100438375,-2.03230542773261e-09,0.00601120270638985}}};


    p.set("integrals_binary", maquis::serialize(integrals));
    p.set("L", 4);
    p.set("LATTICE", "preBO lattice");
    p.set("MODEL", "PreBO");
    p.set("max_bond_dimension", 1000);
    p.set("PreBO_ParticleTypeVector",        "2"   );
    p.set("PreBO_FermionOrBosonVector",      "1"   );
    p.set("PreBO_OrbitalVector",             "4"   );
    p.set("PreBO_InitialStateVector",        "1 1" );
    p.set("nsweeps",4);
    p.set("symmetry", "nu1");

    std::vector<std::string> optimizer;
    optimizer.push_back("singlesite");
    optimizer.push_back("twosite");

    // Attention: const guess in test, so that the results are deterministic.
    p.set("init_type", "const");
    // Measure RDMs
    p.set("MEASURE[1rdm]","1");
    p.set("MEASURE[mutinf]","1");

    for (auto&& o: optimizer)
    {
        p.set("optimization", o);

        maquis::cout << "Running Pre-BO test for symmetry nu1 with optimization: " << o << std::endl;

        maquis::DMRGInterface<double> interface(p);
        interface.optimize();

        // test energy
        BOOST_CHECK_CLOSE(interface.energy(), -0.997454827563674, 1e-7);

        //// test 1-RDM
        const typename maquis::DMRGInterface<double>::meas_with_results_type& meas1 = interface.onerdm();
        for (auto i=0; i<meas1.first.size(); ++i) {
            auto ref = RDM[meas1.first[i][0]][meas1.first[i][1]][meas1.first[i][2]];
            auto val = meas1.second[i];
            auto diff = ref-val;
            BOOST_TEST(diff == 0.0, boost::test_tools::tolerance(1.0E-5));
        }
    }


}