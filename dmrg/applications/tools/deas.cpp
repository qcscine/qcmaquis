/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#include "dmrg/models/measurements/chementropy.h"


int main(int argc, char ** argv)
{
    std::cout.precision(12);

    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << " <result.h5>" << std::endl;
            return 1;
        }

        std::string rfile(argv[1]);

        EntanglementData<matrix> em(rfile);
        std::cout << "single orbital entropy:" << std::endl << em.s1() << std::endl;
        std::cout << "two orbital entropy:" << std::endl << em.s2() << std::endl;
        std::cout << "mutual information:" << std::endl << em.I() << std::endl;

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
