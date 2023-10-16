/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "maquis_dmrg.h"
#include <boost/filesystem/operations.hpp>

namespace test_detail {

    // Compare the (RDM) measurement element by element to the reference.
    // Allow differences in sign (phase) in the whole matrix if allow_phase_difference==true
    template<class Meas>
    void check_measurement_mat(const Meas& meas, const Meas& reference, bool allow_phase_difference=false,
                               double threshold=1.0E-16)
    {
        // check sizes
        BOOST_CHECK_EQUAL(meas.first.size(), reference.first.size());
        BOOST_CHECK_EQUAL(meas.second.size(), reference.second.size());
        // check the matrices themselves
        double phase = 1.;
        for (auto&& it = reference.first.begin(); it != reference.first.end(); it++)
        {
            auto index = std::distance(meas.first.begin(), std::find(meas.first.begin(), meas.first.end(), *it));
            auto index_ref = std::distance(reference.first.begin(), it);
            // the phase can be 1 or -1, so set the phase in the first run
            if (it == reference.first.begin() && allow_phase_difference)
                phase = std::copysign(1., meas.second[index]/reference.second[index_ref]);
            auto difference = std::abs(meas.second[index]-phase*reference.second[index_ref]);
            if (difference > threshold)
                BOOST_CHECK_SMALL(std::abs(meas.second[index]-phase*reference.second[index_ref]), 5e-7);
        }
    }

    // Class that inits a unique directory for checkpoints used in excited state tests and deletes it upon termination
    class TestTmpPath
    {
        public:
            TestTmpPath() : p_(boost::filesystem::unique_path(
                    boost::filesystem::temp_directory_path() / std::string("/qcmaquis_test_%%%%%%%%%%%%/" ) ) ) {
                try {
                    boost::filesystem::create_directories(p_);
                } catch (...) {
                    throw std::runtime_error("Cannot create temporary directory at " + p_.string() + ".\n");
                }
                std::cout << "Initialised temp directory at " << p_.c_str() << std::endl;
            }
            ~TestTmpPath() {  try { boost::filesystem::remove_all(p_); } catch( ... ) { } }
            const boost::filesystem::path& path() const { return p_; }
        private:
            boost::filesystem::path p_;

    };

}
#endif
