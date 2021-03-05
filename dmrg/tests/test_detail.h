/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2021         Leon Freitag <lefreita@ethz.ch>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/
#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include "maquis_dmrg.h"
#include <boost/filesystem/operations.hpp>

namespace test_detail {

    // Compare the (RDM) measurement element by element to the reference.
    // Allow differences in sign (phase) in the whole matrix if allow_phase_difference==true
    template<class Meas>
    void check_measurement_mat(const Meas& meas, const Meas& reference, bool allow_phase_difference=false)
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