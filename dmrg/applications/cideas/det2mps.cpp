/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
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

#include <iostream>
#include <boost/lexical_cast.hpp>

#include "dmrg/sim/matrix_types.h"
#include "dmrg/models/chem/cideas/cideas.hpp"
#include "dmrg/utils/DmrgOptions.h"

#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#endif

//parse determinants
template <class SymmGroup>
std::vector<Determinant<SymmGroup> > dets_from_file(std::string file){
    std::ifstream config_file;
    config_file.open(file.c_str());

    std::vector<Determinant<SymmGroup> > configs;

    for (std::string line; std::getline(config_file, line); ) {
        std::vector<std::string> det_coeff;
        boost::split(det_coeff, line, boost::is_any_of(" "));

        std::string det = det_coeff[0];

        Determinant<SymmGroup> tmp;
        for (std::size_t i = 0; i < det.size(); ++i) {
            int occ = boost::lexical_cast<size_t>(det[i]);
            switch(occ) {
                case 4:
                    tmp.push_back(4); // doubly occ
                    break;
                case 3:
                    tmp.push_back(3); // up
                    break;
                case 2:
                    tmp.push_back(2); // down
                    break;
                case 1:
                    tmp.push_back(1); // empty
                    break;
            }
        }
        if(std::find(configs.begin(), configs.end(), tmp) == configs.end())
            configs.push_back(tmp);
    }
    return configs;
}

int main(int argc, char ** argv){
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << "dmrg-input file" << std::endl;
            return 1;
        }

        DmrgOptions opt(argc, argv);
        DmrgParameters& parms = opt.parms;

        std::string rfile(parms.get<std::string>("resultfile"));
        EntanglementData<matrix> em(rfile);

        matrix s1 = em.s1();

        MPS<matrix,grp> mps = maquis::cideas<matrix, grp>(parms, s1);
        maquis::cout << "MPS created" << std::endl;

        std::string chkp = parms["chkpfile"].str();
        save(chkp, mps);
        storage::archive ar2(chkp+"/props.h5", "w");
        ar2["/parameters"] << parms;
        //TODO obtain version string from cmake
        //ar2["/version"] << DMRG_VERSION_STRING;
        ar2["/status/sweep"] << -1;
        ar2["/status/site"] << -1;

        maquis::cout << "MPS saved" << std::endl;

        } catch (std::exception& e) {
            std::cerr << "Error:" << std::endl << e.what() << std::endl;
            return 1;
        }
}
