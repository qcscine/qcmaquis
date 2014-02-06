/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Bela Bauer <bauerb@comp-phys.org>
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

#ifdef USE_AMBIENT
#include <mpi.h>
#endif
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include <alps/hdf5.hpp>

#include "matrix_selector.hpp" /// define matrix
#include "symm_selector.hpp"   /// define grp

#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/utils/DmrgParameters.h"

#include "utils/timings.h"


int main(int argc, char ** argv)
{
    try {
        if (argc != 2 && argc != 3)
        {
            maquis::cout << "Usage: <parms> [<model_parms>]" << std::endl;
            exit(1);
        }
        
        maquis::cout.precision(10);
        
        /// Load parameters
        std::ifstream param_file(argv[1]);
        if (!param_file)
            throw std::runtime_error("Could not open parameter file.");
        DmrgParameters parms(param_file);
        
        /// Load model parameters from second input (if needed)
        std::string model_file;
        if (parms.is_set("model_file")) model_file = parms["model_file"].str();
        if (argc == 3)                  model_file = std::string(argv[2]);
        if (!model_file.empty()) {
            std::ifstream model_ifs(model_file.c_str());
            if (!model_ifs)
                throw std::runtime_error("Could not open model_parms file.");
            parms << ModelParameters(model_ifs);
        }
        
        /// Timers
#ifdef USE_AMBIENT
        ambient::timer
#else
        Timer
#endif
        tim_model      ("Parsing model"),  tim_load       ("Load MPS"),
        tim_l_boundary ("Left boundary"),  tim_r_boundary ("Right boundary"),
        tim_optim_jcd  ("Optim. JCD"   ),  tim_truncation ("Truncation"),
        tim_ts_obj     ("Two site obj.");
        
        /// Parsing model
        tim_model.begin();
        Lattice lattice(parms);
        Model<matrix, grp> model(lattice, parms);
        
        MPO<matrix, grp> mpo = make_mpo(lattice, model, parms);
        tim_model.end();
        maquis::cout << "Parsing model done!\n";
        
        
        /// Initialize & load MPS
        tim_load.begin();
        int L = lattice.size();
        MPS<matrix, grp> mps;
        load(parms["chkpfile"].str(), mps);
        int _site;
        {
            alps::hdf5::archive ar(parms["chkpfile"].str()+"/props.h5");
            ar["/status/site"] >> _site;
        }
        int site, lr;
        if (_site < L) {
            site = _site;
            lr = 1;
        } else {
            site = 2*L-_site-1;
            lr = -1;
        }
        tim_load.end();
        maquis::cout << "Load MPS done!\n";
        maquis::cout << "Optimization at site " << site << " in " << lr << " direction." << std::endl;
        
        mps.canonize(site);

        /// Compute left boundary
        tim_l_boundary.begin();
        Boundary<matrix, grp> left = mps.left_boundary();
        for (size_t i=0; i<site; ++i)
            left = contraction::overlap_mpo_left_step(mps[i], mps[i], left, mpo[i]);
        {
            std::string fname = "left" + boost::lexical_cast<std::string>(site) + ".h5";
            storage::archive ar(parms["chkpfile"].str()+"/"+fname, "w");
            ar["/props/site"] << site;
            ar["/props/type"] << "left";
            ar["/tensor"]     << left;
        }
        tim_l_boundary.end();
        maquis::cout << "Left boundary done!\n";
        
        /// Compute right boundary
        tim_r_boundary.begin();
        Boundary<matrix, grp> right = mps.right_boundary();
        for (int i=L-1; i>site+1; --i)
            right = contraction::overlap_mpo_right_step(mps[i], mps[i], right, mpo[i]);
        {
            std::string fname = "right" + boost::lexical_cast<std::string>(site+2) + ".h5";
            storage::archive ar(parms["chkpfile"].str()+"/"+fname, "w");
            ar["/props/site"] << site+2;
            ar["/props/type"] << "right";
            ar["/tensor"]     << right;
        }
        tim_r_boundary.end();
        maquis::cout << "Right boundary done!\n";
        
    } catch (std::exception & e) {
        maquis::cerr << "Exception caught:" << std::endl << e.what() << std::endl;
        exit(1);
    }
}

