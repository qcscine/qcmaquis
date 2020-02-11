/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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
#include <cmath>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/models/generate_mpo.hpp"

#include "dmrg/models/chem/su2u1/model.h"

#include <mpi_interface.h>

namespace maquis
{
  Scine::Mpi::MpiInterface* mpi__;
  std::unique_ptr<Scine::Mpi::MpiInterface> mpi;
}


#if defined(USE_TWOU1)
typedef TwoU1 symm;
#elif defined(USE_U1DG)
typedef U1DG symm;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG symm;
#elif defined(USE_SU2U1)
typedef SU2U1 symm;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG symm;
#elif defined(USE_NONE)
typedef TrivialGroup symm;
#elif defined(USE_U1)
typedef U1 symm;
#endif

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

MPO<matrix, symm> rdm1(Model<matrix, symm> const & model, Lattice const & lattice, size_t const & myrank)
{
        typedef operator_selector<matrix, symm>::type op_t;
        typedef OPTable<matrix, symm>::tag_type tag_type;
        typedef matrix::value_type value_type;
        typedef typename Model<matrix, symm>::term_descriptor term_descriptor;
        typename TermMakerSU2<matrix, symm>::OperatorCollection op_collection;

        boost::shared_ptr<TagHandler<matrix, symm> > tag_handler_local = model.operators_table();

        const int nterm = 2;
        int pos_[nterm] = {myrank, myrank};
        std::vector<int> pos(pos_, pos_ + nterm);

        std::vector<tag_type> ident_no_couple;
        std::vector<tag_type> ident_full_no_couple;
        std::vector<tag_type> fill_no_couple;
        for (size_t p = 0; p <= lattice.size(); ++p)
            {
                ident_no_couple.push_back(model.identity_matrix_tag(p));
                ident_full_no_couple.push_back(model.get_operator_tag("ident_full", p));
                fill_no_couple.push_back(model.get_operator_tag("ident_full", p));
            }


        std::vector<term_descriptor> terms;

        if (pos_[0] != pos_[1])
                        // The sqrt(2.) balances the magnitudes of Clebsch coeffs C^{1/2 1/2 0}_{mrm'} which apply at the second spin-1/2 operator
                        terms.push_back(TermMakerSU2<matrix, symm>::positional_two_term(
                            true, op_collection.ident.no_couple, std::sqrt(2.), pos_[0], pos_[1], op_collection.create.couple_down, op_collection.create.fill_couple_up,
                                                              op_collection.destroy.couple_down, op_collection.destroy.fill_couple_up, lattice
                        ));
        else {
                        term_descriptor term;
                        term.coeff = 1.;
                        maquis::cout << "I am creating term " << pos_[0] << pos_[1] << " " << lattice.get_prop<typename symm::subcharge>("type", pos_[0]) << std::endl;
                        term.push_back( boost::make_tuple(pos_[0], model.get_operator_tag("count", lattice.get_prop<typename symm::subcharge>("type", pos_[0]))));
                        terms.push_back(term);
        }

        generate_mpo::TaggedMPOMaker<matrix, symm> mpo_m(lattice, ident_no_couple, ident_full_no_couple,
                                                              fill_no_couple, tag_handler_local, terms);
        MPO<matrix, symm> mpo = mpo_m.create_mpo();


        //maquis::cout << terms << std::endl;
        //maquis::cout << "terms ops: " << std::endl;
        //for (int i=0; i<term.size(); ++i)
        //    maquis::cout << tag_handler_local->get_op(term.operator_tag(i));

        maquis::cout << "MPO generated " << std::endl;

        return mpo;
}

int main(int argc, char ** argv)
{

    // setup MPI interface. It does nothing for serial runs
    if (!maquis::mpi__) {
        maquis::mpi   = std::unique_ptr<Scine::Mpi::MpiInterface>(new Scine::Mpi::MpiInterface(nullptr, 0));
        maquis::mpi__ = maquis::mpi.get();
    }

    DmrgOptions opt(argc, argv);
    if (!opt.valid) return 0;
    DmrgParameters parms = opt.parms;

    maquis::cout.precision(10);

    /// Parsing model
    Lattice lattice = Lattice(parms);
    Model<matrix, symm> model = Model<matrix, symm>(lattice, parms);

    MPS<matrix, symm> mps;

    if(maquis::mpi__->getGlobalRank() == 0){
        // load state
        maquis::cout << "Loading " << parms["chkpfile"] << std::endl;
        std::string wvf = parms["chkpfile"];
        load(wvf, mps);
    } else {
    }

    //MPO<matrix, symm> mpo = rdm2(model, lattice);
    MPO<matrix, symm> mpo = rdm1(model, lattice, maquis::mpi__->getGlobalRank());

    maquis::cout << "enter expval " << std::endl;
    if(maquis::mpi__->getGlobalRank() == 0){
        maquis::traits::real_type<matrix::value_type>::type value = maquis::real(expval(mps, mpo));
        maquis::cout << "Expval is: " << value << std::endl;
    }

    // terminate MPI (does nothing if serial run)
    maquis::mpi.reset(nullptr);
    maquis::mpi__ = nullptr;
}
