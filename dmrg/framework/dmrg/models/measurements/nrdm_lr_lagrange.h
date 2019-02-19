/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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
#ifndef NRDM_LR_LAGRANGE_H
#define NRDM_LR_LAGRANGE_H

#include "dmrg/models/measurements/tagged_nrankrdm.h"
#include "dmrg/mp_tensors/twositetensor.h"

    namespace measurements
    {
        template <class Matrix, class SymmGroup>
        class NRDMLRLagrange : public TaggedNRankRDM<Matrix, SymmGroup>
        {
        typedef TaggedNRankRDM<Matrix, SymmGroup> base;
        typedef typename base::tag_vec tag_vec;
        typedef typename base::scaled_bond_term scaled_bond_term;
        typedef typename base::positions_type positions_type;


        using base::ext_labels;
        public:
            // Specialization to call the correct constructor of the base class
            // TODO: Check if this compiles with disabled SU2U1/SU2U1PG!
            // symm_traits::HasSU2<SU2U1> and symm_traits::HasSU2<SU2U1PG> yields boost::true_type
            NRDMLRLagrange(int lr_site_, boost::true_type, std::string const& filename_, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       typename TermMakerSU2<Matrix, SymmGroup>::OperatorCollection const & op_collection_,
                       positions_type const& positions_ = positions_type(), bool twosite_ = true)
                       : base(name_, lat, tag_handler_, op_collection_, positions_), filename(filename_), lr_site(lr_site_), twosite(twosite_) {};

            // 2U1 and other symmetry groups
            NRDMLRLagrange(int lr_site_, boost::false_type, std::string const& filename_, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       tag_vec const & identities_, tag_vec const & fillings_, std::vector<scaled_bond_term> const& ops_,
                       bool half_only_, positions_type const& positions_ = positions_type(), bool twosite_ = true)
                        : base(name_, lat, tag_handler_, identities_, fillings_, ops_, half_only_, positions_), filename(filename_), lr_site(lr_site_), twosite(twosite_) {};


            virtual void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
            {
                // Wrapper for the TDM evaluation function
                std::function<void(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &)> RDMEvaluator;

                // Assign the appropriate TDM evaluation function to RDMEvaluator based on which derivative should be calculated
                if (this->name() == "onerdmlagrangeR")
                    RDMEvaluator = std::bind(&NRDMLRLagrange<Matrix, SymmGroup>::measure_correlation, this, std::placeholders::_1, std::placeholders::_2);

                if (this->name() == "twordmlagrangeR")
                    RDMEvaluator = std::bind(&NRDMLRLagrange<Matrix, SymmGroup>::measure_2rdm, this, std::placeholders::_1, std::placeholders::_2);

                // for left RDM derivatives, the bra and ket in the RDM evaluation are swapped
                if (this->name() == "onerdmlagrangeL")
                    RDMEvaluator = std::bind(&NRDMLRLagrange<Matrix, SymmGroup>::measure_correlation, this, std::placeholders::_2, std::placeholders::_1);

                if (this->name() == "twordmlagrangeL")
                    RDMEvaluator = std::bind(&NRDMLRLagrange<Matrix, SymmGroup>::measure_2rdm, this, std::placeholders::_2, std::placeholders::_1);

                measure_lagrange_rdm(ket_mps, RDMEvaluator, lr_site);
            }
        protected:
            measurement<Matrix, SymmGroup>* do_clone() const
            {
                return new NRDMLRLagrange(*this);
            }

            // MEASURE_LAGRANGE_RDM
            //
            // Calculate lagrange RDM update (MPS contribution to the Lagrange effective RDM for gradient calculations)

            void measure_lagrange_rdm(MPS<Matrix, SymmGroup> const& mps, std::function<void(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &)> RDMEvaluator, int site = 0)
            {

                typedef typename Matrix::value_type value_type;
                MPS<Matrix, SymmGroup> mps_aux = mps;

                // The MPS needs to be canonized up to the site before the measurement to have the same local basis in all the states
                // However, we canonize all MPS simultaneously, so we don't need to canonize here again.
                // mps_aux.canonize(site);

                // Read auxiliary MPSTensor elements into a vector aux_elements from file written by Yingjin's LR program
                std::vector<value_type> aux_elements;

                maquis::cout << "Reading the auxiliary MPSTensor elements for site " << site << " from file " << filename << std::endl;

                // read and parse the file
                std::ifstream infile(filename);
                if (infile)
                    std::copy(std::istream_iterator<value_type>(infile), std::istream_iterator<value_type>(), std::back_inserter(aux_elements));
                else
                    throw std::runtime_error("File " + filename + " could not be opened!");


                if (twosite) // two-site
                {
                    // forbid the last site because we can't build a two-site tensor
                    if (site > mps.length() - 1)
                        throw std::runtime_error("site > L-1 not allowed in building Lagrange RDM contributions with two-site tensors");

                    // Prepare the two-site tensor from two sites of the MPS
                    TwoSiteTensor<Matrix, SymmGroup> tst(mps_aux[site], mps_aux[site+1]);

                    // To keep the consistency with other measurements (local Hamiltonian), or the yingjin-devel branch
                    // which apparently works with left-paired two site tensors, we introduce left pairing
                    // Note that pairing may introduce additional blocks!
                    tst.make_left_paired();

                    // Fill in the TwoSiteTensor with the values we just read from the file

                    assert(tst.data().num_elements() == aux_elements.size());

                    size_t fileidx = 0;
                    for (size_t i = 0; i < tst.data().n_blocks(); i++)
                    for (size_t j = 0; j < tst.data()[i].num_rows(); j++)
                    for (size_t k = 0; k < tst.data()[i].num_cols(); k++)
                    {
                        parallel::guard::serial guard;
                        tst.data()[i](j,k) = aux_elements[fileidx++];
                    }


                    // Incorporate the modified TwoSiteTensor back into the MPS to form auxiliary state
                    // Note that we provide a m value as large as possible (maxint) and a negative cutoff in order to prevent truncation.

                    truncation_results trunc;
                    std::tie(mps_aux[site], mps_aux[site+1], trunc) = tst.split_mps_l2r(std::numeric_limits<int>::max(), -1.0);

                    // maquis::cout << "Truncated weight in evaluating RDM: " << trunc.truncated_weight << std::endl;

                    // for whatever reason! -- this is to ensure the reproducibility of the results by the yingjin-devel branch
                    if (site < mps_aux.length()-2)
                    {
                        mps_aux[site+2].multiply_from_left(mps_aux[site+1].normalize_left(DefaultSolver()));
                        mps_aux[site+2].make_left_paired();
                    }

                    // Fix the left pairing since manipulation of the MPS messes it up
                    // FIXME: there is the same race condition as in nrdmderivative.h
                    mps_aux[site].make_left_paired();
                    mps_aux[site+1].make_left_paired();

                    // measure the transition RDM <mps_aux|c+...c...|mps>
                    // TODO: Note that we will need both <mps_aux|c+...c...|mps> and <mps|c+...c...|mps_aux> for symmetrised derivatives.
                    // The symmetrisation will be taken care for later
                    RDMEvaluator(mps_aux, mps);
                }
                else // one-site
                {
                    MPSTensor<Matrix, SymmGroup> & mpst = mps_aux[site];
                    mpst.make_left_paired();
                    assert(mpst.data().num_elements() == aux_elements.size());
                    size_t fileidx = 0;
                    for (size_t i = 0; i < mpst.data().n_blocks(); i++)
                    for (size_t j = 0; j < mpst.data()[i].num_rows(); j++)
                    for (size_t k = 0; k < mpst.data()[i].num_cols(); k++)
                    {
                        parallel::guard::serial guard;
                        mpst.data()[i](j,k) = aux_elements[fileidx++];
                    }

                    // for whatever reason! -- this is to ensure the reproducibility of the results by the yingjin-devel branch
                    // TODO: Check if this is really needed for one-site MPS parameters!
                    // if (site < mps_aux.length()-1)
                    // {
                    //     mps_aux[site+1].multiply_from_left(mpst.normalize_left(DefaultSolver()));
                    //     mps_aux[site+1].make_left_paired();
                    // }

                    mpst.make_left_paired();
                    RDMEvaluator(mps_aux, mps);
                }
            }
        private:
            const std::string& filename;
            int lr_site;
            bool twosite;
        };
    } // measurements

#endif