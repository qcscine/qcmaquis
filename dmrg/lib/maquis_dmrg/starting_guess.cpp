/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
*               2019         Stefan Knecht <stknecht@ethz.ch>
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
#include "starting_guess.h"

#include "dmrg/models/measurements/chementropy.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/sim/matrix_types.h"

// namespace chem
// {
//     namespace detail
//     {
//         template <class V>
//         inline
//         BaseParameters set_su2u1_parameters(int L, int Nel, int spin, const chem::integral_map<V> & integrals, const std::vector<int> & site_types)
//         {
//             BaseParameters parms;

//             parms.set("L", L);
//             parms.set("nelec", Nel);
//             parms.set("spin", spin);

//             // take care of site types
//             std::string site_types_str;
//             assert(site_types.size() == L);
//             for (int i = 0; i < L; i++)
//                 site_types_str += std::to_string(site_types[i]) + ((i < L - 1) ? "," : "") ;
//             parms.set("site_types", site_types_str);

//             // integrals
//             parms.set("integrals_binary", chem::serialize(integrals));

//             return parms;
//         }
//     }
// }

namespace maquis
{
    template <class V>
    class StartingGuess<V>::Impl
    {
        public:
            typedef alps::numeric::matrix<V> Matrix;

            Impl(const DmrgParameters& parms, const std::string& pname, int nstates, bool do_fiedler, bool do_cideas)
                : parms_(parms), nstates_(nstates), pname_(pname), do_fiedler_(do_fiedler), do_cideas_(do_cideas)
            {

                parms_.erase_substring("MEASURE");

                if (do_fiedler_)
                {
                    if (parms_.is_set("orbital_order"))
                    {
                        // reset to default orbital order if some order is present
                        // if this isn't done, there're strange side-effects
                        std::string default_order;
                        int L = parms_["L"];
                        std::vector<int> v(L);
                        std::iota(v.begin(), v.end(), 1);
                        parms_.set("orbital_order", vector_tostring(v));
                    }
                    // we need mutual information for the Fiedler ordering
                    parms_.set("MEASURE[ChemEntropy]", 1);
                }

                measurements.reserve(nstates);

                // set sweeps and m, same values as in the old python interface
                parms_.set("nsweeps", 4);
                if (parms_.is_set("L"))
                    parms_.set("max_bond_dimension", parms_["L"] > 24 ? 256 : 128);
                else
                    throw std::runtime_error("L not defined for a starting guess calculation!");

                for (int i = 0; i < nstates; i++)
                {
                    // set correct checkpoints and result file names
                    std::string chkpfile = pname + ".checkpoint_state." + std::to_string(i) + ".h5";
                    std::string rfile = pname + ".results_state." + std::to_string(i) + ".h5";
                    parms_.set("chkpfile", chkpfile);
                    parms_.set("rfile", rfile);

                    if (i > 0)
                    {
                        parms_.set("n_ortho_states", i-1);
                        std::string all_ortho_states;
                        for (int j = 0; j < i; j++)
                            all_ortho_states += pname + ".checkpoint_state." + std::to_string(j-1) + ".h5" + ((j < i-1) ? ";" : "") ;
                        parms_.set("ortho_states", all_ortho_states);
                    }

                    interface.reset(new DMRGInterface<V>(parms_));
                    interface->optimize();
                    measurements.emplace_back(std::move(interface->measurements()));
                }
            }

            // calculate Fiedler order
            std::string getFiedlerOrder()
            {
                if (!do_fiedler_)
                    throw std::runtime_error("Cannot obtain Fiedler ordering if do_fiedler_ was not enabled.");
                // TODO: implement also Block fiedler ordering per symmetry

                // Collect mutual information from all the states
                std::vector<Matrix> mutI;
                mutI.reserve(nstates_);

                for (int i = 0; i < nstates_; i++)
                {

                    // get the entropy data
                    EntanglementData<Matrix> em(get_measurements()[i]);
                    // store data - scaled with state-average factor
                    mutI.emplace_back(std::move(em.I()));
                }

                // accumulate mutual information matrix
                Matrix SAmutI(mutI[0].num_rows(), mutI[0].num_cols());
                for (auto& n : mutI)
                    SAmutI += n;
                // get Laplacian
                Matrix L = get_laplacian(SAmutI);

                if (L.num_rows() < 2)
                    throw std::runtime_error("Fiedler vector orbital ordering doesn't work for only one orbital!");

                Matrix evecs(L.num_rows(), L.num_cols());
                std::vector<V> evals(L.num_rows());
                alps::numeric::syev(L,evecs,evals);

                // get the Fiedler vector, i.e. the eigenvector corresponding to the second lowest eigenvalue of the Laplacian
                // The eigenvalues in evecs are assumed to be sorted starting from the highest eigenvalue
                // i.e. the second lowest eigenvalue has an index L-2
                auto fv_col = evecs.col(L.num_rows()-2);
                std::vector<V> fiedler_vector(fv_col.first, fv_col.second);

                // prepare ordering. first create a vector with indices 0..L-1 in ascending order
                std::vector<int> order(fiedler_vector.size());
                std::iota(order.begin(), order.end(), 0);

                // Sort the order vector according to the Fiedler vector
                std::sort(order.begin(), order.end(), [&fiedler_vector](size_t i1, size_t i2) { return fiedler_vector[i1] < fiedler_vector[i2]; });

                // add 1 to each element because in the parameters our counting starts with 1
                for (auto&& n: order) n++;
                // std::transform(order.begin(), order.end(), order.begin(), [](int i){ return i+1; });

                // convert the ordering into a string
                return vector_tostring(order);
            }

        private:
            DmrgParameters parms_;
            std::string pname_;
            int nstates_;

            std::unique_ptr<DMRGInterface<V> > interface;
            std::vector<results_map_type<V> > measurements;

            bool do_fiedler_, do_cideas_;

            // Calculate the Laplacian
            Matrix get_laplacian(const Matrix & mutI)
            {
                // Ported to C++ from fiedler.py

                Matrix laplacian(mutI.num_rows(), mutI.num_cols(), 0.0);

                for (int i = 0; i < mutI.num_rows(); i++)
                    for (int j = 0; j < mutI.num_cols(); j++)
                        laplacian(i,i) += mutI(i,j);

                laplacian -= mutI;
                return laplacian;
            }

            template <class K>
            std::string vector_tostring(const std::vector<K> & v)
            {
                std::string s;

                for (int i = 0; i < v.size(); i++)
                    s += std::to_string(v[i]) + ((i < v.size()-1) ? "," : "") ;
                return s;
            }

            // Obtain measurements
            const std::vector<results_map_type<V> > & get_measurements() { return measurements; };


    };

    template <class V>
    StartingGuess<V>::StartingGuess(const DmrgParameters& parms, int nstates, const std::string & pname, bool do_fiedler, bool do_cideas)
        : impl_(new Impl(parms, pname, nstates, do_fiedler, do_cideas)) {}

    template <class V>
    StartingGuess<V>::~StartingGuess() = default;

    template <class V>
    std::string StartingGuess<V>::getFiedlerOrder() { return impl_-> getFiedlerOrder(); }

    template class StartingGuess<double>;
}



