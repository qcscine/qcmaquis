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
#ifndef LOCAL_HAMILTONIAN_H
#define LOCAL_HAMILTONIAN_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/measurements/local_hamiltonian_initializer.h"
#include "dmrg/utils/time_stopper.h"

namespace measurements
{
    template<class Matrix, class SymmGroup>
    class LocalHamiltonian : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        public:
            LocalHamiltonian(std::string name_, Lattice const & lat_, BaseParameters & parms_) : base(name_), lat(lat_), parms(parms_) {}
            virtual void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
            {
                throw std::runtime_error("Local Hamiltonian must be measured with the Hamiltonian provided as an MPO!");
            }

            virtual void evaluate(MPS<Matrix, SymmGroup> const& mps, MPO<Matrix, SymmGroup> const& mpo)
            {
                this->vector_results.clear();
                this->labels.clear();
                // Right now this is hard-coded to local Hamiltonian at sites 0 and 1
                measure_local_hamiltonian(mps, mpo);

            }
        protected:
            measurement<Matrix, SymmGroup>* do_clone() const
            {
                return new LocalHamiltonian(*this);
            }

            void measure_local_hamiltonian(MPS<Matrix, SymmGroup> const& mps, MPO<Matrix, SymmGroup> const& mpo, int site = 0, bool twosite = true)
            {
                MPS<Matrix, SymmGroup> mps_aux = mps;
                Model<Matrix, SymmGroup> model(lat, parms);


                // Initialise the boundaries using the LocalHamiltonianInitializer class
                time_stopper stop_callback = static_cast<double>(parms["run_seconds"]);
                // Warning: we're casting away the const'ness of mps
                LocalHamiltonianInitializer<Matrix, SymmGroup, storage::disk> init(const_cast<MPS<Matrix, SymmGroup> &>(mps), mpo, parms, stop_callback, site);

                if (twosite)
                {
                    if (site > mps.length() - 1)
                       throw std::runtime_error("site > L-1 not allowed in Local Hamiltonian with two-site tensors");

                    int site1 = site, site2 = site+1;

                    // Prepare the two-site MPO tensor(s)

                    MPO<Matrix, SymmGroup> ts_mpo;
                    make_ts_cache_mpo(mpo, ts_mpo, mps);
                    ts_mpo[site1].placement_l = mpo[site1].placement_l;
                    ts_mpo[site1].placement_r = parallel::get_right_placement(ts_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);

                    // Prepare the two-site tensor from two sites of the MPS
                    TwoSiteTensor<Matrix, SymmGroup> tst(mps_aux[site1], mps_aux[site2]);
                    MPSTensor<Matrix, SymmGroup> twin_mps = tst.make_mps();


                    // Fetch the boundaries from disk
                    init.fetch_boundaries(site, twosite);

                    // get the local Hamiltonian matrix elements

                    // Currently doing it very inefficiently by explicitly calculating each <i|H|j> matrix element
                    // or actually H|j>[i][j] (this is only half-inefficient, as there's one contraction instead of two)

                    // Prepare the basis state from twin_mps: set all elements to 0
                    twin_mps.multiply_by_scalar(0.0);

                    for (int i = 0; i < twin_mps.data().n_blocks(); i++)
                    for (int j = 0; j < num_rows(twin_mps.data()[i]); j++)
                    for (int k = 0; k < num_cols(twin_mps.data()[i]); k++)
                    {
                        // set the corresponding element in the basis MPS to 1
                        twin_mps.data()[i](j,k) = 1.0;

                        // Calculate H|I>
                        // TODO: The Hamiltonian should be better calculated by calculating the tensor
                        // H^{sigma_L sigma'_L}_{a_{l-1}a_la'_{l-1}a'_l} = L_{a_{l-1}a'_{l-1}b_{l-1}} W_{sigma_L sigma'_L}_{b_{l-1}b_l} R_{a_l a'_l b_l}}
                        // TODO: Look into this!
                        MPSTensor<Matrix, SymmGroup> Hij_mps = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(twin_mps, init.left()[site1], init.right()[site2+1], ts_mpo[site1]);
                        twin_mps.make_left_paired(); // Contraction destroys pairing and for some reason there's a race condition somewhere later if we remove this line!

                        block_matrix<Matrix, SymmGroup> & Hij = Hij_mps.data();

                        // TODO: Check if dim(i) == dim(ii), dim(j) == dim(jj) and dim(k) == dim(kk) !!!
                        for (int ii = 0; ii < Hij.n_blocks(); ii++)
                        for (int jj = 0; jj < num_rows(Hij[ii]); jj++)
                        for (int kk = 0; kk < num_cols(Hij[ii]); kk++)
                        {
                            // prepare labels
                            std::vector<Lattice::pos_t> labels = { i, j, k, ii, jj, kk };
                            std::string label_string = label_string_simple(labels);

                            // TODO: use std::vector::reserve() in a reasonable way
                            // Write back matrix elements of MPSTensor H|I> which results in <I|H|J>
                            this->vector_results.push_back(Hij[ii](jj,kk));
                            this->labels.push_back(label_string);
                        }
                        // set the corresponding element in the basis MPS again to 0
                        twin_mps.data()[i](j,k) = 0.0;
                    }
                }
                else
                {
                    throw std::runtime_error("One-site local Hamiltonian not implemented yet! O_o");
                }

            }
        private:
            Lattice const & lat;
            BaseParameters & parms;
    };
} // measurements

#endif