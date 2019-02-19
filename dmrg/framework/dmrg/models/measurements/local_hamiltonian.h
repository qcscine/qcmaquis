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
    // Class to measure local Hamiltonian matrix elements
    template<class Matrix, class SymmGroup>
    class LocalHamiltonian : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        enum H_measurement_type { HamiltonianMatrix, HamiltonianMatrixDiag, SigmaVector };
        public:
            LocalHamiltonian(std::string name_, Lattice const & lat_, BaseParameters & parms_, const std::string & ext_filename_=std::string(""))
            : base(name_), lat(lat_), parms(parms_), ext_filename(ext_filename_), twosite(parms["lrparam_twosite"]) {}
            virtual void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
            {
                throw std::runtime_error("Local Hamiltonian must be measured with the Hamiltonian provided as an MPO!");
            }

            virtual void evaluate(MPS<Matrix, SymmGroup> const& mps, MPO<Matrix, SymmGroup> const& mpo)
            {
                this->vector_results.clear();
                this->labels.clear();
                if (this->name() == "local_hamiltonian")
                    measure(HamiltonianMatrix, mps, mpo, parms["lrparam_site"]);
                else if (this->name() == "local_hamiltonian_diag")
                    measure(HamiltonianMatrixDiag, mps, mpo, parms["lrparam_site"]);
                else if (this->name() == "sigma_vector")
                    measure(SigmaVector, mps, mpo, parms["lrparam_site"]);
                else
                    throw std::runtime_error("Don't know this measurement type for Local Hamiltonian.");
            }
        protected:
            measurement<Matrix, SymmGroup>* do_clone() const
            {
                return new LocalHamiltonian(*this);
            }

            // Measures local Hamiltonian matrix elements H_{ij} if meas_type = HamiltonianMatrix
            // or \sum_j H_ijc_j with c_j = M^sigma_l_{a_{l-1}a_l} (or its twosite equivalent) if meas_type = SigmaVector
            void measure(H_measurement_type meas_type, MPS<Matrix, SymmGroup> const& mps, MPO<Matrix, SymmGroup> const& mpo, int site = 0)
            {
                typedef typename Matrix::value_type value_type;

                // MPSTensor that we're going to work with
                MPSTensor<Matrix, SymmGroup> mpst;
                // Initialise the boundaries using the LocalHamiltonianInitializer class
                time_stopper stop_callback = static_cast<double>(parms["run_seconds"]);

                // LocalHamiltonianInitializer (or actually optimizer_base) requires a vector<MPS> reference to initialise
                // for now, we copy-initialise the vector
                // TODO: think how can it be initialised with move semantics!
                std::vector<MPS<Matrix, SymmGroup> > mpsv(1, mps);

                LocalHamiltonianInitializer<Matrix, SymmGroup, storage::disk> init(mpsv, mpo, parms, stop_callback, site);

                int site1, site2;

                // Two-site MPO tensor
                MPO<Matrix, SymmGroup> ts_mpo;

                // Prepare either two-site tensor or one-site tensor depending on which kind of MPS parameters we work with
                if (twosite)
                {
                    if (site > mps.length() - 1)
                       throw std::runtime_error("site > L-1 not allowed in Local Hamiltonian with two-site tensors");

                    site1 = site;
                    site2 = site+1;
                    // Prepare the two-site MPO tensor(s)
                    make_ts_cache_mpo(mpo, ts_mpo, mps);
                    ts_mpo[site1].placement_l = mpo[site1].placement_l;
                    ts_mpo[site1].placement_r = parallel::get_right_placement(ts_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);

                    // Prepare the two-site tensor from two sites of the MPS
                    TwoSiteTensor<Matrix, SymmGroup> tst(mps[site1], mps[site2]);
                    mpst = tst.make_mps();
                }
                else // one-site
                {
                    mpst = mps[site];
                    site1 = site2 = site;
                }

                const MPO<Matrix, SymmGroup> & mpo_work = twosite ? ts_mpo : mpo;

                mpst.make_left_paired();

                // Fetch the boundaries from disk
                init.fetch_boundaries(site, twosite);

                if (meas_type == HamiltonianMatrix)
                {
                    // Measure Local Hamiltonian matrix elements
                    // Currently doing it very inefficiently by explicitly calculating each <i|H|j> matrix element
                    // or actually H|j>[i][j] (this is only half-inefficient, as there's one contraction instead of two)
                    // This is DEPRECATED and will be removed soon

                    // Prepare the basis state from mpst: set all elements to 0
                    mpst.multiply_by_scalar(0.0);

                    if (meas_type == HamiltonianMatrix)
                        maquis::cout << "Measuring local Hamiltonian matrix at site " << site << std::endl;
                    else if (meas_type == HamiltonianMatrixDiag)
                        maquis::cout << "Measuring local Hamiltonian diagonal at site " << site << std::endl;

                    parallel::scheduler_balanced_iterative scheduler(mpst.data());

                    for (int i = 0; i < mpst.data().n_blocks(); i++)
                    for (int j = 0; j < num_rows(mpst.data()[i]); j++)
                    for (int k = 0; k < num_cols(mpst.data()[i]); k++)
                    {
                        parallel::guard proc(scheduler(i));
                        // set the corresponding element in the basis MPS to 1
                        mpst.data()[i](j,k) = 1.0;

                        // Calculate H|I>
                        // TODO: The Hamiltonian should be better calculated by calculating the tensor
                        // H^{sigma_L sigma'_L}_{a_{l-1}a_la'_{l-1}a'_l} = L_{a_{l-1}a'_{l-1}b_{l-1}} W_{sigma_L sigma'_L}_{b_{l-1}b_l} R_{a_l a'_l b_l}}
                        // TODO: Look into this!
                        MPSTensor<Matrix, SymmGroup> Hij_mps = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(mpst, init.left()[site1], init.right()[site2+1], mpo_work[site1]);
                        mpst.make_left_paired(); // Contraction destroys pairing and for some reason there's a race condition somewhere later if we remove this line!

                        block_matrix<Matrix, SymmGroup> & Hij = Hij_mps.data();

                        // TODO: Check if dim(i) == dim(ii), dim(j) == dim(jj) and dim(k) == dim(kk) !!!
                        // TODO: The Hamiltonian matrix is symmetric. Currently, we're not exploiting the symmetry, but one should do it!

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
                        mpst.data()[i](j,k) = 0.0;
                    }

                }
                else if (meas_type == HamiltonianMatrixDiag)
                {
                    // Use contraction::su2::diagonal_hamiltonian() to obtain the Hamiltonian diagonal
                    maquis::cout << "Measuring local Hamiltonian diagonal at site " << site << std::endl;

                    block_matrix<Matrix, SymmGroup> Hdiag = contraction::diagonal_hamiltonian(init.left()[site1], init.right()[site2+1], mpo_work[site1], mpst);

                    // dump the Hamiltonian diagonal into the result vector
                    for (int ii = 0; ii < Hdiag.n_blocks(); ii++)
                    for (int jj = 0; jj < num_rows(Hdiag[ii]); jj++)
                    for (int kk = 0; kk < num_cols(Hdiag[ii]); kk++)
                    {
                        // prepare labels: double the labels to ensure compatibility with the old python script
                        // TODO: this should not be needed in the future
                        std::vector<Lattice::pos_t> labels = { ii, jj, kk, ii, jj, kk };
                        std::string label_string = label_string_simple(labels);

                        // TODO: use std::vector::reserve() in a reasonable way
                        // Write the matrix elements
                        this->vector_results.push_back(Hdiag[ii](jj,kk));
                        this->labels.push_back(label_string);
                    }
                }
                else if (meas_type == SigmaVector)
                {
                    // Read the MPSTensor from an external file if requested
                    maquis::cout << "Measuring local Hamiltonian sigma vector at site " << site << std::endl;
                    if (!ext_filename.empty())
                    {
                        // Read MPSTensor elements into a vector aux_elements from a text file

                        std::vector<value_type> aux_elements;
                        maquis::cout << "Reading the auxiliary MPSTensor elements from file " << ext_filename << std::endl;

                        // read and parse the file
                        std::ifstream infile(ext_filename);
                        if (infile)
                            std::copy(std::istream_iterator<value_type>(infile), std::istream_iterator<value_type>(), std::back_inserter(aux_elements));
                        else
                            throw std::runtime_error("File " + ext_filename + " could not be opened!");

                        // fill the MPSTensor with elements from file
                        size_t fileidx = 0;
                        for (size_t i = 0; i < mpst.data().n_blocks(); i++)
                        for (size_t j = 0; j < mpst.data()[i].num_rows(); j++)
                        for (size_t k = 0; k < mpst.data()[i].num_cols(); k++)
                        {
                            parallel::guard::serial guard;
                            mpst.data()[i](j,k) = aux_elements[fileidx++];
                        }
                    }

                    // Measure the sigma_vector \sum_j H_ijc_j with c_j = M^sigma_l_{a_{l-1}a_l}
                    // i.e. just dump the output of site_hamil2()
                    MPSTensor<Matrix, SymmGroup> Hij_mps = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(mpst, init.left()[site1], init.right()[site2+1], mpo_work[site1]);
                    mpst.make_left_paired(); // Contraction destroys pairing and for some reason there's a race condition somewhere later if we remove this line!

                    // save the result
                    for (int i = 0; i < Hij_mps.data().n_blocks(); i++)
                    for (int j = 0; j < Hij_mps.data()[i].num_rows(); j++)
                    for (int k = 0; k < Hij_mps.data()[i].num_cols(); k++)
                    {
                        this->vector_results.push_back(Hij_mps.data()[i](j,k));
                        // Labels are dumped as 'i, j, k'
                        this->labels.push_back(label_string_simple({i, j, k}));
                    }
                }
                else
                    throw std::runtime_error("Don't know this measurement type for Local Hamiltonian.");
            }
        private:
            Lattice const & lat;
            std::string ext_filename; // Filename from where the external MPSTensor should be read if a sigma vector with external MPSTensor is requested
            BaseParameters & parms;
            bool twosite;

    };
} // measurements

#endif