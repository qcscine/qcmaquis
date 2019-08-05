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
#include "dmrg/mp_tensors/xvector.h"
#include "dmrg/models/measurements/local_hamiltonian_initializer.h"
#include "dmrg/utils/time_stopper.h"

namespace measurements
{
    // Class to measure local Hamiltonian matrix elements
    template<class Matrix, class SymmGroup>
    class LocalHamiltonian : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        enum H_measurement_type { HamiltonianDiag, SigmaVector };
        public:
            LocalHamiltonian(std::string name_, Lattice const & lat_, BaseParameters & parms_,
            // parameters for loading and saving XVectors
            const std::string& xvec_filename, // path to the xvector H5 file
            const std::string& xvec_aux_filename, // path to the output text file (for MOLCAS)
            const std::string& xvec_aux_input // optional -- path to the input text file from MOLCAS
            )
            : base(name_), lat(lat_), parms(parms_),
              xvec_filename_(xvec_filename),
              xvec_aux_filename_(xvec_aux_filename),
              xvec_aux_input_(xvec_aux_input)
            {}
            virtual void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
            {
                throw std::runtime_error("Local Hamiltonian must be measured with the Hamiltonian provided as an MPO!");
            }

            virtual void evaluate(MPS<Matrix, SymmGroup> const& mps, MPO<Matrix, SymmGroup> const& mpo)
            {
                this->vector_results.clear();
                this->labels.clear();
                if (this->name() == "local_hamiltonian_diag")
                    measure(HamiltonianDiag, mps, mpo);
                else if (this->name() == "sigma_vector")
                    measure(SigmaVector, mps, mpo);
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
            void measure(H_measurement_type meas_type, MPS<Matrix, SymmGroup> const& mps, MPO<Matrix, SymmGroup> const& mpo)
            {
                typedef typename Matrix::value_type value_type;

                // Flag if we should employ the MPS variation in the calculation of the sigma vectors
                bool use_variation = !xvec_aux_input_.empty();

                // Initialise the boundaries using the LocalHamiltonianInitializer class
                time_stopper stop_callback = static_cast<double>(parms["run_seconds"]);

                // LocalHamiltonianInitializer (or actually optimizer_base) requires a vector<MPS> reference to initialise
                // for now, we copy-initialise the vector

                std::vector<MPS<Matrix, SymmGroup> > mpsv(1, mps);
                MPS<Matrix, SymmGroup> & mps_ = mpsv[0];

                // prepare the x vector
                // small lambda function to construct the X vector and generate the MPS variation -- so that in case we don't need it at all, we avoid the initialisation
                auto generate_mps_variation = [&](const std::string & xvec_filename, const std::string & xvec_aux_input, MPS<Matrix, SymmGroup>& mps) -> MPS<Matrix,SymmGroup>
                {
                    lr::XVector<Matrix, SymmGroup> xvec;
                    xvec.load(xvec_filename);
                    xvec.assign_mps(mps);
                    xvec.update_from_textfile(xvec_aux_input);
                    return xvec.GenerateMPSVariation();
                };

                // prepare the bra MPS. If we don't use the variation, it's identical to ket MPS, otherwise it is the variation
                // obtained from the x vector
                const MPS<Matrix, SymmGroup> & bra_mps = (use_variation) ? generate_mps_variation(xvec_filename_, xvec_aux_input_, mps_) : mps_;

                LocalHamiltonianInitializer<Matrix, SymmGroup, storage::disk> init(mpsv, mpo, parms, stop_callback, mps_);

                mps_.normalize_right();

                // Output sigma vector/ Hamiltonian diagonal in the form of a MPS.
                // Initialise with the same dimensions/indices as the ket MPS. The easiest way to do this is to copy-initialise and then replace the data
                MPS<Matrix, SymmGroup> sigma_vector_out = mps_;

                for (std::size_t site = 0; site < mps_.length(); site++)
                {
                    // init.fetch_boundaries(site); // <-- this will be required if we save the boundaries to disk. to be considered later

                    MPSTensor<Matrix, SymmGroup>& ret_site = sigma_vector_out[site]; // output MPSTensor

                    if (meas_type == HamiltonianDiag)
                    {
                        maquis::cout << "Measuring local Hamiltonian diagonal at site " << site << std::endl;
                        block_matrix<Matrix, SymmGroup> && Hdiag = contraction::diagonal_hamiltonian(init.left()[site], init.right()[site+1], mpo[site], mps_[site]);
                        ret_site.replace_left_paired(Hdiag); // TODO: check the pairing that comes out of diagonal_hamiltonian()
                        // TODO2: implement move semantics
                    }
                    else if (meas_type == SigmaVector)
                    {
                        maquis::cout << "Measuring sigma vector at site " << site << std::endl;
                        ret_site = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(mps_[site], init.left()[site], init.right()[site+1], mpo[site]);
                    }
                    else
                        throw std::runtime_error("Don't know this measurement type for Local Hamiltonian.");

                }

                // Now we have the sigma vector/diagonal Hamiltonian as the MPSTensor. Now convert it into an XVector
                lr::XVector<Matrix, SymmGroup> xvec_out(sigma_vector_out, mps_);
                xvec_out.save(xvec_filename_);
                xvec_out.dump_to_textfile(xvec_aux_filename_);

                // Test: transforming the xvector to MPS must yield sigma_vector_out

                MPS<Matrix, SymmGroup> && test_sigmavec = xvec_out.transformXtoB();

                // dump to HDF5 to make measurements happy, because labels and vector_results may not be empty
                for (int site = 0; site < mps_.length(); site++)
                {
                    parallel::scheduler_balanced_iterative scheduler(xvec_out[site]);

                    for (int i = 0; i < xvec_out[site].n_blocks(); i++)
                    for (int j = 0; j < xvec_out[site][i].num_rows(); j++)
                    for (int k = 0; k < xvec_out[site][i].num_cols(); k++)
                    {
                        parallel::guard proc(scheduler(i));
                        this->vector_results.push_back(xvec_out[site][i](j,k));
                        // Labels are dumped as 'site, i, j, k'
                        this->labels.push_back(label_string_simple({site, i, j, k}));
                    }
                }
            }
        private:
            Lattice const & lat;
            BaseParameters & parms;
            std::string xvec_filename_; // path to the xvector H5 file
            std::string xvec_aux_filename_; // path to the output text file (for MOLCAS)
            std::string xvec_aux_input_; // optional -- path to the input text file from MOLCAS

    };
} // measurements

#endif