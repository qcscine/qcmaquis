/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include <mpi.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <alps/hdf5.hpp>

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> amatrix;

#ifdef USE_AMBIENT
#include "dmrg/block_matrix/detail/ambient.hpp"
typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > pmatrix;
typedef pmatrix matrix;
#else
typedef amatrix matrix;
#endif

#include "dmrg/block_matrix/symmetry.h"
#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_U1)
typedef U1 grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_NU1)
typedef NU1 grp;
#else
#error "No valid symmetry defined."
#endif

#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/twositetensor.h"

#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/optimize/ietl_jacobi_davidson.h"

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

#include "utils/timings.h"


#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

class triag_range {
public:
    typedef std::size_t size_type;
    
    
    triag_range(size_type length, size_type nprocs=1, size_type rank=0)
    : extent_(length)
    {
        double half_chunk = length / 2. / nprocs;
        
        size_type start = half_chunk*rank;
        size_type end   = std::min(size_type(half_chunk*(rank+1)), length / 2);
        
        size_type odd_shift = length % 2;
        indexes_.reserve(2*half_chunk);
        for (size_type i=start; i<end; ++i) {
            indexes_.push_back(i);
            indexes_.push_back(length-1 - i - odd_shift);
        }
        
        if (rank == 0 && length % 2 != 0)
            indexes_.push_back(length-1);
    }
    
    inline size_type extent() const { return extent_; }
    inline size_type size() const { return indexes_.size(); }
    inline size_type operator[](size_type i) const {return indexes_[i]; }

private:
    size_type extent_;
    std::vector<size_type> indexes_;
};


struct corr_measurement {
    typedef block_matrix<matrix, grp> op_t;
    
    std::string name;
    std::vector<std::vector<op_t> > ops;
    bool fermionic;
};

template<class Matrix, class SymmGroup, class Range>
void measure_correlation(Range const& range,
                         MPS<Matrix, SymmGroup> const & bra,
                         MPS<Matrix, SymmGroup> const & ket,
                         Lattice const & lattice,
                         std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
                         std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
                         corr_measurement const & meas,
                         std::vector<Boundary<Matrix,SymmGroup> > const & left,
                         std::vector<Boundary<Matrix,SymmGroup> > const & right,
                         std::vector<typename Matrix::value_type>& dc,
                         std::vector<std::string>& labels)
{
    std::size_t L = lattice.size();
    assert( range.extent() == L-1 );
    
    int d = 1;
    
    MPOTensor<Matrix, SymmGroup> mpo_fill(1,1), mpo_ident(1,1);
    
    std::size_t nops = meas.ops.size();
    int ntypes = lattice.maximum_vertex_type()+1;
    std::vector<std::vector<MPOTensor<Matrix, SymmGroup> > > mpo_ops(nops), mpo_signed_ops(nops);
    for (std::size_t i=0; i<nops; ++i) {
        mpo_ops[i].resize(ntypes);
        mpo_signed_ops[i].resize(ntypes);
        for (int type=0; type<ntypes; ++type) {
            mpo_ops[i][type] = MPOTensor<Matrix, SymmGroup>(1,1);
            block_matrix<Matrix, SymmGroup> tmp;
            gemm(fillings[type], meas.ops[i][type], tmp);
            mpo_ops[i][type].set(0, 0, meas.ops[i][type]);
            mpo_signed_ops[i][type].set(0, 0, tmp);
        }
    }
    
    for (size_t i = 0; i < range.size(); ++i) {        
        const size_t p = range[i];
        if (p >= L-d-1) continue;
        
        Timer tim("measure p="+boost::lexical_cast<std::string>(p)); tim.begin();
        
        std::vector<std::vector<typename Lattice::pos_t> > numeric_labels;
        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
        
        int type_p = lattice.get_prop<int>("type", p);
        Boundary<Matrix, SymmGroup> current = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(bra[p], ket[p], left[p], mpo_signed_ops[0][type_p]);
        
        
        {
            current = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(bra[p+d], ket[p+d], current, mpo_ops[1][lattice.get_prop<int>("type", p+d)]);
            
            for (int q = p+d+1; q < L-d; ++q) {
                int type_q = lattice.get_prop<int>("type", q);
                
                Boundary<Matrix, SymmGroup> current2 = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(bra[q], ket[q], current, mpo_signed_ops[2][type_q]);
                for (int r = q+1; r < q+d; ++r) {
                    int type_r = lattice.get_prop<int>("type", r);
                    mpo_fill.set(0,0, fillings[type_r]);
                    current2 = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(bra[r], ket[r], current2, mpo_fill);
                }
                current2 = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(bra[q+d], ket[q+d], current2, mpo_ops[3][lattice.get_prop<int>("type", q+d)]);
                block_matrix<Matrix, SymmGroup> const& vec = current2[0];
                
                double obs = 0.;
                for (int k=0; k<vec.n_blocks(); ++k) {
                    std::size_t matched_block = right[q+d+1][0].find_block(vec.right_basis()[k].first, vec.left_basis()[k].first);
                    if ( matched_block <  right[q+d+1][0].n_blocks() ) {
                        Matrix m = transpose( right[q+d+1][0][matched_block] );
                        obs += overlap(m, vec[k]);
                    }
                }
                
                std::vector<typename Lattice::pos_t> lab(4);
                lab[0] = p; lab[1] = p+d;
                lab[2] = q; lab[3] = q+d;
                
                dct.push_back( obs );
                numeric_labels.push_back(lab);
                
                {
                    mpo_fill.set(0,0, fillings[type_q]);
                    mpo_ident.set(0,0, identities[type_q]);
                    current = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(bra[q], ket[q], current, mpo_ident);
                }
            }
        }
        
        
        std::copy(dct.begin(), dct.end(), std::back_inserter(dc));
        
        std::vector<std::string> lbt = label_strings(lattice, numeric_labels);
        std::copy(lbt.begin(), lbt.end(), std::back_inserter(labels));
        tim.end();
    }
}

int main(int argc, char ** argv)
{
    try {
        mpi::environment env(argc, argv);
        mpi::communicator comm;
        
        DmrgOptions opt(argc, argv);
        if (!opt.valid) return 0;
        DmrgParameters parms = opt.parms;

        maquis::cout.precision(10);
        
        typedef matrix::value_type value_type;
        typedef block_matrix<matrix, grp> op_t;
        typedef Model<matrix, grp>::table_ptr table_ptr;
        typedef Model<matrix, grp>::tag_type tag_type;
        
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

        /// Parsing model
        Lattice lattice(parms);
        int L = lattice.size();
        int ntypes = lattice.maximum_vertex_type()+1;
        Model<matrix, grp> model(lattice, parms);
        table_ptr tag_handler = model.operators_table();
        
        /// identities and fillings
        std::vector<op_t> identities(ntypes), fillings(ntypes);
        for (size_t type = 0; type < ntypes; ++type) {
            identities[type] = model.identity_matrix(type);
            fillings[type]   = model.filling_matrix(type);
        }
        
        /// parsing measurements
        std::vector<corr_measurement> measurements;
        {
            boost::regex expression("^MEASURE_HALF_NN_CORRELATIONS\\[(.*)]$");
            boost::smatch what;
            for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
                std::string lhs = it->key();
                
                if (boost::regex_match(lhs, what, expression)) {
                    std::string value = it->value();

                    corr_measurement meas;
                    meas.name = what.str(1);
                    
                    
                    int f_ops = 0;
                    std::vector<std::string> corr_tokens;
                    boost::split(corr_tokens, value, boost::is_any_of(":"));
                    
                    if (corr_tokens.size() != 4)
                        throw std::runtime_error("`pcorr_nn` can compute only two-body correlation functions.");
                    
                    meas.ops.resize(corr_tokens.size());
                    
                    int ii;
                    std::vector<std::string>::const_iterator it2;
                    for (it2=corr_tokens.begin(), ii=0; it2 != corr_tokens.end(); it2++, ++ii)
                    {
                        // int ii = std::distance(corr_tokens.begin(), it2);
                        meas.ops[ii].resize(ntypes);
                        
                        enum {uknown, bosonic, fermionic} kind = uknown;
                        
                        for (int type=0; type<ntypes; ++type) {
                            tag_type tag = model.get_operator_tag(*it2, type);
                            meas.ops[ii][type] = tag_handler->get_op(tag);
                            
                            bool is_ferm = tag_handler->is_fermionic(tag);
                            if (kind == uknown)
                                kind = is_ferm ? fermionic : bosonic;
                            else if ((is_ferm && kind==bosonic) || (!is_ferm && kind==fermionic))
                                throw std::runtime_error("Model is inconsitent. On some site the operator " + *it2 + "fermionic, on others is bosonic.");
                        }
                        
                        if (kind == fermionic) ++f_ops;
                    }
                    
                    if (f_ops % 2 != 0)
                        throw std::runtime_error("Number of fermionic operators has to be even.");
                    
                    meas.fermionic = (f_ops > 0);
                    measurements.push_back(meas);
                }
            }
        }
        
        /// Load MPS
        MPS<matrix, grp> mps;
        load(parms["chkpfile"].str(), mps);

        unsigned nprocs = comm.size();
        unsigned rank   = comm.rank();
        triag_range myrange(L-1, nprocs, rank);
        
        maquis::cout << "[" << rank << "] will do " << myrange.size() << " steps." << std::endl;
        
        /// compute left / right boundaries
        /// note: this is the current bottleneck. this part has to run sequentially!
        std::vector<Boundary<matrix,grp> > right(L+1, Boundary<matrix,grp>(Index<grp>(), Index<grp>(), 1));
        std::vector<Boundary<matrix,grp> > left (L+1, Boundary<matrix,grp>(Index<grp>(), Index<grp>(), 1));
        
        if (rank == 0) {
            right[L] = make_right_boundary(mps, mps);
            for (int p=L-1; p>=static_cast<int>(0); --p)
                right[p][0] = contraction::Engine<matrix, matrix, grp>::overlap_right_step(mps[p], MPSTensor<matrix, grp>(mps[p]), right[p+1][0]);
        }
        
        if (rank == nprocs-1) {
            MPOTensor<matrix, grp> mpo_ident(1,1);
            left[0] = make_left_boundary(mps, mps);
            for (int p=0; p<L; ++p) {
                mpo_ident.set(0,0, identity_matrix<matrix>(mps[p].site_dim()) );
                left[p+1] = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(mps[p], mps[p], left[p], mpo_ident);
            }
        }
        
        if (L/2 > 96) {
            for (int p=0; p<L+1; ++p) {
                mpi::broadcast(comm, right[p], 0       );
                mpi::broadcast(comm, left[p],  nprocs-1);
            }
        } else {
            mpi::broadcast(comm, right, 0       );
            mpi::broadcast(comm, left,  nprocs-1);
        }
        
        /// Compute all measurements in local range
        for (std::vector<corr_measurement>::const_iterator it = measurements.begin();
             it != measurements.end(); ++it) {
            
            if (rank == 0) maquis::cout << "Measure " << it->name << "." << std::endl;
            
            std::vector<value_type> values;
            std::vector<std::string> labels;
            measure_correlation(myrange, mps, mps, lattice, identities, (it->fermionic) ? fillings : identities, *it, left, right, values, labels);

            
            if (rank != 0) {
                comm.send(0, 10, values);
                comm.send(0, 11, labels);
            } else {
                
                for (int r=1; r<nprocs; ++r) {
                    std::vector<value_type>  tmp1;
                    std::vector<std::string> tmp2;
                    
                    comm.recv(r, 10, tmp1);
                    comm.recv(r, 11, tmp2);
                    
                    values.reserve(values.size() + tmp1.size());
                    labels.reserve(labels.size() + tmp2.size());
                    std::copy(tmp1.begin(), tmp1.end(), back_inserter(values));
                    std::copy(tmp2.begin(), tmp2.end(), back_inserter(labels));
                }
                
                storage::archive ar(parms["resultfile"].str(), "w");
                ar["/spectrum/results/" + it->name + "/labels"    ] << labels;
                ar["/spectrum/results/" + it->name + "/mean/value"] << std::vector<std::vector<value_type> >(1, values);
            }
        }
        
        
    } catch (std::exception & e) {
        maquis::cerr << "Exception caught:" << std::endl << e.what() << std::endl;
        throw;
    }
}

