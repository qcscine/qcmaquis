/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_LOCAL_AT_H
#define MEASUREMENTS_LOCAL_AT_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/super_mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

namespace measurements {

template <class Matrix, class SymmGroup>
class local_at : public measurement<Matrix, SymmGroup> {
public:
    using base = measurement<Matrix, SymmGroup>;
    using op_t = typename base::op_t;
    using pos_t = Lattice::pos_t;
    using op_vec = std::vector<op_t>;
    using positions_type = std::vector<std::vector<pos_t> >;

    local_at(std::string const& name_, const Lattice & lat, positions_type const& positions_,
             op_vec const & identities_, op_vec const & fillings_, std::vector<std::pair<op_vec, bool> > const& ops_)
    : base(name_), lattice(lat), positions(positions_), identities(identities_), fillings(fillings_), ops(ops_)
    {
        this->labels = label_strings(lattice, positions);
        this->labels_num = positions;
        this->cast_to_real = false;
    }

    void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
    {
        this->vector_results.clear();
        this->vector_results.reserve(positions.size());
        for (pos_t p = 0; p < positions.size(); ++p)
        {
            assert( positions[p].size() == ops.size() );
            for (pos_t i=1; i<ops.size(); ++i)
                if (positions[p][i-1] >= positions[p][i])
                    throw std::runtime_error("measure_local_at requires i1<i2<...<in.");
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lattice, identities, fillings);
            generate_mpo::OperatorTerm<Matrix, SymmGroup> hterm;
            bool with_sign = false;
            for (std::size_t i=0; i<ops.size(); ++i) 
            {
                pos_t pos = positions[p][i];
                op_t const& fill  = fillings[lattice.get_prop<int>("type", pos)];
                op_t const& op    = ops[i].first[lattice.get_prop<int>("type", pos)];
                op_t tmp;
                if (!with_sign && ops[i].second) gemm(fill, op, tmp);
                else                             tmp = op;
                hterm.operators.push_back( std::make_pair(pos, tmp) );
                pos++;
                with_sign = (ops[i].second) ? !with_sign : with_sign;
                if (i != ops.size()-1)
                    for (; pos<positions[p][i+1]; ++pos) {
                        op_t const& fill  = fillings[lattice.get_prop<int>("type", pos)];
                        op_t const& ident = identities[lattice.get_prop<int>("type", pos)];
                        hterm.operators.push_back( std::make_pair(pos, (with_sign) ? fill : ident) );
                    }
            }
            hterm.with_sign = false; // filling already taken care above
            mpom.add_term(hterm);
            MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
            if (!this->is_super_meas){
                typename MPS<Matrix, SymmGroup>::scalar_type val = expval(mps, mpo);
                this->vector_results.push_back(val);
            } else {
                typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, this->phys_psi);
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, this->phys_psi);
                // static_cast needed for icpc 12.x
                typedef typename MPS<Matrix, SymmGroup>::scalar_type (*overlap_func)(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &);
                typename MPS<Matrix, SymmGroup>::scalar_type val = ::overlap(super_mpo, mps);
                this->vector_results.push_back(val/nn);
            }
        }
    }

protected:
    measurement<Matrix, SymmGroup>* do_clone() const { return new local_at(*this); }

private:
    Lattice lattice;
    positions_type positions;
    op_vec identities, fillings;
    std::vector<std::pair<op_vec, bool> > ops;
};

}

#endif
