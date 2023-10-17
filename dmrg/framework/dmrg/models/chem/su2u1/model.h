/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef QC_SU2_H
#define QC_SU2_H

#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <regex>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/chem/util.h"
#include "dmrg/models/chem/parse_integrals.h"
#include "dmrg/models/chem/pg_util.h"
#include "dmrg/models/chem/su2u1/chem_helper.h"
#include "dmrg/models/chem/su2u1/term_maker.h"

template<class Matrix, class SymmGroup>
class qc_su2 : public model_impl<Matrix, SymmGroup>
{
    typedef model_impl<Matrix, SymmGroup> base;

    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;

    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

    typedef typename Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;
    typedef typename alps::numeric::associated_one_matrix<Matrix>::type one_matrix;

public:

    qc_su2(Lattice const & lat_, BaseParameters & parms_);

    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }

    void create_terms();

    // For this model: site_type == point group irrep
    Index<SymmGroup> const & phys_dim(size_t type) const
    {
        return phys_indices[type];
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ops.ident[type];
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return ops.fill[type];
    }

    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms_) const
    {
        return chem::detail::qn_helper<SymmGroup>().total_qn(parms_);
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "ident_full")
            return ops.ident_full[type];

        else if (name == "create")
            return ops.create[type];
        else if (name == "create_fill")
            return ops.create_fill[type];
        else if (name == "create_couple_up")
            return ops.create_couple_up[type];
        else if (name == "create_fill_couple_down")
            return ops.create_fill_couple_down[type];

        else if (name == "destroy")
            return ops.destroy[type];
        else if (name == "destroy_fill")
            return ops.destroy_fill[type];
        else if (name == "destroy_couple_up")
            return ops.destroy_couple_up[type];
        else if (name == "destroy_fill_couple_down")
            return ops.destroy_fill_couple_down[type];

        else if (name == "count")
            return ops.count[type];
        else if (name == "count_fill")
            return ops.count_fill[type];

        else if (name == "docc")
            return ops.docc[type];
        else
            throw std::runtime_error("Operator not valid for this model: " + name);
        return 0;
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }

    measurements_type measurements () const
    {
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

        measurements_type meas;

        typedef std::vector<tag_type> tag_vec;
        typedef std::vector<tag_vec> bond_element;

        std::regex expression_oneptdm("^MEASURE\\[1rdm\\]");
        std::regex expression_twoptdm("^MEASURE\\[2rdm\\]");
        std::regex expression_transition_oneptdm("^MEASURE\\[trans1rdm\\]");
        std::regex expression_transition_twoptdm("^MEASURE\\[trans2rdm\\]");
        std::smatch what;
        std::vector<pos_t> positions;
        std::string bra_ckp("");

        for (auto&& it : parms.get_range()) {
            bool expr_rdm = false;

            std::string lhs = it.first;

            std::string name, value;
            if (std::regex_match(lhs, what, expression_twoptdm)) {

                name = "twoptdm";
                expr_rdm = true;
            }

            if (std::regex_match(lhs, what, expression_transition_twoptdm)) {

                name = "transition_twoptdm";
                bra_ckp = it.second;
                expr_rdm = true;
            }

            if (std::regex_match(lhs, what, expression_oneptdm))
            {

                name = "oneptdm";
                expr_rdm = true;
            }

            if (std::regex_match(lhs, what, expression_transition_oneptdm)) {

                name = "transition_oneptdm";
                bra_ckp = it.second;
                expr_rdm = true;
            }

            if (expr_rdm)
                meas.push_back( new measurements::TaggedNRankRDM<Matrix, SymmGroup>(
                                name, lat, tag_handler, op_collection, positions, bra_ckp));
        }

        return meas;
    }

private:
    Lattice const & lat;
    BaseParameters & parms;
    std::vector<Index<SymmGroup> > phys_indices;

    std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    typename TermMakerSU2<Matrix, SymmGroup>::Operators ops;
    typename TermMakerSU2<Matrix, SymmGroup>::OperatorCollection op_collection;

    typename SymmGroup::subcharge max_irrep;

};

#include "dmrg/models/chem/su2u1/model.hpp"

#endif
